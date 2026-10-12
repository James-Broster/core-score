"""Prepare consistent ligand-site inputs; preserve predicted atomic coordinates."""

from pathlib import Path
from core_score.stages.utils.site_algorithms import assign_sites
from core_score.stages.utils.sites_utils import align_copy_centres, prepare_site_pose
from core_score.results import diagnostic, fingerprint, write_json


def prepare_site_groups(groups, output, *, ost="ost", timeout=None):
    """Prepare one ensemble per ligand site after PoseBusters filtering."""

    expanded, reports = [], []
    for group in groups:
        try:
            if group.get("pb_filter_error"):
                raise ValueError(group["pb_filter_error"])
            if any("prepared" not in p for p in group["poses"]):
                raise ValueError(
                    "One or more poses failed preparation; inspect prepared_manifest.json"
                )
            if not group["poses"]:
                raise ValueError("No retained poses")
            counts = {p["prepared"]["copies"] for p in group["poses"]}
            if len(counts) != 1:
                raise ValueError(
                    "Inconsistent ligand copy counts; no poses silently excluded"
                )
            copies = counts.pop()
            if copies == 1:
                expanded.append(dict(group, site="1", n_input_copies=1))
                continue
            poses = sorted(
                group["poses"], key=lambda p: (p["seed"], p["sample"], p["pose_id"])
            )
            points, transforms = align_copy_centres(poses, copies, ost, timeout)
            assignments, centres = assign_sites(points, copies)
            mapping = {a.pose_id: a.chain_to_site for a in assignments}
            reports.append(
                dict(
                    ensemble=group["ensemble"],
                    anchor_pose=poses[0]["pose_id"],
                    transforms=transforms,
                    centres=centres,
                    assignments=mapping,
                    costs={a.pose_id: a.chain_to_cost for a in assignments},
                    policy="historical_kmeans_hungarian_seed0",
                    biological_correspondence_verified=False,
                )
            )
            for site in range(copies):
                selected = []
                for pose in poses:
                    copy = next(
                        c
                        for c in pose["prepared"]["copy_records"]
                        if mapping[pose["pose_id"]][c["copy_id"]] == site
                    )
                    dest = (
                        Path(output)
                        / "prepared"
                        / fingerprint((group["ensemble"], pose["pose_id"]))[:20]
                        / f"site{site+1}.cif"
                    )
                    selected.append(prepare_site_pose(pose, copy, dest))
                expanded.append(
                    dict(group, poses=selected, site=str(site + 1), n_input_copies=1)
                )
        except Exception as exc:
            report = diagnostic("site_assignment", exc, ensemble=group["ensemble"])
            reports.append(report)
            expanded.append(dict(group, site="unresolved", stage_error=report))
    write_json(Path(output) / "sites.json", reports)
    return expanded



