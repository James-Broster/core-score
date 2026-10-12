"""Split ensembles with several identical ligand copies into one ensemble per site."""

from core_score.results import diagnostic, write_json
from core_score.stages.utils.sites import split_sites


def prepare_site_groups(groups, output, *, timeout=None):
    """One ensemble per ligand site. Ensembles with one copy pass through as site 1."""
    split, reports = [], []
    for group in groups:
        try:
            failed = [pose["pose_id"] for pose in group["poses"] if pose.get("error")]
            if failed:
                raise ValueError(f"Poses failed preparation: {failed}; see prepared_manifest.json")
            counts = {pose["prepared"]["copies"] for pose in group["poses"]}
            if len(counts) != 1:
                # Copies can only be matched to sites when every pose has the same number.
                raise ValueError(f"Poses have different numbers of ligand copies: {sorted(counts)}")
            copies = counts.pop()
            if copies == 1:
                split.append(dict(group, site="1"))
                continue
            sites, report = split_sites(group, copies, output, timeout)
            split += sites
            reports.append(report)
        except Exception as exc:
            report = diagnostic("site_assignment", exc, ensemble=group["ensemble"])
            reports.append(report)
            split.append(dict(group, site="unresolved", error=report))
    write_json(output / "sites.json", reports)
    return split
