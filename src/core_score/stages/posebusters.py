"""PoseBusters dock checks for every pose. The results are saved for use downstream, such as
choosing the top-ranked valid pose. CORE uses every pose whatever its result. Unknown never
counts as a pass."""

from core_score.results import clean, diagnostic, write_json


def run_posebusters(groups, output):
    from posebusters import PoseBusters

    buster = PoseBusters(config="dock", max_workers=0)
    checks = dock_checks(buster)
    records = [check_pose(buster, checks, pose) for group in groups for pose in group["poses"]]
    write_json(output / "posebusters.json", records)
    return records


def dock_checks(buster):
    """Names of the binary checks in PoseBusters' dock configuration."""
    return [
        normalized(module.get("rename_outputs", {}).get(name, name))
        for module in buster.config["modules"]
        for name in module.get("chosen_binary_test_output", [])
    ]


def normalized(name):
    return name.lower().replace(" ", "_")


def aggregate(values):
    """False if any check failed, None if any is unknown (or there are none), otherwise True."""
    values = list(values)
    if any(v is False for v in values):
        return False
    if not values or any(v is None for v in values):
        return None
    return True


def check_pose(buster, checks, pose):
    """A pose is valid only if every ligand copy passes every dock check. Each check's result
    for the pose combines its copies the same way."""
    record = dict(ensemble=pose["ensemble"], pose_id=pose["pose_id"], pb_valid=None,
                  checks=dict.fromkeys(checks))
    if pose.get("error"):
        return dict(record, error="Pose was not prepared; see details/prepared_manifest.json")
    prepared = pose["prepared"]
    copies = [check_copy(buster, checks, copy, prepared["receptor"]) for copy in prepared["copy_records"]]
    combined = {name: aggregate(copy.get("checks", {}).get(name) for copy in copies) for name in checks}
    return dict(record, pb_valid=aggregate(copy["pb_valid"] for copy in copies), checks=combined, copies=copies)


def check_copy(buster, checks, copy, receptor):
    record = dict(copy_id=copy["copy_id"], preparation=copy["provenance"], pb_valid=None)
    if copy["pb_error"]:
        return dict(record, error=copy["pb_error"])
    try:
        report = buster.bust(copy["pb_sdf"], mol_cond=receptor, full_report=True)
    except Exception as exc:
        return dict(record, **diagnostic("posebusters", exc))
    tests = {normalized(k): clean(v) for k, v in report.iloc[0].items()}
    results = {name: tests.get(name) if isinstance(tests.get(name), bool) else None for name in checks}
    return dict(record, pb_valid=aggregate(results.values()), checks=results, tests=tests)


def summary_row(record):
    """One row of posebusters.csv: validity, the failed checks by name, every check's result
    (blank when unknown), and the first error if any."""
    copies = record.get("copies", [])
    return dict(
        ensemble=record["ensemble"],
        pose_id=record["pose_id"],
        pb_valid=record["pb_valid"],
        failed_checks="; ".join(name for name, passed in record["checks"].items() if passed is False),
        **record["checks"],
        error=record.get("error") or next((copy["error"] for copy in copies if copy.get("error")), None),
    )
