"""PoseBusters dock checks for each prepared ligand copy. Unknown never counts as a pass."""

from core_score.results import clean, diagnostic


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
    """A pose is valid only if every ligand copy passes every dock check."""
    record = dict(ensemble=pose["ensemble"], pose_id=pose["pose_id"], pb_valid=None)
    if pose.get("error"):
        return dict(record, error="Pose was not prepared; see prepared_manifest.json")
    prepared = pose["prepared"]
    copies = [check_copy(buster, checks, copy, prepared["receptor"]) for copy in prepared["copy_records"]]
    return dict(record, pb_valid=aggregate(copy["pb_valid"] for copy in copies), copies=copies)


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
    return dict(
        record,
        pb_valid=aggregate(results.values()),
        failed_checks=[name for name, passed in results.items() if passed is False],
        unknown_checks=[name for name, passed in results.items() if passed is None],
        tests=tests,
    )


def summary_row(record):
    """One CSV row per pose: validity, every failed check, and the first error if any."""
    copies = record.get("copies", [])
    return dict(
        ensemble=record["ensemble"],
        pose_id=record["pose_id"],
        pb_valid=record["pb_valid"],
        failed_checks=sorted({name for copy in copies for name in copy.get("failed_checks", [])}),
        error=record.get("error") or next((copy["error"] for copy in copies if copy.get("error")), None),
    )
