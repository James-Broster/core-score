"""Run the bundled dock configuration; preserve every check and unknown result."""

import csv
from importlib.metadata import version
from core_score.results import clean, diagnostic

PAPER_CHECKS = [
    "sanitization",
    "inchi_convertible",
    "all_atoms_connected",
    "bond_lengths",
    "bond_angles",
    "internal_steric_clash",
    "aromatic_ring_flatness",
    "non-aromatic_ring_non-flatness",
    "double_bond_flatness",
    "internal_energy",
    "protein-ligand_maximum_distance",
    "minimum_distance_to_protein",
    "volume_overlap_with_protein",
]


def normalized(name):
    return name.lower().replace(" ", "_")


def boolean(value):
    if value is True or value is False:
        return value
    if hasattr(value, "item"):
        return boolean(value.item())
    if isinstance(value, str):
        if value.lower() in ("true", "1", "1.0"):
            return True
        if value.lower() in ("false", "0", "0.0"):
            return False
    return None


def aggregate(values):
    values = list(values)
    if any(v is False for v in values):
        return False
    if not values or any(v is None for v in values):
        return None
    return True


def run_pose(pose, prepared):
    from posebusters import PoseBusters

    buster = PoseBusters(config="dock", max_workers=0)
    check_names = [
        normalized(module.get("rename_outputs", {}).get(name, name))
        for module in buster.config["modules"]
        for name in module.get("chosen_binary_test_output", [])
    ]
    records = []
    for copy in prepared.get("copy_records", []):
        record = dict(
            ensemble=pose["ensemble"],
            pose_id=pose["pose_id"],
            copy_id=copy["copy_id"],
            config="dock",
            version=version("posebusters"),
            source="computed",
            preparation=copy.get("provenance"),
            pb_valid=None,
            status="unknown",
        )
        try:
            if copy.get("pb_error") or not copy.get("pb_sdf"):
                raise ValueError(
                    copy.get("pb_error") or "No chemically validated ligand SDF for PB"
                )
            report = buster.bust(
                copy["pb_sdf"], mol_cond=prepared["receptor"], full_report=True
            )
            if len(report) != 1:
                raise ValueError(
                    "PoseBusters did not return exactly one result for the ligand copy"
                )
            tests = {normalized(k): clean(v) for k, v in report.iloc[0].items()}
            checks = {k: boolean(tests.get(k)) for k in check_names}
            valid = aggregate(checks.values())
            record.update(
                tests=tests,
                checks=checks,
                pb_valid=valid,
                pb_valid_13_checks=aggregate(
                    boolean(tests.get(k)) for k in PAPER_CHECKS
                ),
                failed_checks=[k for k, v in checks.items() if v is False],
                missing_checks=[k for k, v in checks.items() if v is None],
                status=(
                    "complete"
                    if all(v is not None for v in checks.values())
                    else "incomplete_checks"
                ),
            )
        except Exception as exc:
            record.update(diagnostic("posebusters", exc, structure=pose["path"]))
            record["status"] = "preparation_or_execution_error"
        records.append(record)
    return dict(
        ensemble=pose["ensemble"],
        pose_id=pose["pose_id"],
        pb_valid=aggregate(r["pb_valid"] for r in records),
        copies=records,
        status=(
            "complete"
            if records and all(r["status"] == "complete" for r in records)
            else "incomplete"
        ),
    )


def load_results(path, expected):
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle))
    indexed = {}
    for row in rows:
        key = (row.get("ensemble"), row.get("pose_id"))
        if key in indexed:
            raise ValueError(f"Duplicate PoseBusters result: {key}")
        if row.get("config") != "dock":
            raise ValueError("Supplied PB CSV must declare config=dock")
        row.update(pb_valid=boolean(row.get("pb_valid")), source="supplied")
        indexed[key] = row
    missing = set(expected) - indexed.keys()
    if missing:
        raise ValueError(f"Missing PoseBusters results for {sorted(missing)}")
    return indexed


def pose_result(pose, ensemble, *, skip=False, supplied=None):
    if skip:
        record = dict(ensemble=ensemble, pose_id=pose["pose_id"], pb_valid=None, status="skipped")
    elif supplied is not None:
        record = supplied[ensemble, pose["pose_id"]]
    elif "prepared" not in pose:
        record = dict(ensemble=ensemble, pose_id=pose["pose_id"], pb_valid=None,
                      status="preparation_error", error=pose.get("preparation_error"))
    else:
        record = run_pose(pose, pose["prepared"])
    record.setdefault("config", "dock")
    return record


def filter_group(group, records, filter_invalid):
    poses, exclusions = [], []
    block = None
    for pose, record in zip(group["poses"], records):
        if filter_invalid and record["pb_valid"] is False:
            exclusions.append(dict(pose_id=pose["pose_id"], reason="posebusters_invalid", result=record))
        else:
            poses.append(pose)
            if filter_invalid and record["pb_valid"] is None:
                block = "PoseBusters filtering requested but one or more poses have unknown PB status"
    return dict(group, poses=poses, exclusions=exclusions,
                filter_pb_invalid=filter_invalid, pb_filter_error=block)
