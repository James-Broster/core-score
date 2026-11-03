"""Optional comparison of every prediction against its ensemble's ground-truth structure.

A pose is correct at threshold T when its ligand RMSD to the reference is below T and its
lDDT-PLI is above 0.8, the definition used in the paper."""

import math
from core_score.outputs.files import fingerprint, write_json
from core_score.analysis.pairs import run_pair_tasks
from core_score.inputs.prepare import prepare

LDDT_PLI_CUTOFF = 0.8


def evaluate_groups(groups, output, thresholds=(2.0,), *, workers=1, timeout=None, radius_fallback=True):
    """Ligand RMSD and lDDT-PLI of each prepared prediction to its reference, labelled correct
    or not at each threshold. Ensembles without a reference are skipped. Reference results
    never enter CORE."""
    groups = [group for group in groups if "reference" in group]
    tasks, records = [], []
    for group in groups:
        destination = output / "references" / fingerprint(group["ensemble"])[:20] / "reference.cif"
        reference = prepare(group["reference"], destination, timeout=timeout)
        for pose in group["poses"]:
            error = reference.get("error") or pose.get("error")
            if error:
                records.append(dict(ensemble=group["ensemble"], pose_id=pose["pose_id"], rmsd=None, error=error))
            else:
                tasks.append((dict(group, site="reference"), pose, reference))
    for row, detail in run_pair_tasks(tasks, workers, timeout=timeout, radius_fallback=radius_fallback, lddt_pli=True):
        records.append(dict(row, pose_id=row["pose_i"], detail=detail))
    write_json(output / "reference_details.json", records)
    keys = ((group["ensemble"], pose["pose_id"]) for group in groups for pose in group["poses"])
    order = {key: i for i, key in enumerate(keys)}
    labels = label_reference_results(records, thresholds)
    return sorted(labels, key=lambda label: order[label["ensemble"], label["pose_id"]])


def correct_column(threshold):
    """The column saying whether a pose is correct at `threshold`."""
    return f"correct_{threshold:g}A"


def label_reference_results(records, thresholds):
    """One row per pose. A pose is correct at T when RMSD < T and lDDT-PLI > 0.8. A failed or
    ambiguous comparison, or one missing either score, is unknown (blank scores and labels),
    not incorrect."""
    labels = []
    for record in records:
        rmsd, lddt_pli = record.get("rmsd"), record.get("lddt_pli")
        known = not record.get("error") and all(v is not None and math.isfinite(v) for v in (rmsd, lddt_pli))
        label = dict(ensemble=record["ensemble"], pose_id=record["pose_id"],
                     rmsd=rmsd if known else None, lddt_pli=lddt_pli if known else None)
        for threshold in thresholds:
            label[correct_column(threshold)] = bool(rmsd < threshold and lddt_pli > LDDT_PLI_CUTOFF) if known else None
        label["radius_fallback_used"] = record.get("radius_fallback_used", False)
        error = record.get("error")
        label["error"] = (error["error"] if isinstance(error, dict) else error) or None  # message only; details in JSON
        labels.append(label)
    return labels
