"""Compare pairs of prepared poses with OST and summarise its scores."""

import itertools
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from core_score.results import diagnostic
from core_score.stages.utils import ost
from core_score.stages.utils.radius import score_with_radius_fallback


def pair_tasks(groups):
    """Every pair of poses within each site ensemble; earlier poses are the OST model."""
    return [
        (group, pose_i, pose_j)
        for group in groups
        if not group.get("error")
        for pose_i, pose_j in itertools.combinations(group["poses"], 2)
    ]


def run_pair_tasks(tasks, workers=1, **options):
    """(row, detail) for each task, in task order, using up to `workers` processes."""
    job = partial(pair_job, **options)
    if workers == 1:
        return [job(task) for task in tasks]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        return list(pool.map(job, tasks))


def pair_job(task, *, timeout=None, radius_fallback=True):
    group, pose_i, pose_j = task
    ids = dict(ensemble=group["ensemble"], site=group["site"],
               pose_i=pose_i["pose_id"], pose_j=pose_j["pose_id"])
    try:
        row, detail = compare_pair(pose_i, pose_j, timeout=timeout, radius_fallback=radius_fallback)
    except Exception as exc:
        detail = diagnostic("comparison", exc, inputs=[pose_i["path"], pose_j["path"]])
        row = dict(n_assigned=0, error=detail["error"])
    return dict(ids, **row), dict(ids, **detail)


def compare_pair(pose_i, pose_j, *, timeout=None, radius_fallback=True):
    model, reference = pose_i["prepared"], pose_j["prepared"]
    attempts = []

    def evaluate(radius):
        result = ost.compare_ligands(model, reference, timeout, radius)
        attempts.append(result)
        if result["status"] == "FAILURE":
            return dict(score=None, reason=result["exception"])
        assigned = result.get("rmsd", {}).get("assigned_scores", [])
        return dict(score=assigned[0].get("score") if len(assigned) == 1 else None, reason=result.get("rmsd", {}))

    # Recovery needs one ligand copy per side; ensembles are split by site before comparison.
    single = model["copies"] == reference["copies"] == 1
    _, radius = score_with_radius_fallback(evaluate, (6.0, 8.0) if radius_fallback and single else ())
    result = dict(attempts[0])
    if radius["radius_fallback_used"]:
        # Only the binding-site RMSD metrics use the enlarged site; lDDT-PLI stays from 4 A.
        result["rmsd"] = attempts[-1]["rmsd"]
    row = summarize(result)
    if row["n_assigned"] == 0 and not row["error"]:
        row["error"] = "OST assigned no ligand pair; see the saved OST results"
    row.update(
        standard_rmsd=radius["standard_rmsd"],
        binding_site_radius_angstrom=radius["binding_site_radius_angstrom"],
        radius_fallback_used=radius["radius_fallback_used"],
    )
    detail = dict(error=row["error"], source_pose_i=pose_i["path"], source_pose_j=pose_j["path"],
                  radius_attempts=radius["attempts"], ost_results=attempts)
    return row, detail


def summarize(result):
    """OST's scores for the pair, blank unless exactly one ligand pair was assigned."""
    failed = result["status"] == "FAILURE"
    assigned = [] if failed else [
        s for s in result.get("rmsd", {}).get("assigned_scores", []) if s.get("score") is not None
    ]
    contacts = [] if failed else [
        s["score"] for s in result.get("lddt_pli", {}).get("assigned_scores", []) if s.get("score") is not None
    ]
    single = assigned[0] if len(assigned) == 1 else {}
    return dict(
        n_assigned=len(assigned),
        rmsd=single.get("score"),
        bb_rmsd=single.get("bb_rmsd"),
        lddt_lp=single.get("lddt_lp"),
        lddt_pli=contacts[0] if len(contacts) == 1 else None,
        error=str(result.get("exception") or result.get("traceback") or "OST failed") if failed else "",
    )
