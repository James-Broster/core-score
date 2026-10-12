"""Load prediction ensembles from a JSON manifest or from structure files given directly."""

import json
from pathlib import Path
from core_score.stages.utils.preparation.chemistry import template_mol

MANIFEST_ERROR = (
    "A CORE Score manifest needs an 'ensembles' list. Each ensemble needs a unique "
    "'ensemble' name and a non-empty 'poses' list of objects with a 'path'. "
    "See examples/manifest.json."
)


def load_inputs(inputs, *, smiles=None, ligand_chains=None, reference=None, evaluate_reference=False):
    """Ensemble records from one JSON manifest, or from structure files forming one ensemble.

    Settings fall back from pose to ensemble to the command line. References are loaded
    only for reference evaluation. Poses are kept in (seed, sample, pose_id) order, which
    fixes the site-assignment anchor and which pose of each pair OST treats as the model.
    """
    inputs = [Path(path) for path in inputs]
    if len(inputs) == 1 and inputs[0].suffix.lower() == ".json":
        base = inputs[0].resolve().parent
        ensembles = read_manifest(inputs[0])
    else:
        base = Path.cwd()
        ensembles = [dict(ensemble="ensemble", poses=[dict(path=str(p), pose_id=p.stem) for p in inputs])]
    defaults = {k: v for k, v in dict(smiles=smiles, ligand_chains=ligand_chains).items() if v}
    reference = reference and str(Path(reference).resolve())
    groups = []
    for ensemble in ensembles:
        group = make_group(ensemble, base, defaults)
        if evaluate_reference:
            group["reference"] = make_reference(ensemble.get("reference") or reference, group, base)
        groups.append(group)
    for group in groups:
        for pose in group["poses"]:
            try:
                template_mol(pose.get("smiles") or "")
            except ValueError as exc:
                raise ValueError(
                    f"{group['ensemble']} / {pose['pose_id']}: {exc}. Set smiles in the JSON "
                    "ensemble or pose, or pass --smiles with structure files."
                ) from exc
    return groups


def read_manifest(path):
    data = json.loads(Path(path).read_text())
    ensembles = data.get("ensembles") if isinstance(data, dict) else None
    valid = isinstance(ensembles, list) and ensembles and all(
        isinstance(e, dict) and e.get("ensemble") and isinstance(e.get("poses"), list) and e["poses"]
        and all(isinstance(p, dict) and p.get("path") for p in e["poses"])
        for e in ensembles
    )
    if not valid or len({str(e["ensemble"]) for e in ensembles}) != len(ensembles):
        raise ValueError(MANIFEST_ERROR)
    return ensembles


def make_group(ensemble, base, defaults):
    name = str(ensemble["ensemble"])
    shared = dict(defaults, **{k: ensemble[k] for k in ("smiles", "ligand_chains", "ligand_residue") if k in ensemble})
    poses = [make_pose(record, i, base, shared, name) for i, record in enumerate(ensemble["poses"])]
    if len({pose["pose_id"] for pose in poses}) != len(poses):
        raise ValueError(f"{name}: pose IDs must be unique")
    if len({"confidence" in pose for pose in poses}) != 1:
        raise ValueError(f"{name}: give a confidence for every pose or for none")
    poses.sort(key=lambda pose: (pose["seed"], pose["sample"], pose["pose_id"]))
    return dict(ensemble=name, poses=poses)


def make_pose(record, index, base, shared, ensemble):
    pose = dict(shared, **{k: v for k, v in record.items() if v not in ("", None)})
    pose.update(
        ensemble=ensemble,
        path=resolve(base, pose["path"]),
        pose_id=str(pose.get("pose_id", f"p{index:04d}")),
        seed=int(pose.get("seed", 0)),
        sample=int(pose.get("sample", index)),
    )
    if pose.get("ligand_sdfs"):
        pose["ligand_sdfs"] = [resolve(base, path) for path in pose["ligand_sdfs"]]
    if "confidence" in pose:
        pose["confidence"] = float(pose["confidence"])
    return pose


def make_reference(reference, group, base):
    if not reference:
        raise ValueError(
            f"{group['ensemble']}: reference comparison requires a ground-truth structure. "
            "Set reference in JSON or supply --reference."
        )
    if isinstance(reference, str):
        reference = dict(path=reference, smiles=group["poses"][0].get("smiles"))
    return make_pose(dict(reference, pose_id="reference"), 0, base, {}, group["ensemble"])


def resolve(base, path):
    return str((base / path).resolve(strict=True))
