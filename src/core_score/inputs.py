"""Direct structure paths and JSON manifests become the same ensemble records."""

import json
from pathlib import Path


def resolve_pose(record, base, index=0, defaults=None):
    pose = dict(defaults or {})
    pose.update({k: v for k, v in record.items() if v not in ("", None)})
    path = pose.get("path") or pose.get("structure")
    if not path:
        raise ValueError("Each pose requires a structure/path")
    pose["path"] = str((base / path).resolve(strict=True))
    pose["pose_id"] = str(pose.get("pose_id", f"p{index:04d}"))
    pose["seed"] = int(pose.get("seed", 0))
    pose["sample"] = int(pose.get("sample", index))
    pose["ensemble"] = str(pose.get("ensemble", "ensemble"))
    if pose.get("ligand_chain"):
        pose["ligand_chains"] = str(pose["ligand_chain"]).split(";")
    if pose.get("ligand_sdf"):
        pose["ligand_sdfs"] = str(pose["ligand_sdf"]).split(";")
    if pose.get("ligand_sdfs"):
        pose["ligand_sdfs"] = [
            str((base / p).resolve(strict=True)) for p in pose["ligand_sdfs"]
        ]
    return pose


def load_inputs(paths=(), manifest=None, *, load_references=True, **defaults):
    if bool(paths) == bool(manifest):
        raise ValueError("Supply direct structure paths OR --manifest, not both")
    groups = []
    if manifest:
        manifest = Path(manifest)
        if manifest.suffix.lower() != ".json":
            raise ValueError(
                "Input manifests must be JSON (.json); CSV manifests are no longer supported. "
                "Use an ensembles array containing an ensemble name and poses with file paths."
            )
        manifest = manifest.resolve(strict=True)
        base = manifest.parent
        raw = json.loads(manifest.read_text())
        if isinstance(raw, dict):
            raw = raw.get("ensembles", [raw])
        if not isinstance(raw, list) or not raw:
            raise ValueError("A CORE Score manifest must contain a non-empty ensembles array")
        for group in raw:
            if (
                not isinstance(group, dict)
                or not isinstance(group.get("poses"), list)
                or not group["poses"]
                or any(not isinstance(pose, dict) for pose in group["poses"])
            ):
                raise ValueError(
                    "Each CORE Score ensemble needs a non-empty poses array of objects with "
                    "structure file paths. Model-specific input JSON is not imported automatically."
                )
    else:
        base = Path.cwd()
        raw = [
            dict(
                ensemble="ensemble",
                poses=[dict(path=str(p), pose_id=Path(p).stem) for p in paths],
            )
        ]
    identities = set()
    for g in raw:
        ensemble = str(
            g.get("ensemble")
            or "|".join(str(g.get(k, "")) for k in ("method", "target", "ligand_key"))
        )
        if not ensemble or ensemble in identities:
            raise ValueError(f"Duplicate or empty ensemble: {ensemble}")
        identities.add(ensemble)
        group_defaults = dict(defaults, ensemble=ensemble)
        group_defaults.update(
            {k: g[k] for k in ("smiles", "ligand_chains", "ligand_residue") if k in g}
        )
        poses = [
            resolve_pose(p, base, i, group_defaults) for i, p in enumerate(g["poses"])
        ]
        if not poses or len({p["pose_id"] for p in poses}) != len(poses):
            raise ValueError(f"Empty ensemble or duplicate pose IDs: {ensemble}")
        if len({(p["seed"], p["sample"]) for p in poses}) != len(poses):
            raise ValueError(f"Duplicate seed/sample identifiers: {ensemble}")
        group = dict(
            g,
            ensemble=ensemble,
            method=g.get("method", "input"),
            target=g.get("target", ensemble),
            ligand_key=g.get("ligand_key", "selected"),
            poses=poses,
        )
        if load_references:
            ref = g.get("reference") or defaults.get("reference")
            pose_refs = {p["reference"] for p in poses if p.get("reference")}
            if not ref and pose_refs:
                if len(pose_refs) != 1:
                    raise ValueError(f"Conflicting references in {ensemble}")
                ref = pose_refs.pop()
            if ref:
                if isinstance(ref, str):
                    ref = dict(path=ref)
                    first = poses[0]
                    for key in ("smiles",):
                        if first.get(key):
                            ref[key] = first[key]
                    for name in ("ligand_chain", "ligand_residue", "ligand_sdf"):
                        val = first.get("reference_" + name) or defaults.get(
                            "reference_" + name
                        )
                        if val:
                            ref[name] = val
                group["reference"] = resolve_pose(
                    ref, base, defaults=dict(ensemble=ensemble, pose_id="reference")
                )
        groups.append(group)
    return groups


def load_workflow_inputs(args):
    """Load CLI inputs and require the intended ligand SMILES for each prediction.

    JSON carries ensemble/pose settings. CLI settings are optional fallbacks for
    direct-file usage; JSON values take precedence. Ligand extraction and chemical
    matching remain in preparation, so they are not duplicated here.
    """
    from core_score.stages.utils.preparation.chemistry import template_mol

    input_settings = {
        "smiles": args.smiles,
        "ligand_chains": args.ligand_chains,
        "ligand_residue": args.ligand_residue,
        "ligand_sdfs": args.ligand_sdfs,
    }
    evaluate_reference = args.evaluate_reference
    reference = getattr(args, "reference", None)
    if evaluate_reference and reference is not None:
        # CLI paths are relative to the working directory, not the JSON file.
        input_settings["reference"] = str(reference.resolve())
        input_settings["reference_ligand_chain"] = args.reference_ligand_chain
        input_settings["reference_ligand_residue"] = args.reference_ligand_residue
        input_settings["reference_ligand_sdf"] = args.reference_ligand_sdf

    supplied = {name: value for name, value in input_settings.items() if value is not None}
    structures, manifest = args.structures, args.manifest
    if len(structures) == 1 and structures[0].suffix.lower() in (".json", ".csv"):
        if manifest is not None:
            raise ValueError("Supply one JSON manifest, not both INPUT and --manifest")
        manifest, structures = structures[0], []
    groups = load_inputs(
        structures, manifest, load_references=evaluate_reference, **supplied
    )
    if evaluate_reference:
        for group in groups:
            reference = group.get("reference")
            if not reference or not reference.get("path"):
                raise ValueError(
                    f"{group['ensemble']}: reference comparison requires a ground-truth "
                    "structure. Set reference in JSON or supply --reference."
                )
            if not Path(reference["path"]).is_file():
                raise ValueError(
                    f"{group['ensemble']}: ground-truth structure is not a file: "
                    f"{reference['path']}"
                )
    for group in groups:
        for pose in group["poses"]:
            label = f"{group['ensemble']} / {pose['pose_id']}"
            smiles = pose.get("smiles")
            if not isinstance(smiles, str) or not smiles.strip():
                raise ValueError(
                    f"{label}: ligand SMILES is required. Set smiles in the JSON "
                    "ensemble/pose or supply --smiles with direct structure paths."
                )
            try:
                molecule = template_mol(smiles)
                if molecule.GetNumAtoms() == 0:
                    raise ValueError("Ligand SMILES has no heavy atoms")
            except ValueError as exc:
                raise ValueError(f"{label}: {exc}") from exc
    return groups
