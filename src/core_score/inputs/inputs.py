"""Load prediction ensembles from a CSV file, or from structure files given directly.

The CSV needs three columns. `poses` is one structure file, or a wildcard matching many:

    ensemble,smiles,poses
    7ly6,CC(=O)Nc1ccc(O)cc1,predictions/7ly6/*.cif

Optional columns are confidence (one pose per row), ligand_chain (space separated when the
ligand has several copies), reference (the ground-truth structure) and ligand_sdf (the
ligand in the structure's frame, for docking outputs). Paths are relative to the CSV.
"""

import csv
import glob
import os
from collections import defaultdict
from pathlib import Path
from core_score.analysis.core import NEIGHBOURS
from core_score.inputs.structures import template_mol

REQUIRED = ("ensemble", "smiles", "poses")
OPTIONAL = ("confidence", "ligand_chain", "reference", "ligand_sdf")
MIN_POSES = NEIGHBOURS + 1


def load_inputs(inputs, *, smiles=None, ligand_chains=None, reference=None):
    """Ensemble records from one CSV, or from structure files forming one ensemble.

    Poses are sorted by path, never kept in listing order. That order fixes the
    site-assignment anchor and which pose of each pair OST treats as the model, so
    reordering the input never changes a score. Pose IDs are paths relative to the CSV.
    """
    inputs = [Path(path) for path in inputs]
    if len(inputs) == 1 and inputs[0].suffix.lower() == ".csv":
        base, rows = inputs[0].resolve().parent, read_csv(inputs[0])
    else:
        base = Path.cwd()
        rows = [dict(ensemble="ensemble", smiles=smiles or "", poses=str(path),
                     ligand_chain=" ".join(ligand_chains or []), reference=str(reference or ""))
                for path in inputs]
    ensembles = defaultdict(list)
    for row in rows:
        ensembles[row["ensemble"]].append({k: v.strip() for k, v in row.items() if v and v.strip()})
    groups = [make_group(name, members, base) for name, members in ensembles.items()]
    small = [f"{group['ensemble']} ({len(group['poses'])})" for group in groups if len(group["poses"]) < MIN_POSES]
    if small:
        raise ValueError(f"CORE needs at least {MIN_POSES} poses per ensemble, since each pose is scored "
                         f"against its {NEIGHBOURS} closest others. Too few: {', '.join(small)}")
    return groups


def read_csv(path):
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle)
        columns = reader.fieldnames or []
        missing = [c for c in REQUIRED if c not in columns]
        unknown = [c for c in columns if c not in REQUIRED + OPTIONAL]
        if missing or unknown:
            raise ValueError(f"{path}: missing columns {missing}, unknown columns {unknown}. "
                             f"The CSV needs {', '.join(REQUIRED)} and may add {', '.join(OPTIONAL)}.")
        rows = list(reader)
    blank = [i + 2 for i, row in enumerate(rows) if not all((row[c] or "").strip() for c in REQUIRED)]
    if blank:
        raise ValueError(f"{path}: lines {blank} leave ensemble, smiles or poses empty")
    return rows


def make_group(name, rows, base):
    """One ensemble from its CSV rows. Every row must give the same SMILES and reference."""
    for column in ("smiles", "reference"):
        if len({row.get(column) for row in rows}) > 1:
            raise ValueError(f"{name}: every row of an ensemble needs the same {column}")
    try:
        template_mol(rows[0].get("smiles", ""))
    except ValueError as exc:
        raise ValueError(f"{name}: {exc}. Give the ligand SMILES in the CSV, or --smiles with structure files.") from exc
    poses = [pose for row in rows for pose in make_poses(name, row, base)]
    paths = [pose["path"] for pose in poses]
    if len(set(paths)) != len(paths):
        raise ValueError(f"{name}: a structure is listed more than once")
    if len({"confidence" in pose for pose in poses}) != 1:
        raise ValueError(f"{name}: give a confidence for every pose or for none")
    group = dict(ensemble=name, poses=sorted(poses, key=lambda pose: pose["path"]))
    if rows[0].get("reference"):
        path = resolve(base, rows[0]["reference"])
        group["reference"] = dict(ensemble=name, pose_id="reference", path=path, smiles=rows[0]["smiles"])
    return group


def make_poses(name, row, base):
    """The poses one CSV row names: one file, or every file its wildcard matches."""
    paths = sorted(glob.glob(str(base / row["poses"]), recursive=True))
    if not paths:
        raise FileNotFoundError(f"{name}: no structure file matches {row['poses']}")
    single = {k: row[k] for k in ("confidence", "ligand_sdf") if k in row}
    if single and len(paths) > 1:
        raise ValueError(f"{name}: {row['poses']} matches {len(paths)} files, so give "
                         f"{' and '.join(single)} on a row of its own for each pose")
    poses = []
    for path in paths:
        path = str(Path(path).resolve())
        pose = dict(ensemble=name, pose_id=Path(os.path.relpath(path, base)).as_posix(), path=path,
                    smiles=row["smiles"], ligand_chains=row.get("ligand_chain", "").split())
        if "confidence" in row:
            pose["confidence"] = float(row["confidence"])
        if "ligand_sdf" in row:
            pose["ligand_sdfs"] = [resolve(base, row["ligand_sdf"])]
        poses.append(pose)
    return poses


def resolve(base, path):
    return str((base / path).resolve(strict=True))
