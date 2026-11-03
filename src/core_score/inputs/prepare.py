"""Prepare every pose: select its ligand copies, give identical chains one entity, write
SDFs for OST and PoseBusters and a receptor for PoseBusters, then check OST can read the
ligand. Predicted coordinates are never changed."""

from pathlib import Path
import numpy as np
from rdkit import Chem
from core_score.analysis import ost
from core_score.outputs.files import diagnostic, pose_directory, write_json
from core_score.inputs.structures import (
    choose_ligands, filter_records, heavy_atoms, ligand_residues, make_ligand, merge_identical_entities, read_document,
    read_original, write_molecule, write_receptor,
)


def prepare_groups(groups, output, *, timeout=None):
    cache = {}
    prepared = []
    for group in groups:
        poses = [prepare(pose, pose_directory(output, pose) / "pose.cif", cache, timeout) for pose in group["poses"]]
        # A ligand graph that failed can borrow one from a peer prepared later in the same ensemble.
        for i, pose in enumerate(poses):
            if "prepared" in pose and any(copy["pb_error"] for copy in pose["prepared"]["copy_records"]):
                retry = prepare(group["poses"][i], pose_directory(output, pose) / "pose.cif", cache, timeout)
                if "prepared" in retry:
                    poses[i] = retry
        prepared.append(dict(group, poses=poses))
    write_json(output / "prepared_manifest.json", dict(groups=prepared))
    return prepared


def prepare(pose, destination, cache=None, timeout=None):
    """The pose with its prepared record, or with the error that stopped preparation."""
    try:
        prepared = prepare_pose(pose, destination, cache)
        if not prepared["ligand_sdfs"] and not ost.reads_ligand(prepared, timeout):
            # OST cannot read this ligand from the structure, or connects it differently,
            # so it gets the package's SDFs, with bonds from the coordinates, instead.
            prepared.update(ligand_sdfs=[copy["ost_sdf"] for copy in prepared["copy_records"]],
                            route="explicit_connectivity_sdf")
        return dict(pose, prepared=prepared)
    except Exception as exc:
        return dict(pose, error=diagnostic("preparation", exc, structure=pose["path"]))


def prepare_pose(pose, destination, cache=None):
    """Prepared record for one pose; its files are written beside `destination`."""
    source = Path(pose["path"])
    destination.parent.mkdir(parents=True, exist_ok=True)
    doc = read_document(source)
    merge_identical_entities(doc)
    if pose.get("ligand_sdfs"):
        prepared, ligands = from_original_sdfs(pose, doc, destination)
    else:
        prepared, ligands = from_structure(pose, doc, destination, cache)
    receptor = destination.with_suffix(".receptor.pdb")
    write_receptor(source, ligands, receptor)
    return dict(prepared, receptor=str(receptor))


def from_original_sdfs(pose, doc, destination):
    """Supplied SDFs, already in the receptor's frame, are the ligand copies."""
    doc.write_file(str(destination))
    path = str(destination)
    residues = ligand_residues(doc)
    copies, matched = [], set()
    for i, sdf in enumerate(pose["ligand_sdfs"]):
        mol = Chem.RemoveAllHs(read_original(sdf, pose["smiles"]))
        # Any structure residue holding the same ligand is left out of the receptor.
        matched |= {key for key, rows in residues.items() if same_atoms(rows, mol)}
        copies.append(dict(copy_id=str(i), ost_sdf=sdf, pb_sdf=sdf, pb_error=None, provenance="original_sdf",
                           centre=mol.GetConformer().GetPositions().mean(0).tolist()))
    prepared = dict(path=path, ligand_sdfs=list(pose["ligand_sdfs"]), copies=len(copies),
                    copy_records=copies, route="original_sdf")
    return prepared, matched


def from_structure(pose, doc, destination, cache):
    """Ligand residues selected in the structure, with chemistry built from their atoms."""
    selected = choose_ligands(doc, pose)
    copies = []
    for i, (key, rows) in enumerate(selected.items()):
        heavy = heavy_atoms(rows)
        ost_mol, pb_mol, provenance, pb_error = make_ligand(doc, rows, pose["smiles"], cache, pose["ensemble"])
        ost_sdf = destination.with_name(f"{destination.stem}.copy{i}.ost.sdf")
        write_molecule(ost_mol, ost_sdf)
        pb_sdf = None
        if pb_mol is not None:
            pb_sdf = destination.with_name(f"{destination.stem}.copy{i}.pb.sdf")
            write_molecule(pb_mol, pb_sdf)
            pb_sdf = str(pb_sdf)
        centre = np.mean([[float(r["Cartn_" + c]) for c in "xyz"] for r in heavy], axis=0)
        copies.append(dict(copy_id=str(i), residue=list(key), ost_sdf=str(ost_sdf), pb_sdf=pb_sdf,
                           pb_error=pb_error, provenance=provenance, centre=centre.tolist()))
    filter_records(doc, set(selected), destination)
    prepared = dict(path=str(destination), ligand_sdfs=[], copies=len(copies), copy_records=copies, route="native_cif")
    return prepared, set(selected)


def same_atoms(rows, mol):
    """Whether a structure residue holds the molecule's heavy atoms at the same coordinates."""
    observed = sorted((r["type_symbol"].capitalize(), *(float(r["Cartn_" + c]) for c in "xyz"))
                      for r in heavy_atoms(rows))
    expected = sorted((a.GetSymbol(), *xyz) for a, xyz in zip(mol.GetAtoms(), mol.GetConformer().GetPositions()))
    # SDF coordinates have four decimal places.
    return len(observed) == len(expected) and all(
        a[0] == b[0] and np.allclose(a[1:], b[1:], atol=0.00011, rtol=0) for a, b in zip(observed, expected)
    )
