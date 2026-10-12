"""Prepare one structure: select its ligand copies, write SDFs for OST and PoseBusters, and
write the receptor PoseBusters checks against. Predicted coordinates are never changed."""

from pathlib import Path
import numpy as np
from rdkit import Chem
from core_score.results import fingerprint
from core_score.stages.utils.preparation.chemistry import choose_ligands, heavy_atoms
from core_score.stages.utils.preparation.cif import filter_records, ligand_residues, read_document, write_receptor
from core_score.stages.utils.preparation.sdf import make_ligand, read_original, write_molecule


def pose_directory(output, pose):
    return Path(output) / "prepared" / fingerprint((pose["ensemble"], pose["pose_id"]))[:20]


def prepare_pose(pose, destination, cache=None):
    """Prepared record for one pose; its files are written beside `destination`."""
    source = Path(pose["path"])
    destination.parent.mkdir(parents=True, exist_ok=True)
    doc = read_document(source)
    if pose.get("ligand_sdfs"):
        prepared, ligands = from_original_sdfs(pose, doc, source, destination)
    else:
        prepared, ligands = from_structure(pose, doc, destination, cache)
    receptor = destination.with_suffix(".receptor.pdb")
    write_receptor(source, ligands, receptor)
    return dict(prepared, receptor=str(receptor))


def from_original_sdfs(pose, doc, source, destination):
    """Supplied SDFs, already in the receptor's frame, are the ligand copies."""
    path = str(source)
    if source.suffix.lower() in (".pdb", ".ent"):
        doc.write_file(str(destination))
        path = str(destination)
    residues = ligand_residues(doc)
    copies, matched = [], set()
    for i, sdf in enumerate(pose["ligand_sdfs"]):
        mol = Chem.RemoveAllHs(read_original(sdf, pose.get("smiles")))
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
        cache_key = (pose.get("ensemble"), pose.get("smiles"), tuple(sorted(r.get("label_atom_id", "") for r in heavy)))
        ost_mol, pb_mol, provenance, pb_error = make_ligand(doc, rows, pose.get("smiles"), cache, cache_key)
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
