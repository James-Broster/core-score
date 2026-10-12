"""Prepare each input once; preserve original coordinates and record all routes."""

from pathlib import Path
import numpy as np
from rdkit import Chem
from core_score.stages.utils.preparation.cif import (
    read_document,
    atom_rows,
    filter_records,
    residue_key,
    write_receptor,
    ligand_residues,
)
from core_score.stages.utils.preparation.chemistry import choose_ligands, heavy_atoms
from core_score.stages.utils.preparation.sdf import make_ligand, write_molecule, read_original


def prepare_pose(pose, destination, cache=None):
    source = Path(pose["path"]).resolve(strict=True)
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    doc = read_document(source)
    atom_rows(doc)
    if pose.get("ligand_sdfs"):
        paths = [str(Path(p).resolve(strict=True)) for p in pose["ligand_sdfs"]]
        native = str(source)
        if source.suffix.lower() in (".pdb", ".ent"):
            doc.write_file(str(destination))
            native = str(destination.resolve())
        copies = []
        excluded = set()
        for i, path in enumerate(paths):
            mol = read_original(path, pose.get("smiles"))
            heavy = Chem.RemoveAllHs(mol)
            sdf_atoms = sorted(
                (a.GetSymbol(), *heavy.GetConformer().GetPositions()[j])
                for j, a in enumerate(heavy.GetAtoms())
            )
            for key, records in ligand_residues(doc).items():
                atoms = sorted(
                    (
                        r["type_symbol"].capitalize(),
                        *[float(r["Cartn_" + c]) for c in "xyz"],
                    )
                    for r in heavy_atoms(records)
                )
                if len(atoms) == len(sdf_atoms) and all(
                    a[0] == b[0] and np.allclose(a[1:], b[1:], atol=0.00011, rtol=0)
                    for a, b in zip(atoms, sdf_atoms)
                ):
                    excluded.add(key)
            copies.append(
                dict(
                    copy_id=str(i),
                    pb_sdf=path,
                    ost_sdf=path,
                    centre=Chem.RemoveAllHs(mol)
                    .GetConformer()
                    .GetPositions()
                    .mean(0)
                    .tolist(),
                    provenance="original_sdf",
                    pb_error=None,
                )
            )
        receptor_doc = read_document(source)
        table = receptor_doc.sole_block().find_mmcif_category("_atom_site.")
        for i, row in reversed(list(enumerate(atom_rows(receptor_doc)))):
            if row["group_PDB"] == "HETATM" and residue_key(row) in excluded:
                table.remove_row(i)
        receptor = destination.with_suffix(".receptor.pdb")
        write_receptor(receptor_doc, receptor)
        return dict(
            path=native,
            ligand_sdfs=paths,
            copies=len(paths),
            copy_records=copies,
            source=str(source),
            receptor=str(receptor.resolve()),
            route="original_sdf",
        )
    selected = choose_ligands(doc, pose)
    copy_records = []
    # Keep the full predicted environment for PB, excluding the selected ligand copies.
    receptor_doc = read_document(source)
    table = receptor_doc.sole_block().find_mmcif_category("_atom_site.")
    for i, row in reversed(list(enumerate(atom_rows(receptor_doc)))):
        if row["group_PDB"] == "HETATM" and residue_key(row) in selected:
            table.remove_row(i)
    receptor_path = destination.with_suffix(".receptor.pdb")
    write_receptor(receptor_doc, receptor_path)
    for i, (key, rows) in enumerate(selected.items()):
        cache_key = (
            pose.get("ensemble"),
            pose.get("smiles"),
            tuple(sorted(r.get("label_atom_id", "") for r in heavy_atoms(rows))),
        )
        ost_mol, pb_mol, provenance, pb_error = make_ligand(
            doc, rows, pose.get("smiles"), cache, cache_key
        )
        ost_path = destination.with_name(destination.stem + f".copy{i}.ost.sdf")
        write_molecule(ost_mol, ost_path)
        pb_path = destination.with_name(destination.stem + f".copy{i}.pb.sdf")
        if pb_mol is not None:
            write_molecule(pb_mol, pb_path)
        xyz = [[float(r["Cartn_" + c]) for c in "xyz"] for r in heavy_atoms(rows)]
        copy_records.append(
            dict(
                copy_id=str(i),
                residue=list(key),
                ost_sdf=str(ost_path.resolve()),
                pb_sdf=str(pb_path.resolve()) if pb_mol is not None else None,
                pb_error=pb_error,
                provenance=provenance,
                centre=np.mean(xyz, axis=0).tolist(),
            )
        )
    filter_records(doc, set(selected), destination)
    return dict(
        path=str(destination.resolve()),
        ligand_sdfs=[],
        copies=len(selected),
        source=str(source),
        copy_records=copy_records,
        receptor=str(receptor_path.resolve()),
        route="native_cif",
        ligand_residues=[list(k) for k in selected],
    )


