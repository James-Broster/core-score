"""Read real CIF records; convert PDB without changing source files."""

from collections import defaultdict
from pathlib import Path
import gemmi


def read_document(path):
    path = Path(path)
    if path.suffix.lower() in (".pdb", ".ent"):
        structure = gemmi.read_structure(str(path))
        if len(structure) != 1:
            raise ValueError("Supply one structural model per pose file")
        structure.setup_entities()
        return structure.make_mmcif_document()
    doc = gemmi.cif.read_file(str(path))
    doc.sole_block()
    return doc


def atom_rows(doc):
    table = doc.sole_block().find_mmcif_category("_atom_site.")
    if not len(table):
        raise ValueError("No atom records in structure")
    tags = [str(t).split(".", 1)[1] for t in table.tags]
    rows = [dict(zip(tags, (gemmi.cif.as_string(v) for v in row))) for row in table]
    models = {r.get("pdbx_PDB_model_num", "1") for r in rows}
    if len(models) != 1:
        raise ValueError("Supply one structural model per pose file")
    return rows


def residue_key(row):
    return (
        row.get("auth_asym_id") or row.get("label_asym_id", ""),
        row.get("auth_seq_id") or row.get("label_seq_id", ""),
        row.get("pdbx_PDB_ins_code", ""),
        row.get("label_comp_id") or row.get("auth_comp_id", ""),
    )


def ligand_residues(doc):
    block = doc.sole_block()
    polymer = {
        gemmi.cif.as_string(v) for v in block.find_values("_entity_poly.entity_id")
    }
    groups = defaultdict(list)
    for row in atom_rows(doc):
        if row["group_PDB"] == "HETATM" and row.get("label_entity_id") not in polymer:
            groups[residue_key(row)].append(row)
    return dict(groups)


def filter_records(doc, selected, destination):
    """Keep all polymer records and selected ligand residues, byte-value exact in CIF."""
    rows = atom_rows(doc)
    table = doc.sole_block().find_mmcif_category("_atom_site.")
    tags = list(table.tags)
    nonpoly = set(ligand_residues(doc))
    keep = [
        r["group_PDB"] != "HETATM"
        or residue_key(r) not in nonpoly
        or residue_key(r) in selected
        for r in rows
    ]
    expected = [list(row) for row, retained in zip(table, keep) if retained]
    for i in range(len(rows) - 1, -1, -1):
        if not keep[i]:
            table.remove_row(i)
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    doc.write_file(str(destination))
    actual = (
        gemmi.cif.read_file(str(destination))
        .sole_block()
        .find_mmcif_category("_atom_site.")
    )
    if list(actual.tags) != tags or [list(row) for row in actual] != expected:
        raise ValueError("CIF serialization changed a retained atom record")


def write_receptor(source, ligands, destination):
    """The predicted environment without the selected ligand residues, as PDB for PoseBusters.
    Chain IDs are kept whole."""
    doc = read_document(source)
    table = doc.sole_block().find_mmcif_category("_atom_site.")
    for i, row in reversed(list(enumerate(atom_rows(doc)))):
        if row["group_PDB"] == "HETATM" and residue_key(row) in ligands:
            table.remove_row(i)
    gemmi.make_structure_from_block(doc.sole_block()).write_pdb(str(destination))
