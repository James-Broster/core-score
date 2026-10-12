"""Ligand selection and coordinate-preserving graph construction."""

from collections import Counter
import numpy as np
from rdkit import Chem
from rdkit.Geometry import Point3D
from core_score.stages.utils.preparation.cif import ligand_residues


def template_mol(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or not mol.GetNumHeavyAtoms():
        raise ValueError("ligand SMILES is missing or invalid")
    return Chem.RemoveAllHs(mol)


def heavy_atoms(rows):
    return [r for r in rows if r["type_symbol"].capitalize() not in ("H", "D")]


def atom_name(row):
    return row.get("label_atom_id", row.get("auth_atom_id"))


def composition(rows):
    """Heavy-atom element counts of a residue's atom records."""
    return Counter(r["type_symbol"].capitalize() for r in heavy_atoms(rows))


def choose_ligands(doc, pose):
    """Ligand residues picked out by the pose's chains, residue number and SMILES composition.
    With no selector at all, the structure must hold exactly one ligand residue."""
    candidates = ligand_residues(doc)
    selected = candidates
    if pose.get("ligand_chains"):
        missing = set(pose["ligand_chains"]) - {key[0] for key in candidates}
        if missing:
            raise ValueError(f"Ligand chains absent: {sorted(missing)}; candidates: {list(candidates)}")
        selected = {k: v for k, v in selected.items() if k[0] in pose["ligand_chains"]}
    if pose.get("ligand_residue"):
        selected = {k: v for k, v in selected.items() if k[1] == str(pose["ligand_residue"])}
    if pose.get("smiles"):
        expected = Counter(a.GetSymbol() for a in template_mol(pose["smiles"]).GetAtoms())
        selected = {k: v for k, v in selected.items() if composition(v) == expected}
    if not any(pose.get(k) for k in ("ligand_chains", "ligand_residue", "smiles")) and len(selected) != 1:
        raise ValueError(
            f"Ambiguous ligand selection; candidates: {list(candidates)}. Supply ligand selection or SMILES."
        )
    if not selected:
        raise ValueError(f"Selected ligand absent; candidates: {list(candidates)}")
    if len({tuple(sorted(composition(rows).items())) for rows in selected.values()}) > 1:
        raise ValueError(
            "Selected residues contain different ligand compositions; give each ligand identity its own ensemble"
        )
    for key, rows in selected.items():
        names = [atom_name(r) for r in heavy_atoms(rows)]
        if len(names) != len(set(names)):
            raise ValueError(f"Alternate conformers or duplicate ligand atom names in {key}; supply one conformer")
    return selected


def coordinate_mol(rows):
    """Construct directly from atoms, retaining the old safe-XYZ fix without monkeypatching."""
    rows = heavy_atoms(rows)
    rw = Chem.RWMol()
    conf = Chem.Conformer(len(rows))
    for i, row in enumerate(rows):
        atom = Chem.Atom(row["type_symbol"].capitalize())
        atom.SetNoImplicit(True)
        atom.SetProp(
            "_CIF_ATOM_NAME", row.get("label_atom_id", row.get("auth_atom_id", str(i)))
        )
        rw.AddAtom(atom)
        xyz = [float(row["Cartn_" + c]) for c in "xyz"]
        if not np.isfinite(xyz).all():
            raise ValueError("Non-finite ligand coordinates")
        conf.SetAtomPosition(i, Point3D(*xyz))
    mol = rw.GetMol()
    mol.AddConformer(conf)
    return mol


def identity(mol):
    mol = Chem.RemoveAllHs(Chem.Mol(mol))
    Chem.RemoveStereochemistry(mol)
    return Chem.MolToSmiles(mol)


def native_graph(doc, rows):
    """Use CIF bond tables or AF3's native SMILES/name bijection when available."""
    names = [atom_name(r) for r in heavy_atoms(rows)]
    comp = rows[0].get("label_comp_id", rows[0].get("auth_comp_id"))
    block = doc.sole_block()
    table = block.find_mmcif_category("_chem_comp_bond.")
    if len(table):
        tags = [str(t).split(".", 1)[1] for t in table.tags]
        bonds = [dict(zip(tags, (str(v).strip("\"'") for v in row))) for row in table]
        bonds = [b for b in bonds if b.get("comp_id") == comp]
        if bonds:
            mol = Chem.RWMol(coordinate_mol(rows))
            orders = {
                "SING": Chem.BondType.SINGLE,
                "DOUB": Chem.BondType.DOUBLE,
                "TRIP": Chem.BondType.TRIPLE,
                "AROM": Chem.BondType.AROMATIC,
            }
            for b in bonds:
                a, c = b["atom_id_1"], b["atom_id_2"]
                if a in names and c in names:
                    order = b.get("value_order", "").upper()
                    if order not in orders:
                        raise ValueError(f"Unsupported native bond order {order}")
                    mol.AddBond(names.index(a), names.index(c), orders[order])
            return mol.GetMol(), "native_cif_bonds"
    table = block.find_mmcif_category("_chem_comp.")
    tags = [str(t).split(".", 1)[1] for t in table.tags]
    for record in table:
        row = dict(zip(tags, (str(v).strip("\"'") for v in record)))
        if row.get("id") == comp and row.get("pdbx_smiles") not in (None, "?", "."):
            mol = template_mol(row["pdbx_smiles"])
            counts, expected = Counter(), []
            for atom in mol.GetAtoms():
                symbol = atom.GetSymbol()
                counts[symbol] += 1
                expected.append(symbol + str(counts[symbol]))
            if len(names) == len(expected) and set(names) == set(expected):
                mol = Chem.RenumberAtoms(mol, [expected.index(n) for n in names])
                mol.AddConformer(Chem.Conformer(coordinate_mol(rows).GetConformer()))
                return mol, "native_smiles_atom_names"
    return None, None
