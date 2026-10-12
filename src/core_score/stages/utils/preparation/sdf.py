"""Separate OST connectivity SDFs from chemically validated PoseBusters inputs."""

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
from core_score.stages.utils.preparation.chemistry import atom_name, coordinate_mol, heavy_atoms, identity, native_graph, template_mol


def generic_graph(mol):
    result = Chem.RWMol(Chem.Mol(mol))
    Chem.RemoveStereochemistry(result)
    for atom in result.GetAtoms():
        atom.SetIsAromatic(False)
        atom.SetFormalCharge(0)
        atom.SetNumExplicitHs(0)
        atom.SetNumRadicalElectrons(0)
        atom.SetNoImplicit(True)
    for bond in result.GetBonds():
        bond.SetIsAromatic(False)
        bond.SetBondType(Chem.BondType.SINGLE)
    mol = result.GetMol()
    mol.UpdatePropertyCache(strict=False)
    return mol


def matched_template(template, observed):
    """Transfer the exact template graph only after an element-labelled graph bijection.

    Retains dative-bond direction, charges, radical state and hydrogen counts;
    never changes coordinates or adds a missing observed bond.
    """
    left, right = generic_graph(template), generic_graph(observed)
    if (left.GetNumAtoms(), left.GetNumBonds()) != (
        right.GetNumAtoms(),
        right.GetNumBonds(),
    ):
        raise ValueError("Chemical identity mismatch: atom/bond counts differ")
    mapping = right.GetSubstructMatch(left)
    if len(mapping) != left.GetNumAtoms():
        raise ValueError("Chemical identity mismatch: non-isomorphic heavy-atom graph")
    inverse = [mapping.index(i) for i in range(len(mapping))]
    mol = Chem.RenumberAtoms(template, inverse)
    mol.RemoveAllConformers()
    mol.AddConformer(Chem.Conformer(observed.GetConformer()))
    return mol


def make_ligand(doc, rows, smiles=None, cache=None, ensemble=None):
    """OST and PoseBusters molecules for one ligand copy, with the route its graph came from.

    The OST molecule keeps connectivity from coordinates only. The PoseBusters molecule
    takes the SMILES template's chemistry, matched atom by atom onto a graph from the
    structure, another pose in the ensemble, or the coordinates. If that fails it is None
    and the reason is returned instead.
    """
    cache = {} if cache is None else cache
    coords = coordinate_mol(rows)
    names = [atom_name(r) for r in heavy_atoms(rows)]
    key = (ensemble, smiles, tuple(sorted(names)))
    template = template_mol(smiles) if smiles else None
    graph, provenance = native_graph(doc, rows)
    if provenance == "native_smiles_atom_names" and template is not None and identity(graph) != identity(template):
        raise ValueError("Chemical identity mismatch between native and supplied SMILES")
    if graph is None and key in cache:
        peer_names, peer_graph = cache[key]
        graph = Chem.RenumberAtoms(peer_graph, [peer_names.index(n) for n in names])
        graph.RemoveAllConformers()
        graph.AddConformer(Chem.Conformer(coords.GetConformer()))
        provenance = "same_ensemble_atom_name_graph"
    if graph is None:
        graph = Chem.Mol(coords)
        rdDetermineBonds.DetermineConnectivity(graph)
        provenance = "coordinate_connectivity"
    ost_mol = Chem.Mol(coords)
    rdDetermineBonds.DetermineConnectivity(ost_mol)
    try:
        if template is not None:
            pb_mol = matched_template(template, graph)
        elif provenance == "native_smiles_atom_names":
            pb_mol = Chem.Mol(graph)
        else:
            raise ValueError(
                "PoseBusters chemistry requires ligand SMILES or an original ligand SDF; bond connectivity alone is insufficient"
            )
        Chem.SanitizeMol(pb_mol)
        Chem.AssignStereochemistryFrom3D(pb_mol, replaceExistingTags=True)
        if template is not None and identity(pb_mol) != identity(template):
            raise ValueError("Preparation changed chemical identity; PB not evaluated")
        if not np.array_equal(pb_mol.GetConformer().GetPositions(), coords.GetConformer().GetPositions()):
            raise ValueError("Preparation changed predicted coordinates")
        if [a.GetSymbol() for a in pb_mol.GetAtoms()] != [a.GetSymbol() for a in coords.GetAtoms()]:
            raise ValueError("Preparation changed atom correspondence")
    except (ValueError, RuntimeError) as exc:
        return ost_mol, None, provenance, str(exc)
    cache[key] = (names, Chem.Mol(pb_mol))
    return ost_mol, pb_mol, provenance, None


def write_molecule(mol, path):
    Chem.MolToMolFile(mol, str(path))
    reread = Chem.MolFromMolFile(str(path), sanitize=False, removeHs=False)
    if reread is None or reread.GetNumAtoms() != mol.GetNumAtoms():
        raise ValueError("SDF serialization failed")
    if not np.allclose(
        mol.GetConformer().GetPositions(),
        reread.GetConformer().GetPositions(),
        atol=0.00011,
        rtol=0,
    ):
        raise ValueError(
            "SDF serialization changed coordinates beyond format precision"
        )


def read_original(path, smiles=None):
    molecules = list(Chem.SDMolSupplier(str(path), removeHs=False, sanitize=True))
    if len(molecules) != 1 or molecules[0] is None:
        raise ValueError(f"Expected one readable molecule per original SDF: {path}")
    mol = molecules[0]
    if not mol.GetNumConformers():
        raise ValueError("Original SDF has no coordinates")
    if smiles and identity(mol) != identity(template_mol(smiles)):
        raise ValueError(
            "Chemical identity mismatch between original SDF and supplied SMILES"
        )
    return mol
