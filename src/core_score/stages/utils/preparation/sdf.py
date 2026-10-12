"""Separate OST connectivity SDFs from chemically validated PoseBusters inputs."""

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
from core_score.stages.utils.preparation.chemistry import coordinate_mol, native_graph, template_mol, identity, heavy_atoms


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


def make_ligand(doc, rows, smiles=None, cache=None, cache_key=None):
    coords = coordinate_mol(rows)
    names = [r.get("label_atom_id", r.get("auth_atom_id")) for r in heavy_atoms(rows)]
    template = template_mol(smiles) if smiles else None
    graph, provenance = native_graph(doc, rows)
    if (
        graph is not None
        and provenance == "native_smiles_atom_names"
        and template is not None
    ):
        if identity(graph) != identity(template):
            raise ValueError(
                "Chemical identity mismatch between native and supplied SMILES"
            )
    if graph is None and cache is not None and cache_key in cache:
        old_names, old_graph = cache[cache_key]
        if set(names) == set(old_names) and len(names) == len(old_names):
            graph = Chem.RenumberAtoms(old_graph, [old_names.index(n) for n in names])
            graph.RemoveAllConformers()
            graph.AddConformer(Chem.Conformer(coords.GetConformer()))
            provenance = "same_ensemble_atom_name_graph"
    if graph is None:
        graph = Chem.Mol(coords)
        rdDetermineBonds.DetermineConnectivity(graph)
        provenance = "coordinate_connectivity"
    # The OST fallback deliberately retains historical connectivity-only chemistry.
    ost_mol = Chem.Mol(coords)
    rdDetermineBonds.DetermineConnectivity(ost_mol)
    pb_mol = None
    pb_error = None
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
        if not np.array_equal(
            pb_mol.GetConformer().GetPositions(), coords.GetConformer().GetPositions()
        ):
            raise ValueError("Preparation changed predicted coordinates")
        if [a.GetSymbol() for a in pb_mol.GetAtoms()] != [
            a.GetSymbol() for a in coords.GetAtoms()
        ]:
            raise ValueError("Preparation changed atom correspondence")
        if cache is not None and cache_key is not None:
            cache[cache_key] = (names, Chem.Mol(pb_mol))
    except (ValueError, RuntimeError) as exc:
        pb_mol = None
        pb_error = str(exc)
    return ost_mol, pb_mol, provenance, pb_error


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
