import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem
from core_score.stages.utils.preparation.sdf import (
    matched_template,
    generic_graph,
    write_molecule,
    read_original,
)
from core_score.stages.utils.preparation.chemistry import identity


def conformer(smiles):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol, randomSeed=9)
    return Chem.RemoveAllHs(mol)


def test_transfer_preserves_chemistry_coordinates_and_atom_order():
    template = Chem.MolFromSmiles("CC(=O)[O-]")
    original = conformer("CC(=O)[O-]")
    observed = Chem.RenumberAtoms(generic_graph(original), [3, 1, 2, 0])
    result = matched_template(template, observed)
    Chem.SanitizeMol(result)
    assert identity(result) == identity(template)
    np.testing.assert_array_equal(
        result.GetConformer().GetPositions(), observed.GetConformer().GetPositions()
    )
    assert [a.GetSymbol() for a in result.GetAtoms()] == [
        a.GetSymbol() for a in observed.GetAtoms()
    ]


def test_connectivity_mismatch_is_not_repaired():
    observed = conformer("CCC")
    with pytest.raises(ValueError, match="identity"):
        matched_template(Chem.MolFromSmiles("C1CC1"), observed)


def test_original_sdf_identity_and_coordinate_precision(tmp_path):
    mol = conformer("CCO")
    path = tmp_path / "ligand.sdf"
    write_molecule(mol, path)
    assert identity(read_original(path, "CCO")) == "CCO"
    with pytest.raises(ValueError, match="identity"):
        read_original(path, "CCN")


def test_failed_chemistry_validation_cannot_return_pb_input(monkeypatch):
    import gemmi
    from core_score.stages.utils.preparation import sdf

    # An otherwise valid molecule must be discarded if a later validation fails.
    observed = conformer("CCO")
    rows = [dict(type_symbol=a.GetSymbol(), label_atom_id=f"a{i}",
                 Cartn_x=str(x), Cartn_y=str(y), Cartn_z=str(z))
            for i, (a, (x, y, z)) in enumerate(zip(
                observed.GetAtoms(), observed.GetConformer().GetPositions()))]
    monkeypatch.setattr(sdf, "native_graph", lambda *args: (observed, "native_cif_bonds"))
    def reject(mol):
        raise ValueError("Synthetic sanitization failure")
    monkeypatch.setattr(sdf.Chem, "SanitizeMol", reject)
    cache = {}
    ost, pb, provenance, error = sdf.make_ligand(
        gemmi.cif.Document(), rows, "CCO", cache, ("ensemble", "CCO"))
    assert ost is not None
    assert pb is None
    assert error == "Synthetic sanitization failure"
    assert cache == {}


def test_pb_does_not_run_when_preparation_error_has_a_stale_sdf(monkeypatch):
    from unittest.mock import Mock
    import posebusters
    from core_score.stages.utils.posebusters_utils import run_pose

    buster = Mock(config={"modules": []})
    monkeypatch.setattr(posebusters, "PoseBusters", lambda **kwargs: buster)
    result = run_pose(
        dict(ensemble="e", pose_id="p", path="input.cif"),
        dict(receptor="receptor.pdb", copy_records=[dict(
            copy_id="0", pb_sdf="stale.sdf", pb_error="Identity mismatch")]))
    buster.bust.assert_not_called()
    assert result["pb_valid"] is None
    assert result["copies"][0]["status"] == "preparation_or_execution_error"


def test_peer_graph_recovery_keeps_recipient_coordinates_and_cache_scope(monkeypatch):
    import gemmi
    from core_score.stages.utils.preparation import sdf

    donor = conformer("CCO")
    positions = donor.GetConformer().GetPositions().copy()
    positions[2] += 20  # Distorted pose: distance inference loses the C-O bond.
    rows = [dict(type_symbol=a.GetSymbol(), label_atom_id=f"a{i}",
                 Cartn_x=str(x), Cartn_y=str(y), Cartn_z=str(z))
            for i, (a, (x, y, z)) in enumerate(zip(donor.GetAtoms(), positions))]
    key = ("ensemble", "CCO", ("a0", "a1", "a2"))
    cache = {}
    doc = gemmi.cif.Document()
    monkeypatch.setattr(sdf, "native_graph", lambda *args: (None, None))
    assert sdf.make_ligand(doc, rows, "CCO", cache, key)[3]

    donor_rows = [dict(row, Cartn_x=str(x), Cartn_y=str(y), Cartn_z=str(z))
                  for row, (x, y, z) in zip(rows, donor.GetConformer().GetPositions())]
    monkeypatch.setattr(sdf, "native_graph", lambda *args: (donor, "native_cif_bonds"))
    assert sdf.make_ligand(doc, donor_rows, "CCO", cache, key)[3] is None
    monkeypatch.setattr(sdf, "native_graph", lambda *args: (None, None))
    _, molecule, provenance, error = sdf.make_ligand(doc, rows, "CCO", cache, key)
    assert error is None
    assert provenance == "same_ensemble_atom_name_graph"
    np.testing.assert_array_equal(molecule.GetConformer().GetPositions(), positions)
    # A different ensemble or different declared SMILES cannot borrow this entry.
    assert sdf.make_ligand(doc, rows, "CCO", cache, ("other", *key[1:]))[3]
    assert sdf.make_ligand(doc, rows, "COC", cache, (key[0], "COC", key[2]))[3]
