import json
import pytest
from core_score.cli import parser
from core_score.inputs import load_inputs


def load(words):
    args = parser().parse_args(words)
    return load_inputs(args.inputs, smiles=args.smiles, ligand_chains=args.ligand_chains,
                       reference=args.reference, evaluate_reference=args.evaluate_reference)


def test_json_requires_no_chain_or_residue_selector(tmp_path):
    source = tmp_path / 'pose.cif'
    source.write_text('structure')
    manifest = tmp_path / 'inputs.json'
    manifest.write_text(json.dumps({'ensembles': [
        {'ensemble': 'target', 'smiles': 'CCO', 'poses': [{'path': 'pose.cif'}]}
    ]}))
    groups = load([str(manifest), '--out', 'results'])
    pose = groups[0]['poses'][0]
    assert pose['path'] == str(source)
    assert pose['smiles'] == 'CCO'
    assert not pose.get('ligand_chains')


@pytest.mark.parametrize('smiles', [None, ' ', 'not-a-smiles'])
def test_missing_or_invalid_smiles_has_pose_context(tmp_path, smiles):
    source = tmp_path / 'pose.cif'
    source.write_text('structure')
    words = [str(source), '--out', 'results']
    if smiles is not None:
        words += ['--smiles', smiles]
    with pytest.raises(ValueError, match='ensemble / pose:'):
        load(words)


def test_direct_inputs_and_reference_use_working_directory(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / 'pose.cif').write_text('structure')
    (tmp_path / 'ref.cif').write_text('reference')
    groups = load([
        'pose.cif', '--smiles', 'CCO', '--reference', 'ref.cif', '--evaluate-reference', '--out', 'results'
    ])
    assert groups[0]['poses'][0]['smiles'] == 'CCO'
    assert groups[0]['reference']['path'] == str(tmp_path / 'ref.cif')


def test_missing_structure_is_rejected(tmp_path):
    with pytest.raises(FileNotFoundError):
        load([
            str(tmp_path / 'missing.cif'), '--smiles', 'CCO', '--out', 'results'
        ])


@pytest.mark.parametrize("requested", [False, True])
def test_unused_missing_json_reference_does_not_block_core(tmp_path, requested):
    (tmp_path / "pose.cif").write_text("structure")
    manifest = tmp_path / "inputs.json"
    manifest.write_text(json.dumps({"ensembles": [{
        "ensemble": "target", "smiles": "CCO",
        "poses": [{"path": "pose.cif"}],
        "reference": {"path": "missing-reference.cif"}
    }]}))
    words = [str(manifest), "--out", "results"]
    if requested:
        words.append("--evaluate-reference")
        with pytest.raises(FileNotFoundError, match="missing-reference"):
            load(words)
    else:
        assert load(words)[0]["ensemble"] == "target"


def test_requested_reference_must_be_supplied(tmp_path):
    source = tmp_path / "pose.cif"
    source.write_text("structure")
    words = [str(source), "--smiles", "CCO", "--out", "results", "--evaluate-reference"]
    with pytest.raises(ValueError, match="requires a ground-truth structure"):
        load(words)


def test_unused_cli_reference_does_not_block_core(tmp_path):
    source = tmp_path / "pose.cif"
    source.write_text("structure")
    groups = load([
        str(source), "--smiles", "CCO", "--reference",
        str(tmp_path / "missing.cif"), "--out", "results"
    ])
    assert "reference" not in groups[0]
