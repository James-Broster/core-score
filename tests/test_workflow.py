import json
from unittest.mock import Mock, patch
import numpy as np
import pytest
from core_score.stages.scoring import core_score
from core_score.stages.utils.score_utils import build_matrix
from core_score.cli import parser
from core_score.stages.utils.posebusters_utils import aggregate
from core_score.stages.posebusters import run_pb
from core_score.inputs import load_inputs


def test_core_known_two_pose_score_and_diagonal_independence():
    for diag in (0.0, float("nan"), 99.0):
        result = core_score([[diag, 1.0], [1.0, diag]])
        assert result["S"] == 0.5
        assert result["group_size"] == 2


def test_core_small_tight_group_wins():
    result = core_score([[0, 0.1, 10], [0.1, 0, 10], [10, 10, 0]])
    assert result["group_members"] == [0, 1]
    assert result["D"] == pytest.approx(0.1 / np.sqrt(2 / 3))


@pytest.mark.parametrize(
    "matrix",
    [
        [[0, 1]],
        [[0, float("nan")], [float("nan"), 0]],
        [[0, -1], [-1, 0]],
        [[0, 1], [2, 0]],
        [[0]],
    ],
)
def test_core_rejects_invalid_matrices(matrix):
    with pytest.raises(ValueError):
        core_score(matrix)


def test_matrix_missing_duplicate_and_failed_pairs():
    row = dict(pose_i="a", pose_j="b", rmsd=1, n_assigned=1, error="")
    np.testing.assert_array_equal(build_matrix(["a", "b"], [row]), [[0, 1], [1, 0]])
    for rows in ([], [dict(row, error="OST failed")], [dict(row, n_assigned=2)]):
        with pytest.raises(ValueError):
            build_matrix(["a", "b"], rows)


def test_cli_agreed_defaults():
    args = parser().parse_args(["one.cif", "two.cif", "--out", "results"])
    assert args.workers == 1 and args.timeout is None
    assert not args.skip_posebusters and not hasattr(args, "filter_pb_invalid")
    assert not args.no_radius_fallback and not hasattr(args, "fallback_radii")
    assert args.rmsd_thresholds == [2]
    assert not hasattr(args, "site_mode")


def test_json_groups_shared_settings_overrides_and_relative_paths(tmp_path):
    for name in ("a.cif", "b.cif", "ref.cif", "ligand.sdf"):
        (tmp_path / name).write_text("x")
    manifest = tmp_path / "poses.json"
    manifest.write_text(json.dumps({"ensembles": [
        {"ensemble": "boltz", "smiles": "CCO", "ligand_chains": ["B"],
         "poses": [{"pose_id": "a", "path": "a.cif"},
                   {"pose_id": "b", "path": "b.cif", "ligand_chains": ["C"],
                    "ligand_sdfs": ["ligand.sdf"]}],
         "reference": {"path": "ref.cif", "ligand_chains": ["L"]}},
        {"ensemble": "chai", "smiles": "CCN", "reference": "ref.cif",
         "poses": [{"pose_id": "a", "path": "a.cif"}]}
    ]}))
    groups = load_inputs([manifest], evaluate_reference=True)
    assert [g["ensemble"] for g in groups] == ["boltz", "chai"]
    assert len(groups[0]["poses"]) == 2 and len(groups[1]["poses"]) == 1
    first, second = groups[0]["poses"]
    assert first["path"] == str(tmp_path / "a.cif")
    assert first["ligand_chains"] == ["B"] and first["smiles"] == "CCO"
    assert second["ligand_chains"] == ["C"] and second["smiles"] == "CCO"
    assert second["ligand_sdfs"] == [str(tmp_path / "ligand.sdf")]
    assert groups[0]["reference"]["path"] == str(tmp_path / "ref.cif")
    assert groups[0]["reference"]["ligand_chains"] == ["L"]
    assert groups[1]["reference"]["smiles"] == "CCN"


@pytest.mark.parametrize("value", [None, [], {"ensembles": {}},
    {"name": "model-input", "sequences": []},
    {"ensembles": [{"ensemble": "x", "poses": ["a.cif"]}]}])
def test_invalid_json_manifest_has_actionable_error(tmp_path, value):
    manifest = tmp_path / "input.json"
    manifest.write_text(json.dumps(value))
    with pytest.raises(ValueError, match="ensemble|poses"):
        load_inputs([manifest])


def test_pb_unknown_never_passes():
    assert aggregate([True, None]) is None
    assert aggregate([True, False, None]) is False
    assert aggregate([]) is None


def pb_patches(answers):
    buster = patch("posebusters.PoseBusters", return_value=Mock(config={"modules": []}))
    return buster, patch("core_score.stages.posebusters.check_pose", side_effect=answers)


def test_pb_runs_and_invalid_poses_are_filtered(tmp_path):
    groups = [dict(ensemble="g", poses=[dict(ensemble="g", pose_id="a"), dict(ensemble="g", pose_id="b")])]
    answers = [dict(ensemble="g", pose_id="a", pb_valid=True), dict(ensemble="g", pose_id="b", pb_valid=False)]
    buster, check = pb_patches(answers)
    with buster, check as pb:
        retained = run_pb(groups, tmp_path)
    assert pb.call_count == 2
    assert [p["pose_id"] for p in retained[0]["poses"]] == ["a"]


def test_unknown_pb_blocks_filtered_ensemble(tmp_path):
    groups = [dict(ensemble="g", poses=[dict(ensemble="g", pose_id="a")])]
    buster, check = pb_patches([dict(ensemble="g", pose_id="a", pb_valid=None)])
    with buster, check:
        result = run_pb(groups, tmp_path)
    assert result[0]["error"]
    assert len(result[0]["poses"]) == 1


def test_skip_pb_retains_every_pose_without_running_checks(tmp_path):
    groups = [dict(ensemble="g", poses=[dict(ensemble="g", pose_id="a"), dict(ensemble="g", pose_id="b")])]
    with patch("core_score.stages.posebusters.check_pose") as pb:
        retained = run_pb(groups, tmp_path, skip=True)
    pb.assert_not_called()
    assert retained == groups


def test_pipeline_passes_the_posebusters_setting(tmp_path):
    from core_score.pipeline import execute

    stages = ["prepare_groups", "prepare_site_groups", "compare_groups", "score_groups"]
    patches = [patch(f"core_score.pipeline.{name}", return_value=[]) for name in stages]
    for p in patches:
        p.start()
    try:
        with patch("core_score.pipeline.run_pb", return_value=[]) as pb:
            execute([], tmp_path / "default")
            assert pb.call_args.kwargs["skip"] is False
            execute([], tmp_path / "skip", skip_posebusters=True)
            assert pb.call_args.kwargs["skip"] is True
    finally:
        for p in patches:
            p.stop()


def test_reference_threshold_boundary_is_strict(tmp_path):
    from core_score.stages.evaluation import evaluate_groups

    pose = dict(ensemble="g", pose_id="p", path="pose.cif", prepared=dict(copies=1))
    reference = dict(ensemble="g", pose_id="reference", path="ref.cif")
    group = dict(ensemble="g", poses=[pose], reference=reference)
    row = dict(ensemble="g", site="reference", pose_i="p", pose_j="reference", rmsd=1.0, error="")
    prepared_reference = dict(reference, prepared=dict(copies=1))
    with patch("core_score.stages.evaluation.prepare", return_value=prepared_reference), \
            patch("core_score.stages.utils.comparison_utils.pair_job", return_value=(row, {})) as compare:
        labels = evaluate_groups([group], tmp_path, [0.5, 1.0, 1.5, 2.0])
    assert compare.call_count == 1
    assert [r["success"] for r in labels] == [False, False, True, True]
