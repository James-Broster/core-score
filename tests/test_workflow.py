import csv
import json
from pathlib import Path
from unittest.mock import patch
import numpy as np
import pytest
from core_score.stages.scoring import core_score
from core_score.stages.utils.score_utils import build_matrix
from core_score.cli import parser, main
from core_score.stages.utils.posebusters_utils import aggregate, load_results
from core_score.stages.posebusters import run_pb
from core_score.stages.comparison import compare_groups
from core_score.stages.scoring import score_groups
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
    for rows in (
        [],
        [row, row],
        [dict(row, error="OST failed")],
        [dict(row, rmsd=None)],
        [dict(row, n_assigned=2)],
    ):
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
        {"ensemble": "chai", "poses": [{"pose_id": "a", "path": "a.cif"}]}
    ]}))
    groups = load_inputs(manifest=manifest)
    assert [g["ensemble"] for g in groups] == ["boltz", "chai"]
    assert len(groups[0]["poses"]) == 2 and len(groups[1]["poses"]) == 1
    first, second = groups[0]["poses"]
    assert first["path"] == str(tmp_path / "a.cif")
    assert first["ligand_chains"] == ["B"] and first["smiles"] == "CCO"
    assert second["ligand_chains"] == ["C"] and second["smiles"] == "CCO"
    assert second["ligand_sdfs"] == [str(tmp_path / "ligand.sdf")]
    assert groups[0]["reference"]["path"] == str(tmp_path / "ref.cif")
    assert groups[0]["reference"]["ligand_chains"] == ["L"]


@pytest.mark.parametrize("value", [None, [], {"ensembles": {}},
    {"name": "model-input", "sequences": []},
    {"ensembles": [{"ensemble": "x", "poses": ["a.cif"]}]}])
def test_invalid_json_manifest_has_actionable_error(tmp_path, value):
    manifest = tmp_path / "input.json"
    manifest.write_text(json.dumps(value))
    with pytest.raises(ValueError, match="ensemble|poses"):
        load_inputs(manifest=manifest)


@pytest.mark.parametrize("explicit", [False, True])
def test_csv_manifest_cli_reports_json_migration(tmp_path, capsys, explicit):
    manifest = tmp_path / "input.csv"
    manifest.write_text("ensemble,pose_id,structure\ng,a,a.cif\n")
    args = (["--manifest"] if explicit else []) + [str(manifest), "--out", str(tmp_path / "out")]
    assert main(args) == 2
    assert "CSV manifests are no longer supported" in capsys.readouterr().err


def test_pb_unknown_never_passes():
    assert aggregate([True, None]) is None
    assert aggregate([True, False, None]) is False
    assert aggregate([]) is None


def test_pb_runs_and_invalid_poses_are_filtered(tmp_path):
    source = tmp_path / "protein.pdb"
    source.write_text("protein")
    prepared = dict(receptor=str(source))
    groups = [
        dict(
            ensemble="g",
            poses=[
                dict(ensemble="g", pose_id="a", prepared=prepared),
                dict(ensemble="g", pose_id="b", prepared=prepared),
            ],
        )
    ]
    answers = [
        dict(ensemble="g", pose_id="a", pb_valid=True),
        dict(ensemble="g", pose_id="b", pb_valid=False),
    ]
    with patch("core_score.stages.utils.posebusters_utils.run_pose", side_effect=answers) as pb:
        retained = run_pb(groups, tmp_path, filter_invalid=True)
        assert pb.call_count == 2
        assert [p["pose_id"] for p in retained[0]["poses"]] == ["a"]
        assert retained[0]["exclusions"][0]["pose_id"] == "b"


def test_unknown_pb_blocks_filtered_ensemble(tmp_path):
    source = tmp_path / "protein.pdb"
    source.write_text("protein")
    groups = [
        dict(
            ensemble="g",
            poses=[
                dict(ensemble="g", pose_id="a", prepared=dict(receptor=str(source)))
            ],
        )
    ]
    with patch(
        "core_score.stages.utils.posebusters_utils.run_pose",
        return_value=dict(ensemble="g", pose_id="a", pb_valid=None),
    ):
        result = run_pb(groups, tmp_path, filter_invalid=True)
    assert result[0]["pb_filter_error"]
    assert len(result[0]["poses"]) == 1


def test_skip_pb_retains_every_pose_without_running_checks(tmp_path):
    groups = [
        dict(
            ensemble="g",
            poses=[dict(ensemble="g", pose_id="a"), dict(ensemble="g", pose_id="b")],
        )
    ]
    with patch("core_score.stages.utils.posebusters_utils.run_pose") as pb:
        retained = run_pb(groups, tmp_path, skip=True)
    pb.assert_not_called()
    assert len(retained[0]["poses"]) == 2
    assert not retained[0]["pb_filter_error"]


def test_default_pipeline_automatically_filters(tmp_path):
    from core_score.pipeline import execute

    with patch(
        "core_score.pipeline.prepare_groups", return_value=[]
    ), patch("core_score.pipeline.run_pb", return_value=[]) as pb, patch(
        "core_score.pipeline.prepare_site_groups", return_value=[]
    ), patch(
        "core_score.pipeline.compare_groups", return_value=[]
    ), patch(
        "core_score.pipeline.score_groups", return_value=[]
    ):
        execute([], tmp_path / "default")
        assert pb.call_args.kwargs["skip"] is False
        assert pb.call_args.kwargs["filter_invalid"] is True
        execute([], tmp_path / "skip", skip_posebusters=True)
        assert pb.call_args.kwargs["skip"] is True
        assert pb.call_args.kwargs["filter_invalid"] is False


def test_supplied_pb_does_not_run_pb(tmp_path):
    path = tmp_path / "pb.csv"
    path.write_text("ensemble,pose_id,config,pb_valid\ng,a,dock,True\n")
    groups = [dict(ensemble="g", poses=[dict(ensemble="g", pose_id="a", prepared={})])]
    with patch("core_score.stages.utils.posebusters_utils.run_pose") as pb:
        run_pb(groups, tmp_path, supplied=path)
    pb.assert_not_called()
    with pytest.raises(ValueError):
        load_results(path, [("g", "b")])


def test_comparisons_are_computed_fresh_without_checkpoints(tmp_path):
    source = tmp_path / "source.cif"
    source.write_text("first")
    poses = [
        dict(
            pose_id=p,
            path=str(source),
            seed=0,
            sample=i,
            prepared=dict(path=str(source), copies=1, ligand_sdfs=[]),
        )
        for i, p in enumerate(["a", "b"])
    ]
    groups = [
        dict(
            ensemble="g", site="1", method="m", target="t", ligand_key="l", poses=poses
        )
    ]
    row = dict(
        ensemble="g", site="1", pose_i="a", pose_j="b", rmsd=1, n_assigned=1, error=""
    )
    with patch("core_score.stages.utils.comparison_utils.probe_prepared", side_effect=lambda p, **kw: p), patch(
        "core_score.stages.utils.comparison_utils.pair_job", return_value=(row, {})
    ) as compare:
        compare_groups(groups, tmp_path, 1)
        compare_groups(groups, tmp_path, 1)
        assert compare.call_count == 2
        assert not (tmp_path / "checkpoints").exists()
        source.write_text("changed")
        compare_groups(groups, tmp_path, 1)
        assert compare.call_count == 3


def test_reference_threshold_boundary_is_strict_and_does_not_change_score(tmp_path):
    from core_score.stages.evaluation import evaluate_groups

    pose = dict(
        ensemble="g",
        pose_id="p",
        path="pose.cif",
        seed=0,
        sample=0,
        prepared=dict(copies=1),
    )
    ref = dict(path="ref.cif", pose_id="reference", seed=0, sample=0)
    group = dict(
        ensemble="g",
        method="m",
        target="t",
        ligand_key="B",
        poses=[pose],
        reference=ref,
    )
    row = dict(ensemble="g", pose_i="p", pose_j="reference", rmsd=1.0, error="")
    with patch("core_score.stages.utils.evaluation_utils.prepare_pose", return_value=dict(copies=1)), patch(
        "core_score.stages.utils.evaluation_utils.probe_prepared", side_effect=lambda p, **kw: p
    ), patch("core_score.stages.utils.comparison_utils.pair_job", return_value=(row, {})) as compare:
        labels = evaluate_groups([group], tmp_path, [0.5, 1.0, 1.5, 2.0])
    assert compare.call_count == 1
    assert [r["success"] for r in labels] == [False, False, True, True]
