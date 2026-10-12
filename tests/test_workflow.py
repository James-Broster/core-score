import json
from unittest.mock import Mock, patch
import numpy as np
import pytest
from core_score import core_score
from core_score.stages.scoring import score_groups
from core_score.stages.utils.score_utils import build_matrix
from core_score.cli import parser
from core_score.stages.utils.posebusters_utils import aggregate
from core_score.stages.posebusters import run_posebusters
from core_score.inputs import load_inputs


def test_core_known_two_pose_score_and_diagonal_independence():
    for diag in (0.0, float("nan"), 99.0):
        assert core_score([[diag, 1.0], [1.0, diag]]) == 0.5


def test_core_uses_each_pose_closest_other_pose_capped_at_threshold():
    matrix = [[0, 0.1, 10], [0.1, 0, 10], [10, 10, 0]]
    # Penalties 0.1/T, 0.1/T and 1 (no other pose within T).
    assert core_score(matrix) == pytest.approx(1 - (0.05 + 0.05 + 1) / 3)
    assert core_score(matrix, threshold=1.0) == pytest.approx(1 - (0.1 + 0.1 + 1) / 3)


def test_core_endpoints():
    assert core_score(np.zeros((4, 4))) == 1.0
    assert core_score(np.full((3, 3), 2.0)) == 0.0  # Exactly T apart counts as unsupported.
    with pytest.raises(ValueError):
        core_score([[0]])


def test_scores_at_each_threshold_and_never_from_incomplete_matrices(tmp_path):
    poses = [dict(pose_id=p) for p in "abc"]
    rows = [dict(ensemble="g", site="1", pose_i=a, pose_j=b, rmsd=r, n_assigned=1, error="")
            for a, b, r in [("a", "b", 0.5), ("a", "c", 3.0), ("b", "c", 1.0)]]
    groups = [dict(ensemble="g", site="1", poses=poses), dict(ensemble="h", site="1", poses=poses),
              dict(ensemble="u", site="unresolved", poses=poses, error=dict(error="ValueError: no sites"))]
    scores = score_groups(groups, rows, tmp_path, thresholds=(1.0, 2.0))
    assert [(s["ensemble"], s["threshold_angstrom"]) for s in scores] == [
        ("g", 1.0), ("g", 2.0), ("h", 1.0), ("h", 2.0), ("u", 1.0), ("u", 2.0)]
    assert scores[1]["core_score"] == pytest.approx(1 - (0.25 + 0.25 + 0.5) / 3)
    assert all(s["core_score"] is None and s["error"] for s in scores[2:])


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


def test_posebusters_records_every_pose(tmp_path):
    groups = [dict(ensemble="g", poses=[dict(ensemble="g", pose_id="a"), dict(ensemble="g", pose_id="b")])]
    answers = [dict(ensemble="g", pose_id="a", pb_valid=True), dict(ensemble="g", pose_id="b", pb_valid=False)]
    with patch("posebusters.PoseBusters", return_value=Mock(config={"modules": []})), \
            patch("core_score.stages.posebusters.check_pose", side_effect=answers) as pb:
        records = run_posebusters(groups, tmp_path)
    assert pb.call_count == 2
    assert [r["pb_valid"] for r in records] == [True, False]
    assert json.loads((tmp_path / "posebusters.json").read_text()) == answers


def test_pipeline_scores_every_pose_and_can_skip_posebusters(tmp_path):
    from core_score.pipeline import execute

    stages = ["prepare_groups", "prepare_site_groups", "compare_groups", "score_groups"]
    patches = [patch(f"core_score.pipeline.{name}", return_value=[]) for name in stages]
    for p in patches:
        p.start()
    try:
        with patch("core_score.pipeline.run_posebusters") as pb:
            execute([], tmp_path / "default")
            assert pb.call_count == 1
            execute([], tmp_path / "skip", skip_posebusters=True)
            assert pb.call_count == 1
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
        labels = evaluate_groups([group], [], tmp_path, [0.5, 1.0, 1.5, 2.0])
    assert compare.call_count == 1
    assert [r["success"] for r in labels] == [False, False, True, True]


def test_top_pose_is_the_valid_pose_with_highest_confidence():
    from core_score.stages.utils.decision_utils import top_pose

    group = dict(poses=[dict(pose_id="a", confidence=0.9), dict(pose_id="b", confidence=0.7),
                        dict(pose_id="c", confidence=0.7), dict(pose_id="d", confidence=0.95)])
    validity = dict(a=False, b=True, c=True, d=None)
    assert top_pose(group, validity)[0]["pose_id"] == "b"  # a invalid, d unknown, b wins the tie
    assert top_pose(group, None) == (None, "PoseBusters was skipped")
    assert top_pose(group, dict(a=False))[1] == "no PoseBusters-valid pose"
    assert top_pose(dict(poses=[dict(pose_id="a")]), dict(a=True))[1] == "no confidences given"


def test_decisions_per_cutoff_and_summary():
    from core_score.stages.utils.decision_utils import decide, summarize_cutoffs

    groups = [dict(ensemble="g", poses=[dict(pose_id="a", confidence=0.8), dict(pose_id="b", confidence=0.6)]),
              dict(ensemble="h", poses=[dict(pose_id="a", confidence=0.9)]),
              dict(ensemble="u", poses=[dict(pose_id="a", confidence=0.9)])]
    scores = [dict(ensemble="g", site="1", threshold_angstrom=2.0, core_score=0.9, error=None),
              dict(ensemble="h", site="1", threshold_angstrom=2.0, core_score=0.6, error=None),
              dict(ensemble="u", site="unresolved", threshold_angstrom=2.0, core_score=None, error="ValueError: x")]
    validity = dict(g=dict(a=True, b=True), h=dict(a=True), u=dict(a=True))
    labels = [dict(ensemble="g", pose_id="a", threshold_angstrom=2.0, success=True),
              dict(ensemble="h", pose_id="a", threshold_angstrom=2.0, success=False)]
    rows = decide(groups, scores, validity, labels, [0.5, 0.7])
    assert [(r["ensemble"], r["core_cutoff"], r["likely_correct"]) for r in rows] == [
        ("g", 0.5, True), ("g", 0.7, True), ("h", 0.5, True), ("h", 0.7, False), ("u", 0.5, None), ("u", 0.7, None)]
    assert rows[0]["top_pose"] == "a" and rows[0]["top_pose_correct"] is True
    summary = summarize_cutoffs(rows)
    assert [(s["core_cutoff"], s["scored"], s["accepted"], s["accuracy"]) for s in summary] == [
        (0.5, 2, 2, 0.5), (0.7, 2, 1, 1.0)]


def test_confidence_is_given_for_every_pose_or_none(tmp_path):
    for name in ("a.cif", "b.cif"):
        (tmp_path / name).write_text("x")
    manifest = tmp_path / "poses.json"
    manifest.write_text(json.dumps({"ensembles": [{"ensemble": "g", "smiles": "CCO", "poses": [
        {"path": "a.cif", "confidence": 0.8}, {"path": "b.cif"}]}]}))
    with pytest.raises(ValueError, match="confidence"):
        load_inputs([manifest])


def test_core_cutoffs_default_and_range():
    assert parser().parse_args(["a.cif", "--out", "o"]).core_cutoffs == [0.5, 0.6, 0.7, 0.8, 0.9]
    with pytest.raises(SystemExit):
        parser().parse_args(["a.cif", "--out", "o", "--core-cutoffs", "1.5"])
