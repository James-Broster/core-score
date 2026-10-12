import json
from unittest.mock import Mock, patch
import numpy as np
import pytest
from core_score import core_score
from core_score.stages.scoring import score_groups
from core_score.stages.scoring import build_matrix
from core_score.cli import parser
from core_score.stages.posebusters import aggregate
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
    scores = score_groups(groups, rows, thresholds=(1.0, 2.0))
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


def test_each_check_combines_the_copies_and_unknown_never_passes():
    from core_score.stages.posebusters import check_pose

    pose = dict(ensemble="g", pose_id="a", prepared=dict(receptor="r.pdb", copy_records=["c0", "c1"]))
    copies = [dict(copy_id="0", pb_valid=False, checks=dict(bond_lengths=True, internal_energy=False)),
              dict(copy_id="1", pb_valid=None, error="ValueError: x")]
    with patch("core_score.stages.posebusters.check_copy", side_effect=copies):
        record = check_pose(None, ["bond_lengths", "internal_energy"], pose)
    assert record["pb_valid"] is False
    assert record["checks"] == dict(bond_lengths=None, internal_energy=False)
    unprepared = check_pose(None, ["bond_lengths"], dict(ensemble="g", pose_id="b", error="bad input"))
    assert unprepared["checks"] == dict(bond_lengths=None)


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
    from core_score.stages.evaluation import correct_column, evaluate_groups

    pose = dict(ensemble="g", pose_id="p", path="pose.cif", prepared=dict(copies=1))
    reference = dict(ensemble="g", pose_id="reference", path="ref.cif")
    group = dict(ensemble="g", poses=[pose], reference=reference)
    row = dict(ensemble="g", site="reference", pose_i="p", pose_j="reference", rmsd=1.0, error="")
    prepared_reference = dict(reference, prepared=dict(copies=1))
    with patch("core_score.stages.evaluation.prepare", return_value=prepared_reference), \
            patch("core_score.stages.comparison.pair_job", return_value=(row, {})) as compare:
        labels = evaluate_groups([group], tmp_path, [0.5, 1.0, 1.5, 2.0])
    assert compare.call_count == 1
    assert len(labels) == 1  # one row per pose, one column per threshold
    assert [labels[0][correct_column(t)] for t in (0.5, 1.0, 1.5, 2.0)] == [False, False, True, True]


def test_failed_reference_comparison_is_unknown_not_wrong():
    from core_score.stages.evaluation import label_reference_results

    [label] = label_reference_results([dict(ensemble="g", pose_id="p", rmsd=0.5, error="ambiguous")], [2.0])
    assert (label["rmsd"], label["rmsd_below_2A"], label["error"]) == (None, None, "ambiguous")
    [label] = label_reference_results([dict(ensemble="g", pose_id="p", error=dict(error="ValueError: x"))], [2.0])
    assert label["error"] == "ValueError: x"


def test_top_pose_is_the_valid_pose_with_highest_confidence():
    from core_score.stages.report import top_pose

    poses = [dict(ensemble="g", pose_id=p, confidence=c) for p, c in zip("abcd", (0.9, 0.7, 0.7, 0.95))]
    group = dict(ensemble="g", poses=poses)
    validity = {("g", "a"): False, ("g", "b"): True, ("g", "c"): True, ("g", "d"): None}
    assert top_pose(group, validity)[0]["pose_id"] == "b"  # a invalid, d unknown, b wins the tie
    assert top_pose(group, None) == (None, "PoseBusters was skipped")
    assert top_pose(group, {("g", "a"): False})[1] == "no PoseBusters-valid pose"
    assert top_pose(dict(poses=[dict(ensemble="g", pose_id="a")]), {("g", "a"): True})[1] == "no confidences given"


def report_inputs():
    groups = [dict(ensemble=e, poses=[dict(ensemble=e, pose_id=p, confidence=c) for p, c in poses])
              for e, poses in (("h", [("a", 0.9)]), ("g", [("a", 0.8), ("b", 0.6)]), ("u", [("a", 0.9)]))]
    scores = [dict(ensemble="h", site="1", pose_count=1, threshold_angstrom=2.0, core_score=0.6, error=None),
              dict(ensemble="g", site="1", pose_count=2, threshold_angstrom=2.0, core_score=0.9, error=None),
              dict(ensemble="u", site="unresolved", pose_count=1, threshold_angstrom=2.0, core_score=None,
                   error="ValueError: x")]
    validity = {("g", "a"): True, ("g", "b"): False, ("h", "a"): True, ("u", "a"): True}
    labels = [dict(ensemble="g", pose_id="a", rmsd_below_2A=True),
              dict(ensemble="g", pose_id="b", rmsd_below_2A=False),
              dict(ensemble="h", pose_id="a", rmsd_below_2A=False)]
    return groups, scores, validity, labels


def test_core_table_has_one_row_per_complex_best_first():
    from core_score.stages.report import core_table

    groups, scores, validity, labels = report_inputs()
    rows = core_table(groups, scores, validity, labels, [0.5, 0.7], ranked=True)
    assert [(r["ensemble"], r["above_0.5"], r["above_0.7"]) for r in rows] == [
        ("g", True, True), ("h", True, False), ("u", None, None)]
    g = rows[0]
    assert (g["valid_pose_count"], g["top_pose"], g["top_pose_correct"], g["fraction_correct"]) == (1, "a", True, 0.5)
    assert rows[2]["note"] == "ValueError: x"


def test_report_files_only_have_the_columns_that_apply(tmp_path):
    import csv
    from core_score.stages.report import write_report

    groups, scores, _, _ = report_inputs()
    for group in groups:
        for pose in group["poses"]:
            del pose["confidence"]
    write_report(groups, scores, None, [], tmp_path, cutoffs=(0.5,))
    header = next(csv.reader((tmp_path / "core_scores.csv").open()))
    assert header == ["ensemble", "site", "threshold_angstrom", "core_score", "above_0.5", "pose_count", "note"]
    assert next(csv.reader((tmp_path / "poses.csv").open())) == ["ensemble", "pose_id", "note"]
    assert sorted(path.name for path in tmp_path.iterdir()) == ["core_scores.csv", "poses.csv"]


def test_poses_table_joins_confidence_validity_and_reference():
    from core_score.stages.report import pose_table

    groups, _, _, labels = report_inputs()
    labels = [dict(label, rmsd=0.4, lddt_pli=0.9, radius_fallback_used=True) for label in labels]
    pb_records = [dict(ensemble="g", pose_id="a", pb_valid=True, checks=dict(bond_lengths=True)),
                  dict(ensemble="g", pose_id="b", pb_valid=False,
                       checks=dict(bond_lengths=False, internal_energy=None))]
    groups[2]["poses"][0]["error"] = dict(error="ValueError: no ligand")
    rows = pose_table(groups, pb_records, labels, ranked=True)
    assert [(r["ensemble"], r["pose_id"]) for r in rows] == [
        ("h", "a"), ("g", "a"), ("g", "b"), ("u", "a")]
    assert rows[1] == dict(ensemble="g", pose_id="a", confidence=0.8, pb_valid=True, failed_checks="",
                           rmsd=0.4, lddt_pli=0.9, rmsd_below_2A=True, note=None)
    assert rows[2]["failed_checks"] == "bond_lengths"  # unknown checks are not failures
    assert rows[0]["pb_valid"] is None  # PoseBusters has no record for this pose
    assert rows[3]["note"] == "ValueError: no ligand"


def test_plot_titles_count_the_complexes_above_each_cutoff():
    from core_score.stages.plotting import title

    assert title([0.95, 0.85, 0.4], 2.0, (0.8, 0.9)) == (
        "n = 3 complexes \u00b7 threshold 2 \u00c5\n2 at or above 0.8  \u00b7  1 at or above 0.9")
    assert title([0.5], 1.0, ()) == "n = 1 complexes \u00b7 threshold 1 \u00c5"


def test_plot_writes_one_png_per_threshold(tmp_path):
    pytest.importorskip("matplotlib")
    from core_score.stages.plotting import plot_core_scores

    rows = [dict(threshold_angstrom=t, core_score=c)
            for t in (1.0, 2.0) for c in (0.4, 0.9, None)]
    plot_core_scores(rows, tmp_path, cutoffs=(0.8,))
    assert sorted(path.name for path in tmp_path.iterdir()) == ["core_scores_1A.png", "core_scores_2A.png"]


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
