import json
from unittest.mock import Mock, patch
import numpy as np
import pytest
from core_score import score as core_score
from core_score.analysis.core import score_groups
from core_score.analysis.core import build_matrix
from core_score.cli import parser
from core_score.analysis.posebusters import aggregate
from core_score.analysis.posebusters import run_posebusters


def test_core_known_score_and_diagonal_independence():
    for diag in (0.0, float("nan"), 99.0):
        matrix = np.full((4, 4), 1.0)
        np.fill_diagonal(matrix, diag)
        assert core_score(matrix) == 0.5


def test_core_averages_each_pose_three_closest_other_poses_capped_at_threshold():
    # Four poses within 0.2 A of each other and one far away.
    matrix = np.full((5, 5), 10.0)
    matrix[:4, :4] = 0.2
    np.fill_diagonal(matrix, 0)
    # The four agree: three neighbours at 0.2 A, penalty 0.1 each. The fifth has none within T.
    assert core_score(matrix) == pytest.approx(1 - (4 * 0.1 + 1) / 5)
    assert core_score(matrix, threshold=1.0) == pytest.approx(1 - (4 * 0.2 + 1) / 5)
    # A pose with only one close partner is still penalised for the other two neighbours.
    pair = [[0, 0.1, 10, 10], [0.1, 0, 10, 10], [10, 10, 0, 10], [10, 10, 10, 0]]
    assert core_score(pair) == pytest.approx(1 - (2 * (0.05 + 1 + 1) / 3 + 2) / 4)
    assert core_score(pair, neighbours=1) == pytest.approx(1 - (0.05 + 0.05 + 1 + 1) / 4)


def test_core_needs_four_poses_and_never_uses_fewer_neighbours():
    with pytest.raises(ValueError, match="at least 4 poses"):
        core_score([[0, 0.1, 10], [0.1, 0, 10], [10, 10, 0]])


def test_core_endpoints():
    assert core_score(np.zeros((4, 4))) == 1.0
    assert core_score(np.full((4, 4), 2.0)) == 0.0  # Exactly T apart counts as unsupported.
    with pytest.raises(ValueError):
        core_score([[0]])


def test_scores_at_each_threshold_and_never_from_incomplete_matrices(tmp_path):
    poses = [dict(pose_id=p) for p in "abcd"]
    rows = [dict(ensemble="g", site="1", pose_i=a, pose_j=b, rmsd=r, n_assigned=1, error="")
            for a, b, r in [("a", "b", 0.5), ("a", "c", 3.0), ("b", "c", 1.0),
                            ("a", "d", 3.0), ("b", "d", 3.0), ("c", "d", 3.0)]]
    groups = [dict(ensemble="g", site="1", name="g", poses=poses),
              dict(ensemble="h", site="1", name="h", poses=poses),
              dict(ensemble="u", site="unresolved", name="u", poses=poses, error=dict(error="ValueError: no sites"))]
    scores = score_groups(groups, rows, thresholds=(1.0, 2.0))
    assert [(s["ensemble"], s["threshold_angstrom"]) for s in scores] == [
        ("g", 1.0), ("g", 2.0), ("h", 1.0), ("h", 2.0), ("u", 1.0), ("u", 2.0)]
    # T = 2 A, three other poses each: a (0.5, 3, 3), b (0.5, 1, 3), c (1, 3, 3), d (3, 3, 3).
    penalties = [(0.25 + 1 + 1) / 3, (0.25 + 0.5 + 1) / 3, (0.5 + 1 + 1) / 3, 1]
    assert scores[1]["core_score"] == pytest.approx(1 - sum(penalties) / 4)
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


def test_pb_unknown_never_passes():
    assert aggregate([True, None]) is None
    assert aggregate([True, False, None]) is False
    assert aggregate([]) is None


def test_posebusters_records_every_pose(tmp_path):
    groups = [dict(ensemble="g", poses=[dict(ensemble="g", pose_id="a"), dict(ensemble="g", pose_id="b")])]
    answers = [dict(ensemble="g", pose_id="a", pb_valid=True), dict(ensemble="g", pose_id="b", pb_valid=False)]
    with patch("posebusters.PoseBusters", return_value=Mock(config={"modules": []})), \
            patch("core_score.analysis.posebusters.check_pose", side_effect=answers) as pb:
        records = run_posebusters(groups, tmp_path)
    assert pb.call_count == 2
    assert [r["pb_valid"] for r in records] == [True, False]
    assert json.loads((tmp_path / "posebusters.json").read_text()) == answers


def test_each_check_combines_the_copies_and_unknown_never_passes():
    from core_score.analysis.posebusters import check_pose

    pose = dict(ensemble="g", pose_id="a", prepared=dict(receptor="r.pdb", copy_records=["c0", "c1"]))
    copies = [dict(copy_id="0", pb_valid=False, checks=dict(bond_lengths=True, internal_energy=False)),
              dict(copy_id="1", pb_valid=None, error="ValueError: x")]
    with patch("core_score.analysis.posebusters.check_copy", side_effect=copies):
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
    from core_score.analysis.reference import correct_column, evaluate_groups

    pose = dict(ensemble="g", pose_id="p", path="pose.cif", prepared=dict(copies=1))
    reference = dict(ensemble="g", pose_id="reference", path="ref.cif")
    group = dict(ensemble="g", poses=[pose], reference=reference)
    row = dict(ensemble="g", site="reference", pose_i="p", pose_j="reference", rmsd=1.0, lddt_pli=0.9, error="")
    prepared_reference = dict(reference, prepared=dict(copies=1))
    with patch("core_score.analysis.reference.prepare", return_value=prepared_reference), \
            patch("core_score.analysis.pairs.pair_job", return_value=(row, {})) as compare:
        labels = evaluate_groups([group], tmp_path, [0.5, 1.0, 1.5, 2.0])
    assert compare.call_count == 1
    assert len(labels) == 1  # one row per pose, one column per threshold
    assert [labels[0][correct_column(t)] for t in (0.5, 1.0, 1.5, 2.0)] == [False, False, True, True]


def test_correct_needs_rmsd_below_threshold_and_lddt_pli_above_0_8():
    from core_score.analysis.reference import label_reference_results

    records = [dict(ensemble="g", pose_id=p, rmsd=r, lddt_pli=l, error="")
               for p, r, l in [("a", 1.0, 0.9), ("b", 1.0, 0.7), ("c", 1.0, 0.8), ("d", 3.0, 0.9), ("e", 1.0, None)]]
    labels = label_reference_results(records, [2.0])
    assert [label["correct_2A"] for label in labels] == [True, False, False, False, None]


def test_failed_reference_comparison_is_unknown_not_wrong():
    from core_score.analysis.reference import label_reference_results

    [label] = label_reference_results([dict(ensemble="g", pose_id="p", rmsd=0.5, error="ambiguous")], [2.0])
    assert (label["rmsd"], label["correct_2A"], label["error"]) == (None, None, "ambiguous")
    [label] = label_reference_results([dict(ensemble="g", pose_id="p", error=dict(error="ValueError: x"))], [2.0])
    assert label["error"] == "ValueError: x"


def test_top_pose_is_the_valid_pose_with_highest_confidence():
    from core_score.outputs.report import top_pose

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
    scores = [dict(ensemble="h", site="1", name="h", pose_count=1, threshold_angstrom=2.0, core_score=0.6, error=None),
              dict(ensemble="g", site="1", name="g", pose_count=2, threshold_angstrom=2.0, core_score=0.9, error=None),
              dict(ensemble="u", site="unresolved", name="u", pose_count=1, threshold_angstrom=2.0, core_score=None,
                   error="ValueError: x")]
    validity = {("g", "a"): True, ("g", "b"): False, ("h", "a"): True, ("u", "a"): True}
    labels = [dict(ensemble="g", pose_id="a", correct_2A=True),
              dict(ensemble="g", pose_id="b", correct_2A=False),
              dict(ensemble="h", pose_id="a", correct_2A=False)]
    return groups, scores, validity, labels


def test_core_table_has_one_row_per_complex_best_first():
    from core_score.outputs.report import core_table

    groups, scores, validity, labels = report_inputs()
    rows = core_table(groups, scores, validity, labels, [0.5, 0.7], ranked=True)
    assert [(r["ensemble"], r["above_0.5"], r["above_0.7"]) for r in rows] == [
        ("g", True, True), ("h", True, False), ("u", None, None)]
    g = rows[0]
    assert (g["valid_pose_count"], g["top_pose"], g["top_pose_correct"], g["fraction_correct"]) == (1, "a", True, 0.5)
    assert rows[2]["note"] == "ValueError: x"


def test_core_table_names_each_site_of_a_multicopy_ligand():
    from core_score.outputs.report import core_table

    groups = [dict(ensemble="g", poses=[dict(ensemble="g", pose_id="a")])]
    scores = [dict(ensemble="g", site=site, name=f"g_site_{site}", pose_count=4, threshold_angstrom=2.0,
                   core_score=core, error=None) for site, core in (("1", 0.9), ("2", 0.4))]
    rows = core_table(groups, scores, None, [], [0.5], ranked=False)
    assert [(r["ensemble"], r["core_score"]) for r in rows] == [("g_site_1", 0.9), ("g_site_2", 0.4)]


def test_report_files_only_have_the_columns_that_apply(tmp_path):
    import csv
    from core_score.outputs.report import write_report

    groups, scores, _, _ = report_inputs()
    for group in groups:
        for pose in group["poses"]:
            del pose["confidence"]
    write_report(groups, scores, None, [], tmp_path, cutoffs=(0.5,))
    header = next(csv.reader((tmp_path / "core_scores.csv").open()))
    assert header == ["ensemble", "threshold_angstrom", "core_score", "above_0.5", "pose_count", "note"]
    assert next(csv.reader((tmp_path / "poses.csv").open())) == ["ensemble", "pose_id", "note"]
    assert sorted(path.name for path in tmp_path.iterdir()) == ["core_scores.csv", "poses.csv"]


def test_poses_table_joins_confidence_validity_and_reference():
    from core_score.outputs.report import pose_table

    groups, _, _, labels = report_inputs()
    labels = [dict(label, rmsd=0.4, lddt_pli=0.9, radius_fallback_used=True) for label in labels]
    pb_records = [dict(ensemble="g", pose_id="a", pb_valid=True, checks=dict(bond_lengths=True)),
                  dict(ensemble="g", pose_id="b", pb_valid=False,
                       checks=dict(bond_lengths=False, internal_energy=None))]
    groups[2]["poses"][0]["error"] = dict(error="ValueError: no ligand")
    rows = pose_table(groups, pb_records, labels, ranked=True)
    assert [(r["ensemble"], r["pose_id"]) for r in rows] == [
        ("h", "a"), ("g", "a"), ("g", "b"), ("u", "a")]
    assert rows[1] == dict(ensemble="g", pose_id="a", confidence=0.8, pb_valid=True,
                           rmsd=0.4, lddt_pli=0.9, correct_2A=True, note=None)
    assert rows[2]["pb_valid"] is False
    assert rows[0]["pb_valid"] is None  # PoseBusters has no record for this pose
    assert rows[3]["note"] == "ValueError: no ligand"


def test_plot_dots_shrink_as_the_benchmark_grows():
    from core_score.outputs.plot import dot_size

    assert dot_size(6) > dot_size(900) > dot_size(7000)


def test_plot_writes_one_png_per_threshold(tmp_path):
    pytest.importorskip("matplotlib")
    from core_score.outputs.plot import plot_core_scores

    rows = [dict(threshold_angstrom=t, core_score=c)
            for t in (1.0, 2.0) for c in (0.4, 0.9, None)]
    plot_core_scores(rows, tmp_path, cutoffs=(0.8,))
    assert sorted(path.name for path in tmp_path.iterdir()) == ["core_scores_1A.png", "core_scores_2A.png"]


def test_core_cutoffs_default_and_range():
    assert parser().parse_args(["a.cif", "--out", "o"]).core_cutoffs == [0.5, 0.6, 0.7, 0.8, 0.9]
    with pytest.raises(SystemExit):
        parser().parse_args(["a.cif", "--out", "o", "--core-cutoffs", "1.5"])
