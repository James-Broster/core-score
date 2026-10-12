import unittest
from unittest.mock import patch
from core_score.stages.utils.radius import score_with_radius_fallback


class RadiusPolicyTests(unittest.TestCase):
    def run_case(self, answers, radii=(6.0, 8.0)):
        calls = []

        def evaluate(r):
            calls.append(r)
            return answers[r]

        result, metadata = score_with_radius_fallback(evaluate, radii)
        return result, metadata, calls

    def test_numeric_standard_unchanged(self):
        r, m, c = self.run_case({4.0: dict(score=1.234567890123, reason="OK")})
        self.assertEqual(c, [4.0])
        self.assertEqual(r["score"], 1.234567890123)
        self.assertEqual(m["standard_rmsd"], r["score"])
        self.assertFalse(m["radius_fallback_used"])

    def test_first_success_stops_at_six(self):
        r, m, c = self.run_case(
            {
                4.0: dict(score=None, reason=["target_binding_site", "empty"]),
                6.0: dict(score=3.0, reason="OK"),
            }
        )
        self.assertEqual(c, [4.0, 6.0])
        self.assertEqual(m["binding_site_radius_angstrom"], 6.0)
        self.assertTrue(m["requires_review"])
        self.assertEqual(m["standard_rmsd"], "")

    def test_eight_only_after_six_failed(self):
        bad = dict(score=None, target_report=[{"short desc": "target_binding_site"}])
        r, m, c = self.run_case({4.0: bad, 6.0: bad, 8.0: dict(score=2.0, reason="OK")})
        self.assertEqual(c, [4.0, 6.0, 8.0])
        self.assertEqual(m["binding_site_radius_angstrom"], 8.0)

    def test_graph_failure_not_retried(self):
        r, m, c = self.run_case(
            {4.0: dict(score=None, reason=["identity", "graph mismatch"])}
        )
        self.assertEqual(c, [4.0])
        self.assertFalse(m["radius_fallback_used"])

    def test_failed_retries_stay_missing(self):
        bad = dict(score=None, reason="target_binding_site")
        r, m, c = self.run_case({4.0: bad, 6.0: bad, 8.0: bad})
        self.assertIsNone(r["score"])
        self.assertEqual(m["binding_site_radius_angstrom"], "")
        self.assertFalse(m["radius_fallback_used"])

    def test_disabled(self):
        r, m, c = self.run_case(
            {4.0: dict(score=None, reason="target_binding_site")}, ()
        )
        self.assertEqual(c, [4.0])

    def test_bad_configuration(self):
        for radii in ((8, 7), (7, 7), (3,), (float("nan"),)):
            with self.assertRaises(ValueError):
                score_with_radius_fallback(lambda r: None, radii)

    def test_public_comparison_retains_original_contact_metric(self):
        from core_score.stages.utils.comparison_utils import compare_pair

        baseline = dict(
            status="SUCCESS",
            rmsd=dict(
                assigned_scores=[],
                model_ligand_unassigned_reason={"B": ["binding_site", "empty"]},
            ),
            lddt_pli=dict(assigned_scores=[dict(score=0.5)]),
        )
        retry = dict(
            status="SUCCESS",
            rmsd=dict(assigned_scores=[dict(score=2.0, bb_rmsd=1.0, lddt_lp=0.8)]),
            lddt_pli=dict(assigned_scores=[dict(score=0.9)]),
        )
        model = dict(path="m.cif", copies=1)
        reference = dict(path="r.cif", copies=1)
        with patch(
            "core_score.stages.utils.comparison_utils.run_ost", side_effect=[baseline, retry]
        ) as run:
            row, detail = compare_pair(
                (
                    dict(method="m", target="t", ligand_key="B", n_input_copies=1),
                    dict(seed=1, sample=0, path="m.cif"),
                    dict(seed=1, sample=1, path="r.cif"),
                    model,
                    reference,
                )
            )
        self.assertEqual(row["rmsd"], 2.0)
        self.assertEqual(row["lddt_pli"], 0.5)
        self.assertEqual([c.kwargs["radius"] for c in run.call_args_list], [4.0, 6.0])
        self.assertEqual(len(detail["ost_results"]), 2)


if __name__ == "__main__":
    unittest.main()
