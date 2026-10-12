"""Known site identities and negative controls; no experimental reference needed."""

import itertools
import unittest
import numpy as np
from core_score.stages.utils.site_algorithms import (
    protein_transform,
    assign_centres,
    group_sites,
    same_partition,
    group_sites_multistart,
)


class SiteTests(unittest.TestCase):
    def test_historical_assignment_does_not_accept_nonconvergence(self):
        from core_score.stages.utils.site_algorithms import assign_sites

        with self.assertRaisesRegex(ValueError, "converge"):
            assign_sites(
                {
                    "a": {
                        "B": np.array([0.0, 0.0, 0.0]),
                        "C": np.array([10.0, 0.0, 0.0]),
                    }
                },
                2,
                max_iter=1,
            )

    def test_collinear_fit_rejected(self):
        with self.assertRaises(ValueError):
            protein_transform(
                [[0, 0, 0], [1, 0, 0], [2, 0, 0]], [[0, 0, 0], [1, 0, 0], [2, 0, 0]]
            )

    def test_assignment_and_gap_match_exhaustive_oracle(self):
        rng = np.random.default_rng(77)
        for n, k in [(2, 2), (3, 3), (2, 4), (4, 5)]:
            points = rng.normal(size=(n, 3))
            centres = rng.normal(size=(k, 3))
            costs = sorted(
                sum(np.sum((points[i] - centres[j]) ** 2) for i, j in enumerate(p))
                for p in itertools.permutations(range(k), n)
            )
            result = assign_centres(points, centres)
            self.assertAlmostEqual(sum(result.distances**2), costs[0])
            self.assertAlmostEqual(result.alternative_gap, costs[1] - costs[0])

    def test_more_than_five_copies(self):
        centres = np.array([[10.0 * i, 0.0, 0.0] for i in range(12)])
        order = np.array([8, 2, 0, 11, 6, 5, 9, 4, 10, 3, 7, 1])
        result = assign_centres(centres[order] + [0.1, 0, 0], centres)
        np.testing.assert_array_equal(result.sites, order)

    def test_missing_copy_not_manufactured(self):
        result = assign_centres([[9.0, 0.0, 0.0]], [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]])
        self.assertEqual(result.sites.tolist(), [1])
        self.assertEqual(len(result.distances), 1)

    def test_ambiguous_and_distant_copies_are_visible(self):
        result = assign_centres(
            [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]], [[-10.0, 0.0, 0.0], [10.0, 0.0, 0.0]]
        )
        self.assertEqual(result.alternative_gap, 0.0)
        far = assign_centres(
            [[0.0, 0.0, 0.0], [100.0, 0.0, 0.0]], [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]]
        )
        self.assertEqual(far.distances.tolist(), [0.0, 90.0])

    def test_known_sites_survive_copy_and_pose_reordering(self):
        poses = {
            "a": np.array([[-10.0, 0, 0], [10.0, 0, 0]]),
            "b": np.array([[11.0, 0, 0], [-11.0, 0, 0]]),
            "c": np.array([[-9.0, 0, 0], [9.0, 0, 0]]),
        }
        original = {key: value.copy() for key, value in poses.items()}
        for anchor in poses:
            result, _ = group_sites(poses, poses[anchor])
            labels = {
                (p, i): int(a.sites[i]) for p, a in result.items() for i in range(2)
            }
            expected = {
                ("a", 0): 0,
                ("a", 1): 1,
                ("b", 0): 1,
                ("b", 1): 0,
                ("c", 0): 0,
                ("c", 1): 1,
            }
            self.assertTrue(same_partition(labels, expected))
        reordered = {key: value[::-1] for key, value in reversed(list(poses.items()))}
        result, _ = group_sites(reordered, reordered["b"])
        labels = {
            (p, 1 - i): int(a.sites[i]) for p, a in result.items() for i in range(2)
        }
        self.assertTrue(same_partition(labels, expected))
        for key in poses:
            np.testing.assert_array_equal(poses[key], original[key])

    def test_local_site_swap_is_not_global_renaming(self):
        first = {"p1a": 0, "p1b": 1, "p2a": 0, "p2b": 1}
        self.assertTrue(same_partition(first, {k: 1 - v for k, v in first.items()}))
        self.assertFalse(same_partition(first, dict(first, p2a=1, p2b=0)))

    def test_multistart_keeps_lowest_tested_objective(self):
        poses = {
            "a": np.array([[-10.0, 0, 0], [10.0, 0, 0]]),
            "b": np.array([[11.0, 0, 0], [-11.0, 0, 0]]),
        }
        starts = [poses["a"], poses["b"]]
        result = group_sites_multistart(poses, starts)
        self.assertAlmostEqual(result.objective, 1.0)
        self.assertEqual(result.starts, 2)
        self.assertTrue(np.isinf(result.alternative_objective_gap))
        reverse = group_sites_multistart(poses, starts[::-1])
        self.assertAlmostEqual(result.objective, reverse.objective)

    def test_multistart_empty_starts_rejected(self):
        with self.assertRaises(ValueError):
            group_sites_multistart({"a": np.zeros((2, 3))}, [])

    def test_features_need_matching_dimensions(self):
        with self.assertRaises(ValueError):
            assign_centres([[0.0, 1.0]], [[0.0, 1.0, 2.0]])

    def test_distance_features_invariant_to_rigid_motion(self):
        protein = np.array(
            [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 4.0]]
        )
        ligands = np.array([[-10.0, 0, 0], [10.0, 0, 0]])
        rotation = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
        features = np.linalg.norm(ligands[:, None] - protein[None], axis=2)
        moved_protein = protein @ rotation + [12.0, -4.0, 6.0]
        moved_ligands = ligands[::-1] @ rotation + [12.0, -4.0, 6.0]
        moved_features = np.linalg.norm(
            moved_ligands[:, None] - moved_protein[None], axis=2
        )
        np.testing.assert_array_equal(
            assign_centres(moved_features, features).sites, [1, 0]
        )


if __name__ == "__main__":
    unittest.main()
