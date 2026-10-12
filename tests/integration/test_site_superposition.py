import unittest
import numpy as np
from core_score.stages.utils.site_algorithms import protein_transform


class SiteSuperpositionTests(unittest.TestCase):
    def test_rigid_motion_recovered_without_mutating_input(self):
        protein = np.array(
            [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 4.0]]
        )
        rotation = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
        moving = protein @ rotation + [12.0, -4.0, 6.0]
        original = moving.copy()
        r, t, error = protein_transform(moving, protein)
        ligand = np.array([[8.0, 2.0, 1.0], [-8.0, -2.0, 1.0]])
        # OST's transform has residuals of a few millionths of an Angstrom
        # on this fixture; this tests geometry, not NumPy SVD precision.
        np.testing.assert_allclose(
            (ligand @ rotation + [12.0, -4.0, 6.0]) @ r + t, ligand, rtol=0, atol=1e-5
        )
        np.testing.assert_array_equal(moving, original)
        self.assertLess(error, 1e-5)
        self.assertAlmostEqual(np.linalg.det(r), 1.0, delta=1e-5)

    def test_reflection_is_not_a_rotation(self):
        protein = np.array(
            [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 4.0]]
        )
        r, t, error = protein_transform(protein * [-1, 1, 1], protein)
        self.assertAlmostEqual(np.linalg.det(r), 1.0, delta=1e-5)
        self.assertGreater(error, 0.1)

    def test_fit_is_exactly_native_ost_result(self):
        from ost import geom
        from ost.mol.alg import SuperposeSVD

        moving = np.array(
            [[1.0, 2.0, 3.0], [3.0, 2.0, 3.0], [1.0, 5.0, 3.0], [1.0, 2.0, 7.0]]
        )
        fixed = moving + np.array(
            [[0.1, 0, 0], [0, 0.2, 0], [0, 0, 0.3], [0.2, 0.1, 0]]
        )
        a, b = geom.Vec3List(), geom.Vec3List()
        for x, y in zip(moving, fixed):
            a.append(geom.Vec3(*map(float, x)))
            b.append(geom.Vec3(*map(float, y)))
        native = SuperposeSVD(a, b)
        r, t, error = protein_transform(moving, fixed)
        np.testing.assert_array_equal(
            r, [[native.transformation[j, i] for j in range(3)] for i in range(3)]
        )
        np.testing.assert_array_equal(
            t, [native.transformation[i, 3] for i in range(3)]
        )
        self.assertEqual(error, native.rmsd)
