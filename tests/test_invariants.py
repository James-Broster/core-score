"""Small, dependency-light tests of the scientific output contract."""

import copy
import csv
import io
import math
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import gemmi
from core_score.stages.utils.comparison_utils import summarize, run_ost
from core_score.stages.preparation import prepare_pose


class SummaryTests(unittest.TestCase):
    def summary(self, payload, copies=1):
        return summarize(
            dict(method="test", target="target", ligand_key="B", n_input_copies=copies),
            dict(seed=1, sample=0),
            dict(seed=1, sample=1),
            dict(copies=copies),
            dict(copies=copies),
            payload,
        )

    def test_assigned_values_are_not_replaced_by_a_better_unassigned_match(self):
        payload = dict(
            status="SUCCESS",
            rmsd=dict(
                assigned_scores=[
                    dict(score=1.23456789012345, bb_rmsd=0.2, lddt_lp=0.7)
                ],
                full_results=[dict(score=0.01)],
            ),
            lddt_pli=dict(assigned_scores=[dict(score=0.81234567890123)]),
        )
        original = copy.deepcopy(payload)
        row = self.summary(payload)
        self.assertEqual(row["rmsd"], 1.23456789012345)
        self.assertEqual(row["lddt_pli"], 0.81234567890123)
        self.assertEqual(payload, original)
        stream = io.StringIO()
        writer = csv.DictWriter(stream, fieldnames=list(row))
        writer.writeheader()
        writer.writerow(row)
        stream.seek(0)
        decoded = next(csv.DictReader(stream))
        for metric in ("rmsd", "lddt_pli", "bb_rmsd_mean", "lddt_lp_mean"):
            self.assertEqual(float(decoded[metric]), row[metric])

    def test_missing_assignment_is_not_zero(self):
        row = self.summary(dict(status="SUCCESS", rmsd=dict(assigned_scores=[])))
        self.assertEqual(row["n_assigned"], 0)
        self.assertTrue(math.isnan(row["rmsd"]))
        self.assertTrue(math.isnan(row["rmsd_mean"]))

    def test_repeated_copies_keep_all_assigned_scores(self):
        row = self.summary(
            dict(
                status="SUCCESS",
                rmsd=dict(assigned_scores=[dict(score=1.0), dict(score=3.0)]),
            ),
            copies=2,
        )
        self.assertEqual(row["n_assigned"], 2)
        self.assertTrue(math.isnan(row["rmsd"]))
        self.assertEqual(
            (row["rmsd_min"], row["rmsd_mean"], row["rmsd_max"]), (1.0, 2.0, 3.0)
        )

    def test_failure_cannot_leak_partial_values_as_success(self):
        row = self.summary(
            dict(
                status="FAILURE",
                exception="mapping failed",
                rmsd=dict(assigned_scores=[dict(score=0.1)]),
            )
        )
        self.assertTrue(math.isnan(row["rmsd"]))
        self.assertEqual(row["error"], "mapping failed")

    def test_timeout_is_explicit(self):
        import subprocess

        with patch(
            "core_score.stages.utils.comparison_utils.subprocess.run",
            side_effect=subprocess.TimeoutExpired("ost", 1),
        ):
            result = run_ost(dict(path="model.cif"), dict(path="reference.cif"))
        self.assertEqual(result["status"], "FAILURE")
        self.assertIn("TimeoutExpired", result["exception"])


class PreparationTests(unittest.TestCase):
    CIF = """data_example
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.auth_asym_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
ATOM 1 C CA A 1.23456789 -2.00000001 3.0 0.50
HETATM 2 C "C 1"
B 4.00000001 5.0 6.0 1.00
HETATM 3 N N1 C 7.0 8.0 9.0 1.00
HETATM 4 O O1 X 10.0 11.0 12.0 1.00
"""

    def test_wrapped_quoted_records_preserve_every_retained_field_and_copy(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "source.cif"
            source.write_text(self.CIF)
            original = source.read_bytes()
            destination = Path(tmp) / "prepared.cif"
            result = prepare_pose(
                dict(path=str(source), ligand_chains=["B", "C"]), destination
            )
            before = (
                gemmi.cif.read_file(str(source.resolve()))
                .sole_block()
                .find_mmcif_category("_atom_site.")
            )
            after = (
                gemmi.cif.read_file(str(destination))
                .sole_block()
                .find_mmcif_category("_atom_site.")
            )
            self.assertEqual(
                [list(r) for r in after], [list(r) for r in before if r[4] != "X"]
            )
            self.assertEqual(source.read_bytes(), original)
            self.assertEqual(result["copies"], 2)

    def test_missing_chain_is_not_silently_dropped(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "source.cif"
            source.write_text(self.CIF)
            with self.assertRaisesRegex(ValueError, "absent"):
                prepare_pose(
                    dict(path=str(source), ligand_chains=["Z"]), Path(tmp) / "out.cif"
                )

    def test_explicit_sdfs_are_passed_through_without_rewriting(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "source.cif"
            source.write_text(self.CIF)
            sdf = Path(tmp) / "ligand.sdf"
            sdf.write_bytes(b"original SDF bytes\n")
            result = prepare_pose(
                dict(path=str(source), ligand_sdfs=[str(sdf.resolve())]),
                Path(tmp) / "unused.cif",
            )
            self.assertEqual(result["path"], str(source.resolve()))
            self.assertEqual(result["ligand_sdfs"], [str(sdf.resolve())])
            self.assertEqual(sdf.read_bytes(), b"original SDF bytes\n")


if __name__ == "__main__":
    unittest.main()
