"""Real OST integration tests. Missing regression inputs are an error, never a skip."""

import csv
import hashlib
import itertools
import json
import math
import os
import subprocess
import tempfile
import unittest
from pathlib import Path

from core_score.stages.utils.preparation.prepare import prepare_pose
from core_score.stages.utils.comparison_utils import summarize
from core_score.stages.utils.ost import run_ost


def run(groups, output, workers, options):
    """Prepare and compare every pair without radius recovery, as the frozen fixtures expect."""
    from core_score.stages.comparison import compare_groups
    from core_score.stages.preparation import prepare_groups

    for group in groups:
        group.setdefault("ensemble", group["target"])
        group.setdefault("site", "1")
        for i, p in enumerate(group["poses"]):
            p.setdefault("pose_id", str(i))
            p.setdefault("ensemble", group["ensemble"])
    prepared = prepare_groups(groups, output)
    rows = compare_groups(
        prepared, output, workers, ost=options.get("ost", "ost"),
        chemlib=options.get("chemlib"), timeout=options.get("timeout"), radius_fallback=False
    )
    return len(rows)


def same_number(test, actual, expected):
    actual = float(actual)
    expected = float(expected)
    if math.isnan(expected):
        test.assertTrue(math.isnan(actual))
    else:
        test.assertEqual(actual, expected)  # Deliberately no tolerance or rounding.


class OSTRegressionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = os.environ.get("COFOLD_REGRESSION_MANIFEST")
        if not path:
            raise RuntimeError(
                "Set COFOLD_REGRESSION_MANIFEST to the frozen real-structure fixture manifest"
            )
        cls.fixture = json.loads(Path(path).read_text())
        cls.groups = cls.fixture["groups"]
        cls.options = cls.fixture["options"]
        if not cls.groups:
            raise RuntimeError("Regression fixture set is empty")
        for file, digest in cls.fixture["sha256"].items():
            if hashlib.sha256(Path(file).read_bytes()).hexdigest() != digest:
                raise RuntimeError("Regression input changed: " + file)

    def direct_ost(self, model, reference, destination):
        # Independent baseline: original input files, bypassing our preparation,
        # invocation helper, summarizer, caching and multiprocessing.
        command = [
            self.options["ost"],
            "compare-ligand-structures",
            "-m",
            model["path"],
            "-r",
            reference["path"],
        ]
        if model.get("ligand_sdfs"):
            command += ["-ml", *model["ligand_sdfs"]]
        if reference.get("ligand_sdfs"):
            command += ["-rl", *reference["ligand_sdfs"]]
        command += [
            "--fault-tolerant",
            "--rmsd",
            "--full-bs-search",
            "--full-results",
            "--lddt-pli",
            "--lddt-pli-add-mdl-contacts",
            "-o",
            str(destination),
        ]
        env = dict(
            os.environ,
            OST_COMPOUNDS_CHEMLIB=self.options["chemlib"],
            OMP_NUM_THREADS="1",
            OPENBLAS_NUM_THREADS="1",
            MKL_NUM_THREADS="1",
        )
        subprocess.run(
            command, env=env, check=True, capture_output=True, text=True, timeout=900
        )
        return json.loads(destination.read_text())

    def test_raw_ost_and_prepared_pipeline_preserve_each_ligand_pair_score(self):
        cases = (
            (group, list(poses))
            for group in self.groups
            for poses in itertools.combinations(group["poses"], 2)
        )
        for group, poses in cases:
            with self.subTest(
                method=group["method"], target=group["target"]
            ), tempfile.TemporaryDirectory() as tmp:
                raw = self.direct_ost(*poses, Path(tmp) / "raw.json")
                prepared = [
                    prepare_pose(p, Path(tmp) / f"{i}.cif") for i, p in enumerate(poses)
                ]
                wrapped = run_ost(*prepared, **self.options)
                self.assertEqual(raw["status"], "SUCCESS")
                self.assertEqual(wrapped["status"], "SUCCESS")
                for metric in ("rmsd", "lddt_pli"):
                    baseline = {
                        (r["model_ligand"], r["reference_ligand"]): r
                        for r in raw[metric]["full_results"]
                    }
                    for result in wrapped[metric]["full_results"]:
                        key = (result["model_ligand"], result["reference_ligand"])
                        self.assertIn(key, baseline)
                        for field in ("score", "bb_rmsd", "lddt_lp", "coverage"):
                            self.assertEqual(
                                result.get(field),
                                baseline[key].get(field),
                                (group["method"], metric, key, field),
                            )
                row = summarize(group, *poses, *prepared, wrapped)
                if poses == group["poses"][:2]:
                    for field, value in group["expected"].items():
                        if isinstance(value, (float, int)) and not isinstance(
                            value, bool
                        ):
                            same_number(self, row[field], value)
                        else:
                            self.assertEqual(row[field], value)

    def test_serial_parallel_and_csv_exports_are_identical(self):
        with tempfile.TemporaryDirectory() as tmp:
            outputs = []
            for workers in (1, 2):
                folder = Path(tmp) / str(workers)
                run(self.groups, folder, workers, self.options)
                with open(folder / "pairs.csv") as handle:
                    rows = list(csv.DictReader(handle))
                self.assertEqual(
                    len(rows),
                    sum(
                        len(g["poses"]) * (len(g["poses"]) - 1) // 2
                        for g in self.groups
                    ),
                )
                outputs.append(rows)
                with open(folder / "ligand_details.jsonl") as handle:
                    details = [json.loads(line) for line in handle]
                for row, detail in zip(rows, details):
                    values = [
                        r["score"]
                        for r in detail["ost_results"][0]["rmsd"]["assigned_scores"]
                        if r.get("score") is not None
                    ]
                    self.assertEqual(int(row["n_assigned"]), len(values))
                    if len(values) == 1:
                        same_number(self, row["rmsd"], values[0])
            self.assertEqual(outputs[0], outputs[1])

    def test_source_files_still_match_frozen_hashes(self):
        for path, digest in self.fixture["sha256"].items():
            self.assertEqual(
                hashlib.sha256(Path(path).read_bytes()).hexdigest(), digest, path
            )


if __name__ == "__main__":
    unittest.main()
