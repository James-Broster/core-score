"""Small real-structure regression independent of the external seven-model archive.

Run explicitly with OST on PATH, or CORE_SCORE_TEST_OST set to its executable.
Missing OST or changed frozen inputs are errors, not skipped tests.

OST calculates in single precision, so different builds of the same version give answers
a few float32 steps apart. Frozen values are compared to a relative 1e-5.
"""

import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import pytest
from core_score.pipeline import execute
from core_score.inputs.inputs import load_inputs

FIXTURES = Path(__file__).parents[1] / "fixtures"
TOLERANCE = 1e-5


def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def inputs():
    manifest = json.loads((FIXTURES / "manifest.json").read_text())
    for name, checksum in manifest["structures"].items():
        assert file_hash(FIXTURES / name) == checksum
    return manifest


def options(monkeypatch):
    # Pin the OST version the values were frozen with; production has no preflight check.
    import subprocess

    executable = shutil.which(os.environ.get("CORE_SCORE_TEST_OST", "ost"))
    assert executable, "OST is required for the real-structure regression"
    version = subprocess.check_output([executable, "--version"], text=True)
    assert inputs()["ost_version"] in version
    # The package runs `ost` from PATH.
    monkeypatch.setenv("PATH", str(Path(executable).parent) + os.pathsep + os.environ["PATH"])


def groups(folder):
    """Four poses, the two frozen structures twice each, since CORE needs at least four."""
    folder.mkdir()
    for name, source in zip("abcd", ("7wl4_pose0.cif", "7wl4_pose1.cif") * 2):
        shutil.copy(FIXTURES / source, folder / f"{name}.cif")
    return load_inputs([folder / f"{name}.cif" for name in "abcd"], ligand_chains=["B"], ligand=inputs()["smiles"])


def test_original_ost_values_and_parallel_exports(tmp_path, monkeypatch):
    # Exercise OST's own library discovery, without a package-supplied path.
    monkeypatch.delenv("OST_COMPOUNDS_CHEMLIB", raising=False)
    monkeypatch.delenv("OST_ROOT", raising=False)
    options(monkeypatch)
    expected = inputs()["expected_single_copy"]
    outputs = []
    poses = groups(tmp_path / "poses")
    for workers in (1, 2):
        summary = execute(
            poses,
            tmp_path / str(workers),
            workers=workers,
            skip_posebusters=True,
        )
        with (tmp_path / str(workers) / "core_scores.csv").open() as handle:
            scores = list(csv.DictReader(handle))
        with (tmp_path / str(workers) / "details" / "pairs.csv").open() as handle:
            rows = list(csv.DictReader(handle))
        assert summary["score_count"] == 1
        assert summary["completed_comparisons"] == 6
        # Each pose's three closest others are its copy (0 A) and the two other poses (RMSD each).
        assert float(scores[0]["core_score"]) == pytest.approx(1 - (2 * expected["rmsd"] / 2) / 3, rel=TOLERANCE)
        for field, value in expected.items():
            if isinstance(value, (int, float)):
                assert float(rows[0][field]) == pytest.approx(value, rel=TOLERANCE), field
            else:
                assert rows[0][field] == value, field
        outputs.append(rows)
    assert outputs[0] == outputs[1]


def test_pb_records_invalid_real_poses_and_core_still_uses_them(tmp_path, monkeypatch):
    options(monkeypatch)
    summary = execute(groups(tmp_path / "poses"), tmp_path / "pb")
    pb = json.loads((tmp_path / "pb/details/posebusters.json").read_text())
    assert len(pb) == 4
    assert all(p["pb_valid"] is False for p in pb)
    with (tmp_path / "pb/core_scores.csv").open() as handle:
        scores = list(csv.DictReader(handle))
    assert summary["comparison_count"] == 6
    assert summary["score_count"] == 1
    assert scores[0]["prediction_count"] == "4"
    assert scores[0]["core_score"] is not None
