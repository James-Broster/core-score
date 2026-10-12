"""Small real-structure regression independent of the external seven-model archive.

Run explicitly with OST on PATH, or CORE_SCORE_TEST_OST set to its executable.
Missing OST or changed frozen inputs are errors, not skipped tests.
"""

import csv
import hashlib
import json
import os
from pathlib import Path
import pytest
from core_score.pipeline import execute
from core_score.inputs import load_inputs

FIXTURES = Path(__file__).parents[1] / "fixtures"


def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def inputs():
    manifest = json.loads((FIXTURES / "manifest.json").read_text())
    for name, checksum in manifest["structures"].items():
        assert file_hash(FIXTURES / name) == checksum
    return manifest


def options(monkeypatch):
    # Pin the integration fixture's toolchain; production has no preflight check.
    import shutil
    import subprocess

    executable = shutil.which(os.environ.get("CORE_SCORE_TEST_OST", "ost"))
    assert executable, "OST is required for the real-structure regression"
    library = Path(executable).resolve().parent.parent / "share/openstructure/compounds.chemlib"
    env = dict(
        ost=executable,
        ost_version=subprocess.check_output([executable, "--version"], text=True),
        chemlib_sha256=file_hash(library),
    )
    manifest = inputs()
    assert manifest["ost_version"] in env["ost_version"]
    assert env["chemlib_sha256"] == manifest["chem_library_sha256"]
    # The package runs `ost` from PATH.
    monkeypatch.setenv("PATH", str(Path(executable).parent) + os.pathsep + os.environ["PATH"])
    return env


def groups():
    return load_inputs(
        [FIXTURES / "7wl4_pose0.cif", FIXTURES / "7wl4_pose1.cif"],
        ligand_chains=["B"],
        smiles=inputs()["smiles"],
    )


def test_original_ost_values_and_parallel_exports(tmp_path, monkeypatch):
    # Exercise OST's own library discovery, without a package-supplied path.
    monkeypatch.delenv("OST_COMPOUNDS_CHEMLIB", raising=False)
    monkeypatch.delenv("OST_ROOT", raising=False)
    options(monkeypatch)
    expected = inputs()["expected_single_copy"]
    outputs = []
    for workers in (1, 2):
        summary = execute(
            groups(),
            tmp_path / str(workers),
            workers=workers,
            skip_posebusters=True,
        )
        with (tmp_path / str(workers) / "core_scores.csv").open() as handle:
            scores = list(csv.DictReader(handle))
        with (tmp_path / str(workers) / "details" / "pairs.csv").open() as handle:
            rows = list(csv.DictReader(handle))
        assert summary["score_count"] == 1
        assert summary["completed_comparisons"] == 1
        # Two poses: each one's closest other pose is the other, so CORE = 1 - RMSD / T.
        assert float(scores[0]["core_score"]) == pytest.approx(1 - expected["rmsd"] / 2)
        for field, value in expected.items():
            assert (float(rows[0][field]) if isinstance(value, (int, float)) else rows[0][field]) == value
        outputs.append(rows)
    assert outputs[0] == outputs[1]


def test_pb_records_invalid_real_poses_and_core_still_uses_them(tmp_path, monkeypatch):
    options(monkeypatch)
    summary = execute(groups(), tmp_path / "pb")
    pb = json.loads((tmp_path / "pb/details/posebusters.json").read_text())
    assert len(pb) == 2
    assert all(p["pb_valid"] is False for p in pb)
    with (tmp_path / "pb/core_scores.csv").open() as handle:
        scores = list(csv.DictReader(handle))
    assert summary["comparison_count"] == 1
    assert summary["score_count"] == 1
    assert scores[0]["pose_count"] == 2
    assert scores[0]["core_score"] is not None
