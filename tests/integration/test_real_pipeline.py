"""Small real-structure regression independent of the external seven-model archive.

Run explicitly with OST on PATH, or CORE_SCORE_TEST_OST set to its executable.
Missing OST or changed frozen inputs are errors, not skipped tests.
"""

import csv
import json
import os
from pathlib import Path
import pytest
from core_score.pipeline import execute
from core_score.inputs import load_inputs
from core_score.results import file_hash

FIXTURES = Path(__file__).parents[1] / "fixtures"


def inputs():
    manifest = json.loads((FIXTURES / "manifest.json").read_text())
    for name, checksum in manifest["structures"].items():
        assert file_hash(FIXTURES / name) == checksum
    return manifest


def options():
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
    env = options()
    expected = inputs()["expected_single_copy"]
    outputs = []
    for workers in (1, 2):
        summary = execute(
            groups(),
            tmp_path / str(workers),
            workers=workers,
            ost=env["ost"],
            skip_posebusters=True,
        )
        scores = json.loads((tmp_path / str(workers) / "scores.json").read_text())
        with (tmp_path / str(workers) / "pairs.csv").open() as handle:
            rows = list(csv.DictReader(handle))
        assert summary["score_count"] == 1
        assert summary["completed_comparisons"] == 1
        assert scores[0]["status"] == "complete"
        for field, value in expected.items():
            assert (float(rows[0][field]) if isinstance(value, (int, float)) else rows[0][field]) == value
        outputs.append(rows)
    assert outputs[0] == outputs[1]


def test_pb_runs_and_excludes_invalid_real_poses(tmp_path):
    env = options()
    summary = execute(groups(), tmp_path / "pb", ost=env["ost"])
    pb = json.loads((tmp_path / "pb/posebusters.json").read_text())
    assert len(pb) == 2
    assert all(p["copies"][0]["config"] == "dock" for p in pb)
    assert all(p["pb_valid"] is False for p in pb)
    scores = json.loads((tmp_path / "pb/scores.json").read_text())
    assert summary["comparison_count"] == 0
    assert summary["score_count"] == 0
    assert scores[0]["status"] == "incomplete"
    assert scores[0]["retained_pose_count"] == 0
    assert scores[0]["original_pose_count"] == 2
