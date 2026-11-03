"""Run OpenStructure's command-line actions. `ost` is found on PATH, and its chemical
library through OST's own discovery or OST_COMPOUNDS_CHEMLIB."""

import json
import os
import subprocess
import tempfile
from pathlib import Path
import numpy as np


class OSTError(RuntimeError):
    def __init__(self, message, evidence):
        super().__init__(message)
        self.evidence = evidence


def execute(arguments, timeout=None):
    """Run `ost <arguments> -o <file>` and return its JSON output with the execution record."""
    # One numerical-library thread per process, so --workers alone sets the parallelism.
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    with tempfile.TemporaryDirectory() as directory:
        output = Path(directory) / "ost.json"
        command = ["ost", *map(str, arguments), "-o", str(output)]
        process = subprocess.run(command, env=env, capture_output=True, text=True, timeout=timeout)
        execution = dict(command=command, returncode=process.returncode,
                         stdout=process.stdout, stderr=process.stderr)
        if not output.exists():
            raise OSTError(process.stderr.strip() or f"OST exited with code {process.returncode}", execution)
        return json.loads(output.read_text()), execution


def compare_ligands(model, reference, timeout=None, radius=4.0, lddt_pli=False):
    """`ost compare-ligand-structures` result for one model/reference pair of prepared structures.
    lDDT-PLI is only computed when asked for, since CORE uses RMSD alone."""
    arguments = ["compare-ligand-structures", "-m", model["path"]]
    if model["ligand_sdfs"]:
        arguments += ["-ml", *model["ligand_sdfs"]]
    arguments += ["-r", reference["path"]]
    if reference["ligand_sdfs"]:
        arguments += ["-rl", *reference["ligand_sdfs"]]
    arguments += ["--fault-tolerant", "--rmsd", "--full-bs-search", "--full-results"]
    if lddt_pli:
        arguments += ["--lddt-pli", "--lddt-pli-add-mdl-contacts"]
    if radius != 4.0:
        arguments += ["--radius", radius]
    result, execution = execute(arguments, timeout)
    if execution["returncode"] or result.get("status") != "SUCCESS":
        result["status"] = "FAILURE"
        result.setdefault("exception", execution["stderr"][-2000:] or "OST did not report success")
    result["_execution"] = execution
    return result


def reads_ligand(prepared, timeout=None):
    """Whether OST reads each ligand copy from the structure itself, connected the same way
    as in the package's own SDFs (bonds from the coordinates). A ligand OST cannot find,
    cannot read, or connects differently returns False."""
    sdfs = [copy["ost_sdf"] for copy in prepared["copy_records"]]
    try:
        result = compare_ligands(prepared, dict(prepared, ligand_sdfs=sdfs), timeout)
    except OSTError:
        return False
    assigned = [] if result["status"] == "FAILURE" else result.get("rmsd", {}).get("assigned_scores", [])
    return len(assigned) == len(sdfs)


def superpose(model_path, reference_path, timeout=None):
    """OST rigid superposition of matched protein C-alpha atoms, model onto reference."""
    result, execution = execute(
        ["compare-structures", "-m", model_path, "-r", reference_path, "--rigid-scores"], timeout
    )
    if execution["returncode"] or result.get("status") != "SUCCESS" or result.get("transform") is None:
        raise OSTError(
            f"OST rigid superposition failed (status {result.get('status')}); chains may not map",
            dict(execution, ost_results=result),
        )
    return dict(transform=np.asarray(result["transform"], dtype=float),
                rmsd=result.get("rmsd"), chain_mapping=result.get("chain_mapping") or {})
