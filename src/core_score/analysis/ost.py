"""Run OpenStructure's command-line actions. `ost` is found on PATH, and its chemical
library through OST's own discovery or OST_COMPOUNDS_CHEMLIB."""

import json
import shutil
import subprocess
import tempfile
from pathlib import Path
import numpy as np
from core_score.analysis import transport


INSTALL = "conda install -c conda-forge -c bioconda openstructure"


def require_ost():
    """Stop before any work when OpenStructure is not installed, saying how to install it."""
    if shutil.which("ost") is None:
        raise RuntimeError(f"OpenStructure's ost command was not found. Install it with: {INSTALL}")


class OSTError(RuntimeError):
    def __init__(self, message, evidence):
        super().__init__(message)
        self.evidence = evidence


def execute(arguments, timeout=None):
    """Run `ost <arguments> -o <file>` and return its JSON output with the execution record."""
    retries = []
    for attempt in range(transport.MAX_TIMEOUT_RETRIES + 1):
        # Fresh output path on every attempt: never read partial/stale results after timeout.
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "ost.json"
            try:
                execution = transport.run(arguments, output, timeout)
            except subprocess.TimeoutExpired as exc:
                event = dict(event="timeout", arguments=list(map(str, arguments)), timeout=timeout,
                             attempt=attempt + 1, retry=attempt < transport.MAX_TIMEOUT_RETRIES,
                             error=str(exc))
                retries.append(event)
                transport.audit(event)
                if attempt == transport.MAX_TIMEOUT_RETRIES:
                    raise OSTError(f"TimeoutExpired: {exc}", dict(timeout_retries=retries)) from exc
                continue
            execution["timeout_retries"] = retries
            if not output.exists():
                raise OSTError(execution["stderr"].strip() or
                               f"OST exited with code {execution['returncode']}", execution)
            try:
                result = json.loads(output.read_text())
            except (ValueError, OSError) as exc:
                raise OSTError(f"Invalid OST output: {exc}", execution) from exc
            return result, execution


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
