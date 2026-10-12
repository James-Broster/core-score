"""Run OpenStructure's command-line actions and return their JSON results."""

import dataclasses
import json
import os
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from core_score.results import diagnostic


def run_ost(model, reference, ost="ost", timeout=None, chemlib=None, radius=4.0):
    env = dict(
        os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1"
    )
    if chemlib:
        env["OST_COMPOUNDS_CHEMLIB"] = str(Path(chemlib).resolve())
    with tempfile.TemporaryDirectory() as directory:
        output = str(Path(directory) / "ost.json")
        command = [ost, "compare-ligand-structures", "-m", model["path"]]
        if model.get("ligand_sdfs"):
            command += ["-ml", *model["ligand_sdfs"]]
        command += ["-r", reference["path"]]
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
            output,
        ]
        if radius != 4.0:
            command += ["--radius", str(radius)]
        evidence = dict(
            command=command,
            timeout=timeout,
            radius=radius,
            chemlib=env.get("OST_COMPOUNDS_CHEMLIB"),
            model=model,
            reference=reference,
        )
        try:
            process = subprocess.run(
                command, env=env, capture_output=True, text=True, timeout=timeout
            )
            evidence.update(
                returncode=process.returncode,
                stdout=process.stdout,
                stderr=process.stderr,
            )
            evidence["raw_output"] = (
                Path(output).read_text() if Path(output).exists() else None
            )
            result = json.loads(Path(output).read_text())
            if process.returncode and result.get("status") != "FAILURE":
                result = dict(
                    result, status="FAILURE", exception=process.stderr[-2000:]
                )
            if result.get("status") not in ("SUCCESS", "FAILURE"):
                result.update(
                    status="FAILURE", exception="OST did not report a recognised status"
                )
            result["_execution"] = evidence
            return result
        except (OSError, ValueError, subprocess.TimeoutExpired) as exc:
            if isinstance(exc, subprocess.TimeoutExpired):
                evidence.update(
                    stdout=(
                        exc.stdout.decode(errors="replace")
                        if isinstance(exc.stdout, bytes)
                        else exc.stdout
                    ),
                    stderr=(
                        exc.stderr.decode(errors="replace")
                        if isinstance(exc.stderr, bytes)
                        else exc.stderr
                    ),
                )
            report = diagnostic("comparison", exc, **evidence)
            if isinstance(exc, subprocess.TimeoutExpired):
                report["suggestion"] = (
                    "Retry with a larger --timeout or omit --timeout for no limit."
                )
            return dict(
                status="FAILURE",
                exception=f"{type(exc).__name__}: {exc}",
                _execution=report,
            )


class RigidSuperpositionError(RuntimeError):
    """Raised when `ost compare-structures --rigid-scores` cannot superpose two poses."""

    def __init__(self, message, evidence=None):
        super().__init__(message)
        self.evidence = evidence or {}


@dataclasses.dataclass(frozen=True)
class RigidResult:
    transform: (
        np.ndarray
    )  # (4, 4) homogeneous transform: superposes `model` onto `reference`
    rmsd: float  # protein backbone RMSD after superposition (informational)
    chain_mapping: dict  # OST's own reference-chain -> model-chain mapping


def rigid_transform(
    ost_bin: str, model_cif: Path, reference_cif: Path, timeout: float | None = None
) -> RigidResult:
    """Superpose `model_cif` onto `reference_cif` using OST's own rigid Kabsch superposition
    of matched protein C-alpha atoms (`ost compare-structures --rigid-scores`).

    Raises RigidSuperpositionError on any failure (process error, non-SUCCESS status, or a
    missing transform) rather than returning a sentinel -- callers must handle this
    explicitly, not have it silently absorbed.
    """
    with tempfile.TemporaryDirectory() as tmp:
        out_json = f"{tmp}/rigid.json"
        cmd = [
            ost_bin,
            "compare-structures",
            "-m",
            str(model_cif),
            "-r",
            str(reference_cif),
            "--rigid-scores",
            "-o",
            out_json,
        ]
        evidence = dict(
            command=cmd,
            timeout=timeout,
            model=str(model_cif),
            reference=str(reference_cif),
        )
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        except subprocess.TimeoutExpired as exc:
            raise RigidSuperpositionError(
                f"ost compare-structures timed out after {timeout}s", evidence
            ) from exc
        evidence.update(
            returncode=proc.returncode, stdout=proc.stdout, stderr=proc.stderr
        )
        if proc.returncode != 0:
            raise RigidSuperpositionError(
                proc.stderr.strip()
                or proc.stdout.strip()
                or f"exit code {proc.returncode}",
                evidence,
            )
        data = json.loads(Path(out_json).read_text())
        if data.get("status") != "SUCCESS":
            raise RigidSuperpositionError(
                f"OST status={data.get('status')}", dict(evidence, ost_results=data)
            )
        transform = data.get("transform")
        if transform is None:
            raise RigidSuperpositionError(
                "OST returned no transform (chains likely could not be mapped)",
                dict(evidence, ost_results=data),
            )
        return RigidResult(
            transform=np.asarray(transform, dtype=float),
            rmsd=data.get("rmsd"),
            chain_mapping=data.get("chain_mapping") or {},
        )
