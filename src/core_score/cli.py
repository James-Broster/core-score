"""One command: load inputs, run the analysis, then optionally plot results."""

import argparse
import math
import sys
from pathlib import Path
from core_score.inputs import load_workflow_inputs
from core_score.pipeline import execute


def positive(value):
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise argparse.ArgumentTypeError("must be finite and positive")
    return number


def worker_count(value):
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError("must be at least 1")
    return number


def parser():
    p = argparse.ArgumentParser(
        prog="core-score",
        description="Prepare structures, filter with PoseBusters, compare poses and calculate CORE.",
        epilog="Example: core-score inputs.json --out results/ --workers 4 --plots",
    )
    p.add_argument("structures", metavar="INPUT", nargs="*", type=Path,
                   help="JSON manifest, or structure paths for one ensemble")
    p.add_argument("--manifest", type=Path, help=argparse.SUPPRESS)
    p.add_argument("--out", type=Path, required=True, help="Results directory")
    p.add_argument("--smiles", help="Ligand SMILES for direct structure inputs; otherwise set in JSON")
    p.add_argument("--workers", type=worker_count, default=1,
                   help="Simultaneous OST comparisons (default: 1)")
    p.add_argument("--timeout", type=positive, default=None,
                   help="Seconds per OST call (default: no timeout)")
    p.add_argument("--no-radius-fallback", action="store_true",
                   help="Disable eligible binding-site retries at 6 and 8 angstroms")
    pb = p.add_mutually_exclusive_group()
    pb.add_argument("--skip-posebusters", action="store_true",
                    help="Skip both PB checking and filtering")
    pb.add_argument("--pb-results", type=Path,
                    help="Advanced: use existing dock results instead of running PB")
    p.add_argument("--evaluate-reference", action="store_true",
                   help="Also compare predictions against supplied ground truth")
    p.add_argument("--reference", type=Path, help="Optional ground-truth structure; may instead be set in JSON")
    p.add_argument("--rmsd-thresholds", type=positive, nargs="+", default=[2.0],
                   help="Reference success thresholds in angstroms; uses RMSD < threshold (default: 2)")
    p.add_argument("--plots", action="store_true", help="Generate plots after the analysis")
    # Optional overrides for unusual inputs/installations, not required workflow stages.
    p.add_argument("--ligand-chain", dest="ligand_chains", action="append", help=argparse.SUPPRESS)
    p.add_argument("--ligand-residue", help=argparse.SUPPRESS)
    p.add_argument("--ligand-sdf", dest="ligand_sdfs", action="append", help=argparse.SUPPRESS)
    p.add_argument("--reference-ligand-chain", help=argparse.SUPPRESS)
    p.add_argument("--reference-ligand-residue", help=argparse.SUPPRESS)
    p.add_argument("--reference-ligand-sdf", help=argparse.SUPPRESS)
    p.add_argument("--ost-bin", default="ost", help=argparse.SUPPRESS)
    p.add_argument("--chemlib", help=argparse.SUPPRESS)
    return p


def main(argv=None):
    words = sys.argv[1:] if argv is None else argv
    args = parser().parse_args(words or ["--help"])
    try:
        groups = load_workflow_inputs(args)
        summary = execute(
            groups,
            args.out,
            workers=args.workers,
            timeout=args.timeout,
            ost=args.ost_bin,
            chemlib=args.chemlib,
            radius_fallback=not args.no_radius_fallback,
            skip_posebusters=args.skip_posebusters,
            pb_results=args.pb_results,
            evaluate_reference=args.evaluate_reference,
            thresholds=args.rmsd_thresholds,
            plots=args.plots,
        )
        print(
            f"Completed {summary['completed_comparisons']}/{summary['comparison_count']} "
            f"pair comparisons; calculated {summary['score_count']} CORE scores."
        )
        print(f"Results saved to {args.out}")
        return int(any(summary["status"].values()))
    except (ValueError, OSError, RuntimeError) as exc:
        print(f"core-score: {exc}", file=sys.stderr)
        return 2
