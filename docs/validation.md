# Validation

## Version 0.3.0 (12 October 2026)

Version 0.3.0 streamlines the code and computes CORE as defined in the paper. It has
**not yet been run against real OpenStructure**. Run the real-structure tests below
before relying on it.

What has been checked.

- **57 unit tests** pass with PoseBusters 0.4.4, RDKit and Gemmi.
- **Equivalence with 0.2.0.** Both versions were run on the same ensembles with a
  stand-in `ost` executable that computes ligand RMSD from coordinates. The ensembles
  covered single-copy and two-copy predictions, an OST-undetectable ligand, original
  SDF inputs, PDB inputs, radius recovery and a pose that fails preparation. All 63
  prepared files were byte-identical. Every pair measurement, site assignment,
  PoseBusters result and reference label matched, as did the old-formula CORE values
  before the formula changed.
- **New CORE values** were recomputed independently from `pairs.csv` with the
  paper's equations, and satisfy S = f (1 − ν/T).
- Generated heatmaps and the reference plot were inspected visually.

The stand-in checks plumbing, not OST's scoring. The real-structure tests check that.

## Reproduce

From an installed development environment with OpenStructure 2.12 on PATH:

```bash
python -m pytest
python -m pytest -o addopts='' tests/integration/test_real_pipeline.py
```

`CORE_SCORE_TEST_OST` can point the real-structure test at an `ost` executable that
is not on PATH. The fixture manifest pins source-file and chemical-library hashes
and the expected single-copy measurements. A changed library is an explicit failure,
not an automatic rewrite of expected results. For the two fixture poses, CORE is
expected to equal 1 − RMSD/2.

## Version 0.2.0 (5 October 2026)

Validated locally with OpenStructure 2.12.0 and PoseBusters 0.4.4. Two real 7WL4
predictions reproduced the original cofold-ost single-copy measurements exactly
(ligand RMSD 0.2812016010284424 Å, lDDT-PLI 1.0, pocket backbone RMSD
0.07564989477396011 Å, pocket lDDT 1.0), with identical serial and parallel results.
A two-copy ensemble produced separate site matrices, and PDB inputs and the original
SDF route completed real OST runs. A reference RMSD of 0.610021710395813 Å failed at
0.5 Å and succeeded at 1, 1.5 and 2 Å, so labels use strict `<`.

## Historical coverage still outstanding

The original seven-model fixture archive was not included in this package. Its
regression suite remains in `tests/integration/test_ost_regression.py`.

```bash
COFOLD_REGRESSION_MANIFEST=/path/to/original/manifest.json \
  python -m pytest -o addopts='' tests/integration/test_ost_regression.py
```

It has not been run since it was adapted to 0.3.0. Columns removed from `pairs.csv`
are skipped when comparing its frozen expected values, and renamed ones are mapped.
Native AF3, Protenix, metal and dative chemistry, multicopy ambiguity and the
historical repair cases need their original fixtures before scientific equivalence
across the full dataset can be asserted.
