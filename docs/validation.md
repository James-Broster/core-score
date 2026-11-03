# Validation

## Version 0.4.0, CSV inputs and core-score-csv (3 November 2026)

The same 12 ensembles were laid out exactly as each model writes its output, with the
real structure and confidence files. AlphaFold3 used both its current file names and
the 3.0.0 names, with the copy of the top model beside the seed folders.

- `core-score-csv` found all 25 poses of every ensemble and left out each AlphaFold3
  top-level copy.
- Its confidence equalled the Runs N' Poses ligand ipTM for all 250 AlphaFold3, Boltz-1,
  Boltz-1x, Chai and Protenix poses, with a largest difference of 0. Runs N' Poses has no
  OpenFold3 values. Its OpenFold3 values were checked against the chain pairs by hand.
- `core-score` on those CSVs gave the same CORE as the run below for all 12 ensembles,
  so ordering poses by path changed nothing.
- AlphaFold3 writes a separate entity for each chain listed in its input, so
  OpenStructure 2.12 never mapped one 7LY6 chain onto the other. A pose compared with an
  identical copy of itself as reference had RMSD 0 but lDDT-PLI 0.59. With identical
  chains given one entity, that lDDT-PLI is 1.0 and AlphaFold3 7LY6 scores 0.910
  against 0.911 in the paper, instead of 0.898. The other models already write one
  entity for identical chains, so their values do not change.
- The final run used `core-score-csv` with the PLINDER system structures as references,
  then `core-score` with every default, one ensemble at a time. Against the paper's
  valid-poses table, CORE agreed within 0.002 for all 12 ensembles, the fraction of
  correct poses agreed exactly for all 12, and whether the top-ranked pose is correct
  agreed for all 12. PoseBusters validity differed for the same two poses as before. One
  of them was the paper's AlphaFold3 7LY6 top pose, so the top pose moved to the next
  valid one, confidence 0.925 against 0.93, and both are correct.
- PLINDER chain names such as 1.A are too long for PDB, which stopped every reference
  from preparing. The PoseBusters receptor file now shortens them.
- 81 unit tests and both real-structure tests pass. The real-structure pair RMSD is
  0.2812022 here against 0.2812016 frozen earlier with the same OST version. Raw OST on
  the untouched files gives 0.2812022 too, and cleans no residues with either chemical
  library, so the gap is OST's single-precision arithmetic on a different build, 19
  float32 steps. Frozen values are now compared to a relative 1e-5.

## Real co-folding outputs (3 November 2026)

The full pipeline was run with OpenStructure 2.12.0, PoseBusters 0.6.5, pandas 3.1,
RDKit 2026.09 and Gemmi 0.7.5 on the raw Runs N' Poses predictions for 7LY6 (FMN) and
7UY4 from AlphaFold3, Boltz-1, Boltz-1x, Chai, Protenix and OpenFold3, 25 poses each,
with no ligand chains given.

- All 12 ensembles were scored, with all 300 pairs compared in each.
- CORE matched the paper's values to within 0.002 for 11 of 12 ensembles. AlphaFold3
  7LY6 gave 0.898 against 0.911. 7LY6 is a homodimer, and OpenStructure 2.12 matches
  poses across its two identical binding sites differently from the 2.8.0 used for
  the paper.
- PoseBusters validity agreed with the paper for 298 of 300 poses. The two that differ
  fail the minimum distance to organic cofactors check in PoseBusters 0.6.5.
- Protenix poses, which OpenStructure cannot read ligands from, and one AlphaFold3 pose
  whose ligand OpenStructure connected differently, are compared using the package's
  ligand SDFs, with bonds from the coordinates.

## Version 0.3.0 (12 October 2026)

Version 0.3.0 streamlines the code and computes CORE as defined in the paper. It has
**not yet been run against real OpenStructure**. Run the real-structure tests below
before relying on it.

What has been checked.

- **69 unit tests** pass with PoseBusters 0.4.4, RDKit and Gemmi.
- **The built wheel** installs into a clean environment, and `core-score` runs the whole
  pipeline from that install.
- **Equivalence with 0.2.0.** Both versions were run on the same ensembles with a
  stand-in `ost` executable that computes ligand RMSD from coordinates. The ensembles
  covered single-copy and two-copy predictions, an OST-undetectable ligand, original
  SDF inputs, PDB inputs, radius recovery and a pose that fails preparation. All 63
  prepared files were byte-identical. Every pair measurement, site assignment,
  PoseBusters result and reference label matched, as did the old-formula CORE values
  before the formula changed.
- **New CORE values** were recomputed independently from `pairs.csv` with the
  paper's equations, and satisfy S = f (1 − ν/T).
- **The result files** were checked the same way. Every CORE value, pair measurement,
  PoseBusters result and reference label in them matches the run before they were
  redesigned.

The stand-in checks plumbing, not OST's scoring. The real-structure tests check that.

## Reproduce

From an installed development environment with OpenStructure 2.12 on PATH:

```bash
python -m pytest
python -m pytest -o addopts='' tests/integration/test_real_pipeline.py
```

`CORE_SCORE_TEST_OST` can point the real-structure test at an `ost` executable that
is not on PATH. The fixture manifest pins the source-file hashes, the OST version and
the expected pair measurements. The test scores four poses, the two fixture poses twice
each, so CORE is expected to equal 1 − RMSD/3 at T = 2 Å.

## Version 0.2.0 (5 October 2026)

Validated locally with OpenStructure 2.12.0 and PoseBusters 0.4.4. Two real 7WL4
predictions reproduced the original cofold-ost single-copy measurements exactly
(ligand RMSD 0.2812016010284424 Å, lDDT-PLI 1.0, pocket backbone RMSD
0.07564989477396011 Å, pocket lDDT 1.0), with identical serial and parallel results.
A two-copy ensemble produced separate site matrices, and PDB inputs and the original
SDF route completed real OST runs. A reference RMSD of 0.610021710395813 Å failed at
0.5 Å and succeeded at 1, 1.5 and 2 Å, so labels use strict `<`.
