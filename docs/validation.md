# Validation of the initial refactor

Validated locally on 5 October 2026. This is a tested initial implementation, not
a claim that every historical model/chemistry edge case has been revalidated.

## Automated checks

- The original refactor passed **51 unit tests**. Subsequent CLI simplification
  and input validation changes expanded the suite; the latest run passed **74
  non-integration tests**. Coverage includes CIF records, site assignment, radius
  recovery, chemistry/coordinate guards, complete matrices, CORE, JSON inputs,
  conditional reference validation, PB filtering and fresh comparison execution.
  Checkpoint/resume and cross-host merge functionality has been removed.
- **2 real-structure integration tests passed** against OpenStructure 2.12.0 and
  PoseBusters 0.4.4: exact original single-copy OST measurements, serial/parallel
  agreement, and automatic dock checking/exclusion on two real 7WL4 predictions.
- **3 OST Python-binding superposition tests passed** in the local OST environment.
- Built and installed a wheel into an isolated temporary environment; the CLI
  is available outside the source checkout.
- The conda environment recipe resolved successfully in a **dry run on macOS ARM**.
  This did not create a full new conda environment and does not validate every OS.

## Real-input smoke checks

- Two 7WL4 predicted structures through preparation, PB, comparison, CORE and plots.
- The single-copy ligand RMSD, lDDT-PLI, pocket backbone RMSD and pocket lDDT matched
  the **original copied cofold-ost implementation exactly** for this pair:
  ligand RMSD 0.2812016010284424 Å, lDDT-PLI 1.0, pocket backbone RMSD
  0.07564989477396011 Å, pocket lDDT 1.0.
- Both selected single-copy poses failed dock checks. The revised default excluded
  both, ran zero pair comparisons, and reported CORE unavailable with original and
  retained counts. With `--skip-posebusters`, both remained and CORE was computed.
- A two-copy ensemble produced separate site matrices and scores with two workers.
- PDB inputs and the explicit ligand-SDF/receptor route both completed real OST runs.
- Reference evaluation produced an RMSD of 0.610021710395813 Å for a selected pose:
  failure at 0.5 Å and success at 1, 1.5 and 2 Å. Labels use strict `<`.
- Generated heatmaps were inspected visually.

## Reproduce

From an installed development environment:

```bash
python -m pytest
python -m pytest -o addopts='' tests/integration/test_real_pipeline.py
python -m unittest discover -s tests/integration -p test_site_superposition.py -v
```

`CORE_SCORE_TEST_OST` can point the real-structure test at a custom executable.
The bundled fixture manifest pins source-file and chemical-library hashes and
expected measurements. A changed library is an explicit failure, not an automatic
rewrite of expected results. The superposition suite requires OST Python bindings.

## Historical coverage still outstanding

The original seven-model fixture archive/manifest was not included in the copied
bundle. Its regression suite remains in `tests/integration/test_ost_regression.py`:

```bash
COFOLD_REGRESSION_MANIFEST=/path/to/original/manifest.json \
  python -m pytest -o addopts='' tests/integration/test_ost_regression.py
```

Missing or changed fixtures are errors when that suite is requested. It has **not**
been claimed as passed here. Broader native-AF3, Protenix, metal/dative chemistry,
multicopy ambiguity and historical repair cases need their original fixtures before
scientific equivalence across the full dataset can be asserted.

The new default filters on all PoseBusters dock checks; historical paper summaries
used a named 13-check subset. Both results are exported, and this difference is
documented in [migration notes](migration.md).
