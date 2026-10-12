# What moved and what stayed outside the package

The source copies under `~/Downloads/pairwise-code/` remain the historical reference.
The new runtime imports none of those scripts.

| Original source | New responsibility |
| --- | --- |
| `cofold-ost/prepare.py`, `production/cif_record_fix.py`, analysis workspace `filter_cif_records.py` | `stages/utils/preparation/cif.py`, `stages/utils/prepare_utils.py`: Gemmi parsing and retained-field verification, always enabled. |
| `original-scoring/ligand_identity.py` | `stages/utils/preparation/chemistry.py`, `stages/utils/preparation/sdf.py`: formula selection and explicit connectivity SDFs. |
| `production/xyz_fix.py` | Direct atom/conformer construction; no global RDKit monkeypatch. |
| analysis workspace `pb_compute/worker.py`, `pb_every_copy/copy_worker.py` | `stages/utils/preparation/sdf.py`, `stages/utils/posebusters_utils.py`: template graph bijection, native AF3 graph, peer atom-name graph reuse, chemical identity and coordinate checks, every-copy PB. |
| `cofold-ost/compare.py`, `radius.py` | `stages/utils/comparison_utils.py`, `stages/utils/radius.py`: same OST metric flags, assignment extraction and eligible recovery; stronger execution diagnostics. |
| `original-scoring/site_correspondence.py`, `multicopy_site_pairwise.py` | `stages/utils/site_algorithms.py`, `stages/site_assignment.py` and `stages/utils/sites_utils.py`: historical OST rigid mapping and k-means/Hungarian correspondence, followed by one matrix per site. |
| `cofold-ost/sites.py` | Lower-level explicit-correspondence and multistart utilities remain available with tests. |
| `original-scoring/core_score.py` | `stages/utils/score_utils.py`: unchanged score definition, with square-matrix validation added. |
| `production/resume.py`, `pairs.py`, audit helpers | Pair identifiers, coverage validation and export in pipeline/matrix/results; checkpoint/resume machinery removed. |

Cluster paths, submission/continuation logic, node assertions, hardcoded repair
approvals and production logs are not runtime dependencies. They remain in the
historical reference copy. Pilot contact/pocket investigations and reference-choice
sensitivity experiments remain there too; exported native OST pocket RMSD, pocket
lDDT and lDDT-PLI are included in every comparison.

## Deliberate distinctions

- OST connectivity SDFs and chemically validated PB SDFs are separate artifacts.
  PB graph transfer uses an exact heavy-atom graph mapping and retains the template's
  charges, bond directions, radical state and hydrogen counts. It does not move atoms
  or invent missing bonds. This replaces implicit chemistry normalization; regressions
  are needed for unusual chemistry before claiming complete historical equivalence.
- Native AF3 SMILES/atom-name correspondence and same-ensemble graph reuse are
  explicit preparation routes. Their provenance is recorded.
- The new default uses **all dock-config checks** for exclusion. The historical
  paper's 13-check summary is also exported but is not silently substituted for dock.
- CIF receptor environments retain non-selected cofactors/waters for the full dock
  checks. The historical per-copy production worker had filtered some of that
  environment. PB labels can therefore differ even when ligand coordinates agree.
- PB now runs and filters by default, as requested. `--skip-posebusters` skips both.
- Default timeout is unlimited. Default radius recovery is enabled for eligible
  single-copy comparisons after site separation. Raw baseline lDDT-PLI is retained.
- Default site correspondence keeps the historical production algorithm. It is not
  equivalent to proving biological site identity; ambiguous or symmetric systems
  still require inspection of the saved mapping and costs.
- CORE never consumes experimental-reference measurements. Reference threshold
  summaries are evaluation outputs, not fitted confidence-to-accuracy calibration.

## Still requiring broader validation

The copied package's original seven-model frozen fixture manifest and archived
input SDFs were not in the transferred bundle. The original regression tests are
retained as an explicit integration suite. A small real 7WL4 fixture provides a
local baseline, not a replacement claim of seven-model coverage.
