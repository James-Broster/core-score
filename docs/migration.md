# What moved and what stayed outside the package

The source copies under `~/Downloads/pairwise-code/` remain the historical reference.
The package imports none of those scripts.

| Original source | Responsibility in this package |
| --- | --- |
| `cofold-ost/prepare.py`, `production/cif_record_fix.py`, `filter_cif_records.py` | `structure/cif.py`, `stages/preparation.py`. Gemmi parsing with retained-record verification. |
| `original-scoring/ligand_identity.py` | `structure/chemistry.py`, `structure/sdf.py`. Ligand selection and explicit connectivity SDFs. |
| `production/xyz_fix.py` | Direct atom and conformer construction, with no RDKit monkeypatch. |
| `pb_compute/worker.py`, `pb_every_copy/copy_worker.py` | `structure/sdf.py`, `stages/posebusters.py`. Template graph bijection, native AF3 graph, peer atom-name graph reuse, identity and coordinate checks, every-copy PB. |
| `cofold-ost/compare.py`, `radius.py` | `ost.py`, `stages/comparison.py`. The same OST flags, assignment extraction and eligible recovery. |
| `original-scoring/site_correspondence.py`, `multicopy_site_pairwise.py` | `stages/site_assignment.py`. OST rigid mapping and k-means/Hungarian correspondence, then one matrix per site. |
| `original-scoring/core_score.py` | `stages/scoring.py`. Replaced in 0.3.0 by the paper's nearest-neighbour CORE. |
| `production/resume.py`, `pairs.py` | Removed. Every run is fresh. |

Cluster paths, submission logic, node assertions, hardcoded repair approvals and
production logs stay in the historical copy, as do the pilot contact and pocket
investigations and the reference-choice sensitivity experiments.

## Deliberate distinctions

- OST connectivity SDFs and chemically validated PB SDFs are separate files. PB graph
  transfer uses an exact heavy-atom graph mapping and keeps the template's charges,
  bond directions, radical state and hydrogen counts. It never moves atoms or invents
  missing bonds.
- Native AF3 SMILES and atom-name correspondence and same-ensemble graph reuse are
  explicit preparation routes, and their provenance is recorded.
- Validity uses **all dock-config binary checks**. The historical 13-check summary is
  no longer exported.
- CIF receptors keep non-selected cofactors and waters for the full dock checks. The
  historical per-copy worker had filtered some of that environment, so PB labels can
  differ even when ligand coordinates agree.
- PoseBusters records validity but removes nothing. CORE uses the whole ensemble, as in
  the paper (0.2.0 excluded invalid poses before scoring).
- Default timeout is unlimited. Radius recovery applies to single-copy comparisons
  after site separation, and lDDT-PLI always comes from the 4 Å call.
- Site correspondence keeps the historical production algorithm. It does not prove
  biological site identity, so ambiguous or symmetric systems still need the saved
  mapping and costs inspected.
- CORE never uses experimental-reference measurements. Reference summaries are
  evaluation outputs, not a calibration.
