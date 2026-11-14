![CORE Score](https://raw.githubusercontent.com/James-Broster/core-score/main/assets/logo.png)

CORE Score is a way of measuring how certain a co-folding model is about its prediction,
by measuring how closely ligand predictions agree with each other. Higher agreement has
been shown to indicate that predictions are more likely to be correct (paper link
coming). Additionally, we also calculate the PoseBusters checks for each prediction,
and an experimental structure can be provided and used to calculate each prediction's
ligand RMSD and lDDT-PLI against it.

The full description of how CORE Score is calculated can be found in the paper (paper
link coming).

## Install

```bash
conda install -c conda-forge -c bioconda openstructure
pip install core-score
```

OpenStructure is not on PyPI, so it comes from conda.

## Usage

Write a CSV naming each ensemble's ligand and predictions, then run CORE Score on it. For
example:

```text
ensemble,ligand,predictions
target1,CC(=O)Nc1ccc(O)cc1,predictions/target1/**/*.cif
```

```bash
core-score core_score_input.csv --out results/
```

**Required columns**

| Column | Meaning |
| --- | --- |
| `ensemble` | Name of the ensemble. Rows with the same name form one ensemble. |
| `ligand` | The ligand, as a SMILES or an SDF file. |
| `predictions` | A predicted structure, or a wildcard matching several. Every ensemble needs at least four. |

**Optional columns**

| Column | Meaning |
| --- | --- |
| `confidence` | The model's confidence for the prediction, one prediction per row. |
| `reference` | The experimental structure. |
| `ligand_chain` | The ligand's chain, or chains separated by spaces. |

Paths are relative to the CSV. For AlphaFold3, match `predictions/target1/seed-*/*.cif` so the
copy of its top model beside the seed folders is left out.

For AlphaFold3, Boltz, Chai-1, Protenix and OpenFold3 output, `core-score-csv` can write
this CSV for you, with one row per prediction and each prediction's confidence read from
the model's own files.

```bash
core-score-csv --predictions predictions/ --ligands ligands.csv --output core_score_input.csv
```

Here `predictions/` holds one folder per ensemble, named for it, and `ligands.csv` names
each ensemble's ligand.

```text
ensemble,ligand
target1,CC(=O)Nc1ccc(O)cc1
target2,ligands/target2.sdf
```

Add `--references references/`, a folder holding `<ensemble>.cif`, to fill in the
references too.

## Output

`core_scores.csv` has one row per ensemble, with its CORE Score, how many predictions it
has and how many pass PoseBusters.

| ensemble | core_score | prediction_count | valid_prediction_count |
| --- | --- | --- | --- |
| 7uy4 | 0.885 | 25 | 23 |

`predictions.csv` has one row per prediction, and its columns depend on what you give.

| You give | `predictions.csv` |
| --- | --- |
| Ligands and predictions | `ensemble`, `prediction`, `core_score`, `pb_valid` for each prediction |
| \+ confidences | \+ `confidence` for each prediction |
| \+ references | \+ `rmsd`, `lddt_pli` for each prediction |

With ligands and predictions only:

| ensemble | prediction | core_score | pb_valid |
| --- | --- | --- | --- |
| 7uy4 | predictions/7uy4/seed_2519/pred.model_idx_0.cif | 0.885 | True |
| 7uy4 | predictions/7uy4/seed_2519/pred.model_idx_1.cif | 0.885 | True |

With confidences:

| ensemble | prediction | core_score | confidence | pb_valid |
| --- | --- | --- | --- | --- |
| 7uy4 | predictions/7uy4/seed_2519/pred.model_idx_0.cif | 0.885 | 0.853 | True |
| 7uy4 | predictions/7uy4/seed_2519/pred.model_idx_1.cif | 0.885 | 0.852 | True |

With confidences and references:

| ensemble | prediction | core_score | confidence | pb_valid | rmsd | lddt_pli |
| --- | --- | --- | --- | --- | --- | --- |
| 7uy4 | predictions/7uy4/seed_2519/pred.model_idx_0.cif | 0.885 | 0.853 | True | 2.11 | 0.76 |
| 7uy4 | predictions/7uy4/seed_2519/pred.model_idx_1.cif | 0.885 | 0.852 | True | 2.11 | 0.77 |

These come from Chai-1's 25 predictions for 7UY4, and the
[examples](https://github.com/James-Broster/core-score/tree/main/examples) folder has
every input and output file for each case in full.

`pb_valid` is True when a prediction passes every PoseBusters `dock` check. A ligand present
in several copies is scored once per binding site, as `<ensemble>_site_1`,
`<ensemble>_site_2` and so on. A `note` column appears when something could not be
calculated, saying why.

`results/details/` holds the evidence behind both files.

| File | Contents |
| --- | --- |
| `pairs.csv` | The RMSD for every pair of predictions, the numbers CORE is built from. |
| `posebusters.json` | Every individual PoseBusters check for every prediction, so you can see why one failed. |
| `reference_details.json` | OpenStructure's full output against the reference, when you gave one. |
| `run.json` | The settings and package versions used. |
| `ligand_preparation.json` | Verified ligand chemistry sources and coordinate serialization checks. |

Two more appear only when needed: `sites.json`, saying how the copies of a ligand were
matched to binding sites, and `failed_pairs.json`, with OpenStructure's error for any pair
of predictions it could not compare.
If an OpenStructure worker cannot be reused or an attempt times out, `ost_execution/`
records the reason and any identical-command retries.

## Options

| Option | Meaning |
| --- | --- |
| `--workers N` | Run N OpenStructure comparisons at once (default 1). |
| `--plot` | Also plot the scores. Needs `pip install "core-score[plot]"`. |
| `--skip-posebusters` | Do not run PoseBusters. |
| `--rmsd-thresholds 1 2` | Calculate CORE at each threshold T in Å (default 2). |
| `--timeout SECONDS` | Time limit per OpenStructure attempt (default 120); timeouts retry twice. |
| `--no-radius-fallback` | Do not retry failed binding-site comparisons at 6 and 8 Å. |

Predicted structures can also be given directly, as one ensemble.

```bash
core-score predictions/*.cif --ligand "CCO" --out results/
```

## Changes in 0.6.0

OpenStructure now receives the same verified, SMILES-matched ligand connectivity as
PoseBusters, for predictions and experimental references. Atom order and coordinates
are checked; SDF rounding is limited to 0.0000500001 Å per coordinate. Multiple ligand
copies retain their own verified connectivity when split into sites. If chemistry cannot
be verified, the comparison remains unresolved with a diagnostic rather than silently
using a different inferred graph.

The CORE formula is unchanged: each pose's three closest neighbours, capped at T=2 Å
by default, with every pose included regardless of PoseBusters validity. **Scores may
change from 0.5.1 where its coordinate-inferred bonds differed from the verified ligand
graph.** No minimisation or coordinate repair is performed.

Ligand comparisons reuse an OpenStructure interpreter per worker, avoiding repeated
startup. Each comparison still loads its input structures and runs the installed OST
command's scoring code; structures and scores are not cached. The package discovers
OST through `PATH`, including when OST lives in a separate conda environment.
PoseBusters energy calculations use one thread to avoid CPU oversubscription, with
their thresholds and conformer seed unchanged. Unknown checks remain blank, carry a
diagnostic note, and never count as valid.

The default timeout is 120 seconds per attempt. A timeout kills the worker and its
child processes before retrying the identical comparison, at most twice. Other scoring
errors are not retried. Increase `--timeout` for larger structures; Python callers can
pass `timeout=None` to disable the comparison deadline.

Persistent workers are used automatically on supported POSIX OST installations. If
startup is unsupported, the package records the reason and uses fresh CLI processes.
For diagnostics, `CORE_SCORE_OST_MODE=subprocess` forces fresh processes, while
`CORE_SCORE_OST_MODE=persistent` requires reuse and reports an error if unavailable.
Worker startup has a separate limit of up to 30 seconds. Both modes use the same
scoring arguments and chemical library.

## Citing

- Broster JH, et al. Structural Agreement Indicates Prediction Accuracy For
  Protein-Ligand Co-Folding. Paper link coming.

CORE Score runs OpenStructure and PoseBusters, so please cite them too.

- Biasini M, et al. OpenStructure: an integrated software framework for computational
  structural biology. *Acta Crystallographica D* (2013).
- Robin X, et al. Assessment of protein–ligand complexes in CASP15. *Proteins* (2023).
- Buttenschoen M, Morris GM, Deane CM. PoseBusters: AI-based docking methods fail to
  generate physically valid poses or generalise to novel sequences. *Chemical Science*
  (2024).

## Licence

BSD 3-Clause.
