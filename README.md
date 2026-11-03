# CORE Score

A Python package for prediction-only ensemble agreement. It prepares protein–ligand
structures, runs PoseBusters `dock` checks, compares every pair of poses with
OpenStructure, and calculates CORE from the complete pairwise ligand-RMSD matrix.
Experimental references are optional and never enter the CORE calculation.

See [validation](https://github.com/James-Broster/core-score/blob/main/docs/validation.md) for what has been checked.

## CORE Score

For each pose *i* in an ensemble of *N* poses, δ*ᵢ,ₘ* is the ligand RMSD to its *m*th
closest other pose, for its three closest other poses. Each distance is capped at a
threshold *T* and divided by *T*, so a pose reproduced exactly by its three closest poses
has penalty 0 and a pose with no other pose within *T* has penalty 1. CORE is one minus
the mean penalty.

```text
S = 1 − (1/3N) Σᵢ Σₘ min(δᵢ,ₘ, T) / T
```

Equivalently S = f (1 − ν/T), where f is the fraction of closest-pose distances δᵢ,ₘ
within T and ν is the mean of those distances. Using three closest poses means a pose
scores highly only when at least four poses agree closely. The default T is 2 Å, the RMSD
threshold for a correct prediction. CORE uses every pose in the ensemble, including
poses that fail PoseBusters checks, because they still sample the model's output
distribution. Validity filtering belongs downstream, for example when choosing the
top-ranked prediction.

CORE needs at least four poses per ensemble, so that every pose has three others to be
scored against. An ensemble with fewer is rejected before anything runs.

Each pair of poses is compared once with OpenStructure. Poses are ordered by path, and
the earlier pose of each pair is the OST model, so the same poses always give the same
scores whatever order they are listed in.

## Install

CORE Score needs OpenStructure 2.12, which is installed with conda. Then install the
package with pip.

```bash
conda create -n core-score -c conda-forge -c bioconda python=3.12 openstructure=2.12
conda activate core-score
pip install "core-score[plot]"
core-score --help
```

`plot` adds matplotlib for `--plot`. Plain `pip install core-score` installs the
pipeline alone. The package is `core_score`. It installs two commands, `core-score` and
`core-score-csv`.

To work on the package itself, clone the repository and use its environment, which
installs it in editable mode with the tests.

```bash
conda env create -f environment.yml
conda activate core-score
```

The conda solve has been checked on macOS ARM. Other platforms need an available
OpenStructure 2.12 build. Dependencies track current releases, including PoseBusters
0.6 and pandas 3. OpenStructure is an external executable. `ost` is launched from PATH
and finds its chemical library itself or through `OST_COMPOUNDS_CHEMLIB`.

## What you give, what you get

| You give | You get |
| --- | --- |
| At least four predicted structures per ensemble, and the ligand's SMILES | CORE Score per ensemble, and PoseBusters validity per pose |
| \+ a confidence per pose (read automatically by `core-score-csv`) | \+ the top-ranked pose, the PoseBusters-valid pose with the highest confidence |
| \+ a ground-truth structure per ensemble | \+ ligand RMSD, lDDT-PLI and whether each pose is correct, and whether the top-ranked pose is |

An ensemble is a set of predictions of the same protein–ligand system to compare
together, usually every seed and sample one model made for one target. Keep different
models in separate ensembles to score them separately.

## Walkthrough

**1. Put each ensemble's predictions in a folder named for it.** Leave the files exactly
as the model wrote them. Any folders inside are searched.

```text
predictions/
  7ly6/
    seed-1_sample-0/7ly6_seed-1_sample-0_model.cif
    seed-1_sample-1/7ly6_seed-1_sample-1_model.cif
    ...
  7uy4/
    ...
```

**2. Write each ensemble's ligand SMILES in a CSV.**

```text
ensemble,smiles
7ly6,Cc1cc2nc3c(=O)[nH]c(=O)nc-3n(C[C@H](O)[C@H](O)[C@H](O)COP(=O)(O)O)c2cc1C
7uy4,CN[C@@H]1[C@H](O)[C@H](NC)[C@H]2O[C@]3(O)[C@H](O[C@@H]2[C@H]1O)O[C@H](C)CC3(O)O
```

**3. Make the input CSV.**

```bash
core-score-csv --predictions predictions/ --smiles smiles.csv --output inputs.csv
```

This finds every structure, lists a structure found twice with the same coordinates
once (AlphaFold3 writes a copy of its top model beside the seed folders), and reads each
pose's confidence from the file the model wrote beside it. It stops and lists every
problem, such as an ensemble without a SMILES or with fewer than four structures. Add
`--references references/`, a folder holding `<ensemble>.cif` or `<ensemble>.pdb`, to
evaluate the predictions against ground truth.

The confidence is the ligand's ipTM as in Runs N' Poses, the mean chain-pair ipTM between
the ligand chain and each receptor chain, also averaged over the ligand's copies when it
has several. It is read from these files.

| Model | Structure | Confidence file |
| --- | --- | --- |
| AlphaFold3 | `*model.cif` | `*summary_confidences.json` |
| Boltz | `<name>_model_<k>.cif` | `confidence_<name>_model_<k>.json` |
| Chai-1 | `pred.model_idx_<k>.cif` | `scores.model_idx_<k>.npz` |
| Protenix | `<name>_sample_<k>.cif` | `<name>_summary_confidence_sample_<k>.json` |
| OpenFold3 | `*model.cif` | `*confidences_aggregated.json` |

For any other method, the confidence column is left out and you can add your own.

**4. Run CORE Score.**

```bash
core-score inputs.csv --out results/ --workers 8
```

**5. Read `results/core_scores.csv`**, one row per ensemble, highest CORE first, and
`results/poses.csv`, one row per pose. [Outputs](#outputs) describes every column.

## The input CSV

`core-score-csv` writes this file, and you can also write it yourself. It needs three
columns. A row names one structure, or many with a wildcard.

```text
ensemble,smiles,poses
target1,CC(=O)Nc1ccc(O)cc1,predictions/target1/*.cif
```

| Column | Meaning |
| --- | --- |
| `ensemble` | The ensemble's name. Rows with the same name form one ensemble. |
| `smiles` | The ligand's SMILES, the same on every row of an ensemble. |
| `poses` | A structure file, or a wildcard matching several (`**` searches folders inside). PDB, CIF and mmCIF, one model per file. |
| `confidence` | Optional. The model's own score for the pose, higher meaning more confident. One pose per row, and given for every pose of an ensemble or for none. |
| `reference` | Optional. The ground-truth structure, the same on every row of an ensemble. |
| `ligand_chain` | Optional. The ligand's chain, or chains separated by spaces for several identical copies. |
| `ligand_sdf` | Optional. The ligand as an SDF in the structure's frame, for tools such as docking that write it separately. One pose per row. |

Paths are relative to the CSV. Pose IDs are the structure paths relative to the CSV.
See [the example](https://github.com/James-Broster/core-score/blob/main/examples/inputs.csv), which uses every column a CSV usually needs.
Coordinates are never aligned or changed.

The ligand is found in each structure by its SMILES heavy-atom composition, and
`ligand_chain` narrows that down. Missing ligands are reported, not guessed. Composition
only selects residues. Preparation separately verifies the chemistry against the SMILES.

Structure files given directly form one ensemble.

```bash
core-score poses/*.cif --smiles 'CCO' --out results/
```

The same runs from Python:

```python
import core_score

core_score.run("inputs.csv", "results/", workers=8)
```

`core_score.score(matrix, threshold=2.0)` scores a pairwise RMSD matrix you have
computed yourself. `neighbours=3` is the number of closest poses each pose is scored against.

## Options

| Option | Meaning |
| --- | --- |
| `--workers N` | At most N simultaneous OST comparisons (default 1). Each OST process uses one numerical-library thread. |
| `--timeout SECONDS` | Limit per OST call. Default is no limit. |
| `--rmsd-thresholds 0.5 1 2` | Thresholds T in Å. CORE is calculated at each from the same matrix, and a pose is correct at T when its RMSD to the reference is below T and its lDDT-PLI above 0.8. Default 2. |
| `--core-cutoffs 0.5 0.7 0.9` | CORE cutoffs to report against. A complex at or above a cutoff is above it. Default 0.5, 0.6, 0.7, 0.8 and 0.9. |
| `--no-radius-fallback` | Disable the 4 → 6 → 8 Å binding-site recovery. |
| `--skip-posebusters` | Do not run PoseBusters, which runs by default. CORE is unaffected, but there is no top-ranked pose. |
| `--plot` | Also save `core_scores.png`, a dot per complex ordered by CORE Score. Needs matplotlib. |
| `--smiles`, `--ligand-chain`, `--reference` | The ligand's SMILES and chains, and the ground truth, for structure files given directly. A CSV sets its own. |

## PoseBusters

PoseBusters runs once per selected ligand copy against its predicted environment,
using the official API with `config='dock'`. A pose is valid only if every copy passes
every dock binary check. If none fail but any are unknown, the pose is unknown,
never valid. Missing chemical information is unknown, not a pass. Each check's result
for a pose combines its copies the same way. No pose is removed.

## Binding-site recovery

When OST cannot score a pair because of a binding-site failure, the comparison is
retried at 6 Å and then 8 Å. Graph errors and general failures are not retried. The
ligand and pocket RMSD can come from the larger radius, but lDDT-PLI always stays
from the original 4 Å call. The standard RMSD, chosen radius and every attempt are
recorded.

## Multiple identical copies

When each prediction contains several copies of the ligand, every pose is superposed
onto the anchor pose by OST's rigid protein alignment. Constrained k-means then gives
each copy a site label consistent across the ensemble, and each site gets its own
matrix and CORE Score. Confidence and ground truth are not used. Different copy
counts or alignment failures block the ensemble and are reported. Site labels are
arbitrary and do not prove biological correspondence. Transforms, costs and the
anchor are saved in `sites.json`. The sites of ensemble `7ly6` are reported as
`7ly6_site_1`, `7ly6_site_2` and so on.

## Identical chains

OpenStructure treats chains of different mmCIF entities as different molecules when it
maps one structure's chains onto another's. AlphaFold3 writes an entity per chain listed
in its input, so the two chains of a homodimer could never swap places and pockets
between them would be scored wrongly. Before comparison, polymer chains with the same
sequence are given one entity. Only entity IDs change.

## Reference evaluation

```bash
core-score inputs.csv --out results/ --rmsd-thresholds 0.5 1 1.5 2
core-score poses/*.cif --smiles 'CCO' --reference experimental.cif --out results/
```

Predictions are compared with the reference wherever their ensemble has one, set by the
`reference` column or by `--reference`. Each prediction is compared with OST's
`compare-ligand-structures`, giving the ligand RMSD after binding-site superposition and
lDDT-PLI. A pose is **correct** at T when its RMSD is below T and its lDDT-PLI is above
0.8, the definition used in the paper. A multi-copy comparison has no label unless OST
assigns exactly one pair. A failed comparison, or one missing either score, is unknown,
never wrong.

## Top-ranked pose

When confidences are given, the top-ranked pose of each complex is its
PoseBusters-valid pose with the highest confidence, as in the paper. Ties go to the
earlier pose. A complex with no valid pose, no confidences, or PoseBusters skipped has
no top pose, and `note` in `core_scores.csv` says why.

## Running

Every run loads and checks the inputs, prepares structures, runs PoseBusters, assigns
sites where needed, compares every pair, calculates CORE, evaluates references when
given, and writes the report. Each run is fresh and needs a new or empty
output directory. A missing or failed comparison prevents that ensemble's CORE
Score. Nothing is averaged over a subset or filled in. Other ensembles carry on, and
their results and diagnostics are saved. The exit code is 0 when everything
completed, 1 when anything is incomplete and 2 for invalid input.

## Outputs

A complex is one ensemble, or one ligand site of a multi-copy ensemble.

```text
results/
  core_scores.csv     one row per complex
  poses.csv           one row per pose
  core_scores.png     every complex ordered by CORE Score, with --plot
  details/            the evidence behind them
```

**`core_scores.csv`** has one row per complex and T, highest CORE first and unscored
complexes last.

| Column | Meaning |
| --- | --- |
| `ensemble` | The complex, either the ensemble's name or `<ensemble>_site_<n>` for each site of a ligand with several copies. |
| `threshold_angstrom`, `core_score` | T and CORE at T. Blank when CORE could not be calculated. |
| `above_0.5` … `above_0.9` | True when CORE ≥ that cutoff, one column per cutoff. |
| `pose_count`, `valid_pose_count` | Poses in the complex and how many are PoseBusters-valid. |
| `top_pose`, `top_pose_confidence` | The top-ranked pose, when confidences are given. |
| `fraction_correct`, `top_pose_correct` | With references, the fraction of poses correct at T and whether the top pose is. |
| `note` | Why CORE or the top pose is missing. |

**`poses.csv`** has one row per pose, in path order.

| Column | Meaning |
| --- | --- |
| `ensemble`, `pose_id` | The pose. |
| `confidence` | The model's own score, when confidences are given. |
| `pb_valid` | True when PoseBusters passed every dock check on every ligand copy. Which checks failed is in `details/posebusters.json`. |
| `rmsd`, `lddt_pli` | To the reference, when there is one. |
| `correct_2A` | Whether the pose is correct at T, RMSD below T and lDDT-PLI above 0.8, one column per T. |
| `note` | Why a result is missing. |

A blank is unknown, never a pass or a failure.

**`core_scores.png`**, with `--plot`, is a dot per complex, ordered by CORE Score, one PNG
per T, with a dashed line at each cutoff saying how many complexes reach it. It is set in
Optima where that font is installed, otherwise the nearest humanist sans the machine has.
It needs matplotlib, from `pip install '.[plot]'`.

**`details/`** holds the evidence behind those files.

| File | Contents |
| --- | --- |
| `pairs.csv` | Ligand RMSD, pocket backbone RMSD, pocket lDDT and lDDT-PLI per pair, with radius recovery fields |
| `ligand_details.jsonl` | Raw OST results, commands and output for every pair |
| `posebusters.json` | PoseBusters validity, every dock check and every raw result per pose and copy |
| `reference_details.json`, `references/` | Raw reference comparisons and prepared references |
| `prepared_manifest.json`, `prepared/` | Prepared structures, SDFs and receptors |
| `sites.json` | Site assignments for multi-copy ensembles |
| `run.json`, `status.json` | Settings and completion counts |

Every failure records its stage, error and traceback, plus any evidence such as the
OST command and output.

## Code layout

`pipeline.py` runs the steps in order. Each folder holds one phase of the run.

```text
src/core_score/
  cli.py                  the core-score command
  pipeline.py             run() and execute(): every step in order
  inputs/
    make_csv.py           the core-score-csv command, which writes the input CSV
    inputs.py             load and check the CSV or structure inputs
    structures.py         read structure records, select ligands, write ligand SDFs
    prepare.py            prepare every pose for OST and PoseBusters
    sites.py              one ensemble per ligand site
  analysis/
    posebusters.py        dock checks for every pose
    ost.py                every OpenStructure call
    pairs.py              every pair of poses, with binding-site recovery
    core.py               RMSD matrix and CORE per threshold
    reference.py          optional comparison with the ground truth
  outputs/
    report.py             core_scores.csv and poses.csv
    plot.py               core_scores.png, with --plot
    files.py              output writers, failure records, output folders
```

[The code map](https://github.com/James-Broster/core-score/blob/main/docs/code-map.md) draws the whole run, then each stage, function by
function, and lists every check.

## Tests

```bash
python -m pytest
```

Real OST tests are separate because they need OST and structure fixtures. See
[validation](https://github.com/James-Broster/core-score/blob/main/docs/validation.md).
