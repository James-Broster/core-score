# CORE Score

A Python package for prediction-only ensemble agreement. It prepares protein–ligand
structures, runs PoseBusters `dock` checks, compares every pair of poses with
OpenStructure, and calculates CORE from the complete pairwise ligand-RMSD matrix.
Experimental references are optional and never enter the CORE calculation.

See [validation](docs/validation.md) for what has been checked.

## CORE Score

For each pose *i* in an ensemble of *N* poses, δ*ᵢ* is the ligand RMSD to its closest
other pose. Each distance is capped at a threshold *T* and divided by *T*, so a pose
reproduced exactly has penalty 0 and a pose with no other pose within *T* has
penalty 1. CORE is one minus the mean penalty.

```text
S = 1 − (1/N) Σᵢ min(δᵢ, T) / T
```

Equivalently S = f (1 − ν/T), where f is the fraction of poses with another pose
within T and ν is their mean closest-pose distance. The default T is 2 Å, the RMSD
threshold for a correct prediction. CORE uses every pose in the ensemble, including
poses that fail PoseBusters checks, because they still sample the model's output
distribution. Validity filtering belongs downstream, for example when choosing the
top-ranked prediction.

`core_score(matrix, threshold=2.0)` is importable from `core_score` for matrices
computed elsewhere.

## Install

```bash
conda env create -f environment.yml
conda activate core-score
```

The conda solve has been checked on macOS ARM. Other platforms need an available
OpenStructure 2.12 build. Alternatively, activate an environment containing
**OpenStructure 2.12** and install with pip.

```bash
python -m pip install -e '.[plots,test]'
core-score --help
```

PoseBusters is pinned at **0.4.4**, the version used for the existing analysis.
OpenStructure is an external executable. `ost` is launched from PATH and finds its
chemical library itself or through `OST_COMPOUNDS_CHEMLIB`.

## Inputs

```bash
core-score inputs.json --out results/ --workers 8 --plots
core-score poses/*.cif --smiles 'CCO' --out results/
```

A JSON manifest can hold many ensembles. Structure files given directly form one
ensemble. PDB, CIF and mmCIF are accepted, each with one structural model. Every
prediction needs the SMILES of its intended ligand (`CCO` here is only an example).

```json
{
  "ensembles": [
    {
      "ensemble": "target1_boltz",
      "smiles": "CCO",
      "ligand_chains": ["B"],
      "poses": [
        {"pose_id": "p0", "path": "boltz/run_1/model.cif"},
        {"pose_id": "p1", "path": "boltz/run_2/model.cif"}
      ]
    }
  ]
}
```

Each ensemble contains predictions of the same protein–ligand system to compare
together. Keep different methods in separate ensembles to score them separately.
Ensemble names must be unique. Paths are relative to the JSON file. `smiles`,
`ligand_chains` and `ligand_residue` can be set once per ensemble and overridden per
pose. Pose IDs default to list order. Optional pose fields are `seed`, `sample` and
`ligand_sdfs`, a list of original predicted SDFs in the receptor's coordinate frame
for structures whose ligand is supplied separately. Coordinates are never aligned or
changed.

A structure with exactly one nonpolymer ligand residue is selected automatically.
Otherwise the SMILES heavy-atom composition, `ligand_chains` (or `--ligand-chain`,
repeated for identical copies) and `ligand_residue` select it. Ambiguous or missing
ligands are reported, not guessed. Composition only selects residues. Preparation
separately verifies the chemistry against the SMILES.

Poses are ordered by `(seed, sample, pose_id)`. The first is the anchor for site
assignment, and in each pair the earlier pose is OST's model.

See [the input example](examples/manifest.json). This is a CORE Score manifest of
generated structures, not an importer for a prediction model's own input JSON.

## Options

| Option | Meaning |
| --- | --- |
| `--workers N` | At most N simultaneous OST comparisons (default 1). Each OST process uses one numerical-library thread. |
| `--timeout SECONDS` | Limit per OST call. Default is no limit. |
| `--rmsd-thresholds 0.5 1 2` | Thresholds T in Å. CORE is calculated at each from the same matrix, and reference success means RMSD < T. Default 2. |
| `--no-radius-fallback` | Disable the 4 → 6 → 8 Å binding-site recovery. |
| `--skip-posebusters` | Do not run PoseBusters. CORE is unaffected. |
| `--evaluate-reference` | Also compare predictions with supplied reference structures. |
| `--reference FILE` | Reference for structure files given directly (JSON sets its own). |
| `--plots` | Save heatmaps and, with references, a summary plot. |

## PoseBusters

PoseBusters runs once per selected ligand copy against its predicted environment,
using the official API with `config='dock'`. A pose is valid only if every copy passes
every dock binary check. If none fail but any are unknown, the pose is unknown,
never valid. Missing chemical information is unknown, not a pass. Results are saved
in `posebusters.json` (every raw check) and `posebusters.csv` (one row per pose).
No pose is removed.

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
anchor are saved in `sites.json`.

## Reference evaluation

```bash
core-score inputs.json --out results/ --evaluate-reference --rmsd-thresholds 0.5 1 1.5 2
core-score poses/*.cif --smiles 'CCO' --out results/ --evaluate-reference --reference experimental.cif
```

A JSON ensemble can carry a `reference` path or object with its own ligand selection
and chemistry (see [the reference example](examples/manifest-with-reference.json)).
References are loaded and checked only when `--evaluate-reference` is given. Each
prediction is compared with the reference using OST's ligand RMSD. A multi-copy
comparison has no single success label unless OST assigns exactly one pair.
`reference_summary.csv` places each ensemble's CORE Score beside the fraction of its
poses within the same threshold. It is descriptive, not a calibration.

## Running

Every run loads and checks the inputs, prepares structures, runs PoseBusters, assigns
sites where needed, compares every pair, and calculates CORE, then evaluates
references and draws plots when requested. Each run is fresh and needs a new or empty
output directory. A missing or failed comparison prevents that ensemble's CORE
Score. Nothing is averaged over a subset or filled in. Other ensembles carry on, and
their results and diagnostics are saved. The exit code is 0 when everything
completed, 1 when anything is incomplete and 2 for invalid input.

## Outputs

| File | Contents |
| --- | --- |
| `scores.csv`, `scores.json` | CORE Score per ensemble, site and threshold, or the reason it is missing |
| `pairs.csv` | Ligand RMSD, pocket backbone RMSD, pocket lDDT and lDDT-PLI per pair, with radius recovery fields |
| `ligand_details.jsonl` | Raw OST results, commands and output for every pair |
| `posebusters.csv`, `posebusters.json` | Dock check results per pose and copy |
| `prepared_manifest.json`, `prepared/` | Prepared structures, SDFs and receptors |
| `sites.json` | Site assignments for multi-copy ensembles |
| `reference_evaluation.csv`, `reference_summary.csv`, `reference_details.json` | Reference results, when requested |
| `run.json`, `status.json` | Settings and completion counts |

Every failure records its stage, error and traceback, plus any evidence such as the
OST command and output.

## Code layout

Reading order is `cli.py → pipeline.execute() → stages/<stage>.py → stages/utils/`.
Stages run the steps and write results. Utils compute.

```text
src/core_score/
  cli.py                      parse the command
  inputs.py                   load and check JSON or structure inputs
  pipeline.py                 run the stages in order
  results.py                  output writers and diagnostics
  stages/
    preparation.py            prepare structures and ligand chemistry
    posebusters.py            dock checks for every pose
    site_assignment.py        split multi-copy ensembles by site
    comparison.py             compare every pair
    scoring.py                matrix and CORE per threshold
    evaluation.py             optional reference comparison
    plotting.py               figures from saved results
    utils/
      ost.py                  every OpenStructure call
      preparation/
        prepare.py            prepare one structure
        chemistry.py          ligand selection and chemical identity
        cif.py                read and filter structure records
        sdf.py                coordinate-preserving ligand molecules
      posebusters_utils.py    dock checks and validity
      sites.py                protein superposition and site assignment
      comparison_utils.py     pair tasks and OST score summaries
      radius.py               binding-site recovery
      score_utils.py          RMSD matrix and CORE
      evaluation_utils.py     reference labels and summary
```

## Tests

```bash
python -m pytest
```

Real OST tests are separate because they need OST and structure fixtures. See
[validation](docs/validation.md).
