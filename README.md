# CORE Score

A portable Python package for prediction-only ensemble agreement. It prepares
protein–ligand structures, runs PoseBusters `dock` checks and excludes invalid poses, compares predictions
with OpenStructure, and calculates CORE from a complete pairwise ligand-RMSD
matrix. Experimental references are optional and never enter the CORE calculation.

This is the first refactored implementation. See [validation](docs/validation.md)
for what has actually been checked and remaining historical regression coverage.

## Install

From the project directory, an environment recipe is provided:

```bash
conda env create -f environment.yml
conda activate core-score
```

Its conda dependency solve has been checked on macOS ARM; other platforms still
require an available OpenStructure build. The package itself has no scheduler or
machine-specific paths. Alternatively:

Activate an environment containing **OpenStructure 2.12**, then install:

```bash
python -m pip install -e '.[plots,test]'
core-score --help
```

The Python dependencies include PoseBusters **0.4.4**, preserving the production
version used for the existing analysis. OpenStructure is an external executable,
not installed by pip. `ost` is launched through PATH. OST discovers its own chemical library from
its installation or `OST_COMPOUNDS_CHEMLIB`; the pipeline does not run a separate
installation preflight. There are no Oxford paths or Slurm rules.
Custom installations can use advanced `--ost-bin` / `--chemlib` overrides.
The optional protein-superposition Python API requires OST Python bindings;
the default workflow calls OST through its executable.

## Start with files

```bash
core-score inputs.json --out results/ --workers 8 --plots
core-score poses/*.cif --smiles 'CCO' --out results/
```

PDB, CIF and mmCIF are accepted. Each file must contain one structural model.
Use an original coordinate-bearing ligand SDF with `--ligand-sdf ligand.sdf`
when the receptor and ligand are supplied separately. It must be in the receptor's
coordinate frame; we do not align it silently. SMILES describe the **actual**
intended ligand; `CCO` above is an illustrative ethanol input, not a universal default.

When a structure has exactly one nonpolymer ligand residue it can be selected
automatically. Otherwise supply `--ligand-chain` (repeat for identical copies),
optionally `--ligand-residue`, or SMILES to select matching heavy-atom composition.
Composition is only a selection aid; PB preparation separately verifies chemistry.
Ambiguities and missing ligands are reported, not guessed away.

The CLI requires the intended ligand SMILES for each prediction, supplied in JSON
or with `--smiles`. It checks this before preparation. Native chemical information
and original SDFs remain useful during preparation. Missing chemical information is **unknown**, not
a PB pass. Coordinate-derived connectivity is preserved as a separate OST fallback;
it is not treated as a chemically authoritative PB molecule.

Input validation runs automatically: input loading checks structure paths and ligand
SMILES, preparation checks structure readability and ligand selection, and the
pipeline launches OST when comparisons are needed, recording execution failures. Reference files are required
only for requested reference evaluation. CORE needs at least two retained poses.

## JSON inputs

```bash
core-score inputs.json --out results/
```

JSON is the input manifest format. Each ensemble contains predictions of the same
protein–ligand system to compare together. Files can live in separate directories;
ensemble membership is explicit, not inferred from their locations. Keep different
methods in separate ensembles when you want separate scores for them.

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

Paths are relative to the JSON file; absolute paths also work. Shared `smiles`,
`ligand_chains` and `ligand_residue` can be set once per ensemble; individual poses
can override them. Supply `smiles` for the intended ligand; chain and residue selectors are optional.
`CCO` is an illustrative ligand, not a default. Pose IDs may be omitted and are then
generated from list order. Optional pose fields include `seed`, `sample`, and
`ligand_sdfs` (an array of original predicted SDF paths in the receptor coordinate
frame). An SDF is not required if reliable chemistry can be obtained otherwise.

An optional ensemble-level `reference` object has its own `path`, ligand selection
and chemistry; see [the reference example](examples/manifest-with-reference.json).
It is evaluated only when requested with `--evaluate-reference`.

See [the input example](examples/manifest.json). This is a **CORE Score manifest**
pointing at generated structures, not an automatic importer for a prediction
model's own input JSON schema. Direct structure paths remain available for one
ensemble. CSV input manifests are no longer supported; result tables and the advanced
PB-results CSV interface remain available.

Existing single-ensemble and ensemble-list `cofold-ost` JSON remains readable.
The pipeline writes `prepared_manifest.json` recording prepared structures.
Prepared files must remain available at the recorded paths; move the whole result
only with path updates.

## Defaults and controls

| Option | Meaning |
| --- | --- |
| `--workers N` | At most N simultaneous local pair comparisons, default 1. Each OST subprocess caps numerical-library threads at one. Not CPUs per comparison. |
| `--timeout SECONDS` | Optional limit per OST invocation. Default is **no timeout**. |
| `--no-radius-fallback` | Disable the default eligible 4 → 6 → 8 Å binding-site recovery. |
| `--skip-posebusters` | Skip both PB checking and filtering; retain all input poses. |
| `--pb-results FILE` | Advanced override using existing PB results; otherwise PB runs automatically. |
| `--evaluate-reference` | Also evaluate against explicitly supplied reference structures. |
| `--rmsd-thresholds 0.5 1 1.5 2` | Reference success labels use **RMSD < threshold**, in Å. Default 2. Does not change CORE. |
| `--plots` | Export optional heatmaps and, when available, a reference-evaluation plot. |

Recovery radii are fixed at 6 then 8 Å after the initial 4 Å attempt. Only recognized binding-site
failure reasons are eligible. Graph errors and general process failures are not
fixed by enlarging the pocket. Fallback operates on one assigned copy/site at a
time. Pocket/RMSD results can come from recovery, but **lDDT-PLI stays from the
original 4 Å call**. All attempts, original RMSD and chosen radius remain recorded.

By default, PB runs once per selected ligand copy against its predicted environment
and invalid poses are excluded before comparisons and CORE. Use
`--skip-posebusters` to skip both checking and filtering. There is no check-and-retain
CLI mode. All
configured `dock` binary checks determine `pb_valid`. A separate
`pb_valid_13_checks` retains the historical paper's 13-check summary; these are
different definitions. Every raw check is exported. A pose fails if any copy
fails; if none fail but any are unknown, the pose is unknown. Default PB filtering
blocks that ensemble when status is unknown rather than silently treating it as
valid. Exclusions and original/retained ensemble sizes are saved.

The wrapper uses the official [PoseBusters Python API](https://posebusters.readthedocs.io/en/latest/api.html)
with `config='dock'` and full reports.

Supplied PB CSVs require unique `(ensemble, pose_id)` rows, `config=dock` and
`pb_valid` (`True`, `False`, blank for unknown). `posebusters.csv` from this package
can be reused. Supplied labels are user-provided evidence, not revalidated chemistry.

## Multiple identical copies

The default workflow uses the historical production approach: choose the earliest
`(seed, sample, pose_id)` prediction as an anchor, use OST protein rigid alignment,
then constrained centroid assignment to give each copy a consistent site label.
Each site gets its own pairwise matrix and CORE score. Confidence and ground truth
are not used. Copy-count mismatches or alignment failures block the ensemble and
are reported. Site labels are arbitrary and do not prove biological correspondence.
The assignment costs, transforms and anchor are exported in `sites.json`.

The lower-level site utilities and their tests are retained. The historical
site-assignment policy has not been replaced with the experimental multistart
policy. There is no `--site-mode` flag.

## Optional reference evaluation

```bash
core-score inputs-with-references.json --out results/ \
  --evaluate-reference --rmsd-thresholds 0.5 1 1.5 2 --plots

core-score inputs.json --out evaluation/ --evaluate-reference \
  --reference experimental.cif
```

Reference selections are separate from prediction selections. Ground-truth paths
are checked only when reference evaluation is requested; missing requested
references are errors before preparation. Unused reference paths do not block a
normal CORE run. Reference comparisons use OST's ligand RMSD and preserve
the same diagnostics. Threshold sweeps reuse each computed RMSD. Multiple assigned
reference copies have no scalar correctness label unless OST assigns exactly one
pair; site-specific ground-truth correspondence is not inferred from arbitrary site
numbers. `reference_summary.csv` is descriptive evidence, not a trained calibration
or a promise that a given CORE cutoff achieves a particular accuracy.

## Running the workflow

Use `core-score inputs.json --out results/ --plots`. JSON is the primary input;
direct structure paths with `--smiles` remain available for one ensemble.
There are no stage subcommands. Every run loads and validates inputs, prepares
structures, applies the PB policy, assigns sites where needed, compares all retained
pairs, validates matrices, and calculates CORE. Requested reference evaluation is
part of that same run. Requested plots are generated last.
Use `core-score --help` for workflow options.

Use `--workers N` to run up to N pair comparisons simultaneously within one run.
The pipeline schedules every required pair and collects the results automatically.
Each invocation is a fresh run. Use a new or empty output directory; interrupted
runs must be started again. No checkpoint files are created or reused.
There is no pair-subset input or separate result-merging command.

Duplicate, missing, failed or incompatible comparisons prevent CORE calculation.
There is no averaging over successful subsets and no zero-filling failures. The
package never submits a cluster job. Independent ensembles keep running after a
local scientific failure; their results and diagnostics are included in the output.

## Outputs and code layout

`pairs.csv` contains ligand RMSD, pocket backbone RMSD, pocket lDDT, lDDT-PLI and
fallback metadata. `ligand_details.jsonl` retains raw OST results, commands,
stdout/stderr, inputs and attempts. `scores.json/csv` gives CORE and winning-group
pose IDs; `posebusters.json/csv` retains PB checks. Input preparation, site mapping,
run settings and status are separately recorded. Errors identify the stage and
known next action; unknown exceptions retain their original evidence.

The reading order is `cli.py → pipeline.execute() → stages/<stage>.py → stages/utils/`.
The pipeline shows the complete workflow. Each stage shows its own sequence of work;
its helpers contain the implementation details. Preparation-specific chemistry and
file handling are grouped together inside the stage helpers. JSON input loading
and shared result-writing helpers stay at the package root.

```text
src/core_score/
  cli.py                      parse the user's command
  inputs.py                   load JSON and validate required input fields
  pipeline.py                 orchestrate the complete run
  results.py                  shared output and diagnostic helpers
  stages/
    preparation.py            prepare structures and ligand chemistry
    site_assignment.py        establish consistent ligand sites when needed
    posebusters.py            check → filter → save
    comparison.py             prepare pairs → run workers → save
    scoring.py                matrix → CORE → score export
    evaluation.py             reference comparisons → threshold labels → save
    plotting.py               generate figures from saved results
    utils/
      prepare_utils.py        prepare an individual structure and ligand
      sites_utils.py          align centres and select/write ligand copies
      site_algorithms.py      protein alignment and site-assignment algorithms
      posebusters_utils.py    dock calls, copy results and validity handling
      comparison_utils.py     OST calls, probes, tasks and diagnostics
      matrix.py               complete-matrix validation
      radius.py               established binding-site recovery
      score_utils.py          clustering and CORE mathematics
      evaluation_utils.py     reference preparation, thresholds and summaries
      preparation/
        chemistry.py          ligand selection and chemical identity
        cif.py                read/filter structural records
        sdf.py                coordinate-preserving ligand representations
```


No legacy runtime imports, cluster launchers, generated egg metadata or historical
logs are required. The original Downloads copy is retained outside this package.

## Tests

```bash
python -m pytest
```

Real OST regression tests are explicitly separate because they need OST and
structure fixtures. See [validation](docs/validation.md). Missing requested
integration fixtures cause an error; we do not silently skip that validation.
