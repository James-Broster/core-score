# Code map

How a run moves through `core_score`. Every diagram reads the same way:

- It runs top to bottom, one rung per step.
- Each rung lists its functions left to right, in the order they run.
- The rounded boxes between rungs are what's handed to the next step.
- Dashed boxes are files written. Amber functions call OpenStructure through `analysis/ost.py`.

`pipeline.execute()` runs the rungs of the first diagram in order. The result files go
in the output folder and everything behind them in `details/`.

## The whole run

```mermaid
%%{init: {"flowchart": {"curve": "linear", "nodeSpacing": 18, "rankSpacing": 28, "subGraphTitleMargin": {"top": 6, "bottom": 10}}}}%%
flowchart TB
    classDef fn fill:#e7eff7,stroke:#5b7fa3,color:#13293d
    classDef ost fill:#fbedd5,stroke:#b07a1e,color:#3d2a08
    classDef data fill:#ffffff,stroke:#6f8f4e,color:#1f2d12,stroke-width:2px
    classDef file fill:#f3f3f3,stroke:#9a9a9a,color:#444,stroke-dasharray:3 3
    classDef stage fill:#f7f9fb,stroke:#c9d2dc,color:#1b2430

    start(["core-score inputs.csv --out results/<br/>or core_score.run()"]):::data

    subgraph R1["1 · inputs/inputs.py"]
        direction LR
        a1["load_inputs()"]:::fn --> a2["read_csv()"]:::fn --> a3["make_group()"]:::fn --> a4["make_poses()"]:::fn
    end
    d1(["ensembles"]):::data

    subgraph R2["2 · inputs/prepare.py"]
        direction LR
        b1["prepare_groups()"]:::fn --> b2["prepare()"]:::fn --> b3["prepare_pose()"]:::fn --> b4["reads_ligand()"]:::ost --> b5[/"prepared/ · prepared_manifest.json"/]:::file
    end
    d2(["prepared poses"]):::data

    subgraph R3["3 · analysis/posebusters.py"]
        direction LR
        c1["run_posebusters()"]:::fn --> c2["check_pose()"]:::fn --> c3["check_copy()"]:::fn --> c4[/"posebusters.json"/]:::file
    end
    d3(["prepared poses and pb_valid"]):::data

    subgraph R4["4 · inputs/sites.py"]
        direction LR
        e1["prepare_site_groups()"]:::fn --> e2["split_sites()"]:::fn --> e3["align_copy_centres()"]:::ost --> e4["assign_sites()"]:::fn --> e5[/"sites.json"/]:::file
    end
    d4(["one ensemble per site"]):::data

    subgraph R5["5 · analysis/pairs.py"]
        direction LR
        g1["compare_groups()"]:::fn --> g2["pair_tasks()"]:::fn --> g3["pair_job()"]:::fn --> g4["compare_pair()"]:::ost --> g5["summarize()"]:::fn --> g6[/"pairs.csv · ligand_details.jsonl"/]:::file
    end
    d5(["one row per pair"]):::data

    subgraph R6["6 · analysis/core.py"]
        direction LR
        h1["score_groups()"]:::fn --> h2["build_matrix()"]:::fn --> h3["core_score()"]:::fn
    end
    d6(["CORE per complex and T"]):::data

    subgraph R7["7 · analysis/reference.py"]
        direction LR
        i1["evaluate_groups()<br/>ensembles with a reference"]:::fn --> i2["prepare()<br/>the reference"]:::fn --> i3["run_pair_tasks()"]:::ost --> i4["label_reference_results()"]:::fn --> i5[/"reference_details.json"/]:::file
    end
    d7(["CORE and reference labels"]):::data

    subgraph R8["8 · outputs/report.py"]
        direction LR
        j1["write_report()"]:::fn --> j2["core_table()"]:::fn --> j3["top_pose()"]:::fn --> j4["pose_table()"]:::fn --> j5[/"core_scores.csv · poses.csv"/]:::file
    end

    subgraph R9["9 · outputs/plot.py and pipeline.py"]
        direction LR
        k0["plot_core_scores()<br/>with --plot"]:::fn --> k1[/"core_scores.png"/]:::file --> k2["summarize_run()"]:::fn --> k3[/"status.json · exit code"/]:::file
    end

    start --> R1 --> d1 --> R2 --> d2 --> R3 --> d3 --> R4 --> d4 --> R5 --> d5 --> R6 --> d6 --> R7 --> d7 --> R8 --> R9

    class R1,R2,R3,R4,R5,R6,R7,R8,R9 stage
```

## 1 · Inputs

`core-score-csv` writes the CSV beforehand, from a folder of predictions. It is its own
command and not part of a run. [Its section](#core-score-csv) lists its functions.

```mermaid
%%{init: {"flowchart": {"curve": "linear", "nodeSpacing": 18, "rankSpacing": 28, "subGraphTitleMargin": {"top": 6, "bottom": 10}}}}%%
flowchart TB
    classDef fn fill:#e7eff7,stroke:#5b7fa3,color:#13293d
    classDef data fill:#ffffff,stroke:#6f8f4e,color:#1f2d12,stroke-width:2px
    classDef stage fill:#f7f9fb,stroke:#c9d2dc,color:#1b2430

    start(["command line, or core_score.run()"]):::data

    subgraph I1["cli.py"]
        direction LR
        a1["main()"]:::fn --> a2["parser()"]:::fn --> a3["positive() · fraction() · worker_count()<br/>type check each flag"]:::fn
    end
    d1(["paths and settings"]):::data

    subgraph I2["inputs/inputs.py · read"]
        direction LR
        b1["load_inputs()"]:::fn --> b2["read_csv()<br/>or one row per structure file"]:::fn
    end
    d2(["rows, gathered by ensemble"]):::data

    subgraph I3["inputs/inputs.py · build"]
        direction LR
        c1["make_group()<br/>one SMILES and reference,<br/>SMILES parses"]:::fn --> c2["make_poses()<br/>wildcards expanded"]:::fn --> c3["resolve()<br/>every path must exist"]:::fn
    end
    d3(["ensembles, poses sorted by path"]):::data

    subgraph I4["inputs/inputs.py · finish"]
        direction LR
        e1["at least four poses<br/>per ensemble"]:::fn
    end
    d4(["ensembles, handed to pipeline.execute()"]):::data

    start --> I1 --> d1 --> I2 --> d2 --> I3 --> d3 --> I4 --> d4

    class I1,I2,I3,I4 stage
```

## 2 · Preparation

The functions in this stage come from `inputs/prepare.py` and the helpers in
`inputs/structures.py`.

```mermaid
%%{init: {"flowchart": {"curve": "linear", "nodeSpacing": 18, "rankSpacing": 28, "subGraphTitleMargin": {"top": 6, "bottom": 10}}}}%%
flowchart TB
    classDef fn fill:#e7eff7,stroke:#5b7fa3,color:#13293d
    classDef ost fill:#fbedd5,stroke:#b07a1e,color:#3d2a08
    classDef data fill:#ffffff,stroke:#6f8f4e,color:#1f2d12,stroke-width:2px
    classDef file fill:#f3f3f3,stroke:#9a9a9a,color:#444,stroke-dasharray:3 3
    classDef stage fill:#f7f9fb,stroke:#c9d2dc,color:#1b2430

    start(["one ensemble's poses"]):::data

    subgraph P1["inputs/prepare.py"]
        direction LR
        a1["prepare_groups()"]:::fn --> a2["prepare()<br/>each pose"]:::fn --> a3["pose_directory()"]:::fn --> a4["prepare_pose()"]:::fn
    end
    d1(["one pose and its folder"]):::data

    subgraph P2["inputs/structures.py · read"]
        direction LR
        b1["read_document()<br/>CIF, or PDB converted"]:::fn --> b2["merge_identical_entities()<br/>one entity per sequence"]:::fn
    end

    subgraph P3["from_structure() · pick"]
        direction LR
        c1["choose_ligands()"]:::fn --> c2["ligand_residues()"]:::fn --> c3["composition()<br/>matches the SMILES"]:::fn
    end
    d3(["selected ligand residues"]):::data

    subgraph P4["inputs/structures.py · each copy"]
        direction LR
        m1["make_ligand()"]:::fn --> m2["coordinate_mol()"]:::fn --> m3["graph from native_graph(),<br/>else the peer cache,<br/>else bonds from coordinates"]:::fn --> m4["matched_template()<br/>SMILES chemistry onto the atoms"]:::fn --> m5["identity, coordinate<br/>and atom order checks"]:::fn
    end
    d4(["OST and PoseBusters molecules per copy"]):::data

    subgraph P5["from_structure() · write"]
        direction LR
        w1["write_molecule()<br/>ost.sdf and pb.sdf"]:::fn --> w2["filter_records()<br/>pose.cif"]:::fn
    end

    subgraph P6["from_original_sdfs()"]
        direction LR
        o1["read_original()<br/>matches the SMILES"]:::fn --> o2["same_atoms()<br/>finds the residue it replaces"]:::fn --> o3["write pose.cif"]:::fn
    end
    d5(["copy records"]):::data

    subgraph P7["inputs/structures.py · receptor"]
        direction LR
        r1["write_receptor()"]:::fn --> r2[/"pose.receptor.pdb"/]:::file
    end
    d6(["prepared record"]):::data

    subgraph P8["inputs/prepare.py · OST"]
        direction LR
        t1["reads_ligand()<br/>structure route only"]:::ost --> t2["compare_ligands()<br/>with its own SDFs"]:::ost --> t3["if OST can't read it or connects it differently,<br/>use the ost.sdf files"]:::fn
    end
    d7(["prepared pose, or pose with an error"]):::data

    subgraph P9["inputs/prepare.py · retry"]
        direction LR
        f1["retry poses with a pb_error<br/>once the peer cache is full"]:::fn --> f2[/"prepared_manifest.json"/]:::file
    end

    start --> P1 --> d1 --> P2
    P2 -->|"no ligand_sdfs"| P3
    P2 -->|"ligand_sdfs given"| P6
    P3 --> d3 --> P4 --> d4 --> P5 --> d5
    P6 --> d5
    d5 --> P7 --> d6 --> P8 --> d7 --> P9

    class P1,P2,P3,P4,P5,P6,P7,P8,P9 stage
```

## 3 · PoseBusters

```mermaid
%%{init: {"flowchart": {"curve": "linear", "nodeSpacing": 18, "rankSpacing": 28, "subGraphTitleMargin": {"top": 6, "bottom": 10}}}}%%
flowchart TB
    classDef fn fill:#e7eff7,stroke:#5b7fa3,color:#13293d
    classDef data fill:#ffffff,stroke:#6f8f4e,color:#1f2d12,stroke-width:2px
    classDef file fill:#f3f3f3,stroke:#9a9a9a,color:#444,stroke-dasharray:3 3
    classDef stage fill:#f7f9fb,stroke:#c9d2dc,color:#1b2430

    start(["prepared poses"]):::data

    subgraph B1["analysis/posebusters.py · setup"]
        direction LR
        a1["run_posebusters()"]:::fn --> a2["PoseBusters(config='dock')"]:::fn --> a3["dock_checks()"]:::fn --> a4["normalized()"]:::fn
    end
    d1(["names of the 21 dock checks"]):::data

    subgraph B2["analysis/posebusters.py · check"]
        direction LR
        c1["check_pose()<br/>each pose"]:::fn --> c2["check_copy()<br/>each ligand copy"]:::fn --> c3["buster.bust()<br/>pb.sdf against the receptor"]:::fn --> c4["aggregate()<br/>False, then unknown, then True"]:::fn
    end
    d2(["pb_valid and every check, per copy and pose"]):::data

    subgraph B3["analysis/posebusters.py · write"]
        direction LR
        e1["write_json()"]:::fn --> e2[/"posebusters.json"/]:::file
    end
    d3(["records, kept for the report"]):::data

    start --> B1 --> d1 --> B2 --> d2 --> B3 --> d3

    class B1,B2,B3 stage
```

## 4 · Site assignment

Only ensembles with several copies of the same ligand are split. Single-copy ensembles
pass straight through as site 1.

```mermaid
%%{init: {"flowchart": {"curve": "linear", "nodeSpacing": 18, "rankSpacing": 28, "subGraphTitleMargin": {"top": 6, "bottom": 10}}}}%%
flowchart TB
    classDef fn fill:#e7eff7,stroke:#5b7fa3,color:#13293d
    classDef ost fill:#fbedd5,stroke:#b07a1e,color:#3d2a08
    classDef data fill:#ffffff,stroke:#6f8f4e,color:#1f2d12,stroke-width:2px
    classDef file fill:#f3f3f3,stroke:#9a9a9a,color:#444,stroke-dasharray:3 3
    classDef stage fill:#f7f9fb,stroke:#c9d2dc,color:#1b2430

    start(["prepared poses"]):::data

    subgraph S1["site_assignment · check"]
        direction LR
        a1["prepare_site_groups()"]:::fn --> a2["no failed poses,<br/>same copy count"]:::fn --> a3["one copy<br/>passes through, named for the ensemble"]:::fn
    end
    d1(["ensembles with k copies per pose"]):::data

    subgraph S2["site_assignment · align"]
        direction LR
        b1["split_sites()"]:::fn --> b2["align_copy_centres()"]:::fn --> b3["superpose()<br/>each pose onto the anchor"]:::ost
    end
    d2(["copy centres in the anchor's frame"]):::data

    subgraph S3["site_assignment · assign"]
        direction LR
        c1["assign_sites()"]:::fn --> c2["k-means<br/>seeds k centres"]:::fn --> c3["Hungarian match<br/>one copy per site per pose"]:::fn --> c4["repeat until<br/>labels stop changing"]:::fn
    end
    d3(["site label for every copy"]):::data

    subgraph S4["site_assignment · split"]
        direction LR
        e1["site_pose()<br/>one copy per site"]:::fn --> e2["filter_records()<br/>siteN.cif"]:::fn --> e3["named<br/>ensemble_site_N"]:::fn --> e4[/"sites.json"/]:::file
    end
    d4(["one ensemble per site"]):::data

    start --> S1 --> d1 --> S2 --> d2 --> S3 --> d3 --> S4 --> d4

    class S1,S2,S3,S4 stage
```

## 5 · Comparison

```mermaid
%%{init: {"flowchart": {"curve": "linear", "nodeSpacing": 18, "rankSpacing": 28, "subGraphTitleMargin": {"top": 6, "bottom": 10}}}}%%
flowchart TB
    classDef fn fill:#e7eff7,stroke:#5b7fa3,color:#13293d
    classDef ost fill:#fbedd5,stroke:#b07a1e,color:#3d2a08
    classDef data fill:#ffffff,stroke:#6f8f4e,color:#1f2d12,stroke-width:2px
    classDef file fill:#f3f3f3,stroke:#9a9a9a,color:#444,stroke-dasharray:3 3
    classDef stage fill:#f7f9fb,stroke:#c9d2dc,color:#1b2430

    start(["one ensemble per site"]):::data

    subgraph C1["analysis/pairs.py · pairs"]
        direction LR
        a1["compare_groups()"]:::fn --> a2["pair_tasks()<br/>every pair, earlier pose is the model"]:::fn
    end
    d1(["list of pairs"]):::data

    subgraph C2["analysis/pairs.py · run"]
        direction LR
        b1["run_pair_tasks()<br/>serial or --workers processes"]:::fn --> b2["pair_job()<br/>a crash becomes an error row"]:::fn --> b3["compare_pair()"]:::fn
    end
    d2(["one pair"]):::data

    subgraph C3["compare_pair() · radius"]
        direction LR
        c1["compare_ligands()<br/>at 4 Å"]:::ost --> c2["binding_site_failure()<br/>retry only for this"]:::fn --> c3["compare_ligands()<br/>at 6 Å, then 8 Å"]:::ost
    end
    d3(["every OST attempt"]):::data

    subgraph C4["analysis/pairs.py · row"]
        direction LR
        e1["summarize()<br/>values only if one ligand pair assigned,<br/>lDDT-PLI from the 4 Å call"]:::fn --> e2[/"pairs.csv · ligand_details.jsonl"/]:::file
    end
    d4(["one row per pair"]):::data

    start --> C1 --> d1 --> C2 --> d2 --> C3 --> d3 --> C4 --> d4

    class C1,C2,C3,C4 stage
```

## 6 · Scoring

```mermaid
%%{init: {"flowchart": {"curve": "linear", "nodeSpacing": 18, "rankSpacing": 28, "subGraphTitleMargin": {"top": 6, "bottom": 10}}}}%%
flowchart TB
    classDef fn fill:#e7eff7,stroke:#5b7fa3,color:#13293d
    classDef data fill:#ffffff,stroke:#6f8f4e,color:#1f2d12,stroke-width:2px
    classDef stage fill:#f7f9fb,stroke:#c9d2dc,color:#1b2430

    start(["one row per pair"]):::data

    subgraph K1["analysis/core.py"]
        direction LR
        a1["score_groups()"]:::fn --> a2["build_matrix()<br/>every pair usable, or no score"]:::fn --> a3["core_score(matrix, T)<br/>for each T"]:::fn
    end
    d1(["CORE per complex and T,<br/>or the reason it's missing"]):::data

    start --> K1 --> d1

    class K1 stage
```

## 7 · Reference evaluation

Only for ensembles with a reference. Nothing here feeds back into CORE.

```mermaid
%%{init: {"flowchart": {"curve": "linear", "nodeSpacing": 18, "rankSpacing": 28, "subGraphTitleMargin": {"top": 6, "bottom": 10}}}}%%
flowchart TB
    classDef fn fill:#e7eff7,stroke:#5b7fa3,color:#13293d
    classDef ost fill:#fbedd5,stroke:#b07a1e,color:#3d2a08
    classDef data fill:#ffffff,stroke:#6f8f4e,color:#1f2d12,stroke-width:2px
    classDef file fill:#f3f3f3,stroke:#9a9a9a,color:#444,stroke-dasharray:3 3
    classDef stage fill:#f7f9fb,stroke:#c9d2dc,color:#1b2430

    start(["prepared poses"]):::data

    subgraph E1["analysis/reference.py · prepare"]
        direction LR
        a1["evaluate_groups()"]:::fn --> a2["prepare()<br/>same as the predictions"]:::fn
    end
    d1(["prepared reference"]):::data

    subgraph E2["analysis/pairs.py"]
        direction LR
        b1["run_pair_tasks()"]:::fn --> b2["compare_pair()<br/>each pose against the reference"]:::ost
    end
    d2(["RMSD and lDDT-PLI of each pose to the reference"]):::data

    subgraph E3["analysis/reference.py · label"]
        direction LR
        c1["write_json()"]:::fn --> c2[/"reference_details.json"/]:::file --> c3["label_reference_results()<br/>correct means RMSD < T and lDDT-PLI > 0.8"]:::fn
    end
    d3(["one label row per pose"]):::data

    start --> E1 --> d1 --> E2 --> d2 --> E3 --> d3

    class E1,E2,E3 stage
```

## 8 · Report

```mermaid
%%{init: {"flowchart": {"curve": "linear", "nodeSpacing": 18, "rankSpacing": 28, "subGraphTitleMargin": {"top": 6, "bottom": 10}}}}%%
flowchart TB
    classDef fn fill:#e7eff7,stroke:#5b7fa3,color:#13293d
    classDef data fill:#ffffff,stroke:#6f8f4e,color:#1f2d12,stroke-width:2px
    classDef file fill:#f3f3f3,stroke:#9a9a9a,color:#444,stroke-dasharray:3 3
    classDef stage fill:#f7f9fb,stroke:#c9d2dc,color:#1b2430

    start(["prepared poses, PoseBusters records,<br/>CORE scores, reference labels"]):::data

    subgraph D1["outputs/report.py · setup"]
        direction LR
        a1["write_report()"]:::fn --> a2["pb_valid by<br/>ensemble and pose"]:::fn
    end
    d1(["validity lookup"]):::data

    subgraph D2["outputs/report.py · complexes"]
        direction LR
        b1["core_table()"]:::fn --> b2["top_pose()<br/>valid pose with highest confidence"]:::fn --> b3["above each cutoff<br/>CORE ≥ cutoff"]:::fn --> b4[/"core_scores.csv"/]:::file
    end
    d2(["one row per complex and T"]):::data

    subgraph D3["outputs/report.py · poses"]
        direction LR
        c1["pose_table()<br/>confidence, validity and reference<br/>for every pose"]:::fn --> c2[/"poses.csv"/]:::file
    end

    subgraph D4["outputs/plot.py"]
        direction LR
        e1["plot_core_scores()<br/>with --plot"]:::fn --> e2["a dot per complex,<br/>ordered by CORE"]:::fn --> e3[/"core_scores.png"/]:::file
    end

    start --> D1 --> d1 --> D2 --> d2 --> D3 --> D4

    class D1,D2,D3,D4 stage
```

## Every check, and what happens when it fails

Input problems stop the run before anything is computed. After that, a failure is
recorded as an `error` on the pose or ensemble, later stages skip it, and the rest of the
run carries on.

| Where | Check | If it fails |
| --- | --- | --- |
| `cli.py` | Thresholds are finite and positive, cutoffs are between 0 and 1, workers is at least 1 | Usage error |
| `inputs/inputs.py` | The CSV has the columns `ensemble`, `smiles` and `poses`, filled on every row, and no unknown columns | Exit 2 |
| `inputs/inputs.py` | Every row of an ensemble gives the same SMILES and reference, and no structure is listed twice | Exit 2 |
| `inputs/inputs.py` | A row with a `confidence` or `ligand_sdf` names one structure | Exit 2 |
| `inputs/inputs.py` | `confidence` is given for every pose of an ensemble or for none | Exit 2 |
| `inputs/inputs.py` | Every structure, SDF and reference file exists, and every wildcard matches a file | Exit 2 |
| `inputs/inputs.py` | Every ensemble's SMILES parses and has heavy atoms | Exit 2 |
| `inputs/inputs.py` | Every ensemble has at least four poses | Exit 2 |
| `pipeline.py` | The output folder is new or empty | Exit 2 |
| `inputs/structures.py` | One structural model per file, with atom records | Pose error |
| `inputs/structures.py` | Requested chains exist, the ligand is present, selected copies share one composition, no duplicate atom names | Pose error |
| `inputs/structures.py` | Ligand coordinates are finite, native bond orders are supported | Pose error |
| `inputs/structures.py` | A native AF3 SMILES matches the supplied SMILES | Pose error |
| `inputs/structures.py` | The PoseBusters molecule matches the SMILES graph, sanitizes, and keeps identity, coordinates and atom order | `pb_error`, PoseBusters unknown, pose still scored |
| `inputs/structures.py` | Written SDFs reread with the same coordinates. Supplied SDFs hold one molecule with coordinates that matches the SMILES | Pose error |
| `inputs/structures.py` | The filtered structure rereads with every kept record unchanged | Pose error |
| `analysis/ost.py` | OST can read the prepared structure | Pose error |
| `analysis/posebusters.py` | Every dock check passes on every copy. Unknown is never a pass | `pb_valid` False or unknown |
| `inputs/sites.py` | No pose in the ensemble failed preparation | Ensemble blocked |
| `inputs/sites.py` | Every pose has the same number of ligand copies | Ensemble blocked |
| `analysis/ost.py`, `inputs/sites.py` | OST superposes each pose onto the anchor, and site assignment converges | Ensemble blocked |
| `analysis/pairs.py` | OST succeeds and assigns exactly one ligand pair | Pair error |
| `analysis/pairs.py` | Only binding-site failures are retried at a larger radius | Other failures stay failed |
| `analysis/core.py` | Every pair has a usable RMSD, and there are at least four poses | No CORE, reason saved |
| `analysis/reference.py` | The reference and the prediction both prepared, and the RMSD and lDDT-PLI are known | Label unknown (blank) |
| `outputs/report.py` | Confidences given, PoseBusters run, and at least one valid pose | No top pose, `note` says why |
| `outputs/plot.py` | matplotlib is installed, only with `--plot` | ImportError, after the files are written |
| `pipeline.py` | Nothing above left anything incomplete | Exit 1 |

## Every function

### `cli.py`
| Function | Does |
| --- | --- |
| `main(argv)` | Parses the command, calls `run()`, prints counts. Returns 0, 1 if anything is incomplete, 2 for bad input. |
| `parser()` | Defines every flag. |
| `positive`, `fraction`, `worker_count` | Argument types for thresholds, cutoffs and workers. |

### `inputs/inputs.py`
| Function | Does |
| --- | --- |
| `load_inputs(inputs, ...)` | Turns one CSV, or a list of structures, into ensembles of at least four poses. |
| `read_csv(path)` | Reads the CSV and checks its columns and required cells, with one error message. |
| `make_group(name, rows, base)` | One ensemble from its rows. Checks one SMILES and reference, that the SMILES parses, no structure twice and all-or-none confidences. Sorts poses by path. |
| `make_poses(name, row, base)` | The poses one row names, expanding a wildcard. Pose IDs are paths relative to the CSV. |
| `resolve(base, path)` | Absolute path that must exist. |

### `inputs/make_csv.py`, the `core-score-csv` command <a id="core-score-csv"></a>
| Function | Does |
| --- | --- |
| `main(argv)` | Parses the command, writes the CSV, prints what it found. Returns 2 if anything is wrong, listing every problem. |
| `build_rows(predictions, smiles_csv, output, references)` | One row per pose of every ensemble folder, with confidence and reference where found. |
| `find_poses(folder)` | Structure files in the folder tree, each set of coordinates once, the deepest copy kept. |
| `confidence_file(structure)` | The model and confidence file written beside a structure, by `CONFIDENCE_FILES`. |
| `ligand_iptm(structure, smiles, rows, scores)` | Mean chain-pair ipTM between the ligand's chains and the receptor chains. |
| `af3_or_protenix`, `boltz`, `chai`, `openfold3` | Read each model's chain-pair ipTM. |
| `read_smiles(path)`, `relative(path, output)` | The ensemble,smiles CSV, and paths relative to the output CSV. |

### `pipeline.py`
| Function | Does |
| --- | --- |
| `run(inputs, out, ...)` | The Python entry point: `load_inputs()` then `execute()`, exactly as the command does. |
| `execute(groups, output, ...)` | Runs every stage in order and writes `details/run.json` and `details/status.json`. |
| `summarize_run(...)` | Counts incomplete scores, failed pairs, blocked ensembles and unknown reference labels. |

### `outputs/files.py`
| Function | Does |
| --- | --- |
| `clean(value)` | Makes values JSON-safe. NaN and pandas NA become None. |
| `write_json`, `write_jsonl`, `write_csv` | Write output files. |
| `fingerprint(value)` | Stable hash, used for folder names. |
| `diagnostic(stage, exc, ...)` | The failure record: stage, error, traceback, evidence. |
| `pose_directory(output, pose)` | `prepared/<hash>/` for a pose. |

### `analysis/ost.py`
| Function | Does |
| --- | --- |
| `execute(arguments, timeout)` | Runs `ost ... -o file` and returns its JSON and the command, return code and output. |
| `compare_ligands(model, reference, timeout, radius, lddt_pli)` | `compare-ligand-structures` for one pair, with lDDT-PLI only for references. |
| `reads_ligand(prepared, timeout)` | Does OST read every ligand copy from the structure, connected as in the package's SDFs? |
| `superpose(model_path, reference_path, timeout)` | `compare-structures --rigid-scores` transform. |
| `OSTError` | Error carrying the OST evidence. |

### `inputs/prepare.py`
| Function | Does |
| --- | --- |
| `prepare_groups(groups, output)` | Prepares every pose, retries poses whose chemistry failed once the peer cache is full, writes `prepared_manifest.json`. |
| `prepare(pose, destination, ...)` | Prepares one pose and asks OST whether it sees the ligand. Returns the pose with `prepared`, or with `error`. Also used for references. |
| `prepare_pose(pose, destination, cache)` | Picks a route, then writes the receptor. Returns the prepared record. |
| `from_structure(...)` | Selects ligand residues, builds OST and PoseBusters SDFs per copy, writes the filtered structure. |
| `from_original_sdfs(...)` | Uses supplied SDFs, finds which structure residues they replace and writes the structure. |
| `same_atoms(rows, mol)` | Whether a residue holds the molecule's atoms at the same coordinates. |

### `analysis/posebusters.py`
| Function | Does |
| --- | --- |
| `run_posebusters(groups, output)` | Checks every pose and writes `posebusters.json`. Removes nothing. |
| `dock_checks(buster)` | Names of the dock config's binary checks. |
| `normalized(name)` | Check names in lower case with underscores. |
| `aggregate(values)` | False if any False, None if any unknown, else True. |
| `check_pose(buster, checks, pose)` | Validity of one pose, and each check's result, across its copies. |
| `check_copy(buster, checks, copy, receptor)` | Runs PoseBusters on one copy. |

### `inputs/sites.py`
| Function | Does |
| --- | --- |
| `prepare_site_groups(groups, output)` | Blocks ensembles with failed poses or mismatched copy counts, passes one-copy ensembles through, splits the rest and names each site `<ensemble>_site_<n>`. |
| `split_sites(group, copies, output)` | One ensemble per site, and the report for `sites.json`. |
| `align_copy_centres(poses)` | Copy centres in the anchor pose's frame. |
| `assign_sites(points, k)` | Site label per copy, one copy per site per pose. |
| `site_pose(pose, copy, destination)` | The pose cut to one copy. |

### `analysis/pairs.py`
| Function | Does |
| --- | --- |
| `compare_groups(groups, output, ...)` | Runs every pair and writes `pairs.csv` and `ligand_details.jsonl`. |
| `pair_tasks(groups)` | Every pair within each site ensemble. Skips groups with an error. |
| `run_pair_tasks(tasks, workers, ...)` | Runs tasks serially or in parallel, in order. |
| `pair_job(task, ...)` | One pair, with identity fields added and any crash turned into an error row. |
| `compare_pair(pose_i, pose_j, ...)` | OST at 4 Å, then 6 and 8 Å after a binding-site failure. Returns the row and its details. |
| `summarize(result)` | Metrics only when exactly one ligand pair was assigned. |
| `usable(score)`, `binding_site_failure(result)` | A finite non-negative score, and a binding-site reason. |

### `analysis/core.py`
| Function | Does |
| --- | --- |
| `score_groups(groups, rows, thresholds)` | CORE per site ensemble and T, or the reason it's missing. |
| `build_matrix(pose_ids, rows)` | Pairwise RMSD matrix. Raises unless every pair is usable. |
| `core_score(matrix, threshold, neighbours)` | 1 − mean of min(RMSD, T) / T over each pose's three closest other poses. Raises with fewer than four poses. |

### `analysis/reference.py`
| Function | Does |
| --- | --- |
| `evaluate_groups(groups, output, thresholds, ...)` | For ensembles with a reference, prepares it, compares every pose with it, writes `reference_details.json`. Returns one label row per pose, in path order, for `poses.csv`. |
| `label_reference_results(records, thresholds)` | Correct means RMSD < T and lDDT-PLI > 0.8. Failed, ambiguous or missing either score is unknown. |
| `correct_column(threshold)` | The `poses.csv` column name, such as `correct_2A`. |

### `outputs/report.py`
| Function | Does |
| --- | --- |
| `write_report(groups, scores, pb_records, labels, output, cutoffs)` | Writes `core_scores.csv` and `poses.csv`. |
| `top_pose(group, validity)` | The PoseBusters-valid pose with the highest confidence, or the reason there is none. |
| `core_table(...)` | One row per complex and T: CORE, above each cutoff, pose counts, the top pose, and reference results. |
| `pose_table(groups, pb_records, labels, ranked)` | One row per pose: confidence, PoseBusters validity, and reference RMSD, lDDT-PLI and correctness. |
| `message(error)` | The one-line message of a diagnostic record. |

### `outputs/plot.py`
| Function | Does |
| --- | --- |
| `plot_core_scores(rows, output, cutoffs)` | `core_scores.png` per T: a dot per complex, ordered by CORE Score, with a line at each cutoff. Imports matplotlib only when called. |
| `style(matplotlib)` | The figure's type and hairline axes, using the first of `FONTS` installed. |
| `draw(axes, scores, threshold, cutoffs)` | The dots, the cutoff lines and their labels. |
| `dot_size(count)` | Smaller dots for bigger benchmarks. |

### `inputs/structures.py`
| Function | Does |
| --- | --- |
| `read_document(path)` | Reads CIF, or converts PDB, requiring one model. |
| `merge_identical_entities(doc)` | Gives polymer chains with the same sequence one entity, so OST can swap identical chains. Only entity IDs change. |
| `atom_rows(doc)` | Atom records as dictionaries. |
| `residue_key(row)` | (chain, number, insertion code, residue name). |
| `ligand_residues(doc)` | Nonpolymer residues and their atoms. |
| `filter_records(doc, selected, destination)` | Writes the structure without unselected ligands and checks every kept record is unchanged. |
| `write_receptor(source, ligands, destination)` | Writes the PoseBusters receptor without the selected ligands. |


| Function | Does |
| --- | --- |
| `template_mol(smiles)` | Heavy-atom molecule from SMILES. Raises if missing, invalid or without heavy atoms. |
| `heavy_atoms`, `atom_name`, `composition` | Small helpers on atom records. |
| `choose_ligands(doc, pose)` | Ligand residues picked by chain and SMILES composition. Raises on absence, mixed compositions or duplicate atom names. |
| `coordinate_mol(rows)` | Molecule with atoms and coordinates only. |
| `canonical_smiles(mol)` | Canonical SMILES without stereo, for comparing chemistry. |
| `native_graph(doc, rows)` | Bonds from the structure's bond table or AF3's SMILES and atom names, if present. |


| Function | Does |
| --- | --- |
| `make_ligand(doc, rows, smiles, cache, ensemble)` | The OST molecule and the PoseBusters molecule for one copy, with provenance, or the PoseBusters error. |
| `matched_template(template, observed)` | Puts the SMILES template's chemistry onto the observed atoms after an exact graph match. |
| `generic_graph(mol)` | Element-and-connectivity-only copy, for matching. |
| `write_molecule(mol, path)` | Writes an SDF and rereads it to check coordinates. |
| `read_original(path, smiles)` | Reads a supplied SDF and checks it matches the SMILES. |
