# Code map

How information moves through `core_score`, file by file and then function by function.
Arrows are calls. Arrow labels name the data passed along. Every file writes its outputs
through `results.py` (`write_json`, `write_csv`, `write_jsonl`, `diagnostic`), which is
left off the diagrams to keep them readable.

Read the code in this order: `cli.py`, `pipeline.execute()`, then each file in `stages/`
in the order `execute()` calls them, then the helpers each one uses in `stages/utils/`.

## 1. Files

```mermaid
flowchart TD
    main["__main__.py"] --> cli["cli.py"]
    cli -->|"paths and flags"| inputs["inputs.py"]
    inputs -->|"SMILES check"| chem
    inputs -->|"ensembles"| cli
    cli -->|"ensembles and settings"| pipeline["pipeline.py"]

    pipeline --> s_prep["stages/preparation.py"]
    pipeline --> s_pb["stages/posebusters.py"]
    pipeline --> s_sites["stages/site_assignment.py"]
    pipeline --> s_cmp["stages/comparison.py"]
    pipeline --> s_score["stages/scoring.py"]
    pipeline --> s_eval["stages/evaluation.py"]
    pipeline --> s_dec["stages/decisions.py"]
    pipeline --> s_plot["stages/plotting.py"]

    s_prep --> prepare["utils/preparation/prepare.py"]
    s_prep --> ost["utils/ost.py"]
    prepare --> chem["utils/preparation/chemistry.py"]
    prepare --> cif["utils/preparation/cif.py"]
    prepare --> sdf["utils/preparation/sdf.py"]
    sdf --> chem
    chem --> cif

    s_pb --> pbu["utils/posebusters_utils.py"]
    pbu --> pblib[("posebusters library")]

    s_sites --> sites["utils/sites.py"]
    sites --> ost
    sites --> cif
    sites -->|"pose_directory"| prepare

    s_cmp --> cmpu["utils/comparison_utils.py"]
    cmpu --> ost
    cmpu --> radius["utils/radius.py"]

    s_score --> score["utils/score_utils.py"]

    s_eval -->|"prepare()"| s_prep
    s_eval -->|"run_pair_tasks()"| cmpu
    s_eval --> evalu["utils/evaluation_utils.py"]

    s_dec --> decu["utils/decision_utils.py"]

    s_plot -->|"reads saved CSV and JSON"| files[("output folder")]
    ost --> ostexe[("ost executable on PATH")]
```

## 2. Data passed between stages

Every record is a plain dictionary. Each stage adds to it or derives a new list from it.

```mermaid
flowchart LR
    E(["ensembles<br/>name + poses"]) -->|prepare_groups| P(["prepared groups<br/>each pose has prepared or error"])
    P -->|run_posebusters| B(["PoseBusters records<br/>pb_valid per pose"])
    P -->|prepare_site_groups| S(["site groups<br/>one ensemble per site, or error"])
    S -->|compare_groups| R(["pair rows<br/>RMSD per pair"])
    S --> X
    R -->|score_groups| X(["scores<br/>CORE per site and T"])
    P -->|evaluate_groups| L(["reference labels<br/>RMSD < T per pose"])
    X --> L
    P --> D
    B --> D
    X --> D
    L --> D(["decisions<br/>accepted at each cutoff, top pose"])
```

The record for one pose, as it grows:

| After | Adds | Key fields |
| --- | --- | --- |
| `load_inputs` | the pose | `ensemble`, `pose_id`, `path`, `smiles`, `ligand_chains`, `seed`, `sample`, optional `ligand_sdfs` and `confidence` |
| `prepare_groups` | `prepared` or `error` | `prepared.path`, `ligand_sdfs`, `copies`, `copy_records` (one per ligand copy: `residue`, `ost_sdf`, `pb_sdf`, `pb_error`, `provenance`, `centre`), `route`, `receptor` |
| `prepare_site_groups` | `site` on the group | each pose cut to one copy per site |
| `compare_groups` | one row per pair | `pose_i`, `pose_j`, `n_assigned`, `rmsd`, `bb_rmsd`, `lddt_lp`, `lddt_pli`, radius fields, `error` |
| `score_groups` | one row per site and T | `pose_count`, `threshold_angstrom`, `core_score`, `error` |
| `decide_groups` | one row per site, T and cutoff | `core_cutoff`, `likely_correct`, `top_pose`, `top_pose_correct`, `note` |

## 3. Functions, stage by stage

### 3.1 Command and inputs

```mermaid
flowchart TD
    main["cli.main()"] --> parser["cli.parser()"]
    parser -.->|"type checks"| types["positive() / fraction() / worker_count()"]
    main -->|"paths, --smiles, --ligand-chain, --reference"| load["inputs.load_inputs()"]
    load -->|"one .json path"| read["read_manifest()"]
    read -->|"ensembles list, shape checked once"| load
    load -->|"structure paths"| direct["one ensemble named 'ensemble'"]
    load --> group["make_group()"]
    group --> pose["make_pose()"]
    pose --> resolve["resolve()<br/>path must exist"]
    group -->|"unique pose IDs, confidence all or none, sort by seed, sample, pose_id"| load
    load -->|"only with --evaluate-reference"| ref["make_reference()"]
    ref --> pose
    load -->|"every pose's SMILES"| tm["chemistry.template_mol()"]
    load -->|"ensembles"| main
    main -->|"ensembles and settings"| exec["pipeline.execute()"]
```

### 3.2 The run

```mermaid
flowchart TD
    exec["pipeline.execute()"] -->|"output folder must be empty"| run[/"run.json"/]
    exec --> prep["preparation.prepare_groups()"]
    exec -->|"unless --skip-posebusters"| pb["posebusters.run_posebusters()"]
    exec --> sites["site_assignment.prepare_site_groups()"]
    exec --> cmp["comparison.compare_groups()"]
    exec --> score["scoring.score_groups()"]
    exec -->|"only with --evaluate-reference"| ev["evaluation.evaluate_groups()"]
    exec --> dec["decisions.decide_groups()"]
    exec --> sum["summarize_run()"]
    sum --> status[/"status.json and exit code"/]
    exec -->|"only with --plots"| plot["plotting.plot_results()"]
```

### 3.3 Preparation

```mermaid
flowchart TD
    pg["preparation.prepare_groups()"] -->|"each pose, shared peer cache"| pr["preparation.prepare()"]
    pg -.->|"poses with a pb_error, once more"| pr
    pg --> pm[/"prepared_manifest.json"/]
    pr --> pd["prepare.pose_directory()"]
    pr --> pp["prepare.prepare_pose()"]
    pr -->|"structure route only"| det["ost.detects_ligands()"]
    det --> cl["ost.compare_ligands()"]
    cl --> ex["ost.execute()"]
    det -->|"OST cannot see the ligand"| switch["use the ost.sdf files<br/>route explicit_connectivity_sdf"]
    pr -->|"any exception"| err["pose error = diagnostic()"]

    pp --> rd["cif.read_document()"]
    pp -->|"no ligand_sdfs"| fs["prepare.from_structure()"]
    pp -->|"ligand_sdfs given"| fo["prepare.from_original_sdfs()"]
    pp --> wr["cif.write_receptor()"]

    fs --> cho["chemistry.choose_ligands()"]
    cho --> lr["cif.ligand_residues()"]
    lr --> ar["cif.atom_rows()"]
    lr --> rk["cif.residue_key()"]
    cho --> comp["chemistry.composition()"]
    fs --> ml["sdf.make_ligand()"]
    fs --> wm["sdf.write_molecule()<br/>rereads to check coordinates"]
    fs --> fr["cif.filter_records()<br/>rereads to check records"]

    ml --> cm["chemistry.coordinate_mol()"]
    ml --> ng["chemistry.native_graph()<br/>CIF bonds or AF3 SMILES"]
    ml -->|"no native graph"| cache["peer cache<br/>same ensemble, SMILES, atom names"]
    ml -->|"nothing else"| conn["bonds from coordinates"]
    ml --> mt["sdf.matched_template()"]
    mt --> gg["sdf.generic_graph()"]
    ml --> idn["chemistry.identity()"]

    fo --> ro["sdf.read_original()"]
    ro --> idn
    fo --> sa["prepare.same_atoms()"]
```

### 3.4 PoseBusters

```mermaid
flowchart TD
    rp["posebusters.run_posebusters()"] --> bus["PoseBusters(config='dock')"]
    rp --> dc["dock_checks()"]
    dc --> nm["normalized()"]
    rp -->|"every pose"| cp["check_pose()"]
    cp -->|"pose has an error"| unk["pb_valid None"]
    cp -->|"every copy"| cc["check_copy()"]
    cc -->|"copy has a pb_error"| unk
    cc --> bust["buster.bust(pb_sdf, receptor)"]
    cc --> ag["aggregate()<br/>False beats None beats True"]
    cp --> ag
    rp --> sr["summary_row()"]
    rp --> out[/"posebusters.json and .csv"/]
    rp -->|"records, used later for the top pose"| dec["decisions"]
```

### 3.5 Site assignment

```mermaid
flowchart TD
    psg["site_assignment.prepare_site_groups()"] -->|"any pose error, or copy counts differ"| blocked["group error, site unresolved"]
    psg -->|"one copy"| one["site 1, unchanged"]
    psg -->|"k copies"| ss["sites.split_sites()"]
    ss --> acc["sites.align_copy_centres()"]
    acc -->|"each pose onto the anchor"| sp["ost.superpose()"]
    sp --> ex["ost.execute()"]
    ss -->|"copy centres in the anchor frame"| as["sites.assign_sites()<br/>k-means seed, Hungarian, repeat"]
    ss -->|"each site and pose"| spose["sites.site_pose()"]
    spose --> fr["cif.filter_records()"]
    ss --> pd["prepare.pose_directory()"]
    psg --> out[/"sites.json"/]
```

### 3.6 Comparison

```mermaid
flowchart TD
    cg["comparison.compare_groups()"] --> pt["comparison_utils.pair_tasks()<br/>every pair, skips groups with an error"]
    cg --> rt["run_pair_tasks()<br/>1 process or a pool of --workers"]
    rt -->|"each pair"| pj["pair_job()<br/>any crash becomes an error row"]
    pj --> cpair["compare_pair()"]
    cpair --> swr["radius.score_with_radius_fallback()"]
    swr -->|"4, then 6 and 8 A on a binding-site failure"| evl["evaluate(radius)"]
    evl --> cl["ost.compare_ligands()"]
    cl --> ex["ost.execute()"]
    swr --> bsf["radius.binding_site_failure()"]
    swr --> us["radius.usable()"]
    cpair -->|"4 A result, RMSD from the radius used"| sm["summarize()<br/>values only if one pair assigned"]
    cg --> out[/"pairs.csv and ligand_details.jsonl"/]
```

### 3.7 Scoring

```mermaid
flowchart TD
    sg["scoring.score_groups()"] -->|"group error"| blank["core_score blank, error kept"]
    sg -->|"rows for this ensemble and site"| bm["score_utils.build_matrix()<br/>every pair or ValueError"]
    bm -->|"matrix"| cs["score_utils.core_score(matrix, T)<br/>for each T"]
    sg --> out[/"scores.csv and .json"/]
```

### 3.8 Reference evaluation (with `--evaluate-reference`)

```mermaid
flowchart TD
    eg["evaluation.evaluate_groups()"] -->|"the reference"| pr["preparation.prepare()"]
    eg -->|"each pose against the reference"| rt["comparison_utils.run_pair_tasks()"]
    eg --> lab["evaluation_utils.label_reference_results()<br/>success = RMSD < T"]
    eg --> sref["evaluation_utils.summarize_reference()"]
    eg --> out[/"reference_evaluation.csv, reference_summary.csv, reference_details.json"/]
```

### 3.9 Decisions

```mermaid
flowchart TD
    dg["decisions.decide_groups()"] -->|"pb_valid by ensemble and pose"| dd["decision_utils.decide()"]
    dd --> tp["top_pose()<br/>valid pose with highest confidence"]
    dd -->|"CORE >= cutoff"| lc["likely_correct per cutoff"]
    dg --> scs["summarize_cutoffs()<br/>coverage and accuracy"]
    dg --> out[/"decisions.csv and cutoff_summary.csv"/]
```

### 3.10 Plots (with `--plots`)

```mermaid
flowchart TD
    pl["plotting.plot_results()"] --> rc["read_csv()"]
    rc -->|"pairs.csv, reference_summary.csv"| pl
    pl -->|"scores.json"| pl
    pl --> out[/"figures/ heatmaps, reference plot, index.json"/]
```

## 4. Every function

### `cli.py`
| Function | Does |
| --- | --- |
| `main(argv)` | Parses the command, loads inputs, runs `execute()`, prints counts. Returns 0, 1 if anything is incomplete, 2 for bad input. |
| `parser()` | Defines every flag. |
| `positive`, `fraction`, `worker_count` | Argument types for thresholds, cutoffs and workers. |

### `inputs.py`
| Function | Does |
| --- | --- |
| `load_inputs(inputs, ...)` | Turns one JSON path or a list of structures into ensembles. Checks every SMILES. Loads references only when evaluating. |
| `read_manifest(path)` | Reads the JSON and checks its shape once, with one error message. |
| `make_group(ensemble, ...)` | Builds an ensemble's poses, checks unique pose IDs and all-or-none confidences, sorts by `(seed, sample, pose_id)`. |
| `make_pose(record, ...)` | One pose record with settings filled from the ensemble and command line. |
| `make_reference(reference, ...)` | The reference as a pose record with `pose_id` "reference". |
| `resolve(base, path)` | Absolute path that must exist. |

### `pipeline.py`
| Function | Does |
| --- | --- |
| `execute(groups, output, ...)` | Runs every stage in order and writes `run.json` and `status.json`. |
| `summarize_run(...)` | Counts incomplete scores, failed pairs, blocked ensembles and unknown labels. |

### `results.py`
| Function | Does |
| --- | --- |
| `clean(value)` | Makes values JSON-safe. NaN and pandas NA become None. |
| `write_json`, `write_jsonl`, `write_csv` | Write output files. |
| `fingerprint(value)` | Stable hash, used for folder and figure names. |
| `diagnostic(stage, exc, ...)` | The failure record: stage, error, traceback, evidence. |

### `stages/preparation.py`
| Function | Does |
| --- | --- |
| `prepare_groups(groups, output)` | Prepares every pose, retries poses whose chemistry failed once the peer cache is full, writes `prepared_manifest.json`. |
| `prepare(pose, destination, ...)` | Prepares one pose and asks OST whether it sees the ligand. Returns the pose with `prepared`, or with `error`. Also used for references. |

### `stages/utils/preparation/prepare.py`
| Function | Does |
| --- | --- |
| `pose_directory(output, pose)` | `prepared/<hash>/` for a pose. |
| `prepare_pose(pose, destination, cache)` | Picks a route, then writes the receptor. Returns the prepared record. |
| `from_structure(...)` | Selects ligand residues, builds OST and PoseBusters SDFs per copy, writes the filtered structure. |
| `from_original_sdfs(...)` | Uses supplied SDFs and finds which structure residues they replace. |
| `same_atoms(rows, mol)` | Whether a residue holds the molecule's atoms at the same coordinates. |

### `stages/utils/preparation/chemistry.py`
| Function | Does |
| --- | --- |
| `template_mol(smiles)` | Heavy-atom molecule from SMILES. Raises if missing or invalid. |
| `heavy_atoms`, `atom_name`, `composition` | Small helpers on atom records. |
| `choose_ligands(doc, pose)` | Ligand residues picked by chain, residue number and SMILES composition. Raises on ambiguity, absence, mixed compositions or duplicate atom names. |
| `coordinate_mol(rows)` | Molecule with atoms and coordinates only. |
| `identity(mol)` | Canonical SMILES without stereo, for comparing chemistry. |
| `native_graph(doc, rows)` | Bonds from the structure's bond table or AF3's SMILES and atom names, if present. |

### `stages/utils/preparation/sdf.py`
| Function | Does |
| --- | --- |
| `make_ligand(doc, rows, smiles, cache, ensemble)` | The OST molecule and the PoseBusters molecule for one copy, with provenance, or the PoseBusters error. |
| `matched_template(template, observed)` | Puts the SMILES template's chemistry onto the observed atoms after an exact graph match. |
| `generic_graph(mol)` | Element-and-connectivity-only copy, for matching. |
| `write_molecule(mol, path)` | Writes an SDF and rereads it to check coordinates. |
| `read_original(path, smiles)` | Reads a supplied SDF and checks it matches the SMILES. |

### `stages/utils/preparation/cif.py`
| Function | Does |
| --- | --- |
| `read_document(path)` | Reads CIF, or converts PDB, requiring one model. |
| `atom_rows(doc)` | Atom records as dictionaries. |
| `residue_key(row)` | (chain, number, insertion code, residue name). |
| `ligand_residues(doc)` | Nonpolymer residues and their atoms. |
| `filter_records(doc, selected, destination)` | Writes the structure without unselected ligands and checks every kept record is unchanged. |
| `write_receptor(source, ligands, destination)` | Writes the PoseBusters receptor without the selected ligands. |

### `stages/utils/ost.py`
| Function | Does |
| --- | --- |
| `execute(arguments, timeout)` | Runs `ost ... -o file` and returns its JSON and the command, return code and output. |
| `compare_ligands(model, reference, timeout, radius)` | `compare-ligand-structures` for one pair. |
| `detects_ligands(prepared, timeout)` | Self-comparison: does OST find the ligand by itself? |
| `superpose(model_path, reference_path, timeout)` | `compare-structures --rigid-scores` transform. |
| `OSTError` | Error carrying the OST evidence. |

### `stages/posebusters.py` and `stages/utils/posebusters_utils.py`
| Function | Does |
| --- | --- |
| `run_posebusters(groups, output)` | Checks every pose and writes the results. Removes nothing. |
| `dock_checks(buster)` | Names of the dock config's binary checks. |
| `check_pose(buster, checks, pose)` | Validity of one pose across its copies. |
| `check_copy(buster, checks, copy, receptor)` | Runs PoseBusters on one copy. |
| `aggregate(values)` | False if any False, None if any unknown, else True. |
| `normalized(name)`, `summary_row(record)` | Check name format, and the one-row-per-pose CSV. |

### `stages/site_assignment.py` and `stages/utils/sites.py`
| Function | Does |
| --- | --- |
| `prepare_site_groups(groups, output)` | Blocks ensembles with failed poses or mismatched copy counts, passes one-copy ensembles through, splits the rest. |
| `split_sites(group, copies, output)` | One ensemble per site, and the report for `sites.json`. |
| `align_copy_centres(poses)` | Copy centres in the anchor pose's frame. |
| `assign_sites(points, k)` | Site label per copy, one copy per site per pose. |
| `site_pose(pose, copy, destination)` | The pose cut to one copy. |

### `stages/comparison.py` and `stages/utils/comparison_utils.py`, `radius.py`
| Function | Does |
| --- | --- |
| `compare_groups(groups, output, ...)` | Runs every pair and writes `pairs.csv` and `ligand_details.jsonl`. |
| `pair_tasks(groups)` | Every pair within each site ensemble. Skips groups with an error. |
| `run_pair_tasks(tasks, workers, ...)` | Runs tasks serially or in parallel, in order. |
| `pair_job(task, ...)` | One pair, with identity fields added and any crash turned into an error row. |
| `compare_pair(pose_i, pose_j, ...)` | OST at 4 Å with radius recovery, then the row and its details. |
| `summarize(result)` | Metrics only when exactly one ligand pair was assigned. |
| `score_with_radius_fallback(evaluate, radii)` | Tries 4 Å, then larger radii after binding-site failures. |
| `binding_site_failure(result)`, `usable(score)` | Recognise a binding-site reason, and a finite non-negative score. |

### `stages/scoring.py` and `stages/utils/score_utils.py`
| Function | Does |
| --- | --- |
| `score_groups(groups, rows, output, thresholds)` | CORE per site ensemble and T, or the reason it's missing. |
| `build_matrix(pose_ids, rows)` | Pairwise RMSD matrix. Raises unless every pair is usable. |
| `core_score(matrix, threshold)` | 1 − mean of min(closest other pose RMSD, T) / T. |

### `stages/evaluation.py` and `stages/utils/evaluation_utils.py`
| Function | Does |
| --- | --- |
| `evaluate_groups(groups, scores, output, ...)` | Prepares each reference, compares every pose with it, writes the reference files. |
| `label_reference_results(records, thresholds)` | Success means RMSD < T. Failed or ambiguous is unknown. |
| `summarize_reference(scores, labels)` | Each ensemble's CORE beside its fraction of poses correct. |

### `stages/decisions.py` and `stages/utils/decision_utils.py`
| Function | Does |
| --- | --- |
| `decide_groups(groups, scores, pb_records, labels, output, cutoffs)` | Writes `decisions.csv` and `cutoff_summary.csv`. |
| `top_pose(group, validity)` | The PoseBusters-valid pose with the highest confidence, or the reason there is none. |
| `decide(...)` | One row per site ensemble, T and cutoff: likely correct or not, the top pose, and whether it is correct. |
| `summarize_cutoffs(decisions)` | Per cutoff: how many accepted, coverage, and accuracy of the accepted top poses. |

### `stages/plotting.py`
| Function | Does |
| --- | --- |
| `plot_results(results, output)` | Heatmaps per site ensemble and the reference plot, from saved files only. |
| `read_csv(path)` | Rows of a saved CSV. |
