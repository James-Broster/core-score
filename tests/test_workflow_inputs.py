import pytest
from core_score.cli import parser
from core_score.inputs.inputs import load_inputs


def load(words):
    args = parser().parse_args(words)
    return load_inputs(args.inputs, smiles=args.smiles, ligand_chains=args.ligand_chains, reference=args.reference)


def structures(folder, names):
    folder.mkdir(parents=True, exist_ok=True)
    for name in names:
        (folder / name).write_text("structure")


def write_csv(path, text):
    path.write_text(text.strip() + "\n")
    return path


def test_csv_wildcard_row_gives_every_pose_sorted_with_ids_relative_to_the_csv(tmp_path):
    structures(tmp_path / "preds" / "7ly6", ["d.cif", "b.cif", "a.cif", "c.cif"])
    table = write_csv(tmp_path / "inputs.csv", "ensemble,smiles,poses\n7ly6,CCO,preds/7ly6/*.cif")
    [group] = load([str(table), "--out", "results"])
    assert group["ensemble"] == "7ly6" and "reference" not in group
    assert [p["pose_id"] for p in group["poses"]] == [f"preds/7ly6/{n}.cif" for n in "abcd"]
    assert all(p["smiles"] == "CCO" and p["ligand_chains"] == [] for p in group["poses"])


def test_csv_rows_per_pose_carry_confidence_and_ligand_chains(tmp_path):
    structures(tmp_path, [f"{n}.cif" for n in "abcd"])
    rows = "\n".join(f"g,CCO,{n}.cif,{c},B C" for n, c in zip("abcd", (0.1, 0.4, 0.3, 0.2)))
    table = write_csv(tmp_path / "inputs.csv", "ensemble,smiles,poses,confidence,ligand_chain\n" + rows)
    [group] = load([str(table), "--out", "results"])
    assert [p["confidence"] for p in group["poses"]] == [0.1, 0.4, 0.3, 0.2]
    assert group["poses"][0]["ligand_chains"] == ["B", "C"]


def test_pose_order_ignores_listing_order(tmp_path):
    structures(tmp_path, [f"{n}.cif" for n in "abcd"])
    orders = []
    for listing in ("abcd", "dbca"):
        rows = "\n".join(f"g,CCO,{n}.cif" for n in listing)
        table = write_csv(tmp_path / f"{listing}.csv", "ensemble,smiles,poses\n" + rows)
        orders.append([p["pose_id"] for p in load([str(table), "--out", "results"])[0]["poses"]])
    assert orders[0] == orders[1] == [f"{n}.cif" for n in "abcd"]


def test_reference_column_is_loaded_with_the_ensemble_smiles(tmp_path):
    structures(tmp_path, ["a.cif", "b.cif", "c.cif", "d.cif", "ref.cif"])
    table = write_csv(tmp_path / "inputs.csv", "ensemble,smiles,poses,reference\ng,CCO,[abcd].cif,ref.cif")
    [group] = load([str(table), "--out", "results"])
    assert group["reference"] == dict(ensemble="g", pose_id="reference", path=str(tmp_path / "ref.cif"), smiles="CCO")


@pytest.mark.parametrize("text, message", [
    ("ensemble,smiles\ng,CCO", "missing columns \\['poses'\\]"),
    ("ensemble,smiles,poses,confidance\ng,CCO,*.cif,1", "unknown columns \\['confidance'\\]"),
    ("ensemble,smiles,poses\ng,,*.cif", "lines \\[2\\]"),
    ("ensemble,smiles,poses,confidence\ng,CCO,*.cif,0.5", "on a row of its own"),
    ("ensemble,smiles,poses\ng,CCO,[ab].cif\ng,CCN,[cd].cif", "same smiles"),
    ("ensemble,smiles,poses\ng,CCO,*.cif\ng,CCO,a.cif", "more than once"),
    ("ensemble,smiles,poses\ng,CCO,[abc].cif", "at least 4 poses.*g \\(3\\)"),
    ("ensemble,smiles,poses\ng,CCO,missing/*.cif", "no structure file matches"),
    ("ensemble,smiles,poses\ng,not-a-smiles,*.cif", "^g: "),
])
def test_csv_mistakes_have_actionable_errors(tmp_path, text, message):
    structures(tmp_path, [f"{n}.cif" for n in "abcd"])
    with pytest.raises((ValueError, FileNotFoundError), match=message):
        load([str(write_csv(tmp_path / "inputs.csv", text)), "--out", "results"])


def test_confidence_is_given_for_every_pose_or_none(tmp_path):
    structures(tmp_path, [f"{n}.cif" for n in "abcd"])
    rows = "g,CCO,a.cif,0.8\ng,CCO,b.cif,\ng,CCO,c.cif,0.1\ng,CCO,d.cif,0.2"
    with pytest.raises(ValueError, match="confidence for every pose"):
        load([str(write_csv(tmp_path / "inputs.csv", "ensemble,smiles,poses,confidence\n" + rows)), "--out", "o"])


@pytest.mark.parametrize("smiles", [None, " ", "not-a-smiles"])
def test_missing_or_invalid_smiles_names_the_ensemble(tmp_path, smiles):
    structures(tmp_path, [f"{n}.cif" for n in "abcd"])
    words = [str(tmp_path / f"{n}.cif") for n in "abcd"] + ["--out", "results"]
    if smiles is not None:
        words += ["--smiles", smiles]
    with pytest.raises(ValueError, match="^ensemble: "):
        load(words)


def test_direct_inputs_and_reference_use_working_directory(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    structures(tmp_path, ["a.cif", "b.cif", "c.cif", "d.cif", "ref.cif"])
    [group] = load(["a.cif", "b.cif", "c.cif", "d.cif", "--smiles", "CCO", "--ligand-chain", "B",
                    "--reference", "ref.cif", "--out", "results"])
    assert [p["pose_id"] for p in group["poses"]] == ["a.cif", "b.cif", "c.cif", "d.cif"]
    assert group["poses"][0]["ligand_chains"] == ["B"]
    assert group["reference"]["path"] == str(tmp_path / "ref.cif")


def test_missing_structure_or_reference_is_rejected(tmp_path):
    structures(tmp_path, [f"{n}.cif" for n in "abc"])
    paths = [str(tmp_path / f"{n}.cif") for n in "abc"]
    with pytest.raises(FileNotFoundError):
        load(paths + [str(tmp_path / "missing.cif"), "--smiles", "CCO", "--out", "results"])
    structures(tmp_path, ["d.cif"])
    with pytest.raises(FileNotFoundError):
        load(paths + [str(tmp_path / "d.cif"), "--smiles", "CCO", "--reference", "nope.cif", "--out", "results"])
