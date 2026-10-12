from pathlib import Path
import pytest
from core_score.cli import parser, main


@pytest.mark.parametrize("values", [
    ["p0.cif", "p1.pdb", "--out", "results"],
    ["--out", "results", "p0.cif", "p1.pdb"],
])
def test_direct_structures(values):
    args = parser().parse_args(values)
    assert args.inputs == [Path("p0.cif"), Path("p1.pdb")]
    assert not args.skip_posebusters and args.workers == 1


@pytest.mark.parametrize("name", ["poses.json", "poses.JSON"])
def test_json_input(name):
    args = parser().parse_args([name, "--out", "results", "--plots", "--skip-posebusters"])
    assert args.inputs == [Path(name)]
    assert args.plots and args.skip_posebusters


@pytest.mark.parametrize("words", [[], ["--help"]])
def test_help(words, capsys):
    with pytest.raises(SystemExit) as exc:
        main(words)
    assert exc.value.code == 0
    assert "inputs.json" in capsys.readouterr().out


def test_python_run_loads_inputs_then_executes(tmp_path):
    from unittest.mock import patch
    import core_score

    with patch("core_score.pipeline.load_inputs", return_value=["groups"]) as load, \
            patch("core_score.pipeline.execute", return_value="summary") as execute:
        assert core_score.run("inputs.json", tmp_path, smiles="CCO", workers=4) == "summary"
    assert load.call_args.args[0] == ["inputs.json"]
    assert load.call_args.kwargs["smiles"] == "CCO"
    assert execute.call_args.args == (["groups"], tmp_path)
    assert execute.call_args.kwargs == dict(evaluate_reference=False, workers=4)
