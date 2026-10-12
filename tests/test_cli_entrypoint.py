from pathlib import Path
import pytest
from core_score.cli import parser, main


@pytest.mark.parametrize("values", [
    ["p0.cif", "p1.pdb", "--out", "results"],
    ["--out", "results", "p0.cif", "p1.pdb"],
])
def test_direct_structures(values):
    args = parser().parse_args(values)
    assert args.structures == [Path("p0.cif"), Path("p1.pdb")]
    assert not hasattr(args, "command")
    assert not args.skip_posebusters and args.workers == 1


@pytest.mark.parametrize("name", ["poses.json", "poses.JSON"])
def test_json_input(name):
    args = parser().parse_args([name, "--out", "results", "--plots"])
    assert args.structures == [Path(name)]
    assert args.plots


def test_explicit_manifest_and_pb_skip():
    args = parser().parse_args([
        "--manifest", "poses.json", "--out", "results", "--skip-posebusters"
    ])
    assert args.manifest == Path("poses.json") and args.skip_posebusters


@pytest.mark.parametrize("words", [[], ["--help"]])
def test_help(words, capsys):
    with pytest.raises(SystemExit) as exc:
        main(words)
    assert exc.value.code == 0
    assert "inputs.json" in capsys.readouterr().out
