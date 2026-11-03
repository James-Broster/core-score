"""CORE Score: agreement within an ensemble of predicted protein-ligand poses."""

from core_score.pipeline import run
from core_score.analysis.core import core_score as score

__version__ = "0.5.1"

__all__ = ["run", "score"]
