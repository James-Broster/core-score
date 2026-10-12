"""CORE Score: agreement within an ensemble of predicted protein-ligand poses."""

__version__ = "0.3.0"

from core_score.pipeline import run  # noqa: E402  (after __version__, which pipeline imports)
from core_score.stages.scoring import core_score  # noqa: E402

__all__ = ["run", "core_score"]
