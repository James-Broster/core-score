"""Binding-site radius recovery: retry at larger radii only after a binding-site failure."""

import math


def usable(score):
    return score is not None and math.isfinite(float(score)) and float(score) >= 0


def binding_site_failure(result):
    """Whether OST's reasons name a binding-site failure (not a graph or general failure)."""
    if isinstance(result, dict):
        return any(binding_site_failure(v) for v in result.values())
    if isinstance(result, (list, tuple)):
        return any(binding_site_failure(v) for v in result)
    return result in ("target_binding_site", "model_binding_site", "binding_site")


def score_with_radius_fallback(evaluate, fallback_radii=(6.0, 8.0)):
    """Score at 4 A; after a binding-site failure, try each larger radius until one scores.

    `evaluate(radius)` returns dict(score, reason). Returns the selected attempt and a
    record of the standard RMSD, the radius used and every attempt.
    """
    attempts = [dict(radius=4.0, **evaluate(4.0))]
    if not usable(attempts[0]["score"]) and binding_site_failure(attempts[0]):
        for radius in fallback_radii:
            attempts.append(dict(radius=radius, **evaluate(radius)))
            if usable(attempts[-1]["score"]) or not binding_site_failure(attempts[-1]):
                break
    selected = attempts[-1] if usable(attempts[-1]["score"]) else attempts[0]
    scored = usable(selected["score"])
    return selected, dict(
        standard_rmsd=attempts[0]["score"] if usable(attempts[0]["score"]) else None,
        binding_site_radius_angstrom=selected["radius"] if scored else None,
        radius_fallback_used=scored and selected["radius"] != 4.0,
        attempts=attempts,
    )
