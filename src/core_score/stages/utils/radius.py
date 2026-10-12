"""Explicit, auditable binding-site radius fallback around a native scorer."""

import json
import math


def finite_score(value):
    return value is not None and math.isfinite(float(value)) and float(value) >= 0


def binding_site_failure(result):
    # Match native reason identifiers, not generic scoring failures/exceptions.
    def walk(value):
        if isinstance(value, dict):
            return any(walk(v) for v in value.values())
        if isinstance(value, (list, tuple)):
            return any(walk(v) for v in value)
        return (
            value in ("target_binding_site", "model_binding_site", "binding_site")
            if isinstance(value, str)
            else False
        )

    return walk(result)


def score_with_radius_fallback(evaluate, fallback_radii=(6.0, 8.0)):
    radii = tuple(float(r) for r in fallback_radii)
    if (
        any(not math.isfinite(r) or r <= 4 for r in radii)
        or tuple(sorted(set(radii))) != radii
    ):
        raise ValueError(
            "Fallback radii must be finite, strictly increasing and above 4 A"
        )
    baseline = evaluate(4.0)
    attempts = [dict(radius=4.0, **baseline)]
    selected = baseline
    used = 4.0
    if not finite_score(baseline.get("score")) and binding_site_failure(baseline):
        for radius in radii:
            candidate = evaluate(radius)
            attempts.append(dict(radius=radius, **candidate))
            if finite_score(candidate.get("score")):
                selected = candidate
                used = radius
                break
            if not binding_site_failure(candidate):
                break
    ok = finite_score(selected.get("score"))
    encode = lambda v: json.dumps(v, default=lambda x: x.item())
    metadata = dict(
        standard_rmsd=baseline["score"] if finite_score(baseline.get("score")) else "",
        binding_site_radius_angstrom=used if ok else "",
        radius_fallback_used=ok and used != 4,
        requires_review=ok and used != 4,
        standard_failure_reason=(
            "" if finite_score(baseline.get("score")) else encode(baseline)
        ),
        radius_attempts=encode(attempts),
    )
    return selected, metadata
