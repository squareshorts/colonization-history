"""Ordered counterfactual contrasts for colonization-history carriers.

Operational bit order A=1, S=2, E=4, R=8. The six-arm decomposition
requires selective initial-state interventions; its algebra is exact for
nonlinear outcomes, but biological implementability is an unverified premise.
"""
from typing import Mapping

REQUIRED = (0, 2, 3, 4, 12, 15)


def ordered_contrasts(f: Mapping[int, float]) -> dict[str, float]:
    """Return five ordered contrasts and the whole-preparation effect."""
    missing = set(REQUIRED).difference(f)
    if missing:
        raise ValueError(f"Missing counterfactual arms: {sorted(missing)}")
    out = {
        "selection_first": f[2] - f[0],
        "acquired_given_selection": f[3] - f[2],
        "environment_first": f[4] - f[0],
        "residual_given_environment": f[12] - f[4],
        "cell_surface_interaction": f[15] - f[3] - f[12] + f[0],
        "total": f[15] - f[0],
    }
    assert abs(sum(v for k, v in out.items() if k != "total") - out["total"]) < 1e-9
    return out
