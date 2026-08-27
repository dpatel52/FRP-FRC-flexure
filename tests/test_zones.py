"""
Checks the translated zone solvers against the MATLAB originals they came from.

tools/dump_matlab_reference.m evaluates all fourteen zones in MATLAB over twelve
random parameter sets and eight values of beta, and writes the results next to it.
This test replays exactly those inputs through the Python translation and compares.

The translation is operator substitution only, so agreement should be at round-off,
not merely close. The tolerance below is set accordingly.

If the reference CSVs are absent the test skips rather than fails, so the suite still
runs on a machine without MATLAB.
"""

import csv
import math
import os

import pytest

from frp_flexure import zones

HERE = os.path.dirname(os.path.abspath(__file__))
TOOLS = os.path.join(os.path.dirname(HERE), "tools")
REF = os.path.join(TOOLS, "matlab_reference.csv")
PAR = os.path.join(TOOLS, "matlab_params.csv")

pytestmark = pytest.mark.skipif(
    not (os.path.exists(REF) and os.path.exists(PAR)),
    reason="MATLAB reference not generated, run tools/dump_matlab_reference.m",
)

ARG_ORDER = ["L", "b", "h", "alpha", "E", "epcr", "beta_1", "beta_2", "beta_3",
             "eta_1", "eta_2", "eta_3", "xi", "omega", "eta_c", "n", "kappa",
             "eta_s", "rho_c", "rho_t", "iota", "psi", "rho_x"]


def _load():
    params = {}
    with open(PAR, newline="") as fh:
        for row in csv.DictReader(fh):
            params[int(row["set"])] = {k: float(v) for k, v in row.items() if k != "set"}
    ref = []
    with open(REF, newline="") as fh:
        for row in csv.DictReader(fh):
            ref.append((row["zone"], int(row["set"]), float(row["beta"]),
                        float(row["k"]), float(row["M"])))
    return params, ref


PARAMS, REFROWS = _load() if os.path.exists(REF) and os.path.exists(PAR) else ({}, [])
ZONE_NAMES = sorted({r[0] for r in REFROWS})


def _call(zone_name, s, beta_vector):
    """Call a zone with the WHOLE beta vector, exactly as MATLAB was called.

    This matters for zone111, which builds a scalar M_cr from k at the last entry
    of the array. Its moment therefore depends on the whole vector passed to it,
    not only on the state being evaluated, so calling it one beta at a time gives
    a different answer. The other thirteen de-normalise per step and are immune.
    """
    p = PARAMS[s]
    args = [p[a] for a in ARG_ORDER]
    fn = getattr(zones, zone_name)
    if zone_name == "zone111":
        return fn(beta_vector, *args)
    return fn(p["M_cr"], beta_vector, *args)


@pytest.mark.parametrize("zone_name", ZONE_NAMES)
def test_zone_matches_matlab(zone_name):
    """Every zone must reproduce MATLAB to round-off on every sampled state."""
    worst_k = worst_M = 0.0
    worst_at = None
    checked = 0

    # group the reference rows by parameter set, preserving order, so each zone is
    # called once per set with the same beta vector MATLAB used
    by_set = {}
    for name, s, beta, k_ref, M_ref in REFROWS:
        if name == zone_name:
            by_set.setdefault(s, []).append((beta, k_ref, M_ref))

    for s, rows in by_set.items():
        betas = [r[0] for r in rows]
        k_py, M_py = _call(zone_name, s, betas)

        for j, (beta, k_ref, M_ref) in enumerate(rows):
            if not (math.isfinite(k_ref) and math.isfinite(M_ref)):
                continue      # MATLAB itself produced nan or inf, nothing to compare
            checked += 1
            dk = abs(float(k_py[j]) - k_ref) / max(abs(k_ref), 1e-12)
            dM = abs(float(M_py[j]) - M_ref) / max(abs(M_ref), 1e-12)
            if max(dk, dM) > max(worst_k, worst_M):
                worst_at = (s, beta, k_ref, float(k_py[j]), M_ref, float(M_py[j]))
            worst_k = max(worst_k, dk)
            worst_M = max(worst_M, dM)

    assert checked > 0, f"no finite reference rows for {zone_name}"
    assert worst_k < 1e-10, (
        f"{zone_name}: worst relative k error {worst_k:.3e} at {worst_at}")
    assert worst_M < 1e-10, (
        f"{zone_name}: worst relative M error {worst_M:.3e} at {worst_at}")


def test_all_fourteen_zones_are_covered():
    assert len(ZONE_NAMES) == 14, f"expected 14 zones in the reference, got {ZONE_NAMES}"


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-q"]))
