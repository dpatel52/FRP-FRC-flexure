"""
Checks the fourteen zone envelope end to end against MATLAB.

tools/dump_matlab_envelope.m builds the full envelope in MATLAB for three cases,
one without a bonded skin and two with, wired exactly as the driver wires it. This
test replays the same inputs through the Python state machine and compares every
point of the assembled curve, not only the zone solvers in isolation.

Skips if the reference has not been generated, so the suite runs without MATLAB.
"""

import csv
import math
import os

import numpy as np
import pytest

from frp_flexure.stages import assemble_envelope

HERE = os.path.dirname(os.path.abspath(__file__))
TOOLS = os.path.join(os.path.dirname(HERE), "tools")
ENV = os.path.join(TOOLS, "matlab_envelope.csv")
PAR = os.path.join(TOOLS, "matlab_envelope_params.csv")

pytestmark = pytest.mark.skipif(
    not (os.path.exists(ENV) and os.path.exists(PAR)),
    reason="MATLAB envelope reference not generated, run tools/dump_matlab_envelope.m",
)

ARG_ORDER = ["L", "b", "h", "alpha", "E", "epcr", "beta_1", "beta_2", "beta_3",
             "eta_1", "eta_2", "eta_3", "xi", "omega", "eta_c", "n", "kappa",
             "eta_s", "rho_c", "rho_t", "iota", "psi", "rho_x"]


def _load():
    params = {}
    if os.path.exists(PAR):
        with open(PAR, newline="") as fh:
            for row in csv.DictReader(fh):
                params[row["case"]] = {k: float(v) for k, v in row.items() if k != "case"}
    ref = {}
    if os.path.exists(ENV):
        with open(ENV, newline="") as fh:
            for row in csv.DictReader(fh):
                ref.setdefault(row["case"], []).append(
                    (float(row["beta"]), float(row["k"]), float(row["M"])))
    return params, ref


PARAMS, REF = _load()
CASES = sorted(REF.keys())


@pytest.mark.parametrize("case", CASES)
def test_envelope_matches_matlab(case):
    p = PARAMS[case]
    args = [p[a] for a in ARG_ORDER]

    segments = [np.linspace(0.0, 1.0, 200),
                np.linspace(1.0, p["beta_1"], 500),
                np.linspace(p["beta_1"], p["beta_2"], 500),
                np.linspace(p["beta_2"], p["beta_3"], 500)]

    env, beta, stages = assemble_envelope(
        segments, args, p["kappa"], p["omega"], p["epcr"],
        p["beta_1"], p["beta_2"], p["beta_3"], p["alpha"], p["M_cr"])

    rows = REF[case]
    assert len(rows) == env.shape[0], (
        f"{case}: MATLAB produced {len(rows)} points, Python {env.shape[0]}")

    worst_k = worst_M = 0.0
    worst_at = None
    for i, (beta_ref, k_ref, M_ref) in enumerate(rows):
        assert abs(beta[i] - beta_ref) < 1e-9, f"{case}: beta grid differs at {i}"
        if not (math.isfinite(k_ref) and math.isfinite(M_ref)):
            continue
        dk = abs(env[i, 0] - k_ref) / max(abs(k_ref), 1e-12)
        dM = abs(env[i, 1] - M_ref) / max(abs(M_ref), 1e-12)
        if max(dk, dM) > max(worst_k, worst_M):
            worst_at = (i, beta_ref, stages[i], k_ref, env[i, 0], M_ref, env[i, 1])
        worst_k = max(worst_k, dk)
        worst_M = max(worst_M, dM)

    assert worst_k < 1e-9, f"{case}: worst relative k error {worst_k:.3e} at {worst_at}"
    assert worst_M < 1e-9, f"{case}: worst relative M error {worst_M:.3e} at {worst_at}"


def test_three_cases_including_a_bonded_skin():
    assert len(CASES) == 3, f"expected three reference cases, got {CASES}"
    assert any("skin" in c for c in CASES), "no bonded-skin case in the reference"


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-q"]))
