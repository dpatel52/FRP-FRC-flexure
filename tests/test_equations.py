"""
Checks the closed forms against direct force and moment integration of the same
assumed laws. The integrator shares no algebra with the equations, it just walks the
depth and sums, so agreement is evidence and not a tautology.

Run with:  python -m pytest -q
        or python tests/test_equations.py
"""

import math
import numpy as np
import pytest

from frp_flexure import (groups, k_stage1, k_stage2, k_stage3,
                         moment_stage1, moment_stage2, moment_stage3,
                         build_envelope, run_model)

REF = dict(n=1.25, rho=0.010, alpha=0.85, mu=0.35, omega=10.0)


# ------------------------------------------------------- independent integration
def integrate(beta, n, rho, alpha, mu, omega, N=200001):
    """Neutral axis and M/M_cr by quadrature on the normalised section.

    b = h = E = eps_cr = 1. y is measured from the compression face, tension is
    positive below the neutral axis. The axis is found by bisecting on zero axial
    force, then the moment is taken about it.
    """
    y = np.linspace(0.0, 1.0, N)

    def forces(k):
        e = beta * (y - k) / (1.0 - k)
        s = np.where(e >= 0.0,
                     np.where(e <= 1.0, e, mu),
                     -np.where(-e <= omega, -e, omega))
        F = np.trapezoid(s, y)
        M = np.trapezoid(s * (y - k), y)
        e_bar = beta * (alpha - k) / (1.0 - k)
        F_bar = n * rho * e_bar
        return F + F_bar, M + F_bar * (alpha - k)

    lo, hi = 1e-12, 1.0 - 1e-12
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if forces(mid)[0] > 0.0:
            lo = mid
        else:
            hi = mid
    k = 0.5 * (lo + hi)
    return k, forces(k)[1] * 6.0


def closed_form(beta, n, rho, alpha, mu, omega):
    if beta <= 1.0:
        k = k_stage1(n, rho, alpha)
        return k, moment_stage1(beta, k, n, rho, alpha)
    k = k_stage2(beta, n, rho, alpha, mu)
    if beta * k / (1.0 - k) < omega:
        return k, moment_stage2(beta, k, n, rho, alpha, mu)
    k = k_stage3(beta, n, rho, alpha, mu, omega)
    return k, moment_stage3(beta, k, n, rho, alpha, mu, omega)


# --------------------------------------------------------------------------- tests
@pytest.mark.parametrize("beta", [0.5, 1.0, 2.0, 5.0, 20.0, 51.0, 60.0, 80.0, 108.0])
def test_matches_integration_on_reference_section(beta):
    kc, Mc = closed_form(beta, **REF)
    ki, Mi = integrate(beta, **REF)
    # the trapezoid smears the kinks at the crack front and the plastic front, so
    # 1e-5 is the quadrature floor rather than the accuracy of the equations
    assert abs(kc - ki) < 1e-5
    assert abs(Mc - Mi) / max(abs(Mi), 1e-12) < 1e-4


def test_stage1_is_exact():
    """Stage 1 has no kink, so the integrator should hit machine precision."""
    kc, Mc = closed_form(0.5, **REF)
    ki, Mi = integrate(0.5, **REF)
    assert abs(kc - ki) < 1e-12


def test_k1_limits():
    """rho to zero gives the plain elastic section, and a bar on the elastic axis
    contributes no first moment."""
    assert k_stage1(1.25, 0.0, 0.85) == pytest.approx(0.5)
    for nr in (0.001, 0.5, 5.0):
        assert k_stage1(nr, 1.0, 0.5) == pytest.approx(0.5)


def test_minus_branch_is_the_admissible_one():
    """Every root returned must lie in (0,1) across a wide random sweep. The plus
    branch of stage 3 does not, which is what fixes the sign."""
    rng = np.random.default_rng(7)
    for _ in range(400):
        n = rng.uniform(0.1, 2.0)
        rho = rng.uniform(0.002, 0.03)
        alpha = rng.uniform(0.70, 0.95)
        mu = rng.uniform(0.05, 0.95)
        omega = rng.uniform(6.0, 25.0)
        beta = rng.uniform(1.5, 150.0)
        k, _ = closed_form(beta, n, rho, alpha, mu, omega)
        assert 0.0 < k < 1.0


def test_plus_branch_is_inadmissible():
    """The discarded stage 3 branch returns k above 1 on the reference section."""
    n, rho, alpha, mu, omega = (REF[x] for x in ("n", "rho", "alpha", "mu", "omega"))
    _, r, _ = groups(80.0, mu, omega)
    nr = n * rho
    disc = (nr ** 2 * 80.0 ** 2
            + 2.0 * nr * ((1.0 - alpha) * r + 80.0 * omega * (1.0 - 2.0 * alpha))
            + omega ** 2)
    k_plus = (r + omega * 80.0 + nr * 80.0 ** 2 + 80.0 * math.sqrt(disc)) / (r + 2.0 * omega * 80.0)
    assert k_plus > 1.0


def test_stage_continuity():
    """The stage 2 and stage 3 roots must agree where the top fibre reaches omega."""
    n, rho, alpha, mu, omega = (REF[x] for x in ("n", "rho", "alpha", "mu", "omega"))
    lo, hi = 2.0, 200.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        k2 = k_stage2(mid, n, rho, alpha, mu)
        if mid * k2 / (1.0 - k2) < omega:
            lo = mid
        else:
            hi = mid
    b = 0.5 * (lo + hi)
    assert abs(k_stage2(b, n, rho, alpha, mu)
               - k_stage3(b, n, rho, alpha, mu, omega)) < 1e-9


def test_envelope_terminates_on_the_first_limit():
    env = build_envelope(**REF, lambda_cu=23.3, chi_fu=1e9, beta_max=400.0)
    assert env["governing"] == "crushing"
    assert env["lambda_top"][-1] >= 23.3 - 1e-6
    assert env["beta"][-1] < 400.0


def test_worked_example_reproduces_the_paper():
    """Shabani, Asadian and Galal 2025, beam 4#4-F0.5-45, Appendix A of the paper."""
    res = run_model(b=250.0, h=400.0, d=354.5,
                    E=31529.0, sigma_cr=4.025, sigma_res=0.587, f_c=45.0,
                    E_f=61200.0, A_f=507.0,
                    eps_cu=0.005, eps_fu=0.0174,
                    beta_max=300.0, n_points=30000)
    assert res["n"] == pytest.approx(1.941, abs=0.002)
    assert res["rho"] == pytest.approx(0.00507, abs=1e-5)
    assert res["alpha"] == pytest.approx(0.886, abs=0.001)
    assert res["mu"] == pytest.approx(0.146, abs=0.001)
    assert res["omega"] == pytest.approx(11.17, abs=0.02)
    assert res["M_cr"] / 1e6 == pytest.approx(26.83, abs=0.01)
    assert res["governing"] == "bar_rupture"
    assert res["beta"][-1] == pytest.approx(157.7, abs=0.2)
    assert res["M"][-1] / 1e6 == pytest.approx(187.2, abs=0.3)


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-q"]))
