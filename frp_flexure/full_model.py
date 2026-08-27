"""
Fourteen zone model, top level.

Takes a section in physical units with its four material laws, forms the parameters
the zone solvers need, runs the state machine, and returns moment and curvature.

This is the solution the paper uses. The three stage form in equations.py is its
reduction for the case of a bilinear tension law and a bar that never yields.
"""

import numpy as np

from .stages import assemble_envelope
from . import zones as Z

__all__ = ["run_full_model", "derive_parameters"]


def derive_parameters(b, h, cover, E, epsilon_cr,
                      mu_1, mu_2, mu_3, beta_1, beta_2, beta_3,
                      xi, omega, mu_c, eps_cu,
                      E_bar, eps_sy, mu_s, eps_su,
                      rho_t, rho_c=0.0,
                      E_skin=0.0, A_skin=0.0, iota=0.0):
    """Form the dimensionless groups the zone solvers take.

    The tension law is quad-linear through (beta_i, mu_i). The compression law is
    bilinear through omega and mu_c. The bar is bilinear through eps_sy and mu_s,
    and becomes linear elastic to rupture when mu_s = 1. The bonded skin enters
    through psi, rho_x and the activation strain iota.
    """
    alpha = (h - cover) / h
    lambda_cu = eps_cu / epsilon_cr
    kappa = eps_sy / epsilon_cr
    chi_su = eps_su / epsilon_cr
    return {
        "alpha": alpha,
        "eta_1": (mu_1 - 1.0) / (beta_1 - 1.0),
        "eta_2": (mu_2 - mu_1) / (beta_2 - beta_1),
        "eta_3": (mu_3 - mu_2) / (beta_3 - beta_2),
        # each hardening slope carries the stress at the start of its own branch,
        # so eta_c is scaled by omega and eta_s by kappa. The tension slopes need
        # no such factor, since normalisation makes the stress unity at cracking.
        "eta_c": omega * (mu_c - 1.0) / (lambda_cu - omega),
        "n": E_bar / E,
        "kappa": kappa,
        "eta_s": kappa * (mu_s - 1.0) / (chi_su - kappa),
        "lambda_cu": lambda_cu,
        "chi_su": chi_su,
        "psi": (E_skin / E) if E_skin else 0.0,
        "rho_x": (A_skin / (b * h)) if A_skin else 0.0,
        "iota": iota,
    }


def run_full_model(b, h, cover, L, E, epsilon_cr,
                   mu_1, mu_2, mu_3, beta_1, beta_2, beta_3,
                   xi, omega, mu_c, eps_cu,
                   E_bar, eps_sy, mu_s, eps_su,
                   rho_t, rho_c=0.0,
                   E_skin=0.0, A_skin=0.0, iota=0.0,
                   n_points=(200, 500, 500, 500)):
    """Moment-curvature response from the fourteen zone solution.

    Returns a dict with beta, k, curvature, moment, the stage active at each step,
    and the derived parameters.

    xi must differ from unity. The zone expressions divide by (xi - 1), so a value
    of exactly 1 is inadmissible. The usual choice for concrete is 1.001.
    """
    if abs(xi - 1.0) < 1e-12:
        raise ValueError("xi must differ from 1, the zone solvers divide by (xi - 1)")

    p = derive_parameters(b, h, cover, E, epsilon_cr, mu_1, mu_2, mu_3,
                          beta_1, beta_2, beta_3, xi, omega, mu_c, eps_cu,
                          E_bar, eps_sy, mu_s, eps_su, rho_t, rho_c,
                          E_skin, A_skin, iota)

    args = [L, b, h, p["alpha"], E, epsilon_cr, beta_1, beta_2, beta_3,
            p["eta_1"], p["eta_2"], p["eta_3"], xi, omega, p["eta_c"],
            p["n"], p["kappa"], p["eta_s"], rho_c, rho_t,
            p["iota"], p["psi"], p["rho_x"]]

    # M_cr from the cracking neutral axis, as the driver builds it
    kcr, _ = Z.zone111([1.0], *args)
    M_cr = (epsilon_cr * E * b * h ** 2) / (12.0 * (1.0 - float(kcr[-1])))

    n1, n2, n3, n4 = n_points
    segments = [np.linspace(0.0, 1.0, n1),
                np.linspace(1.0, beta_1, n2),
                np.linspace(beta_1, beta_2, n3),
                np.linspace(beta_2, beta_3, n4)]

    env, beta, stages = assemble_envelope(
        segments, args, p["kappa"], omega, epsilon_cr,
        beta_1, beta_2, beta_3, p["alpha"], M_cr)

    k = env[:, 0]
    M = env[:, 1]
    with np.errstate(divide="ignore", invalid="ignore"):
        phi = beta * epsilon_cr / ((1.0 - k) * h)

    return {
        "beta": beta,
        "k": k,
        "c": k * h,
        "M": M,
        "phi": phi,
        "M_cr": M_cr,
        "stage": stages,
        "eps_top": k * beta * epsilon_cr / (1.0 - k),
        "eps_bar": ((-p["alpha"] + k) * beta * epsilon_cr) / (k - 1.0),
        "params": p,
    }
