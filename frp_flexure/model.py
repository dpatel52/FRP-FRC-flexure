"""
Top-level entry point. Takes a section in physical units, forms the dimensionless
groups, runs the envelope, and converts back to moment, curvature and secant stiffness.
"""

import numpy as np

from .envelope import build_envelope

__all__ = ["run_model", "section_parameters"]


def section_parameters(b, h, d, E, sigma_cr, sigma_res, f_c, E_f, A_f):
    """Form the dimensionless groups from a section in physical units.

    Returns a dict carrying n, rho, alpha, mu, omega, eps_cr and M_cr in N-mm.
    """
    eps_cr = sigma_cr / E
    return {
        "n": E_f / E,
        "rho": A_f / (b * h),
        "alpha": d / h,
        "mu": sigma_res / sigma_cr,
        "omega": (f_c / E) / eps_cr,
        "eps_cr": eps_cr,
        "M_cr": sigma_cr * b * h ** 2 / 6.0,
    }


def run_model(b, h, d, E, sigma_cr, sigma_res, f_c, E_f, A_f,
              eps_cu=None, eps_fu=None, beta_max=400.0, n_points=4000):
    """Moment-curvature response of an FRP-reinforced FRC or UHPC section.

    Parameters
    ----------
    b, h, d : float
        Width, overall depth and effective depth, mm.
    E : float
        Matrix modulus, MPa.
    sigma_cr : float
        Matrix cracking stress, MPa. With no direct tension test, 0.60 sqrt(f_c)
        is the usual correlation.
    sigma_res : float
        Residual tensile stress carried across the crack by the fibres, MPa. Take it
        from a direct tension test on the same matrix, or from an ASTM C1609 residual
        strength as f_D150 / 3. Set it near zero for a plain matrix.
    f_c : float
        Compressive strength, MPa.
    E_f, A_f : float
        Bar modulus in MPa and total bar area in mm^2.
    eps_cu, eps_fu : float or None
        Crushing strain and bar rupture strain. Either may be None to leave that
        limit unchecked, but at least one should be given or the curve runs on
        past anything physical.
    beta_max, n_points : float, int
        Sweep range and resolution in beta.

    Returns
    -------
    dict with the dimensionless groups, and arrays for beta, k, curvature (1/mm),
    moment (N-mm), normalised moment, secant stiffness (N-mm^2), and the limit-state
    events with the one that governs.
    """
    p = section_parameters(b, h, d, E, sigma_cr, sigma_res, f_c, E_f, A_f)

    env = build_envelope(
        p["n"], p["rho"], p["alpha"], p["mu"], p["omega"],
        lambda_cu=None if eps_cu is None else eps_cu / p["eps_cr"],
        chi_fu=None if eps_fu is None else eps_fu / p["eps_cr"],
        beta_max=beta_max, n_points=n_points,
    )

    phi = env["beta"] * p["eps_cr"] / ((1.0 - env["k"]) * h)
    M = env["M"] * p["M_cr"]

    with np.errstate(divide="ignore", invalid="ignore"):
        EI_sec = np.where(phi > 0, M / phi, np.nan)

    out = dict(p)
    out.update({
        "beta": env["beta"],
        "k": env["k"],
        "c": env["k"] * h,
        "phi": phi,
        "M": M,
        "M_over_Mcr": env["M"],
        "EI_sec": EI_sec,
        "stage": env["stage"],
        "eps_top": env["lambda_top"] * p["eps_cr"],
        "eps_bar": env["chi_bar"] * p["eps_cr"],
        "events": env["events"],
        "governing": env["governing"],
    })
    return out
