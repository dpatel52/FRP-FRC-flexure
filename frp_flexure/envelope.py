"""
Envelope assembly. Sweeps beta, picks the active stage at each step, and records the
limit states the section passes through.

The stages carry no failure strain of their own, so the response has to be terminated
externally at whichever limit the section reaches first, concrete crushing or bar
rupture. Left unterminated the stage 3 root stays admissible, it just stops meaning
anything, k drifts up toward alpha and the compression block is asked to carry the
section on its own.
"""

import numpy as np

from .equations import (groups, k_stage1, k_stage2, k_stage3,
                        moment_stage1, moment_stage2, moment_stage3)

__all__ = ["build_envelope"]


def _k_and_M(beta, n, rho, alpha, mu, omega):
    """Neutral axis, moment and stage at one value of beta."""
    if beta <= 1.0:
        k = k_stage1(n, rho, alpha)
        return k, moment_stage1(beta, k, n, rho, alpha), 1

    k = k_stage2(beta, n, rho, alpha, mu)
    if beta * k / (1.0 - k) < omega:
        return k, moment_stage2(beta, k, n, rho, alpha, mu), 2

    k = k_stage3(beta, n, rho, alpha, mu, omega)
    return k, moment_stage3(beta, k, n, rho, alpha, mu, omega), 3


def build_envelope(n, rho, alpha, mu, omega,
                   lambda_cu=None, chi_fu=None,
                   beta_max=400.0, n_points=4000):
    """Sweep beta and return the terminated envelope.

    Parameters
    ----------
    n, rho, alpha, mu, omega : float
        The five dimensionless section parameters.
    lambda_cu : float or None
        eps_cu / eps_cr. Crushing limit. None leaves crushing unchecked.
    chi_fu : float or None
        eps_fu / eps_cr. Bar rupture limit. None leaves rupture unchecked.
    beta_max, n_points : float, int
        Sweep range and resolution in beta.

    Returns
    -------
    dict with beta, k, M (as M/M_cr), stage, and the strain histories at the top
    fibre and the bar, all truncated at the governing limit, plus `events`, which
    records the beta and index of each limit state the section reached.
    """
    beta = np.linspace(1e-6, beta_max, int(n_points))
    k = np.empty_like(beta)
    M = np.empty_like(beta)
    stage = np.empty(beta.shape, dtype=int)

    for i, b in enumerate(beta):
        k[i], M[i], stage[i] = _k_and_M(b, n, rho, alpha, mu, omega)

    lam_top = beta * k / (1.0 - k)              # top fibre, in eps_cr units
    chi_bar = beta * (alpha - k) / (1.0 - k)    # bar strain, in eps_cr units

    events = {}
    stop = len(beta) - 1

    i_cr = int(np.argmax(beta >= 1.0))
    if beta[i_cr] >= 1.0:
        events["first_crack"] = (i_cr, beta[i_cr])

    i_cy = np.flatnonzero(stage == 3)
    if i_cy.size:
        events["compression_yield"] = (int(i_cy[0]), beta[i_cy[0]])

    if lambda_cu is not None:
        hit = np.flatnonzero(lam_top >= lambda_cu)
        if hit.size:
            events["crushing"] = (int(hit[0]), beta[hit[0]])
            stop = min(stop, int(hit[0]))

    if chi_fu is not None:
        hit = np.flatnonzero(chi_bar >= chi_fu)
        if hit.size:
            events["bar_rupture"] = (int(hit[0]), beta[hit[0]])
            stop = min(stop, int(hit[0]))

    governing = None
    for name in ("crushing", "bar_rupture"):
        if name in events and events[name][0] == stop:
            governing = name

    sl = slice(0, stop + 1)
    return {
        "beta": beta[sl],
        "k": k[sl],
        "M": M[sl],
        "stage": stage[sl],
        "lambda_top": lam_top[sl],
        "chi_bar": chi_bar[sl],
        "events": events,
        "governing": governing,
    }
