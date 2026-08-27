"""
Closed-form neutral axis and moment for an FRP-reinforced FRC or UHPC section.

Derived in Maple 2022 from first principles, one stage at a time, strain field then
forces then equilibrium then moment. Three stages instead of the four the steel case
needs, because an FRP bar is linear elastic to rupture and never yields.

Material laws
    tension       sigma = E eps            eps <= eps_cr
                  sigma = mu sigma_cr      eps >  eps_cr
    compression   sigma = E eps            |eps| <= omega eps_cr
                  sigma = omega sigma_cr   |eps| >  omega eps_cr
    FRP bar       sigma = E_f eps          throughout

Notation
    beta   = eps_bot / eps_cr        the driver
    k      = c / h                   neutral axis from the compression face
    n      = E_f / E                 bar modular ratio
    rho    = A_f / (b h)             bar ratio, gross section
    alpha  = d / h                   bar depth ratio
    mu     = sigma_res / sigma_cr    residual tension
    omega  = eps_cy / eps_cr         compressive yield
    M_cr   = sigma_cr b h^2 / 6

The bar appears only as the product n rho, so modulus and area are interchangeable
in equilibrium. Three shorthand groups keep every expression to one line.
"""

import numpy as np

__all__ = ["groups", "k_stage1", "k_stage2", "k_stage3",
           "moment_stage1", "moment_stage2", "moment_stage3", "stage_of"]


def groups(beta, mu, omega):
    """The three shorthand groups q, r and s."""
    q = 2.0 * mu * (beta - 1.0) + 1.0
    r = q + omega ** 2
    s = 3.0 * mu * (beta ** 2 - 1.0) + 2.0
    return q, r, s


# --------------------------------------------------------------------- neutral axis
def k_stage1(n, rho, alpha):
    """Uncracked, everything elastic. Explicit and independent of beta."""
    return (1.0 + 2.0 * alpha * n * rho) / (2.0 * (1.0 + n * rho))


def k_stage2(beta, n, rho, alpha, mu):
    """Cracked tension, elastic compression. Minus branch is the root in (0,1)."""
    q, _, _ = groups(beta, mu, 0.0)
    nr = n * rho
    disc = nr ** 2 * beta ** 2 + 2.0 * nr * (alpha * beta ** 2 + (1.0 - alpha) * q) + q
    return (q + nr * beta ** 2 - beta * np.sqrt(disc)) / (q - beta ** 2)


def k_stage3(beta, n, rho, alpha, mu, omega):
    """Cracked tension, plastic compression. Minus branch is the root in (0,1).

    The sign is fixed by admissibility, not by the algebra. Both branches satisfy
    the same quadratic, and the plus branch returns k just above 1, which puts the
    neutral axis below the tension face.
    """
    _, r, _ = groups(beta, mu, omega)
    nr = n * rho
    disc = (nr ** 2 * beta ** 2
            + 2.0 * nr * ((1.0 - alpha) * r + beta * omega * (1.0 - 2.0 * alpha))
            + omega ** 2)
    return (r + omega * beta + nr * beta ** 2 - beta * np.sqrt(disc)) / (r + 2.0 * omega * beta)


# --------------------------------------------------------------------------- moment
def moment_stage1(beta, k, n, rho, alpha):
    """M / M_cr, uncracked."""
    return 2.0 * beta / (1.0 - k) * (k ** 3 + (1.0 - k) ** 3 + 3.0 * n * rho * (alpha - k) ** 2)


def moment_stage2(beta, k, n, rho, alpha, mu):
    """M / M_cr, cracked tension and elastic compression."""
    _, _, s = groups(beta, mu, 0.0)
    return (2.0 * beta / (1.0 - k) * (k ** 3 + 3.0 * n * rho * (alpha - k) ** 2)
            + s * (1.0 - k) ** 2 / beta ** 2)


def moment_stage3(beta, k, n, rho, alpha, mu, omega):
    """M / M_cr, cracked tension and plastic compression.

    The three terms are the compression block, the cracked tension block and the bar,
    in that order.
    """
    _, _, s = groups(beta, mu, omega)
    return (3.0 * omega * k ** 2
            + (s - omega ** 3) * (1.0 - k) ** 2 / beta ** 2
            + 6.0 * beta * n * rho * (alpha - k) ** 2 / (1.0 - k))


# ---------------------------------------------------------------------- stage rule
def stage_of(beta, k, omega):
    """Which stage a state belongs to.

    Stage 1 while uncracked. Stage 3 once the top fibre has passed compressive yield.
    Stage 2 in between. The 2 to 3 boundary has no closed form because k depends on
    beta, so it is found by evaluating the condition rather than solving for it.
    """
    if beta <= 1.0:
        return 1
    return 3 if beta * k / (1.0 - k) >= omega else 2
