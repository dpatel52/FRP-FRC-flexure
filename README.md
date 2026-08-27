# FRP-FRC-Flexure

[![Python Versions](https://img.shields.io/badge/python-3.8%2B-blue)](https://www.python.org)&nbsp;
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

A **closed-form** solution for the flexural response of **FRC** and **UHPC** sections reinforced with **FRP bars**.
Three stages, explicit expressions for the neutral axis and the moment in each, no iteration anywhere.

Cite Correct Reference: *Generalized Solutions for Strain Compatibility-Based Flexural Design of Doubly Reinforced Concrete Beams with FRP Composites*, Patel, Pleesudjai and Mobasher (DOI to follow).

It calculates

* **Moment-curvature** envelopes
* **Neutral axis** depth through the loading history
* **Limit states**, first cracking, compression yield, crushing and bar rupture
* **Secant stiffness** at any load level

---

## Why three stages and not four

The steel derivation needs four because the bar yields. An FRP bar is linear elastic to rupture, so the yielded-bar stage disappears.

| Stage | Tension | Compression | Enters when |
|---|---|---|---|
| 1 | elastic, uncracked | elastic | `beta <= 1` |
| 2 | cracked | elastic | `beta k / (1 - k) < omega` |
| 3 | cracked | plastic | `beta k / (1 - k) >= omega` |

---

## Parameters

Six dimensionless quantities describe any section, and `beta = eps_bot / eps_cr` sets the load level.

```
n     = E_f / E              rho   = A_f / (b h)        alpha = d / h
mu    = sigma_res / sigma_cr omega = eps_cy / eps_cr    M_cr  = sigma_cr b h^2 / 6
```

The bar appears only as the product `n rho`, so modulus and area are interchangeable in equilibrium. Three shorthand groups keep every expression to one line.

```
q = 2 mu (beta - 1) + 1        r = q + omega^2        s = 3 mu (beta^2 - 1) + 2
```

---

## Installation

```bash
git clone https://github.com/dpatel52/FRP-FRC-flexure.git
cd FRP-FRC-flexure
pip install -e .
```

Only `numpy` is required. `matplotlib` is needed for the examples and `pytest` for the tests.

---

## Quick-start

```python
from frp_flexure import run_model

res = run_model(
    # Geometry
    b = 250.0,        # mm
    h = 400.0,        # mm
    d = 354.5,        # mm

    # Matrix
    E         = 31529.0,   # MPa
    sigma_cr  = 4.025,     # MPa, 0.60 sqrt(f_c) where no direct tension test exists
    sigma_res = 0.587,     # MPa, f_D150 / 3 from ASTM C1609, or 0 for a plain matrix
    f_c       = 45.0,      # MPa

    # FRP bars
    E_f = 61200.0,    # MPa
    A_f = 507.0,      # mm2

    # Limits, either may be None to leave that one unchecked
    eps_cu = 0.005,
    eps_fu = 0.0174,
)

res["phi"]        # curvature, 1/mm
res["M"]          # moment, N-mm
res["k"]          # neutral axis ratio
res["governing"]  # 'crushing' or 'bar_rupture'
res["events"]     # index and beta of each limit state reached
```

`examples/Shabani_et_al_(2025).py` runs a measured beam end to end and plots the response.

---

## The equations

`frp_flexure/equations.py` carries them directly. Stage 1 is explicit and independent of `beta`. Stages 2 and 3 are each a single closed form, and in both **the minus branch on the radical is the root in (0,1)**. The sign is fixed by admissibility, not by the algebra, since both branches satisfy the same quadratic. The plus branch of stage 3 returns `k` just above 1, which puts the neutral axis below the tension face.

`frp_flexure/envelope.py` sweeps `beta`, picks the active stage at each step and terminates the response at whichever limit the section reaches first.

The stages carry no failure strain of their own, so termination is external. Left unterminated the stage 3 root stays admissible, it just stops meaning anything, `k` drifts up toward `alpha` and the compression block is asked to carry the section on its own.

---

## Verification

```bash
python -m pytest -q
```

The equations are checked against direct force and moment integration of the same laws. The integrator shares no algebra with them, it walks the depth and sums, so agreement is evidence rather than a tautology.

| Check | Result |
|---|---|
| Stage 1 against integration, no kink in the stress field | `k` to 1e-12 |
| Stages 2 and 3 against integration | at the quadrature floor, `k` under 1e-5 |
| 400 random sections spanning all three stages | every root in (0,1) |
| Stage 2 and stage 3 roots at the transition | agree to 1e-9 |
| `rho` to zero, and a bar on the elastic axis | `k = 1/2` exactly |
| Worked example against the measured beam | reproduces the published values |

---

## Scope

This release is the closed-form sectional solution only. It does not carry bond slip, interfacial debonding, shear, or the externally bonded laminate of the parent formulation. A member that fails by any of those routes is outside what these equations describe.
