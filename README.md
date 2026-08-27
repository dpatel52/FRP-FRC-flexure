# FRP-FRC-Flexure

[![Python Versions](https://img.shields.io/badge/python-3.8%2B-blue)](https://www.python.org)&nbsp;
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

**Closed-form** flexural response of FRC and UHPC sections reinforced with steel or FRP bars, with an optional externally bonded FRP skin that activates after the section has already been strained.

Fourteen zones, each solved once symbolically, so the response at any strain increment follows by substitution rather than by iteration.

Cite Correct Reference: *Generalized Solutions for Strain Compatibility-Based Flexural Design of Doubly Reinforced Concrete Beams with FRP Composites*, Patel, Pleesudjai and Mobasher (DOI to follow).

It calculates

* **Moment-curvature** envelopes across all fourteen zones
* **Neutral axis** depth through the loading history
* **Stage sequence**, which subsystem changed state and where
* **Strain histories** at the top fibre and at the bar

---

## The model

Four subsystems each carry their own state, and the zone is the combination of the four.

| | states | changes at |
|---|---|---|
| `T` tension in the matrix | 1 to 4 | `beta` = 1, `beta_1`, `beta_2` |
| `C` compression in the matrix | 1 to 2 | compressive yield strain |
| `R` tension bar | 1 to 2 | bar yield |
| `RC` compression bar | 1 to 2 | yield, inside `T` = 4 only |

The matrix carries a quad-linear tension law through `(beta_i, mu_i)` and a bilinear compression law through `omega` and `mu_c`. The bar is bilinear through `eps_sy` and `mu_s`, and becomes linear elastic to rupture when `mu_s = 1`, which is how an FRP bar is entered. The bonded skin adds `psi`, `rho_x` and the activation strain `iota`, and creates no new zone.

Each hardening slope carries the stress at the start of its own branch, so `eta_c` is scaled by `omega` and `eta_s` by `kappa`. The tension slopes need no such factor, since normalisation makes the stress exactly unity at first cracking.

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
from frp_flexure import run_full_model

res = run_full_model(
    # Geometry
    b = 250.0, h = 400.0, cover = 45.5, L = 3836.0,     # mm

    # Matrix, quad-linear tension and bilinear compression
    E = 31528.6, epsilon_cr = 1.27658e-4,
    mu_1 = 0.146, mu_2 = 0.146, mu_3 = 0.146,
    beta_1 = 2.0, beta_2 = 50.0, beta_3 = 300.0,
    xi = 1.001, omega = 11.18, mu_c = 1.0, eps_cu = 0.005,

    # Bars. mu_s = 1 makes the bar linear elastic to rupture, i.e. FRP
    E_bar = 61200.0, eps_sy = 0.0174, mu_s = 1.0, eps_su = 0.0175,
    rho_t = 507.0 / (250.0 * 400.0), rho_c = 0.0,

    # Bonded FRP skin, omit for an unstrengthened section
    E_skin = 209000.0, A_skin = 105.6, iota = 10.0,
)

res["phi"]    # curvature, 1/mm
res["M"]      # moment, N-mm
res["k"]      # neutral axis ratio
res["stage"]  # the zone active at each step, e.g. '4221'
```

`xi` must differ from 1. The zone expressions divide by `(xi - 1)`, and 1.001 is the usual choice for concrete.

---

## Three stage reduction

For a bilinear tension law and a bar that never yields, the fourteen zones collapse to three and every result fits on one line. `frp_flexure/equations.py` carries that form, and `run_model` drives it.

```python
from frp_flexure import run_model
res = run_model(b=250.0, h=400.0, d=354.5, E=31529.0, sigma_cr=4.025,
                sigma_res=0.587, f_c=45.0, E_f=61200.0, A_f=507.0,
                eps_cu=0.005, eps_fu=0.0174)
```

Three shorthand groups hold the expressions to one line each.

```
q = 2 mu (beta - 1) + 1        r = q + omega^2        s = 3 mu (beta^2 - 1) + 2
```

In stages 2 and 3 **the minus branch on the radical is the root in (0,1)**. The sign is fixed by admissibility, not by the algebra, since both branches satisfy the same quadratic. The plus branch of stage 3 returns `k` just above 1, which puts the neutral axis below the tension face.

---

## Where the zone equations come from

`frp_flexure/zones.py` is generated from the MATLAB sources by `tools/translate_zones.py`, which does operator substitution only. The expressions are Maple output, thousands of characters each, and are never retyped by hand. `sqrt` is emitted as `np.emath.sqrt`, because MATLAB returns a complex root for a negative argument and these expressions rely on it, taking the real part at the end.

To regenerate after changing the MATLAB:

```bash
python tools/translate_zones.py <matlab_dir> frp_flexure/zones.py
```

---

## Verification

```bash
python -m pytest -q
```

Two independent checks, 35 tests.

**Against MATLAB.** Every zone is evaluated in MATLAB over twelve random parameter sets and eight values of `beta`, and the assembled envelope is built for three further cases, one plain and two carrying a bonded skin. The Python then replays exactly those inputs.

| | worst relative error |
|---|---|
| `k`, all fourteen zones, 1344 states | 3.3e-16 |
| `M`, all fourteen zones, 1344 states | 8.1e-11 |
| full envelope, 3 cases, 5100 points | under 1e-9 |

**Against integration.** The three stage form is checked against direct force and moment integration of the same laws. The integrator shares no algebra with the equations, it walks the section depth and sums, so agreement is evidence rather than a tautology. Stage 1 agrees in `k` to 1e-12, and across 400 random sections every root returned lies in (0,1).

Regenerate the MATLAB reference with

```bash
matlab -batch "addpath('<matlab_dir>'); dump_matlab_reference; dump_matlab_envelope"
```

The reference CSVs are committed, so the tests run without MATLAB installed.

---

## Scope

Sectional response only. There is no bond slip, no interfacial debonding and no shear. A member that fails by any of those routes is outside what these equations describe, and the bonded skin is carried as a linear elastic layer with a rupture strain and no interfacial limit.
