# Shabani_et_al_(2025).py
#
# Beam 4#4-F0.5-45, macro-synthetic FRC with four #4 GFRP bars, four point bending.
# Shabani, Asadian and Galal, Construction and Building Materials 494 (2025) 143217.
#
# Nothing here is fitted to the beam. The residual tensile stress is f_D150 / 3 from
# the ASTM C1609 prisms of that paper, and the cracking stress is the 0.60 sqrt(f_c)
# correlation.

import matplotlib.pyplot as plt

from frp_flexure import run_model

res = run_model(
    # Geometry
    b = 250.0,        # mm
    h = 400.0,        # mm
    d = 354.5,        # mm

    # Matrix
    E         = 31529.0,   # MPa
    sigma_cr  = 4.025,     # MPa, 0.60 sqrt(45)
    sigma_res = 0.587,     # MPa, f_D150 / 3 from ASTM C1609
    f_c       = 45.0,      # MPa

    # GFRP bars, four #4
    E_f = 61200.0,    # MPa
    A_f = 507.0,      # mm2

    # Limits
    eps_cu = 0.005,
    eps_fu = 0.0174,
)

M_kNm = res["M"] / 1e6
phi   = res["phi"] * 1e6          # 1e-6 / mm

print(f"n rho   = {res['n'] * res['rho']:.5f}")
print(f"alpha   = {res['alpha']:.3f}")
print(f"mu      = {res['mu']:.3f}")
print(f"omega   = {res['omega']:.2f}")
print(f"M_cr    = {res['M_cr'] / 1e6:.2f} kN m")
print()
print(f"governing limit = {res['governing']}")
print(f"beta at the end = {res['beta'][-1]:.1f}")
print(f"peak moment     = {M_kNm[-1]:.1f} kN m   (measured 206.0)")

# four point bending, 3836 mm clear span, load points 1000 mm apart
S1 = (3836.0 - 1000.0) / 2.0
print(f"peak total load = {2.0 * res['M'][-1] / S1 / 1000.0:.1f} kN   (measured 290.5)")

fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(phi, M_kNm, "-k", lw=1.6)
for name, (i, b) in res["events"].items():
    if i < len(M_kNm):
        ax.plot(phi[i], M_kNm[i], "o", ms=6, label=name.replace("_", " "))
ax.set_xlabel(r"Curvature, $10^{-6}$/mm")
ax.set_ylabel("Moment, kN m")
ax.set_title("Shabani et al. (2025), beam 4#4-F0.5-45")
ax.legend(frameon=False, fontsize=8)
fig.tight_layout()
plt.show()
