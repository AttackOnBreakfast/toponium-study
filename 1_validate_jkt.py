#!/usr/bin/env python3
"""STAGE 1 -- validate our solver against JKT.  Uses: toponium.JKTPotential, LSSolver.  No swData, no events.

Our solver (LSSolver) with the JKT potential (JKTPotential) against
Jezabek, Kuhn, Teubner, Z. Phys. C 56 (1992) 653 (reference/BF01474740.pdf):
reproduces their Figs. 2-5.

Run from toponium_study/ (~15 s):  python3 1_validate_jkt.py

=====================================================================================
EQUATION MAP  (JKT equation  ->  where it is used here)
=====================================================================================
  JKT (7)    dsigma/d^3p = (24 pi^2 alpha^2 Q_t^2/(s m_t^2)) Gamma_t |G|^2
             (d^3p means d^3p/(2pi)^3), divided by sigma_point = 4 pi alpha^2/(3s):
             R = (8 pi/m_t^2) Gamma_t Int d^3p/(2pi)^3 |G|^2     (Q_t = 2/3)   r_ratio
             free limit: R -> 3 Q_t^2 (3/2) v
  JKT (12)   Int d^3p/(2pi)^3 Gamma_t |G|^2 = -Im G(x=0, x'=0)  (optical theorem):
             R = (8 pi/m_t^2) (-Im G(x=0))                               r_ratio
  JKT (8)    Gamma_t = G_F m_t^3/(8 sqrt2 pi) (1-y)^2 (1+2y) [1 - (2 alpha_s/3pi) f(y)],
             y = m_W^2/m_t^2                                              gamma_t
  JKT (9)    f(y) (QCD correction to the width)                           gamma_t
  JKT Figs. 2-5   G(p), |pG|^2, R(E) for m_t = 120/150/180 GeV            main

Settings as in JKT: alpha_s = 0.125 (their Fig. 5), m_t = 120/150/180 GeV,
Gamma_t from JKT (8)-(9) (0.299/0.807/1.563 GeV; their text rounds to 0.3/0.8/1.5).
JKT plot -G in the sign convention of JKT (6); the plots below use their sign.
The JKT_READ numbers were read by eye from the PDF figures (~2-3% precision).

Output -> results/1_validate_jkt/{jkt_figs.png, summary.txt}
"""

import os

import numpy as np
from scipy.special import spence
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from toponium import JKTPotential, LSSolver


ALPHA = 0.125
OUT = "results/1_validate_jkt"

# (m_t, E, observable) -> value read from the PDF figures
JKT_READ = {
    "fig2 Re(-G)(0), m=120 E=-2.3": -6.3, "fig2 Im(-G)(0), m=120 E=-2.3": 20.6,
    "fig3 Re(-G)(0), m=120 E=0": 3.7, "fig3 Im(-G)(0), m=120 E=0": -2.7,
    "fig3 min Im(-G), m=120 E=0": -7.2, "fig3 max Re(-G), m=120 E=0": 6.5,
    "fig4a peak |pG|^2, E=-2.3": 10800, "fig4a peak |pG|^2, E=0": 4000, "fig4a peak |pG|^2, E=2": 7050,
    "fig4b peak |pG|^2, E=-2.7": 2070, "fig4b peak |pG|^2, E=0": 1490, "fig4b peak |pG|^2, E=2": 2150,
    "fig4c peak |pG|^2, E=-2.9": 740, "fig4c peak |pG|^2, E=0": 725, "fig4c peak |pG|^2, E=2": 930,
    "fig5 R peak (m=120)": 3.78, "fig5 E of R peak": -2.33, "fig5 R min (1S-2S)": 0.66,
    "fig5 R 2S bump": 1.23, "fig5 R(E=0)": 1.09, "fig5 R(E=4)": 1.18,
}


def gamma_t(m, alpha_mz=ALPHA, mW=80.2, GF=1.16637e-5):
    """JKT (8):  Gamma_t = G_F m_t^3/(8 sqrt2 pi) (1-y)^2 (1+2y) [1 - (2 alpha_s/(3 pi)) f(y)],
    y = m_W^2/m_t^2, with JKT (9):
      f(y) = pi^2 + 2 Li2(y) - 2 Li2(1-y)
             + [4y(1-y-2y^2) ln y + 2(1-y)^2 (5+4y) ln(1-y) - (1-y)(5+9y-6y^2)] / [2(1-y)^2 (1+2y)].
    alpha_s at m_t from JKT (14). Check: f(0) = 2 pi^2/3 - 5/2."""

    Li2 = lambda z: spence(1 - z)
    y = mW ** 2 / m ** 2
    f = (np.pi ** 2 + 2 * Li2(y) - 2 * Li2(1 - y)
         + (4 * y * (1 - y - 2 * y ** 2) * np.log(y) + 2 * (1 - y) ** 2 * (5 + 4 * y) * np.log(1 - y)
            - (1 - y) * (5 + 9 * y - 6 * y ** 2)) / (2 * (1 - y) ** 2 * (1 + 2 * y)))
    a = JKTPotential.alpha_2loop(m * m, JKTPotential.lambda_msbar_2loop(alpha_mz))
    return GF * m ** 3 / (8 * np.sqrt(2) * np.pi) * (1 - y) ** 2 * (1 + 2 * y) * (1 - 2 * a / (3 * np.pi) * f)


def r_ratio(solver, E):
    """R(E) two ways (JKT (7) and JKT (12)), from the node solution G(q_i) = G0(q_i) K(q_i) of JKT (21):
      R_7  = (8 pi/m_t^2) Gamma_t (1/(2 pi^2)) sum_i w_i q_i^2 |G(q_i)|^2
      R_12 = (8 pi/m_t^2) (-(1/(2 pi^2)) sum_i w_i q_i^2 Im G(q_i))     (= -Im G(x=0))"""

    q, w, m = solver.q, solver.w, solver.m
    g = solver.G0(E, q) * solver.solve_K(E)
    R7 = 8 * np.pi / m ** 2 * solver.Gam * np.sum(w * q ** 2 * np.abs(g) ** 2) / (2 * np.pi ** 2)
    R12 = 8 * np.pi / m ** 2 * (-np.sum(w * q ** 2 * g.imag) / (2 * np.pi ** 2))
    return R7, R12


def main():
    os.makedirs(OUT, exist_ok=True)
    pot = JKTPotential(ALPHA)
    ours = {}
    solvers = {}
    for m in (120, 150, 180):
        solvers[m] = LSSolver(pot, m_t=m, Gamma_t=gamma_t(m))   # grid breaks 0, 4, 16, 60, m_t

    fig, ax = plt.subplots(2, 3, figsize=(18, 10))

    # --- JKT Figs. 2 and 3: m=120, E=-2.3 and 0 ---
    p = np.linspace(0.02, 40, 400)
    for k, E in enumerate((-2.3, 0.0)):
        g = -solvers[120].G(E, p)  # JKT sign
        a = ax[0, k]
        a.plot(p, g.real, "--", label="Re G [1/GeV]")
        a.plot(p, g.imag, ":", label="Im G")
        a.plot(p, 0.002 * np.abs(p * g) ** 2, "-", label="|pG|² × 0.002")
        a.axhline(0, color="k", lw=0.5)
        a.set(title=f"JKT Fig. {2 + k}: m_t=120, E={E} GeV", xlabel="p [GeV]", ylim=(-8, 22))
        a.legend()
        tag = "fig2" if k == 0 else "fig3"
        ours[f"{tag} Re(-G)(0), m=120 E={E:g}"] = g[0].real
        ours[f"{tag} Im(-G)(0), m=120 E={E:g}"] = g[0].imag
        if k == 1:
            ours["fig3 min Im(-G), m=120 E=0"] = g.imag.min()
            ours["fig3 max Re(-G), m=120 E=0"] = g.real.max()

    # --- JKT Fig. 4a-c ---
    a = ax[0, 2]
    for m, E1, ls, sub in ((120, -2.3, "C0", "a"), (150, -2.7, "C1", "b"), (180, -2.9, "C2", "c")):
        pp = np.linspace(0.05, m / 2, 500)
        for E, style in ((E1, "-"), (0.0, ":"), (2.0, "--")):
            y = np.abs(pp * solvers[m].G(E, pp)) ** 2
            ours[f"fig4{sub} peak |pG|^2, E={E:g}"] = y.max()
            a.plot(pp, y / {120: 1, 150: 0.2, 180: 0.08}[m], style, color=ls,
                   label=f"m_t={m}, E={E:g}" if E == E1 else None)
    a.set(title="JKT Fig. 4: |pG|² (150 ×5, 180 ×12.5 for display)", xlabel="p [GeV]")
    a.legend()

    # --- JKT Fig. 5: R, two ways ---
    Es = np.arange(-5, 5.001, 0.02)
    S = solvers[120]
    R_int, R_opt = [], []
    for E in Es:
        r7, r12 = r_ratio(S, E)
        R_int.append(r7)
        R_opt.append(r12)
    R_int, R_opt = np.array(R_int), np.array(R_opt)
    a = ax[1, 0]
    a.plot(Es, R_int, "-", label="Γ ∫|G|² (Eq. 7)")
    a.plot(Es, R_opt, ":", label="−Im G(x=0) (Eq. 12)")
    a.set(title=f"JKT Fig. 5: R, m_t=120, αs={ALPHA}", xlabel="E [GeV]", ylabel="R", ylim=(0, 4))
    a.legend()
    i = R_int.argmax()
    ours["fig5 R peak (m=120)"], ours["fig5 E of R peak"] = R_int[i], Es[i]
    between = (Es > Es[i]) & (Es < -0.5)
    j = np.where(between)[0][R_int[between].argmin()]
    ours["fig5 R min (1S-2S)"] = R_int[j]
    after = (Es > Es[j]) & (Es < -0.4)
    ours["fig5 R 2S bump"] = R_int[after].max()
    ours["fig5 R(E=0)"] = R_int[np.argmin(abs(Es))]
    ours["fig5 R(E=4)"] = R_int[np.argmin(abs(Es - 4))]

    a = ax[1, 1]
    a.semilogy(Es, np.abs(R_int / R_opt - 1))
    a.set(title="optical theorem: |R(Eq.7)/R(Eq.12) − 1|", xlabel="E [GeV]")

    # --- table ---
    lines = [f"JKTPotential(alpha_s(mZ)={ALPHA}): Lambda_MSbar^(5)={pot.Lam * 1e3:.0f} MeV, "
             f"C={pot.C:+.4f}, q_cut={pot.q_cut * 1e3:.0f} MeV, energy shift (JKT (20))={pot.shift:+.3f} GeV",
             "Gamma_t (JKT Eq. 8): " + ", ".join(f"m={m}: {gamma_t(m):.3f}" for m in (120, 150, 180)),
             f"optical theorem, max |Eq.7/Eq.12 - 1| over E in [-5,5]: {np.max(np.abs(R_int / R_opt - 1)):.1e}",
             "", f"{'observable':36s} {'JKT (read)':>11s} {'ours':>10s} {'diff':>7s}"]
    for key, ref in JKT_READ.items():
        v = ours[key]
        d = f"{v - ref:+.2f}" if "E of" in key else f"{v / ref - 1:+.1%}"
        lines.append(f"{key:36s} {ref:11.2f} {v:10.2f} {d:>7s}")
    ax[1, 2].axis("off")
    ax[1, 2].text(0, 1, "\n".join(lines[3:]), family="monospace", fontsize=8.5, va="top")
    fig.suptitle("Our solver (LSSolver, JKTPotential) reproducing JKT Figs. 2-5")
    fig.tight_layout()
    fig.savefig(f"{OUT}/jkt_figs.png", dpi=120)
    open(f"{OUT}/summary.txt", "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
