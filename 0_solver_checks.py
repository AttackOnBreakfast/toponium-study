#!/usr/bin/env python3
"""STAGE 0 -- checks of our solver (toponium.LSSolver) on its own.  No swData, no events.

Run from toponium_study/ (~3 min):  python3 0_solver_checks.py
Output -> results/0_solver_checks/{summary.txt, kernel_subtraction.png}

=====================================================================================
EQUATION MAP  (check  ->  which equation it tests)
=====================================================================================
  (a) Coulomb, Fuks (12) at fixed alpha = 0.15: JKT (21)-(23) against the exact
      solution (CoulombExact), for N = 5n = 80 .. 320 (JKT: N = 100-400)
  (b) JKT potential (13)-(15), (19)-(20): self-convergence in N
  (c) JKT text below (20): "K(p) is independent from q_cut for q_cut < 200 MeV";
      q_cut = 200, 100, 20, 10 MeV against JKT's 50 MeV
  (d) JKT (22): A(p) on the auxiliary rule against adaptive quadrature
  (e) CoulombExact: alpha -> 0 gives JKT (6) / Fuks (11); vectorised vs adaptive
  (f) JKT (23) grid: JKT's three subintervals vs the default (extra breaks at 4 and
      60 GeV), above threshold at narrow width (m_t = 120, JKT Fig. 4a, E = 2 GeV)
  Figure: the integrands of (17), (22) and (21) near q = p against q - p (symmetric log
      scale), with the grid nodes of (23) and the auxiliary rule of A(p)
"""

import os

import numpy as np
from scipy import integrate
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from toponium import Constants, JKTPotential, LSSolver, CoulombExact


OUT = "results/0_solver_checks"
BLUE, ORANGE, INK, MUTED = "#2a78d6", "#eb6834", "#222222", "#8a8a85"

E = np.round(np.arange(-6.0, 4.01, 0.5), 2)
P = np.round(np.arange(0.5, 30.01, 1.0), 2)
BELOW = E <= 0   # above threshold the largest errors sit on the on-shell ridge p = sqrt(m_t E)
lines = []


def log(s=""):
    print(s, flush=True)
    lines.append(s)


def stats(d):
    return (f"median {np.median(d):.1e}   max {np.max(d):.1e}   "
            f"(max E<=0: {np.max(d[BELOW]):.1e}, E>0: {np.max(d[~BELOW]):.1e})")


def cquad(f, edges):
    """Complex adaptive quadrature of f over consecutive intervals in `edges` (last may be inf)."""

    re = sum(integrate.quad(lambda q: f(q).real, a, b, limit=500, epsabs=1e-13, epsrel=1e-11)[0]
             for a, b in zip(edges[:-1], edges[1:]))
    im = sum(integrate.quad(lambda q: f(q).imag, a, b, limit=500, epsabs=1e-13, epsrel=1e-11)[0]
             for a, b in zip(edges[:-1], edges[1:]))
    return re + 1j * im


def check_A(pot):
    """(d): max |A_aux / A_quad - 1| over a few p and E."""

    S = LSSolver(pot, n=16)
    worst = 0.0
    for p in (0.03, 0.3, 2.0, 10.0, 25.0, 80.0):
        aux = [S._aux(p)]
        f_edges = sorted({0.0, p, 200.0} | {k for k in pot.kinks} |
                         ({p - pot.q_cut, p + pot.q_cut} if pot.q_cut else set()))
        f_edges = [x for x in f_edges if x >= 0] + [np.inf]
        for En in (-3.0, 0.0, 2.0, 4.0):
            f = lambda q: S.kernel(np.array([p]), np.array([q]))[0] * S.G0(En, q) if q != p else 0.0
            worst = max(worst, abs(S.script_A(En, aux)[0] / cquad(f, f_edges) - 1))
    return worst


def figure():
    """JKT (17), (22), (21) integrands at p = 5 GeV, E = -2 GeV (JKT potential, q_cut = 50 MeV),
    against q - p on a symmetric log scale, so the region of the peak (width ~ q_cut) is spread out.
    Left: (17) with the grid nodes of (23). Middle: (22) with the auxiliary rule of A(p), whose panels
    are graded geometrically towards q = p (evenly spaced on this axis). Right: (21) with the grid nodes."""

    S = LSSolver(JKTPotential(0.12))
    p, En, half = 5.0, -2.0, 1.5
    d = np.geomspace(1e-5, half, 3000)
    q = np.concatenate([p - d[::-1], p + d])
    Kq = S.G(En, q) / S.G0(En, q)
    Kp = S.G(En, np.array([p]))[0] / S.G0(En, p)
    VG = S.kernel(np.full_like(q, p), q) * S.G0(En, q)
    near = np.abs(S.q - p) < half
    qn, Kn = S.q[near], S.solve_K(En)[near]
    VGn = S.kernel(np.full_like(qn, p), qn) * S.G0(En, qn)
    Qa, _ = LSSolver.aux_rule(p)
    Qa = Qa[np.abs(Qa - p) < half]
    VGa = S.kernel(np.full_like(Qa, p), Qa) * S.G0(En, Qa)

    plt.rcParams.update({"font.size": 14, "legend.fontsize": 12})
    fig, ax = plt.subplots(1, 3, figsize=(16, 5.2), layout="constrained")
    panels = (
        (ax[0], VG * Kq, qn, VGn * Kn, BLUE, "grid nodes (23)",
         r"JKT (17): $\hat V G_0 K(q)$" "\npeak missed by the grid"),
        (ax[1], VG, Qa, VGa, BLUE, "auxiliary nodes",
         r"JKT (22): $\hat V G_0$ in $\mathcal{A}(p)$" "\nauxiliary rule, graded towards $q=p$"),
        (ax[2], VG * (Kp - Kq), qn, VGn * (Kp - Kn), ORANGE, "grid nodes (23)",
         r"JKT (21): $\hat V G_0 [K(p)-K(q)]$" "\n~1000x smaller, odd about $q=p$"))
    for a, y, xn, yn, c, lab, title in panels:
        a.plot(q - p, y.real, color=c, lw=2)
        a.plot(xn - p, yn.real, "o", ms=4 if len(xn) < 100 else 2, color=INK, label=lab)
        a.axhline(0, color=MUTED, lw=0.6)
        a.axvspan(-S.pot.q_cut, S.pot.q_cut, color=MUTED, alpha=0.15, lw=0, label=r"$|q-p|<q_\mathrm{cut}$")
        a.set_xscale("symlog", linthresh=1e-3)
        a.set_xticks([-1, -1e-2, 0, 1e-2, 1], ["-1", "-0.01", "0", "0.01", "1"])
        a.set(xlim=(-half, half), title=title, ylabel="real part", xlabel="q - p [GeV]")
        a.legend(frameon=True, framealpha=0.9, edgecolor="none", loc="center right" if a is ax[0] else "best")
        a.spines[["top", "right"]].set_visible(False)
    fig.savefig(f"{OUT}/kernel_subtraction.png", dpi=140)
    plt.close(fig)


def main():
    os.makedirs(OUT, exist_ok=True)
    EE, PP = np.meshgrid(E, P, indexing="ij")
    log(f"Test region: E in [{E[0]}, {E[-1]}] GeV (step 0.5), p in [{P[0]}, {P[-1]}] GeV (step 1); m_t = 173, Gamma_t = 1.49")
    log("Grid, JKT (23): n Gauss-Legendre points on [0,4], [4,16], [16,60], [60,173] GeV and t = 1/q on (0, 1/173)")

    log("")
    log("(a) Coulomb alpha = 0.15: |G/G_exact - 1| vs n (N = 5n)")
    coul = JKTPotential.coulomb(0.15)
    Gex = CoulombExact.G(EE, PP, 0.15)
    for n in (16, 24, 32, 48, 64):
        log(f"    n = {n:3d} (N = {5 * n:3d}):   {stats(np.abs(LSSolver(coul, n).table(E, P) / Gex - 1))}")

    log("")
    jkt = JKTPotential(0.12)
    log(f"(b) JKT potential, alpha_s(mZ) = 0.12, q_cut = 50 MeV (energy shift from JKT (20): {jkt.shift:+.3f} GeV):"
        " |G_n/G_96 - 1|")
    Gref = LSSolver(jkt, 96).table(E, P)
    for n in (16, 24, 32, 48, 64):
        log(f"    n = {n:3d} (N = {5 * n:3d}):   {stats(np.abs(LSSolver(jkt, n).table(E, P) / Gref - 1))}")

    log("")
    log("(c) q_cut dependence, n = 48: |G(q_cut)/G(50 MeV) - 1|  (JKT: independent below 200 MeV)")
    G50 = LSSolver(jkt).table(E, P)
    for qc in (0.2, 0.1, 0.02, 0.01):
        pot = JKTPotential(0.12, q_cut=qc)
        d = np.abs(LSSolver(pot).table(E, P) / G50 - 1)
        log(f"    q_cut = {qc * 1e3:3.0f} MeV (shift {pot.shift:+.3f} GeV):   {stats(d)}")

    log("")
    log("(d) JKT (22) A(p) on the auxiliary rule vs adaptive quadrature, p = 0.03..80 GeV, E = -3..4 GeV:")
    for pot in (coul, jkt):
        log(f"    {pot.label:32s} max rel. difference {check_A(pot):.1e}")

    log("")
    log("(e) CoulombExact")
    d0 = max(abs(CoulombExact.G_adaptive(En, p, 0.0) / Constants.G0(En, p) - 1) for En, p in [(-3, 0.5), (0, 10), (2, 20)])
    log(f"    alpha -> 0 vs G0 (JKT (6) / Fuks (11)): max rel. difference {d0:.1e}")
    dv = max(abs(CoulombExact.G(np.array(En), np.array(p), a) / CoulombExact.G_adaptive(En, p, a) - 1)
             for a in (0.12, 0.15, 0.2) for En, p in [(-4., 0.3), (-1.5, 0.3), (0., 5.), (3., 25.)])
    log(f"    G_coulomb_vec vs G_coulomb (adaptive): max rel. difference {dv:.1e}")

    log("")
    log("(f) grid choice, m_t = 120, Gamma_t = 0.299 (JKT Fig. 4a), alpha_s(mZ) = 0.125:"
        " peak of |p G|^2 at E = 2 GeV vs n")
    pot, pp = JKTPotential(0.125), np.linspace(0.05, 60, 500)
    for br, lab in (((0, np.sqrt(120 * 0.299), 120), "JKT (0, a1, a2) "), (None, "default (+4, 60)")):
        peaks = [np.max(np.abs(pp * LSSolver(pot, n, 120, 0.299, br).G(2.0, pp)) ** 2) for n in (48, 96, 160)]
        log(f"    {lab}  n = 48, 96, 160:  " + ", ".join(f"{x:.0f}" for x in peaks) + "   (JKT figure: ~7050)")

    figure()
    open(f"{OUT}/summary.txt", "w").write("\n".join(lines) + "\n")
    print(f"\nSaved {OUT}/summary.txt, {OUT}/kernel_subtraction.png")


if __name__ == "__main__":
    main()
