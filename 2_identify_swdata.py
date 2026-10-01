#!/usr/bin/env python3
"""STAGE 2 -- identify what Fuks's swData contains.  Uses: toponium.FuksTables, CoulombExact, JKTPotential, LSSolver.

Which potential produced Fuks et al.'s swData tables? Competing hypotheses
are fitted to the tables (complex G~, so the phase counts too) and the data
picks the model.

Run from toponium_study/ (~3 min):  python3 2_identify_swdata.py

=====================================================================================
EQUATION MAP  (paper equation  ->  hypothesis)
=====================================================================================
  Fuks (12) / JKT (13) with alpha fixed, tree level      H1 (exact: CoulombExact)
  Fuks (12) / JKT (13) with alpha fixed, NLO bracket     H2
  JKT (13)-(15), (19)-(20) running + Richardson IR       H3 (JKTPotential)
    H3-tree: JKT (13) without its square bracket         (nlo_bracket=False)
    H3-NLO:  JKT (13) as printed (= Fuks (12))            (nlo_bracket=True)
  Fuks (11)  G0; swData stores -G (sign of Fuks (7))

  H1  fixed-alpha tree-level Coulomb, exact (CoulombExact); alpha, m_t,
      Gamma_t and an energy offset dE all free
  H2  fixed alpha = 0.12 with the NLO bracket (paper Eq. 12, alpha frozen)
  H3  JKT'92 running potential, Eqs. 13-15 + 19-20, solved with LSSolver;
      alpha_s(mZ) scanned, dE either JKT-normalized (V(r=1/GeV) = -1/4 GeV) or free.
      The Fuks paper describes its potential in four places, not consistently:
        p.3 text:      "tree-level part of the QCD potential (12) with a running strong
                        coupling constant normalised at alpha_s(mZ) = 0.12"      -> H3-tree
        footnote 1:    Richardson-like regularisation [19] "following the methodology
                        outlined in [16]" ([16] = JKT)                           -> JKT (15), (19)-(20)
        Fig. 1 caption: "tree-level Coulombic potential ... fixing the strong coupling
                        constant at the Z-pole alpha_s(mZ) = 0.12"               -> H1 at alpha = 0.12
        Sec. 4 (after Eq. 23): "running coupling constant at two-loop accuracy";
        Sec. 3: "verified by Green's functions calculated using the QCD potential with
                        NLO corrections"                                         -> H3-NLO

swData stores -G~ in the sign convention of the paper's Eq. 11.
Output -> results/2_identify_swdata/summary.txt
"""

import os

import numpy as np
from scipy.optimize import least_squares, minimize_scalar

from toponium import JKTPotential, LSSolver, CoulombExact, FuksTables


OUT = "results/2_identify_swdata"
os.makedirs(OUT, exist_ok=True)
lines = []


def log(s=""):
    print(s, flush=True)
    lines.append(s)


sw_E, sw_p, sw_G = FuksTables.points()
lookup = FuksTables.lookup()

log("swData: E in [-20,20] GeV (dE 0.1 below -5, 0.02 above), p in [0.01,99.76] GeV (dp 0.25)")

# ================================================================
# H1: exact fixed Coulomb, 4 free parameters
# ================================================================

log("")
log("H1  exact fixed-alpha Coulomb, alpha, m_t, Gamma_t, dE free (complex residual)")
rng = np.random.default_rng(0)
ridge = (sw_E > 5) & (np.abs(sw_p - np.sqrt(173 * np.clip(sw_E, 0, None))) < 12)
for label, sel in [("Fig. 1 window E[-6,4], p<30", (sw_E >= -6) & (sw_E <= 4) & (sw_p <= 30)),
                   ("full table, E>5 on-shell ridge excluded", ~ridge)]:
    idx = rng.choice(np.where(sel)[0], 5000, replace=False)
    E, p, G = sw_E[idx], sw_p[idx], sw_G[idx]

    def res(x):
        r = G / CoulombExact.G(E + x[3], p, x[0], m=x[1], Gam=x[2]) - 1
        return np.concatenate([r.real, r.imag])

    f = least_squares(res, [0.13, 170, 1.4, 0.1], x_scale=[0.01, 1, 0.1, 0.1])
    log(f"    {label}: alpha={f.x[0]:.5f} (C_F*alpha={4 / 3 * f.x[0]:.5f}), m_t={f.x[1]:.2f}, "
        f"Gamma_t={f.x[2]:.4f}, dE={f.x[3] * 1e3:+.1f} MeV, rms {np.sqrt(2 * np.mean(f.fun ** 2)):.1e}")

pts = [(5.0, 29.51), (10.0, 41.51), (15.0, 51.01), (19.9, 59.01)]
log("    swData's own accuracy on the on-shell ridge p = sqrt(m_t E) (vs exact, alpha=0.15): " +
    ", ".join(f"E={E}: {abs(lookup[(E, p)] / CoulombExact.G(np.array(E), np.array(p), 0.15) - 1):.1%}"
              for E, p in pts))

# ================================================================
# All hypotheses on a common grid
# ================================================================

Es = np.round(np.arange(-6, 4.001, 0.2), 2)
ps = np.round(np.arange(0.01, 30, 1.0), 2)
Gsw = np.array([[lookup[(E, p)] for p in ps] for E in Es])
DE_MAX = 3.0  # |dE| range of the free-shift fit below
# the solver table must cover every E + dE the fit asks for; np.interp would silently clamp
Efine = np.arange(Es.min() - DE_MAX, Es.max() + DE_MAX + 0.001, 0.02)


def solver_table(pot):
    return LSSolver(pot).table(Efine, ps)


def shifted(T, dE):
    return np.array([[np.interp(E + dE, Efine, T[:, j].real) + 1j * np.interp(E + dE, Efine, T[:, j].imag)
                      for j in range(len(ps))] for E in Es])


def score(Gm):
    return (np.sqrt(np.mean(np.abs(Gsw / Gm - 1) ** 2)),
            np.sqrt(np.mean(np.log10(np.abs(Gsw / Gm) ** 2) ** 2)))


log("")
log("Comparison on E in [-6,4], p in [0,30] GeV: rms |swData/model - 1| (complex), rms log10 of |G|^2 ratio")
EE, PP = np.meshgrid(Es, ps, indexing="ij")
for al in (0.15, 0.12):
    c, lg = score(CoulombExact.G(EE, PP, al))
    log(f"    H1 Coulomb alpha={al}                 {c:6.3f}   {lg:.3f} dex")
c, lg = score(shifted(solver_table(JKTPotential.coulomb(0.12, nlo_bracket=True)), 0.0))
log(f"    H2 Coulomb alpha=0.12 + NLO bracket    {c:6.3f}   {lg:.3f} dex")

for nlo, lab in ((False, "H3-tree: JKT'92 running potential, tree level (Fuks p.3 text)"),
                 (True, "H3-NLO:  JKT'92 running potential, NLO bracket (JKT (13) as printed)")):
    log(f"    {lab}:")
    log("        alpha_s(mZ)   JKT-normalized        best free dE")
    for amz in np.arange(0.09, 0.161, 0.01):
        T = solver_table(JKTPotential(amz, nlo_bracket=nlo))
        c0, l0 = score(shifted(T, 0.0))
        # coarse scan first (the score has several local minima in dE), then refine
        dE_scan = np.arange(-DE_MAX, DE_MAX + 0.001, 0.1)
        x0 = dE_scan[np.argmin([score(shifted(T, x))[0] for x in dE_scan])]
        f = minimize_scalar(lambda x: score(shifted(T, x))[0],
                            bounds=(max(-DE_MAX, x0 - 0.1), min(DE_MAX, x0 + 0.1)), method="bounded")
        c1, l1 = score(shifted(T, f.x))
        edge = f"  (at the {np.sign(f.x) * DE_MAX:+.0f} GeV scan limit)" if abs(f.x) > DE_MAX - 0.01 else ""
        log(f"        {amz:.2f}        {c0:6.3f}  {l0:.3f} dex    dE={f.x:+.2f} GeV: {c1:.3f}  {l1:.3f} dex{edge}")

log("")
log("alpha_s(mZ)=0.12 runs to 0.15 at ~23 GeV (1-loop) / ~25 GeV (2-loop), near the top Bohr momentum")
log("C_F alpha m_t/2 ~ 17 GeV -- a plausible origin of the 0.15, not stated anywhere by the authors.")
log("(Fuks Eq. 23 and text: with a Coulomb approximation, 1/a0 ~ 20 GeV and alpha_s(1/a0) ~ 0.17.)")

open(f"{OUT}/summary.txt", "w").write("\n".join(lines) + "\n")
print(f"\nSaved {OUT}/summary.txt")
