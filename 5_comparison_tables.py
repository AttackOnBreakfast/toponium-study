#!/usr/bin/env python3
"""STAGE 5 -- tables for comparing our solver with an independent one, number by number.

Run from toponium_study/ (a few minutes):  uv run python 5_comparison_tables.py
Finer grid:  uv run python 5_comparison_tables.py 120   (n = 120 per panel, N = 600; writes results/5_comparison_tables_N600/)
Output -> results/5_comparison_tables/{params.csv, potential.csv, G_small.csv, G_dense.csv, summary.txt}

All tables are for the COLOUR SINGLET (C_F = 4/3, attractive), the only case JKT treat.
Every CSV opens with '#' lines giving the conventions; the other code must match them.

=====================================================================================
EQUATION MAP  (column  ->  JKT equation  ->  toponium.py)
=====================================================================================
  params.csv
    Lambda_MSbar   (14) solved at Q = mZ for Lambda             JKTPotential.lambda_msbar_2loop
    C              (15) continuity at |p| = 5 GeV                JKTPotential.continuity_constant_C
    C1             (19) continuity at |p| = q_cut                JKTPotential.cut_constant_C1
    V0_tilde       (19)-(20) coefficient of d3(p), fixed by
                   V_cut(r = 1 GeV^-1) = -1/4 GeV                JKTPotential.energy_shift x (2pi)^3
    shift          V0_tilde/(2pi)^3: G0 is evaluated at E - shift
    Vcut_r1        (20) check: V_cut(r = 1 GeV^-1), must be -0.25  JKTPotential.V_cut_position + shift
  potential.csv
    alpha_2loop    (14)                                          JKTPotential.alpha_2loop
    alpha_fixed    alpha held fixed (alpha_s(mZ) for JKT rows, alpha for Coulomb rows)
    V_pert         (13), alpha from (14), NLO bracket on         JKTPotential.V_pert
    V_log          (15), |p| < 5 GeV line, at all Q              JKTPotential.V_ir
    V_jkt          (15): V_pert above 5 GeV, V_log below         JKTPotential.V_jkt
    V_cut          (19) without the d3 term                      JKTPotential.V_cut
    V_coulomb      (13) with alpha = alpha_fixed, no bracket: -4 pi C_F alpha / Q^2
  G_small.csv, G_dense.csv
    G              (16)-(23), the full JKT potential (V_pert, V_log, V_cut all on), or Coulomb
                                                                 LSSolver(pot).table
    G0             (6) at the physical E (no shift)              Constants.G0
    G_exact        Coulomb rows only: exact solution             CoulombExact.G
"""

import csv
import os
import sys

import numpy as np

from toponium import Constants, JKTPotential, LSSolver, CoulombExact, C_F, N_F, M_T, GAMMA_T


N_PER_PANEL = int(sys.argv[1]) if len(sys.argv) > 1 else 48
OUT = "results/5_comparison_tables" + ("" if N_PER_PANEL == 48 else f"_N{5 * N_PER_PANEL}")
MZ = 91.1876

# settings: (id, kind, alpha, q_cut)   alpha = alpha_s(mZ) for JKT, the fixed alpha for Coulomb
SETTINGS = ([(f"jkt_a{a:.3f}", "jkt", a, 0.050) for a in (0.112, 0.118, 0.120, 0.125)]
            + [(f"jkt_qc{qc * 1e3:.0f}", "jkt", 0.120, qc) for qc in (0.010, 0.020, 0.100, 0.200)]
            + [(f"coul_a{a:.2f}", "coulomb", a, None) for a in (0.12, 0.15)])

E_SMALL = np.array([-5., -3., -2., -1., 0., 1., 2., 4.])
P_SMALL = np.array([0.5, 1., 2., 5., 10., 20., 30., 50.])
E_DENSE = np.round(np.arange(-6.0, 4.01, 0.5), 2)
P_DENSE = np.round(np.arange(0.5, 30.01, 1.0), 2)

HEADER = [
    "Toponium S-wave Green's function, JKT (Z. Phys. C 56 (1992) 653) Secs. 2-3. Colour singlet.",
    f"Units: GeV. m_t = {M_T}, Gamma_t = {GAMMA_T}, C_F = {C_F:.6f}, n_F = {N_F}, mZ = {MZ}.",
    "JKT potential: Lambda_MSbar from (14) at Q = mZ; (13) with the NLO bracket (31/3 - 10 n_F/9) alpha/4pi;",
    "  (15) with Lambda_R = 0.4 GeV, switch at |p| = 5 GeV, C from continuity there;",
    "  (19) below q_cut with C1 from continuity at q_cut; the d3 term V0_tilde fixed by (20): V_cut(r = 1/GeV) = -1/4 GeV.",
    "Position-space convention: V(r) = Int d^3q/(2pi)^3 e^{-iq.r} V(q), so the d3 term adds V0_tilde/(2pi)^3 = shift to V(r).",
    "Coulomb rows: V = -4 pi C_F alpha/Q^2, alpha fixed, no bracket, no cut, no shift.",
    "Sign of G: JKT (6), G0 = 1/(E + i Gamma_t - p^2/m_t); G -> G0 as alpha -> 0.",
    "G at energy E includes the d3 term, i.e. it is solved with G0 at E - shift. G0 columns are at the physical E.",
    f"Solver: JKT (21) on Gauss-Legendre nodes (23), {N_PER_PANEL} per panel on (0,4,16,60,m_t) + t=1/q tail (N = {5 * N_PER_PANEL}).",
]


def potential(kind, alpha, q_cut):
    return JKTPotential(alpha_mz=alpha, q_cut=q_cut) if kind == "jkt" else JKTPotential.coulomb(alpha)


def write_csv(name, extra_header, columns, rows):
    with open(f"{OUT}/{name}", "w", newline="") as f:
        for line in HEADER + extra_header:
            f.write(f"# {line}\n")
        w = csv.writer(f)
        w.writerow(columns)
        for r in rows:
            w.writerow([x if isinstance(x, str) else f"{x:.10g}" for x in r])


def params_row(sid, kind, alpha, pot):
    """One row of params.csv. Coulomb rows have no Lambda, C, C1 or shift."""

    if kind != "jkt":
        nan = float("nan")
        return [sid, kind, alpha, nan, nan, nan, nan, nan, nan, nan, nan]
    Vr1 = JKTPotential.V_cut_position(pot.V, 1.0, pot.q_cut, pot.kinks) + pot.shift
    a_check = JKTPotential.alpha_2loop(np.array(MZ ** 2), pot.Lam)
    return [sid, kind, alpha, pot.q_cut, pot.Lam, pot.C, pot.C1, pot.shift * (2 * np.pi) ** 3, pot.shift,
            float(a_check), Vr1]


def q_grid(pot):
    """Log-spaced 1 MeV .. 1 TeV, plus the points at and either side of each kink and mZ."""

    special = [MZ] + ([pot.q_cut, pot.p_match] if pot.q_cut else [])
    near = [k * f for k in special for f in (1 - 1e-6, 1.0, 1 + 1e-6)]
    return np.unique(np.concatenate([np.geomspace(1e-3, 1e3, 121), near]))


def potential_rows(sid, kind, alpha, pot):
    """(14), (13), (15), (19) evaluated at each Q. (14) and (13) are nan where ln(Q^2/Lambda^2) <= 1
    (Q <~ 1.65 Lambda): there the two-loop formula has no meaning (and (15) does not use it)."""

    Q = q_grid(pot)
    Q2 = Q ** 2
    V_coul = -4 * np.pi * C_F * alpha / Q2
    if kind != "jkt":
        nan = np.full_like(Q, np.nan)
        cols = (Q, Q2, nan, np.full_like(Q, alpha), nan, nan, nan, nan, V_coul)
    else:
        valid = np.log(Q2 / pot.Lam ** 2) > 1
        a2 = np.where(valid, JKTPotential.alpha_2loop(Q2, pot.Lam), np.nan)
        Vp = np.where(valid, JKTPotential.V_pert(Q2, pot.Lam), np.nan)
        Vlog = JKTPotential.V_ir(Q2, pot.Lambda_R, pot.C)
        Vj = JKTPotential.V_jkt(Q2, pot.Lam, pot.Lambda_R, pot.C, pot.p_match)
        Vc = pot.V(Q2)
        cols = (Q, Q2, a2, np.full_like(Q, alpha), Vp, Vlog, Vj, Vc, V_coul)
    return [[sid] + list(r) for r in zip(*cols)]


def G_rows(sid, kind, alpha, solver, E_vals, p_vals):
    G = solver.table(E_vals, p_vals)
    EE, PP = np.meshgrid(E_vals, p_vals, indexing="ij")
    G0 = Constants.G0(EE, PP)
    Gx = CoulombExact.G(EE, PP, alpha) if kind == "coulomb" else np.full(G.shape, complex(np.nan, np.nan))
    rows = [[sid, e, p, g.real, g.imag, g0.real, g0.imag, abs(g / g0) ** 2, gx.real, gx.imag]
            for e, p, g, g0, gx in zip(EE.ravel(), PP.ravel(), G.ravel(), G0.ravel(), Gx.ravel())]
    return rows, G, Gx


def main():
    os.makedirs(OUT, exist_ok=True)
    lines = ["STAGE 5: tables for comparison with an independent solver (colour singlet)", ""]
    log = lambda s="": (print(s), lines.append(s))

    params, pot_rows, small, dense, G_ref = [], [], [], [], {}
    for sid, kind, alpha, q_cut in SETTINGS:
        print(f"  {sid} ...", flush=True)
        pot = potential(kind, alpha, q_cut)
        solver = LSSolver(pot, n=N_PER_PANEL)
        params.append(params_row(sid, kind, alpha, pot))
        pot_rows += potential_rows(sid, kind, alpha, pot)
        rows, _, _ = G_rows(sid, kind, alpha, solver, E_SMALL, P_SMALL)
        small += rows
        rows, G, Gx = G_rows(sid, kind, alpha, solver, E_DENSE, P_DENSE)
        dense += rows
        G_ref[sid] = (G, Gx)

    p_cols = ["id", "potential", "alpha", "q_cut", "Lambda_MSbar", "C", "C1", "V0_tilde", "shift",
              "alpha_2loop_at_mZ", "Vcut_r1"]
    write_csv("params.csv", ["alpha = alpha_s(mZ) for jkt rows, the fixed alpha for coulomb rows.",
                             "shift = V0_tilde/(2pi)^3. Vcut_r1 = V_cut(r = 1/GeV) including the d3 term: -0.25 by (20)."],
              p_cols, params)
    write_csv("potential.csv", ["V in GeV^-2 (momentum space), Q in GeV. V_cut excludes the d3 term (it is the shift).",
                                "alpha_2loop, V_pert = nan where ln(Q^2/Lambda^2) <= 1 (two-loop formula meaningless).",
                                "V_log = the |p| < 5 GeV line of (15) evaluated at all Q; V_jkt switches at 5 GeV."],
              ["id", "Q", "Q2", "alpha_2loop", "alpha_fixed", "V_pert", "V_log", "V_jkt", "V_cut", "V_coulomb"],
              pot_rows)
    g_cols = ["id", "E", "p", "ReG", "ImG", "ReG0", "ImG0", "absG_over_G0_sq", "ReG_exact", "ImG_exact"]
    g_note = ["G in GeV^-1, E and p in GeV. G_exact: exact Coulomb solution (coulomb rows only)."]
    write_csv("G_small.csv", g_note, g_cols, small)
    write_csv("G_dense.csv", g_note, g_cols, dense)

    log("Parameters (GeV):")
    log(f"  {'id':12s} {'alpha':>6s} {'q_cut':>6s} {'Lambda':>8s} {'C':>9s} {'C1':>10s} {'V0_tilde':>10s}"
        f" {'shift':>8s} {'a(mZ)':>8s} {'Vcut(r=1)':>10s}")
    for r in params:
        if r[1] != "jkt":
            log(f"  {r[0]:12s} {r[2]:6.3f}   (fixed-alpha Coulomb: no Lambda, C, C1, cut or shift)")
            continue
        log(f"  {r[0]:12s} {r[2]:6.3f} {r[3]:6.3f} {r[4]:8.5f} {r[5]:9.5f} {r[6]:10.3e} {r[7]:10.3f}"
            f" {r[8]:8.4f} {r[9]:8.5f} {r[10]:10.6f}")

    log("")
    log("Checks:")
    for sid in ("coul_a0.12", "coul_a0.15"):
        G, Gx = G_ref[sid]
        d = np.abs(G / Gx - 1)
        log(f"  {sid}: solver vs exact on the dense grid, median {np.median(d):.1e}, max {np.max(d):.1e}")
    G50 = G_ref["jkt_a0.120"][0]
    for sid in ("jkt_qc10", "jkt_qc20", "jkt_qc100", "jkt_qc200"):
        log(f"  {sid}: max |G/G(q_cut = 50 MeV) - 1| = {np.max(np.abs(G_ref[sid][0] / G50 - 1)):.1e}")

    log("")
    log("Files: params.csv (one row per setting), potential.csv (V vs Q per setting),")
    log(f"  G_small.csv ({len(E_SMALL)} E x {len(P_SMALL)} p per setting), G_dense.csv ({len(E_DENSE)} E x {len(P_DENSE)} p).")
    open(f"{OUT}/summary.txt", "w").write("\n".join(lines) + "\n")
    print(f"\nSaved {OUT}/")


if __name__ == "__main__":
    main()
