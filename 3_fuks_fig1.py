#!/usr/bin/env python3
"""STAGE 3 -- Fuks Fig. 1.  Uses: toponium.FuksTables (plotted directly); CoulombExact, LSSolver for comparison.  No events.

Reproduce Fuks et al. (arXiv:2411.18962) Figure 1 from their own swData
tables, and compare it with this project's independent solutions.

=====================================================================================
EQUATION MAP  (paper equation  ->  where it lives here)
=====================================================================================
  Fuks Fig. 1 left   |G(E,p)|^2 from swData (columns: E, p, Re, Im of -G)   main
  Fuks (11)          G0(E, p) = 1/(E + i Gamma_t - p^2/m_t): Fig. 1 right   main
  Fuks (10) / JKT (5) solved two ways for comparison:
                     exact fixed-alpha Coulomb (CoulombExact) and
                     our solver (LSSolver)                      main
  JKT (13)-(15), (19)-(20)  running potential with the NLO bracket of JKT (13)
                     (JKTPotential); Fuks cite JKT for their
                     potential (footnote 1)                                 main

Run from toponium_study/:
    python3 3_fuks_fig1.py                   # alpha = 0.15 for the Coulomb comparison
    python3 3_fuks_fig1.py --alpha 0.12      # any other fixed coupling

swData is read by FuksTables, which returns G~ in the sign convention of the
paper's Eq. 11 (the files store -G~); |G~|^2, which Figure 1 shows, is unaffected.

Outputs -> results/3_fuks_fig1/
  fig1_from_swdata.png   Figure 1: |G~|^2 from swData (left), |G~0|^2 from Eq. 11 (right)
  fig1_ours_vs_swdata.png  swData vs exact fixed-alpha Coulomb (CoulombExact), our
                           solver at the same alpha (LSSolver), and the
                           JKT'92 running potential (NLO bracket), which Fuks cite (footnote 1)
"""

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from toponium import Constants, JKTPotential, LSSolver, CoulombExact, FuksTables


E_LO, E_HI, P_HI = -6.0, 4.0, 30.0
OUT = "results/3_fuks_fig1"


def in_window(E, p, G):
    """Restrict one swData block to the Figure 1 window."""

    e_sel = (E >= E_LO - 1e-9) & (E <= E_HI + 1e-9)
    p_sel = p <= P_HI
    return E[e_sel], p[p_sel], G[np.ix_(e_sel, p_sel)]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--alpha", type=float, default=0.15,
                    help="fixed coupling for the Coulomb comparison (swData fits 0.150; "
                         "the paper's caption says 0.12)")
    ap.add_argument("--alpha_mz", type=float, default=0.12, help="alpha_s(mZ) for the JKT'92 potential")
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)

    blocks = [in_window(*b) for b in FuksTables.blocks()]   # E < -5 (dE = 0.1), E >= -5 (dE = 0.02)
    (Ea_k, pa, _), (Eb, pb, _) = blocks
    # explicit cell edges so the two blocks (dE = 0.1 and 0.02) tile without a gap
    edges_b = np.concatenate([[Eb[0] - 0.01], 0.5 * (Eb[1:] + Eb[:-1]), [Eb[-1] + 0.01]])
    edges_a = np.concatenate([[Ea_k[0] - 0.05], 0.5 * (Ea_k[1:] + Ea_k[:-1]), [edges_b[0]]])
    p_edges = np.concatenate([[0.0], 0.5 * (pa[1:] + pa[:-1]), [pa[-1] + 0.125]])
    E_edges = [edges_a, edges_b]

    G2_all = np.concatenate([np.abs(G).ravel() ** 2 for _, _, G in blocks])
    EE_all = np.concatenate([np.meshgrid(E, p, indexing="ij")[0].ravel() for E, p, _ in blocks])
    PP_all = np.concatenate([np.meshgrid(E, p, indexing="ij")[1].ravel() for E, p, _ in blocks])
    G02_all = np.abs(Constants.G0(EE_all, PP_all)) ** 2
    norm = LogNorm(vmin=min(G2_all.min(), G02_all.min()), vmax=max(G2_all.max(), G02_all.max()))

    # ---------------- Figure 1 from swData ----------------
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for (E, p, G), Ee in zip(blocks, E_edges):
        EE, PP = np.meshgrid(E, p, indexing="ij")
        m = axes[0].pcolormesh(Ee, p_edges, (np.abs(G) ** 2).T, shading="flat", cmap="gist_ncar", norm=norm)
        axes[1].pcolormesh(Ee, p_edges, (np.abs(Constants.G0(EE, PP)) ** 2).T, shading="flat", cmap="gist_ncar", norm=norm)
    for ax, title in zip(axes, [r"$|\tilde G(E,p)|^2$ [GeV$^{-2}$]  (swData)",
                                r"$|\tilde G_0(E,p)|^2$ [GeV$^{-2}$]  (Eq. 11)"]):
        ax.set_title(title)
        ax.set_xlabel("E [GeV]")
        ax.set_ylabel("p [GeV]")
        ax.set_xlim(E_LO, E_HI)
        ax.set_ylim(0, P_HI)
        fig.colorbar(m, ax=ax)
    fig.tight_layout()
    fig.savefig(f"{OUT}/fig1_from_swdata.png", dpi=150)
    plt.close(fig)

    # ---------------- ours vs swData ----------------
    # comparison grid: every 5th E of the fine block, every 2nd p (0.5 GeV)
    E = np.concatenate([blocks[0][0], blocks[1][0][::5]])
    Gsw = np.vstack([blocks[0][2], blocks[1][2][::5]])[:, ::2]
    p = pa[::2]
    EE, PP = np.meshgrid(E, p, indexing="ij")

    models = [
        (f"exact fixed Coulomb, α={args.alpha}", CoulombExact.G(EE, PP, args.alpha, m=Constants.M_T, Gam=Constants.GAMMA_T)),
        (f"our LS solver, Coulomb α={args.alpha}", LSSolver(JKTPotential.coulomb(args.alpha)).table(E, p)),
        (f"JKT'92 running + NLO bracket, αs(mZ)={args.alpha_mz}",
         LSSolver(JKTPotential(args.alpha_mz)).table(E, p)),
    ]

    fig, axes = plt.subplots(2, 4, figsize=(24, 10.5), layout="constrained")
    axes[0, 0].pcolormesh(E, p, (np.abs(Gsw) ** 2).T, shading="nearest", cmap="gist_ncar", norm=norm)
    axes[0, 0].set_title("swData |G̃|²")
    axes[1, 0].pcolormesh(E, p, np.angle(Gsw).T, shading="nearest", cmap="twilight")
    axes[1, 0].set_title("swData arg G̃  (phase, sign as Eq. 11)")
    lines = [f"swData vs models over E in [{E_LO},{E_HI}], p in [0,{P_HI}] GeV:"]
    for k, (label, Gm) in enumerate(models, start=1):
        m = axes[0, k].pcolormesh(E, p, (np.abs(Gm) ** 2).T, shading="nearest", cmap="gist_ncar", norm=norm)
        axes[0, k].set_title(label + "\n|G|²")
        dev = np.abs(Gsw / Gm - 1)
        r = axes[1, k].pcolormesh(E, p, dev.T, shading="nearest", cmap="viridis", norm=LogNorm(1e-4, 3))
        axes[1, k].set_title(f"|swData/model − 1|  (median {np.median(dev):.1e})")
        fig.colorbar(r, ax=axes[1, k])
        lines.append(f"  {label:48s} median {np.median(dev):.2e}, max {dev.max():.2e} (complex G)")
    fig.colorbar(m, ax=axes[0, -1], label="|G|² [GeV⁻²]")
    for ax in axes.flat:
        ax.set_xlabel("E [GeV]")
        ax.set_ylabel("p [GeV]")
    fig.savefig(f"{OUT}/fig1_ours_vs_swdata.png", dpi=110)
    open(f"{OUT}/summary.txt", "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"Saved {OUT}/fig1_from_swdata.png, {OUT}/fig1_ours_vs_swdata.png, {OUT}/summary.txt")


if __name__ == "__main__":
    main()
