#!/usr/bin/env python3
"""STAGE 4 -- Fuks Figs. 2-5 from MadGraph events.  Uses: all six classes of toponium.py.

Reproduce Fuks et al. (arXiv:2411.18962v2, reference/2411.18962v2.pdf) Figures 2-5:
the SAME events and the SAME re-weighting (Events.weights, Fuks (15)), with the
Green's function taken from three sources:

  "Fuks: swData via CALCGREEN"   what Fuks et al. did: their tables through their own
                                 interpolation routine (FuksTables)
  "ours: Coulomb exact 0.15"     CoulombExact -- the potential swData turns out
                                 to contain (stage 2)
  "ours: JKT'92 LS"              LSSolver with JKTPotential
                                 (alpha_s(mZ) = 0.12, NLO bracket; Fuks cite JKT, footnote 1)

=====================================================================================
EQUATION MAP  (Fuks equation / code  ->  where it lives)
=====================================================================================
  Fuks (14), (21)  process and p* (Lorentz-invariant form)       Events.kinematics
  Fuks (11)        G0                                            Constants.G0
  Fuks (15)        |M|^2 -> |M|^2 |G/G0|^2                       Events.weights
  Fuks Sec. 3      CALCGREEN on swData                           FuksTables.calcgreen_ratio
  Fuks (22)  <p(E)> = Int d^3p p f / Int d^3p f,  f = d^2 sigma/(p^2 dp dE)
             = dsigma-weighted mean of p* at fixed E                      main ("<p*>(E)")
  Fuks Fig. 2  "d^2 sigma / p*^2 per bin" [fb GeV^-2], E bins 0.4, p* bins 1 GeV
  Fuks Fig. 3  dsigma per 0.5 GeV bin in W = 330..350, p* < 50 (solid) / 100 (dashed)
  Fuks Fig. 4  dsigma/sigma per 1 GeV bin of p* (p* < 50)
  Fuks Fig. 5  m_tL, m_tH for -4 <= E <= 0, p* <= 50 (per GeV; their label says bin^-1)

Events (13 TeV, CT18NLO, m_t = 173, Gamma_t = 1.49, MG5 default lepton cuts), read from
results/4_fuks_figs2to5/_cache_*.npz (the samples in ../output/ were deleted; see README):
  singlet330  colour-singlet (matrix1.f CF = 2/3, Fuks (20)), 330 <= W <= 350 GeV
  full        full colour, 340 <= W <= 350; octet = full - singlet

Run from toponium_study/:  python3 4_fuks_figs2to5.py        Output -> results/4_fuks_figs2to5/
"""

import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.interpolate import RegularGridInterpolator

from toponium import Constants, JKTPotential, LSSolver, CoulombExact, FuksTables, Events


LEPTON_FACTOR = 4.0
SAMPLES = {
    # singlet: 330 <= W <= 350 GeV (Fuks Fig. 3 range; Figs. 2, 5 cut to their own windows)
    "singlet330": "../output/GGtoTTbar_Singlet_Fuks13TeV_W330/Events/run_01/unweighted_events.lhe.gz",
    "full": "../output/GGtoTTbar_Full_Fuks13TeV/Events/run_01/unweighted_events.lhe.gz",
}
OUT = "results/4_fuks_figs2to5"


def solver_source(pot, E_rng=(-16.2, 4.2), p_max=100.0):
    """Our solver as a source G(E, p): tabulated once on the swData grid spacing, then
    interpolated (one linear solve per tabulated energy instead of one per event)."""

    Eg = np.arange(E_rng[0], E_rng[1] + 1e-9, 0.02)
    pg = np.arange(0.01, p_max + 0.26, 0.25)
    T = LSSolver(pot).table(Eg, pg)
    re = RegularGridInterpolator((Eg, pg), T.real)
    im = RegularGridInterpolator((Eg, pg), T.imag)
    return lambda E, p: re(np.column_stack([E, np.clip(p, 0.01, None)])) + 1j * im(
        np.column_stack([E, np.clip(p, 0.01, None)]))


FUKS, COUL, JKT = "Fuks: swData via CALCGREEN", "ours: Coulomb exact 0.15", "ours: JKT'92 LS"


def main():
    os.makedirs(OUT, exist_ok=True)
    ev = {k: Events.kinematics(p, f"{OUT}/_cache_{k}.npz") for k, p in SAMPLES.items()}
    s = ev["singlet330"]
    for k, e in ev.items():
        # run_card event_norm = average: sigma is the MEAN event weight, so each event
        # carries sigma/N; pb -> fb; x4 lepton-flavour combinations
        e["w"] = e["w"] / len(e["w"]) * LEPTON_FACTOR * 1e3

    # The ONLY difference between the curves: where G(E, p) comes from.
    ratio = FuksTables.calcgreen_ratio()
    sources = {
        FUKS: lambda E, p: ratio(E, p) * Constants.G0(E, p),               # CALCGREEN returns G/G0
        COUL: lambda E, p: CoulombExact.G(E, p, 0.15),
        JKT: solver_source(JKTPotential(0.12)),
    }
    rw = {k: Events.weights(G, s["E"], s["pstar"]) for k, G in sources.items()}   # Fuks (15)
    colors = {FUKS: "tab:green", COUL: "k", JKT: "tab:purple"}
    # swData and our Coulomb agree to <1%: draw swData as a wide translucent band under the others
    style = {FUKS: dict(lw=4, alpha=0.45)}

    lines = []

    def log(x=""):
        print(x, flush=True)
        lines.append(x)

    inwin = (s["W"] >= 340) & (s["W"] <= 350)
    log(f"singlet: {len(s['w'])} events in 330<=W<=350, {inwin.mean():.1%} of the rate in 340<=W<=350")
    log(f"full colour: {len(ev['full']['w'])} events in 340<=W<=350")
    sig_s, sig_f = s["w"][inwin].sum(), ev["full"]["w"].sum()
    log(f"sigma(singlet) = {sig_s:.3f} fb, sigma(full colour) = {sig_f:.3f} fb, singlet fraction {sig_s / sig_f:.3f}")
    c50 = inwin & (s["pstar"] < 50)
    log(f"p*<50 keeps {s['w'][c50].sum() / sig_s:.4f} of the singlet rate (340<=W<=350)")
    log("")
    log("Green's-function enhancement sigma_rw/sigma (singlet, W window, p*<50):")
    for k, r in rw.items():
        log(f"   {k:28s} {np.sum(s['w'][c50] * r[c50]) / s['w'][c50].sum():6.3f}")
    d = np.abs(rw[COUL] / rw[FUKS] - 1)[c50]
    log(f"per-event |w_ours(Coulomb 0.15)/w_Fuks - 1|: median {np.median(d):.1e}, max {d.max():.1e}")

    # ---------------- Figure 2 ----------------
    # As in Fuks's figure: "d^2 sigma / p*^2 per bin" [fb GeV^-2] -- the cross section in
    # each (E, p*) bin divided by p*^2 (per event), NOT divided by the bin widths.
    # Binning read off their figure: E in [-6,4] with 0.4 GeV bins, p* in [2,30] with
    # 1 GeV bins; common log colour scale 5e-5 .. 1e-2.
    Eb, pb = np.linspace(-6, 4, 26), np.linspace(2, 30, 29)

    def h2(wt):
        H, _, _ = np.histogram2d(s["E"][c50], s["pstar"][c50], bins=[Eb, pb],
                                 weights=(wt / s["pstar"] ** 2)[c50])
        return H

    panels = [("unweighted", h2(s["w"]))] + [(k, h2(s["w"] * r)) for k, r in rw.items()]
    norm2 = LogNorm(5e-5, 1e-2)
    fig, ax = plt.subplots(1, 4, figsize=(24, 5.6), layout="constrained")
    for a, (k, P) in zip(ax, panels):
        m = a.pcolormesh(Eb, pb, np.where(P > 0, P, np.nan).T, cmap="gist_ncar", norm=norm2)
        a.set(title=("Fig. 2 " + ("without" if k == "unweighted" else "with") + " G re-weighting\n" + k),
              xlabel="E [GeV]", ylabel="p* [GeV]")
    fig.colorbar(m, ax=ax, label="d²σ/p*² per bin  [fb GeV⁻²]")
    fig.savefig(f"{OUT}/fig2.png", dpi=120)
    plt.close(fig)
    log("")
    log("Fig. 2 (d2sigma/p*^2 per bin, fb GeV^-2; Fuks colour scale 5e-5..1e-2):")
    for k, P in panels:
        i0, j0 = np.unravel_index(P.argmax(), P.shape)
        log(f"   {k:28s} max {P.max():.2e} at E={0.5 * (Eb[i0] + Eb[i0 + 1]):+.1f}, "
            f"p*={0.5 * (pb[j0] + pb[j0 + 1]):.1f}; cells < 5e-5: {np.mean(P < 5e-5):.0%}")

    log("")
    log("<p*>(E) with re-weighting (Fuks quote ~20 GeV at E=-2 GeV):")
    band = c50 & (np.abs(s["E"] + 2) < 0.25)
    for k, r in [("unweighted", np.ones_like(s["w"]))] + list(rw.items()):
        wt = s["w"][band] * r[band]
        log(f"   {k:28s} <p*>(E=-2+-0.25) = {np.sum(wt * s['pstar'][band]) / wt.sum():5.1f} GeV")

    # ---------------- Figure 3 ----------------
    # As in the paper: W in [330, 350] GeV, 0.5 GeV bins, d(sigma) per bin [fb], log scale;
    # solid p* < 50, dashed p* < 100 GeV.
    Wb = np.arange(330, 350 + 1e-9, 0.5)
    Wc = 0.5 * (Wb[1:] + Wb[:-1])
    fig, ax = plt.subplots(1, 2, figsize=(16, 6), layout="constrained")
    log("")
    log("Fig. 3, d(sigma) per 0.5 GeV bin [fb] (p*<50 | p*<100):")
    hist3 = {}
    for pcut, ls in ((50, "-"), (100, "--")):
        sel = s["pstar"] < pcut
        hist3[("no re-weighting", pcut)] = np.histogram(s["W"][sel], Wb, weights=s["w"][sel])[0]
        for k, r in rw.items():
            hist3[(k, pcut)] = np.histogram(s["W"][sel], Wb, weights=(s["w"] * r)[sel])[0]
        for k in ["no re-weighting"] + list(rw):
            ax[0].step(Wc, hist3[(k, pcut)], ls, where="mid", color=colors.get(k, "tab:red"),
                       label=f"{k}, p*<{pcut}", **style.get(k, dict(lw=1.2)))
    for W0 in (330.25, 335.25, 340.25, 344.75, 349.75):
        i0 = np.argmin(abs(Wc - W0))
        log(f"   W={Wc[i0]:6.2f}: " + "  ".join(f"{k.split(' via')[0].replace('ours: ', '')}: "
                                          f"{hist3[(k, 50)][i0]:.3g}|{hist3[(k, 100)][i0]:.3g}"
                                          for k in ["no re-weighting"] + list(rw)))
    for k in rw:
        h = hist3[(k, 50)]
        log(f"   peak {k:28s} W = {Wc[h.argmax()]:.2f} GeV (E = {Wc[h.argmax()] - 2 * Constants.M_T:+.2f})")
        ax[1].step(Wc, h / hist3[(FUKS, 50)], where="mid", color=colors[k], label=k,
                   **style.get(k, dict(lw=1.2)))
    ax[0].set(title="Fig. 3: colour-singlet dσ/dW", xlabel="W [GeV]", ylabel="dσ [fb/bin]", yscale="log")
    ax[0].legend(fontsize=8)
    ax[1].axhline(1, color="gray", lw=0.5)
    ax[1].set(title="ratio to Fuks (swData via CALCGREEN), p*<50", xlabel="W [GeV]", ylim=(0, 2.2))
    ax[1].legend()
    fig.savefig(f"{OUT}/fig3.png", dpi=120)
    plt.close(fig)

    # ---------------- Figure 4 ----------------
    # As in the paper: p* < 50 (generator level), 1 GeV bins, dsigma/sigma per bin, log scale.
    # All generated W (330-350 GeV); Fuks do not state the W range of this figure.
    pb4 = np.arange(0, 50 + 1e-9, 1.0)
    fig, a = plt.subplots(figsize=(8, 6), layout="constrained")
    log("")
    log("Fig. 4, dsigma/sigma per 1 GeV bin, p*<50:")
    sel = s["pstar"] < 50
    for k, r in [("no re-weighting", np.ones_like(s["w"]))] + list(rw.items()):
        h = np.histogram(s["pstar"][sel], pb4, weights=(s["w"] * r)[sel])[0]
        h = h / h.sum()
        a.step(0.5 * (pb4[1:] + pb4[:-1]), h, where="mid", color=colors.get(k, "tab:red"), label=k,
               **style.get(k, dict(lw=1.2)))
        log(f"   {k:28s} peak {h.max():.4f} at p*={0.5 + h.argmax():.1f}; p*=4.5: {h[4]:.4f}; "
            f"p*=49.5: {h[-1]:.4f}")
    a.set(title="Fig. 4: top recoil momentum (p* < 50 GeV)", xlabel="p* [GeV]", ylabel="dσ/σ [bin⁻¹]",
          yscale="log", ylim=(5e-4, 1e-1))
    a.legend()
    fig.savefig(f"{OUT}/fig4.png", dpi=120)
    plt.close(fig)

    # ---------------- Figure 5 ----------------
    # As in the paper: -4 <= E <= 0, p* <= 50; 155-180 GeV in 0.25 GeV bins. Their y values
    # (peaks ~0.19 and ~0.54) are 1/sigma dsigma/dm per GeV, although labelled "bin^-1".
    f = ev["full"]
    cut5 = lambda e: (e["E"] >= -4) & (e["E"] <= 0) & (e["pstar"] <= 50)
    mb = np.arange(155, 180 + 1e-9, 0.25)
    mc = 0.5 * (mb[1:] + mb[:-1])
    fig, ax = plt.subplots(1, 2, figsize=(16, 6), layout="constrained")
    log("")
    log("Fig. 5 (-4<=E<=0, p*<=50): 1/sigma dsigma/dm peak [1/GeV] @ m, and mean mass [GeV]:")
    for a, var in zip(ax, ("mtL", "mtH")):
        hs = np.histogram(s[var][cut5(s)], mb, weights=s["w"][cut5(s)])[0]
        hf = np.histogram(f[var][cut5(f)], mb, weights=f["w"][cut5(f)])[0]
        curves = {"octet (full − singlet)": (hf - hs, "tab:blue", {}),
                  "singlet, no re-weighting": (hs, "tab:red", {})}
        for k, r in rw.items():
            curves[f"singlet, {k}"] = (np.histogram(s[var][cut5(s)], mb, weights=(s["w"] * r)[cut5(s)])[0],
                                       colors[k], style.get(k, {}))
        for k, (h, c, st) in curves.items():
            d = h / h.sum() / 0.25
            a.step(mc, d, where="mid", color=c, label=k, **(st or dict(lw=1.2)))
            log(f"   {var} {k:40s} peak {d.max():.3f} @ {mc[d.argmax()]:.2f}   mean {np.sum(h * mc) / h.sum():.2f}")
        a.set(title=f"Fig. 5: {var}  (−4 ≤ E ≤ 0 GeV, p* ≤ 50 GeV)", xlabel=f"{var} [GeV]",
              ylabel="1/σ dσ/dm [GeV⁻¹]", xlim=(155, 180))
        a.legend(fontsize=8)
    fig.savefig(f"{OUT}/fig5.png", dpi=120)
    plt.close(fig)

    open(f"{OUT}/summary.txt", "w").write("\n".join(lines) + "\n")
    print(f"\nSaved {OUT}/fig2..5.png and summary.txt")


if __name__ == "__main__":
    main()
