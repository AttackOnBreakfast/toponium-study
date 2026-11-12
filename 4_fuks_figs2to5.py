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

Second event set, Fuks's process exactly (sample "fuksproc_ct25"):
  generate g g > t t~ > w+ b w- b~, w+ > l+ vl, w- > l- vl~   (l = e, mu: all 4 lepton combinations)
  13 TeV, CT25NNLO (local LHAPDF id, see 6_tfactors.py), MadGraph default cuts, 1M events per sample over
  full phase space (no W window), one seed each: colour singlet (CF = 2/3, Fuks (19)-(20)) and full colour.
  The figures cut the events to their W windows; no lepton factor (all combinations generated).
    uv run --with six python 4_fuks_figs2to5.py --generate             # MadGraph (both samples, ~40 min)
    uv run python 4_fuks_figs2to5.py --sample fuksproc_ct25            # Output -> results/4_fuks_figs2to5/fuksproc_ct25/
    uv run python 4_fuks_figs2to5.py --octet      # Figs. 1-5 with the octet in the signal -> .../fuksproc_ct25_octet/
"""

import argparse
import glob
import importlib
import re
import shutil
import subprocess
import sys

import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.interpolate import RegularGridInterpolator

from toponium import Constants, JKTPotential, LSSolver, CoulombExact, FuksTables, Events


VARIANTS = {
    # original reproduction (Sept 23): e+ mu- only (x4 for l, l' = e, mu), CT18NLO, generated inside the W windows
    "orig": dict(out="results/4_fuks_figs2to5", lepton_factor=4.0, windowed=False, samples={
        # singlet: 330 <= W <= 350 GeV (Fuks Fig. 3 range; Figs. 2, 5 cut to their own windows)
        "singlet330": "../output/GGtoTTbar_Singlet_Fuks13TeV_W330/Events/run_01/unweighted_events.lhe.gz",
        "full": "../output/GGtoTTbar_Full_Fuks13TeV/Events/run_01/unweighted_events.lhe.gz"}),
    # Fuks's exact process (all 4 lepton combinations), CT25NNLO, full phase space: cut to the windows here
    "fuksproc_ct25": dict(out="results/4_fuks_figs2to5/fuksproc_ct25", lepton_factor=1.0, windowed=True, samples={
        "singlet330": "../output/fuks_process_13TeV/TT_fuks_singlet/Events/singlet_s7301/unweighted_events.lhe.gz",
        "full": "../output/fuks_process_13TeV/TT_fuks_full/Events/full_s7302/unweighted_events.lhe.gz"}),
}
LEPTON_FACTOR, SAMPLES, OUT = 4.0, VARIANTS["orig"]["samples"], VARIANTS["orig"]["out"]
GEN_DIR = "../output/fuks_process_13TeV"
GEN = {"singlet": ("TT_fuks_singlet", 7301), "full": ("TT_fuks_full", 7302)}   # sample: (MadGraph dir, seed)
FUKS_PROCESS = "generate g g > t t~ > w+ b w- b~, w+ > l+ vl, w- > l- vl~"     # Fuks Sec. 3, verbatim


def singlet_colour(path):
    """Fuks (19)-(20): replace the 2x2 gg colour matrix CF by the singlet projection (all entries 2/3)."""

    src = open(path).read()
    if "colour-singlet projection" in src:
        return
    shutil.copy(path, path + ".orig")
    new, n = re.subn(r"(      DATA \(CF\(I,  1\),I=  1,  2\) /)5\.333333333333333D\+00,\n(     \$ )-6\.666666666666666D-01/",
                     r"C     colour-singlet projection, Fuks (19)-(20)\n\g<1>6.666666666666666D-01,\n\g<2>6.666666666666666D-01/", src)
    new, m = re.subn(r"(      DATA \(CF\(I,  2\),I=  1,  2\) /)-6\.666666666666666D-01\n(     \$ ),5\.333333333333333D\+00/",
                     r"\g<1>6.666666666666666D-01\n\g<2>,6.666666666666666D-01/", new)
    assert n == 1 and m == 1 and "5.333333333333333D+00" not in new, f"colour matrix not found in {path}"
    open(path, "w").write(new)


def generate():
    """Fuks's process verbatim, 1M events per sample over full phase space, one seed each; per-event kinematics
    cached to the variant's _cache_*.npz, then the LHE files are deleted."""

    tf = importlib.import_module("6_tfactors")
    var = VARIANTS["fuksproc_ct25"]
    os.makedirs(var["out"], exist_ok=True)
    for key, (name, seed) in GEN.items():
        d = os.path.abspath(os.path.join(GEN_DIR, name))
        if not os.path.isdir(d):
            os.makedirs(os.path.dirname(d), exist_ok=True)
            script = os.path.join(os.path.dirname(d), f"_{name}.mg5")
            open(script, "w").write(f"set automatic_html_opening False\nimport model sm\n{FUKS_PROCESS}\n"
                                    f"output {d} --hel_recycling=False\n")
            subprocess.run([sys.executable, "-W", "ignore", os.path.join(tf.MG_ROOT, "bin", "mg5_aMC"), script],
                           cwd=tf.MG_ROOT, check=True, env=tf.mg_env())
            os.remove(script)
        if key == "singlet":
            files = glob.glob(os.path.join(d, "SubProcesses", "P*", "matrix*.f"))
            assert files, "no matrix files"
            for f in files:
                singlet_colour(f)
            print(f"singlet colour matrix in {len(files)} matrix file(s)")
        tf.set_card(os.path.join(d, "Cards", "run_card.dat"), {      # MadGraph default cuts kept
            "pdlabel": "lhapdf", "lhaid": tf.CT25_ID, "use_syst": "False", "event_norm": "average",
            "ebeam1": tf.EBEAM, "ebeam2": tf.EBEAM, "nevents": 1000000, "iseed": seed})
        pc = os.path.join(d, "Cards", "param_card.dat")
        src, n = re.subn(r"^DECAY\s+6\s+\S+\s+# WT", "DECAY   6 1.490000e+00 # WT", open(pc).read(), flags=re.M)
        assert n == 1 and re.search(r"^\s+6 1.730000e\+02 # MT", src, re.M)
        open(pc, "w").write(src)
        tf.set_config(os.path.join(d, "Cards", "me5_configuration.txt"), {
            "run_mode": 2, "nb_core": 10, "lhapdf_py3": os.path.join(tf.LHAPDF_DIR, "lhapdf-config"),
            "automatic_html_opening": "False"})
        run = f"{key}_s{seed}"
        cache = os.path.join(var["out"], f"_cache_{'singlet330' if key == 'singlet' else 'full'}.npz")
        if os.path.exists(cache):
            print(f"{run}: done")
            continue
        print(f"{run}: generating 1000000 events (seed {seed}) ...", flush=True)
        os.makedirs(os.path.join(GEN_DIR, "logs"), exist_ok=True)
        with open(os.path.join(GEN_DIR, "logs", f"{run}.log"), "w") as fl:
            r = subprocess.run([sys.executable, "-W", "ignore", "bin/generate_events", run, "-f"], cwd=d,
                               stdout=fl, stderr=subprocess.STDOUT, env=tf.mg_env())
        lhe = os.path.join(d, "Events", run, "unweighted_events.lhe.gz")
        if r.returncode != 0 or not os.path.exists(lhe):
            sys.exit(f"{run}: MadGraph failed, see {GEN_DIR}/logs/{run}.log")
        ev = Events.kinematics(lhe, cache)
        for f in glob.glob(os.path.join(d, "Events", run, "*.lhe*")):
            os.remove(f)
        print(f"{run}: {len(ev['w'])} events, sigma = {ev['w'].mean():.5g} pb", flush=True)


def solver_source(pot, E_rng=(-16.2, 4.2), p_max=100.0):
    """Our solver as a source G(E, p): tabulated once on the swData grid spacing, then
    interpolated (one linear solve per tabulated energy instead of one per event)."""

    Eg = np.arange(E_rng[0], E_rng[1] + 1e-9, 0.02)
    pg = np.arange(0.01, p_max + 0.26, 0.25)
    T = LSSolver(pot).table(Eg, pg)
    re = RegularGridInterpolator((Eg, pg), T.real)
    im = RegularGridInterpolator((Eg, pg), T.imag)
    # p clamped to the table (p* > p_max only occurs in the full-phase-space samples, and every figure cuts p* <= 100)
    return lambda E, p: re(np.column_stack([E, np.clip(p, 0.01, p_max)])) + 1j * im(
        np.column_stack([E, np.clip(p, 0.01, p_max)]))


FUKS, COUL, JKT = "Fuks: swData via CALCGREEN", "ours: Coulomb exact 0.15", "ours: JKT'92 LS"


def main(variant="orig"):
    var = VARIANTS[variant]
    OUT = var["out"]
    os.makedirs(OUT, exist_ok=True)
    ev = {k: Events.kinematics(p, f"{OUT}/_cache_{k}.npz") for k, p in var["samples"].items()}
    for k, e in ev.items():
        # run_card event_norm = average: sigma is the MEAN event weight, so each event
        # carries sigma/N; pb -> fb; x lepton-flavour factor (4 for e+ mu- only, 1 if all combinations generated)
        e["w"] = e["w"] / len(e["w"]) * var["lepton_factor"] * 1e3
    if var["windowed"]:   # full phase space: keep the W windows of the original samples (after normalising)
        for k, (lo, hi) in (("singlet330", (330, 350)), ("full", (340, 350))):
            sel = (ev[k]["W"] >= lo) & (ev[k]["W"] <= hi)
            print(f"{k}: {sel.sum()} of {len(sel)} events in {lo} <= W <= {hi} GeV")
            ev[k] = {c: v[sel] for c, v in ev[k].items()}
    s = ev["singlet330"]

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


def octet_figures(mode="reweight"):
    """Fuks's Figs. 1-5 as they would look with the colour octet in the toponium signal: all of g g -> t tbar
    (Fuks's process; no qqbar), singlet x |G1/G0|^2 and octet x |G8/G0|^2 (repulsive Coulomb, V8 = +alpha/(6r)).
    Events: the fuksproc_ct25 samples (Fuks's exact process, CT25NNLO, full phase space). The octet is
    (full - singlet), so an octet histogram with weight R8 is  H_full(w R8) - H_singlet(w R8)  (R depends on E, p*).
    Curves: LO; Fuks (singlet x swData, octet LO); A = Fuks + octet at alpha = 0.15 (the coupling his tables
    contain); ours (singlet x our JKT G with CT25 alpha_s(Q), octet LO); B = ours + octet at CT25 alpha_s(25 GeV)."""

    tf = importlib.import_module("6_tfactors")
    var = VARIANTS["fuksproc_ct25"]
    OUT8 = "results/4_fuks_figs2to5/fuksproc_ct25_octet" if mode == "reweight" else OCT_OUT
    os.makedirs(OUT8, exist_ok=True)
    if mode == "reweight":
        ev = {k: Events.kinematics(p, f"{var['out']}/_cache_{k}.npz") for k, p in var["samples"].items()}
    else:   # the five MadGraph samples with the Green's functions in MATRIX1 (generate_octet_mg)
        ev = {k: dict(np.load(os.path.join(OCT_OUT, f"_cache_{k}.npz"))) for k in OCT_SAMPLES}
    for k in ev:
        ev[k]["w"] = ev[k]["w"] / len(ev[k]["w"]) * 1e3                     # fb per event (all lepton combinations)
        sel = (ev[k]["W"] >= 330) & (ev[k]["W"] <= 350)
        ev[k] = {c: v[sel] for c, v in ev[k].items()}
    a8_ours = float(tf.ct25_alphas(tf.OCTET_Q ** 2))
    ratio = FuksTables.calcgreen_ratio()
    ours = solver_source(JKTPotential(alpha_mz=0.118, alpha_fn=tf.ct25_alphas))

    def r_octet(e, a):
        out = np.empty(len(e["E"]))
        for c in range(0, len(out), 20000):
            k = slice(c, c + 20000)
            out[k] = np.abs(CoulombExact.G(e["E"][k], np.clip(e["pstar"][k], 0, None), -a / 8)
                            / Constants.G0(e["E"][k], e["pstar"][k])) ** 2
        return out

    if mode == "mg":
        return _octet_plots(tf, ev, None, OUT8, a8_ours, ratio, ours, r_octet, mode)
    s, f = ev["singlet330"], ev["full"]
    R1 = {"fuks": np.abs(ratio(s["E"], s["pstar"])) ** 2,
          "ours": np.abs(ours(s["E"], s["pstar"]) / Constants.G0(s["E"], s["pstar"])) ** 2}
    R8 = {"A": (r_octet(s, 0.15), r_octet(f, 0.15)), "B": (r_octet(s, a8_ours), r_octet(f, a8_ours))}
    LO, FS, A, OS, B = ("LO (no Green's function)", "Fuks: singlet x swData, octet LO",
                        "A: Fuks + octet x |G8/G0|^2 (alpha 0.15)", "ours: singlet x JKT (CT25 alpha_s), octet LO",
                        f"B: ours + octet x |G8/G0|^2 (alpha {a8_ours:.3f})")
    # each curve: list of (sample, per-event weight) whose histograms are summed
    curves = {
        LO: [(f, f["w"])],
        FS: [(f, f["w"]), (s, s["w"] * (R1["fuks"] - 1))],
        A: [(f, f["w"]), (s, s["w"] * (R1["fuks"] - 1)), (f, f["w"] * (R8["A"][1] - 1)), (s, -s["w"] * (R8["A"][0] - 1))],
        OS: [(f, f["w"]), (s, s["w"] * (R1["ours"] - 1))],
        B: [(f, f["w"]), (s, s["w"] * (R1["ours"] - 1)), (f, f["w"] * (R8["B"][1] - 1)), (s, -s["w"] * (R8["B"][0] - 1))],
    }
    return _octet_plots(tf, ev, curves, OUT8, a8_ours, ratio, ours, r_octet, mode)


def _octet_plots(tf, ev, curves, OUT8, a8_ours, ratio, ours, r_octet, mode):
    """Figs. 1-5 for the octet study; curves = {label: [(events, per-event weight), ...]} (reweight mode), or built
    here from the five MadGraph samples (mg mode: one sample per curve, no reweighting)."""

    LO, FS, A, OS, B = ("LO (no Green's function)", "Fuks: singlet x swData, octet LO",
                        "A: Fuks + octet x |G8/G0|^2 (alpha 0.15)", "ours: singlet x JKT (CT25 alpha_s), octet LO",
                        f"B: ours + octet x |G8/G0|^2 (alpha {a8_ours:.3f})")
    if curves is None:
        curves = {lab: [(ev[k], ev[k]["w"])] for lab, k in ((LO, "LO"), (FS, "A1"), (A, "A18"), (OS, "B1"), (B, "B18"))}
    s = f = None
    color = {LO: "tab:red", FS: "tab:green", A: "k", OS: "tab:purple", B: "tab:orange"}
    style = {LO: dict(lw=1.3), FS: dict(lw=4, alpha=0.45), A: dict(lw=1.4), OS: dict(lw=1.3, ls="--"), B: dict(lw=1.4)}

    def hist(var_, bins, cut, curve, extra=None):
        tot = 0
        for e, wt in curves[curve]:
            c = cut(e)
            x = e[var_][c]
            ww = wt[c] if extra is None else wt[c] * extra(e)[c]
            tot = tot + np.histogram(x, bins, weights=ww)[0]
        return tot

    lines = []

    def log(x=""):
        print(x, flush=True)
        lines.append(x)

    log("Fuks's figures with the colour octet in the toponium signal: g g -> t tbar (Fuks's process), full colour")
    if mode == "reweight":
        log("events: fuksproc_ct25 (CT25NNLO, all e/mu, full phase space), octet = full - singlet, reweighted offline;")
        log(f"  in 330 <= W <= 350 GeV: singlet {len(ev['singlet330']['w'])}, full colour {len(ev['full']['w'])} events")
    else:
        log("events: five MadGraph samples with the Green's functions in MATRIX1 (Fuks's process, CT25NNLO, all e/mu,")
        log("  default cuts, generated in 326 <= W <= 350 GeV, 500k each); in 330 <= W <= 350 GeV:")
        for k in OCT_SAMPLES:
            log(f"  {k:4s} seed {OCT_SAMPLES[k][0]}: {len(ev[k]['w'])} events, sigma = {ev[k]['w'].sum():.3f} fb")
    win = lambda e: (e["W"] >= 340) & (e["W"] <= 350) & (e["pstar"] < 50)
    tot = {k: hist("W", [340, 350], win, k)[0] for k in curves}
    log("")
    log("sigma(340 <= W <= 350, p* < 50) [fb], all gg (singlet + octet):")
    for k, v in tot.items():
        log(f"   {k:48s} {v:8.3f}   ({(v / tot[LO] - 1) * 100:+.1f}% vs LO)")
    for k8, k1 in ((A, FS), (B, OS)):
        log(f"   octet change {k8.split(':')[0]}: {(tot[k8] - tot[k1]) / (tot[k1] - tot[LO]) * 100:+.1f}% of the singlet "
            f"toponium excess")

    # ---------------- Figure 1: the Green's functions ----------------
    Eg, pg = np.linspace(-6, 4, 201), np.linspace(0.05, 30, 121)
    EE, PP = np.meshgrid(Eg, pg, indexing="ij")
    G0 = Constants.G0(EE, PP)
    maps = [("singlet: Fuks swData (CALCGREEN)", np.abs(ratio(EE.ravel(), PP.ravel()).reshape(EE.shape) * G0) ** 2),
            ("singlet: ours, JKT with CT25 $\\alpha_s(Q)$", np.abs(ours(EE.ravel(), PP.ravel()).reshape(EE.shape)) ** 2),
            ("free: $|G_0|^2$", np.abs(G0) ** 2),
            ("octet: repulsive Coulomb, $\\alpha$ = 0.15", np.abs(CoulombExact.G(EE, PP, -0.15 / 8)) ** 2),
            (f"octet: repulsive Coulomb, $\\alpha$ = {a8_ours:.3f}", np.abs(CoulombExact.G(EE, PP, -a8_ours / 8)) ** 2),
            ("octet / free: $|G_8/G_0|^2$, $\\alpha$ = 0.15", np.abs(CoulombExact.G(EE, PP, -0.15 / 8) / G0) ** 2)]
    norm = LogNorm(min(m.min() for _, m in maps[:5]), max(m.max() for _, m in maps[:5]))
    fig, ax = plt.subplots(2, 3, figsize=(20, 11), layout="constrained")
    for a, (title, M) in zip(ax.flat, maps):
        if title.startswith("octet / free"):
            im = a.pcolormesh(Eg, pg, M.T, shading="nearest", cmap="viridis", vmin=0.4, vmax=1.0)
            fig.colorbar(im, ax=a, label=r"$|G_8/G_0|^2$")
        else:
            im = a.pcolormesh(Eg, pg, M.T, shading="nearest", cmap="gist_ncar", norm=norm)
        a.set(title=title, xlabel="E [GeV]", ylabel="p [GeV]")
    fig.colorbar(im if False else ax.flat[0].collections[0], ax=ax[:, :2], label=r"$|\tilde G(E,p)|^2$ [GeV$^{-2}$]")
    fig.suptitle("Fig. 1: Green's functions, colour singlet (attractive) and colour octet (repulsive)", fontsize=13)
    fig.savefig(f"{OUT8}/fig1.png", dpi=110)
    plt.close(fig)
    log("")
    for E0_, p0 in ((-2.0, 5.0), (0.0, 10.0), (2.0, 25.0)):
        r8 = abs(CoulombExact.G(np.array(E0_), np.array(p0), -0.15 / 8) / Constants.G0(E0_, p0)) ** 2
        log(f"Fig. 1: |G8/G0|^2 (alpha 0.15) at E = {E0_:+.0f}, p = {p0:4.1f}: {r8:.3f}")

    # ---------------- Figure 2: d2sigma/p*^2 per bin (E, p*), 340 <= W <= 350, p* < 50 ----------------
    Eb, pb = np.linspace(-6, 4, 26), np.linspace(2, 30, 29)

    def h2(curve):
        tot2 = 0
        for e, wt in curves[curve]:
            c = win(e)
            tot2 = tot2 + np.histogram2d(e["E"][c], e["pstar"][c], bins=[Eb, pb], weights=(wt / e["pstar"] ** 2)[c])[0]
        return tot2

    H2 = {k: h2(k) for k in (LO, FS, A, B)}
    fig, ax = plt.subplots(1, 5, figsize=(29, 5.6), layout="constrained")
    for a, k in zip(ax[:4], H2):
        m = a.pcolormesh(Eb, pb, np.where(H2[k] > 0, H2[k], np.nan).T, cmap="gist_ncar", norm=LogNorm(5e-5, 5e-2))
        a.set(title=f"Fig. 2: {k}", xlabel="E [GeV]", ylabel="p* [GeV]")
    fig.colorbar(m, ax=ax[:4], label="d²σ/p*² per bin [fb GeV⁻²]")
    with np.errstate(divide="ignore", invalid="ignore"):
        rr = np.where(H2[FS] > 0, H2[A] / H2[FS], np.nan)
    m = ax[4].pcolormesh(Eb, pb, rr.T, cmap="viridis", vmin=0.6, vmax=1.05)
    ax[4].set(title="A / Fuks: effect of the octet", xlabel="E [GeV]", ylabel="p* [GeV]")
    fig.colorbar(m, ax=ax[4], label="ratio")
    fig.suptitle("all g g -> t tbar (singlet + octet), 340 <= W <= 350 GeV, p* < 50 GeV", fontsize=12)
    fig.savefig(f"{OUT8}/fig2.png", dpi=110)
    plt.close(fig)

    # ---------------- Figure 3: dsigma/dW, 330 <= W <= 350, 0.5 GeV bins, p* < 50 ----------------
    Wb = np.arange(330, 350 + 1e-9, 0.5)
    Wc = 0.5 * (Wb[1:] + Wb[:-1])
    c50 = lambda e: e["pstar"] < 50
    H3 = {k: hist("W", Wb, c50, k) for k in curves}
    fig, ax = plt.subplots(1, 2, figsize=(16, 6), layout="constrained")
    for k, h in H3.items():
        ax[0].step(Wc, h, where="mid", color=color[k], label=k, **style[k])
        ax[1].step(Wc, h / H3[FS], where="mid", color=color[k], label=k, **style[k])
    ax[0].set(title="Fig. 3: dσ/dW, all g g → t t̄ (singlet + octet), p* < 50 GeV", xlabel="W [GeV]",
              ylabel="dσ [fb/bin]", yscale="log")
    ax[0].legend(fontsize=8)
    ax[1].axhline(1, color="gray", lw=0.5)
    ax[1].set(title="ratio to Fuks (singlet x swData, octet LO)", xlabel="W [GeV]", ylim=(0, 1.6))
    ax[1].legend(fontsize=8)
    fig.savefig(f"{OUT8}/fig3.png", dpi=120)
    plt.close(fig)
    log("")
    log("Fig. 3, dsigma per 0.5 GeV bin [fb], p* < 50:")
    for W0 in (335.25, 340.25, 344.75, 349.75):
        i0 = np.argmin(abs(Wc - W0))
        log(f"   W = {Wc[i0]:6.2f}: " + "  ".join(f"{k.split(':')[0].split(' (')[0]}: {H3[k][i0]:.3g}" for k in H3))
    for k in (FS, A, OS, B):
        log(f"   peak {k:48s} W = {Wc[H3[k].argmax()]:.2f} GeV")

    # ---------------- Figure 4: p* distribution, p* < 50, 1 GeV bins, normalised ----------------
    pb4 = np.arange(0, 50 + 1e-9, 1.0)
    fig, a = plt.subplots(figsize=(8.5, 6), layout="constrained")
    log("")
    log("Fig. 4, dsigma/sigma per 1 GeV bin, p* < 50, 330 <= W <= 350:")
    for k in curves:
        h = hist("pstar", pb4, c50, k)
        h = h / h.sum()
        a.step(0.5 * (pb4[1:] + pb4[:-1]), h, where="mid", color=color[k], label=k, **style[k])
        log(f"   {k:48s} peak {h.max():.4f} at p* = {0.5 + h.argmax():.1f}; mean p* {np.sum(h * (0.5 + np.arange(50))):.2f}")
    a.set(title="Fig. 4: top recoil momentum, all g g → t t̄ (p* < 50 GeV)", xlabel="p* [GeV]", ylabel="dσ/σ [bin⁻¹]",
          yscale="log", ylim=(5e-4, 1e-1))
    a.legend(fontsize=8)
    fig.savefig(f"{OUT8}/fig4.png", dpi=120)
    plt.close(fig)

    # ---------------- Figure 5: m_tL, m_tH, -4 <= E <= 0, p* <= 50 ----------------
    cut5 = lambda e: (e["E"] >= -4) & (e["E"] <= 0) & (e["pstar"] <= 50)
    mb = np.arange(155, 180 + 1e-9, 0.25)
    mc = 0.5 * (mb[1:] + mb[:-1])
    fig, ax = plt.subplots(1, 2, figsize=(16, 6), layout="constrained")
    log("")
    log("Fig. 5 (-4 <= E <= 0, p* <= 50), all gg: 1/sigma dsigma/dm peak [1/GeV] and mean mass [GeV]:")
    for a, v in zip(ax, ("mtL", "mtH")):
        for k in curves:
            h = hist(v, mb, cut5, k)
            d = h / h.sum() / 0.25
            a.step(mc, d, where="mid", color=color[k], label=k, **style[k])
            log(f"   {v} {k:48s} peak {d.max():.3f} @ {mc[d.argmax()]:.2f}   mean {np.sum(h * mc) / h.sum():.2f}")
        a.set(title=f"Fig. 5: {v}, all g g → t t̄ (−4 ≤ E ≤ 0 GeV, p* ≤ 50 GeV)", xlabel=f"{v} [GeV]",
              ylabel="1/σ dσ/dm [GeV⁻¹]", xlim=(155, 180))
        a.legend(fontsize=8)
    fig.savefig(f"{OUT8}/fig5.png", dpi=120)
    plt.close(fig)
    open(f"{OUT8}/summary.txt", "w").write("\n".join(lines) + "\n")
    print(f"\nSaved {OUT8}/fig1..5.png and summary.txt")


# ---------------- the octet study with the Green's functions in MadGraph's MATRIX1 ----------------
OCT_DIR = "../output/fuks_octet_13TeV"
OCT_OUT = "results/4_fuks_figs2to5/fuks_octet_mg"
OCT_SAMPLES = {   # name: (seed, singlet table, octet table or None); all: Fuks's process, full colour, 326 <= W <= 350
    "LO": (8401, None, None),
    "A1": (8402, "A_singlet", None),          # Fuks's prescription: singlet x |G1/G0|^2 (swData), octet LO
    "A18": (8403, "A_singlet", "A_octet"),    # + octet x |G8/G0|^2, alpha = 0.15
    "B1": (8404, "B_singlet", None),          # ours: singlet x JKT with CT25 alpha_s(Q), octet LO
    "B18": (8405, "B_singlet", "B_octet"),    # + octet x |G8/G0|^2, alpha = CT25(25 GeV)
}


def octet_tables():
    """G/G0 tables for MATRIX1 (grid of 6_tfactors): A = Fuks's swData through his CALCGREEN, and the repulsive
    octet Coulomb at alpha = 0.15; B = Stage 6's tables (our JKT with CT25 alpha_s, octet at CT25(25 GeV))."""

    tf = importlib.import_module("6_tfactors")
    os.makedirs(OCT_DIR, exist_ok=True)
    Eg, Pg = tf.E0 + tf.DE * np.arange(tf.NE), tf.P0 + tf.DP * np.arange(tf.NP)
    EE, PP = np.meshgrid(Eg, Pg, indexing="ij")
    paths = {k: os.path.join(OCT_DIR, f"table_{k}.dat") for k in ("A_singlet", "A_octet", "B_singlet", "B_octet")}
    if not os.path.exists(paths["A_singlet"]):
        r = FuksTables.calcgreen_ratio()(EE.ravel(), np.maximum(PP, 0.01).ravel()).reshape(EE.shape)
        tf.write_table(paths["A_singlet"], r)
    if not os.path.exists(paths["A_octet"]):
        tf.write_table(paths["A_octet"], tf.octet_ratio(EE, PP, 0.15))
    shutil.copy(tf.TABLE, paths["B_singlet"])
    shutil.copy(tf.OCTET_TABLE, paths["B_octet"])
    return paths


def generate_octet_mg():
    """Five MadGraph samples of Fuks's process (full colour, CT25NNLO, default cuts, 500k events, one seed each),
    generated only in 326 <= W <= 350 GeV (tau window in genps.f). MATRIX1 (gg) is modified inside W <= 350, p* <= 50,
    E >= -20 GeV with the Stage 6 Fortran: |M|^2 + (|G1/G0|^2 - 1)|M_1|^2 [+ (|G8/G0|^2 - 1)|M_8|^2]."""

    tf = importlib.import_module("6_tfactors")
    tables = octet_tables()
    os.makedirs(os.path.join(OCT_DIR, "logs"), exist_ok=True)
    for name, (seed, t1, t8) in OCT_SAMPLES.items():
        d = os.path.abspath(os.path.join(OCT_DIR, f"TT_{name}"))
        cache = os.path.join(OCT_OUT, f"_cache_{name}.npz")
        os.makedirs(OCT_OUT, exist_ok=True)
        if os.path.exists(cache):
            print(f"{name}: done")
            continue
        if not os.path.isdir(d):
            script = os.path.abspath(os.path.join(OCT_DIR, f"_{name}.mg5"))   # mg5_aMC runs from the MG5 root
            open(script, "w").write(f"set automatic_html_opening False\nimport model sm\n{FUKS_PROCESS}\n"
                                    f"output {d} --hel_recycling=False\n")
            subprocess.run([sys.executable, "-W", "ignore", os.path.join(tf.MG_ROOT, "bin", "mg5_aMC"), script],
                           cwd=tf.MG_ROOT, check=True, env=tf.mg_env())
            os.remove(script)
        sub = os.path.join(d, "SubProcesses")
        genps = os.path.join(sub, "genps.f")
        src = open(genps).read()
        if "toponium window" not in src:
            old = "      TAUMIN = 0d0 !SMIN/S !keep scale fix\n      TAUMAX = 1D0\n"
            assert src.count(old) == 1
            open(genps, "w").write(src.replace(old, "      TAUMIN = 326.0D0**2/S ! toponium window 326 <= W <= 350 GeV\n"
                                                    "      TAUMAX = 350.0D0**2/S\n"))
        if t1:
            tf.CFG = {"octet": t8 is not None, "singlet_only": False}
            tf.GREEN_FILE = os.path.join(sub, "green_table.dat")
            tf.GREEN8_FILE = os.path.join(sub, "green8_table.dat")
            shutil.copy(tables[t1], tf.GREEN_FILE)
            if t8:
                shutil.copy(tables[t8], tf.GREEN8_FILE)
            for f in glob.glob(os.path.join(sub, "P*", "matrix1.f")):
                print(f"{name}: {f.split('SubProcesses/')[1]}: {tf.patch_matrix1(f)}")
        tf.set_card(os.path.join(d, "Cards", "run_card.dat"), {      # MadGraph default cuts kept
            "pdlabel": "lhapdf", "lhaid": tf.CT25_ID, "use_syst": "False", "event_norm": "average",
            "ebeam1": tf.EBEAM, "ebeam2": tf.EBEAM, "nevents": 500000, "iseed": seed})
        pc = os.path.join(d, "Cards", "param_card.dat")
        src, n = re.subn(r"^DECAY\s+6\s+\S+\s+# WT", "DECAY   6 1.490000e+00 # WT", open(pc).read(), flags=re.M)
        assert n == 1
        open(pc, "w").write(src)
        tf.set_config(os.path.join(d, "Cards", "me5_configuration.txt"), {
            "run_mode": 2, "nb_core": 10, "lhapdf_py3": os.path.join(tf.LHAPDF_DIR, "lhapdf-config"),
            "automatic_html_opening": "False"})
        run = f"{name}_s{seed}"
        print(f"{run}: generating 500000 events (seed {seed}) ...", flush=True)
        with open(os.path.join(OCT_DIR, "logs", f"{run}.log"), "w") as fl:
            r = subprocess.run([sys.executable, "-W", "ignore", "bin/generate_events", run, "-f"], cwd=d,
                               stdout=fl, stderr=subprocess.STDOUT, env=tf.mg_env())
        lhe = os.path.join(d, "Events", run, "unweighted_events.lhe.gz")
        if r.returncode != 0 or not os.path.exists(lhe):
            sys.exit(f"{run}: MadGraph failed, see {OCT_DIR}/logs/{run}.log")
        ev = Events.kinematics(lhe, cache)
        for f in glob.glob(os.path.join(d, "Events", run, "*.lhe*")):
            os.remove(f)
        print(f"{run}: {len(ev['w'])} events, sigma = {ev['w'].mean():.5g} pb", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--generate", action="store_true", help="MadGraph: Fuks's exact process (needs uv run --with six)")
    ap.add_argument("--sample", default="orig", choices=list(VARIANTS))
    ap.add_argument("--octet", action="store_true",
                    help="Figs. 1-5 with the colour octet in the toponium signal (fuksproc_ct25 samples, reweighted)")
    ap.add_argument("--octet-mg-generate", action="store_true",
                    help="MadGraph: the five octet-study samples with G in MATRIX1 (needs uv run --with six)")
    ap.add_argument("--octet-mg", action="store_true", help="Figs. 1-5 from the five MATRIX1 samples")
    a = ap.parse_args()
    if a.generate:
        generate()
    elif a.octet_mg_generate:
        generate_octet_mg()
    elif a.octet_mg:
        octet_figures(mode="mg")
    elif a.octet:
        octet_figures()
    else:
        main(a.sample)
