#!/usr/bin/env python3
"""STAGE 6 -- toponium T-factors for binned t-tbar distributions at 13 TeV.

Each RUN is two LO MadGraph samples of the same process at 13 TeV, 10 x 1M events each, one seed per run:
  WITHOUT  the standard matrix element of that run
  WITH     the gg matrix element changed with our Green's function, for W <= 350 GeV, p* <= 50 GeV:
               |M|^2 -> |M|^2 + (|G/G0|^2 - 1) |M_1|^2  =  |M_8|^2 + |G/G0|^2 |M_1|^2   (M_1: colour singlet)
T-factor per bin = sigma_with / sigma_without.
  run1_pp_fullcolour_13TeV  p p > t t~ > b e+ ve b~ mu- vm~: gg singlet + octet, and qqbar (octet)    [archived]
  run2_gg_singlet_13TeV     g g > t t~ > b e+ ve b~ mu- vm~, colour singlet only (CF = 2/3): no octet, no qqbar
  run3_pp_singlet_octet_13TeV  as run 1, but WITH also gives the octet (gg octet and qqbar) |G8/G0|^2, repulsive
                            Coulomb at alpha_s = CT25(25 GeV); WITHOUT = run 1's sample (octet threshold effect)
  run4_pp_singlet_octet_ll_13TeV  as run 3, Fuks's decay-chain syntax with all l+/l- (e, mu) combinations, own WITHOUT
  run5_window_kfactors_13TeV  the student's design (his kfactors_13TeV_W326-350_*.csv), with our G: K = 1 + Delta/SM_full,
                            SM_full = run 4's WITHOUT sample, Delta = (topo_win_gg - SM_win_gg) + (topo_win_qq - SM_win_qq)
                            from runs generated only
                            in W <= 350, p* <= 100 GeV (1M events each, his seeds); |G/G0|^2 for p* < 80 GeV; four
                            octet potentials (freeze10/25/50: alpha frozen below Q0; coulomb_as25)
  run6_final_kfactors_13TeV  the final K-factors: run 5's design and physics with its own SM_full (10 x 1M, full phase
                            space, seeds 7101..7191) and new window seeds (7400..7500); check topo_full = 1M full-phase-space
                            WITH sample (freeze25, seed 7601) to test the window design (window_check)
  run7_mt<m>_13TeV         m_t scan, m = 170, 171.25, 172.5, 173.75, 175 GeV (run7_mt170p00_13TeV, ...): run 6's design
                            (own SM_full 10 x 1M; SM_win, topo_win gg and qqbar, 1M each) at that m_t, Gamma_t scaled at LO
                            from 1.49 GeV at 173 (top_width); window W <= 2 m_t + 4, selection 2 m_t - 20 <= W <= 2 m_t + 4;
                            singlet = full JKT (its two-loop alpha_s, alpha_s(mZ) = 0.12), octet = coulomb_as25 only;
                            tables per mass (table_mt, made by setup); no check sample. MSU HPCC: submit_mtscan.sh
alpha_s everywhere (runs 1-6) = CT25NNLO's alpha_s(Q), as CTEQ use it: MadGraph takes alpha_s from the PDF (LHAPDF), and the
JKT potential behind G uses the same LHAPDF alpha_s(Q) in (13) instead of JKT's two-loop (14).

Run from toponium_study/ with MadGraph's extra dependency (RUN = one of the names above):
    uv run --with six python 6_tfactors.py table           # CT25 alpha_s(Q) from LHAPDF, G/G0 table (all runs)
    uv run --with six python 6_tfactors.py setup RUN       # LHAPDF dir, MadGraph outputs, MATRIX1 patches, cards
    uv run --with six python 6_tfactors.py smoke RUN       # 20k WITHOUT + 20k WITH, sanity check
    uv run --with six python 6_tfactors.py run RUN         # 20 runs of 1M events (resumable)       (~4-5 h)
    uv run --with six python 6_tfactors.py analyse RUN     # T-factors per bin, CSV, plots
  run 5 instead:  table run5 (p <= 80 tables, four octets); setup / smoke / run as above;
    uv run --with six python 6_tfactors.py kfactors run5_window_kfactors_13TeV [HIS_CSV ...]   # the two K-factor CSVs (+ comparison, + plots)
    uv run --with six python 6_tfactors.py plotk run5_window_kfactors_13TeV   # plots only, from the CSVs -> RUN/plots/pstar{100,50}/
  run 6 (the run 5 p80 tables, no table step): setup / smoke / run / kfactors / plotk run6_final_kfactors_13TeV  (~4.5 h)
    uv run --with six python 6_tfactors.py export run6_final_kfactors_13TeV    # easy-to-read K per binning -> RUN/kfactors_final/
    uv run --with six python 6_tfactors.py arbiter run6_final_kfactors_13TeV   # SM_cut (window by dummy_cuts only) vs SM_full, SM_win
    uv run --with six python 6_tfactors.py octetcheck run6_final_kfactors_13TeV   # singlet only (octet table = 1) vs singlet + octet
  run 7 (m_t scan): setup / smoke / one RUN NAME / kfactors / export per run7_mt*_13TeV; then
    uv run --with six python 6_tfactors.py mtscan                  # K vs m_t in the bins with an effect -> run7_mtscan_summary.txt
Output -> results/6_tfactors/{alphas_ct25nnlo.dat, green_jkt_ct25as.dat}           (shared)
          results/6_tfactors/RUN/{tfactors_13TeV.csv, summary.txt, plots/}
MadGraph and events -> ../output/RUN/{<without dir>, <with dir>, events/*.npz} (outside the repo)

=====================================================================================
EQUATION MAP  (paper  ->  where)
=====================================================================================
  Fuks (14)       process, l = e (top side), l' = mu (antitop side)            MG_PROCESSES
  Fuks (16)       |M|^2 -> |M|^2 + (|G/G0|^2 - 1)|M_1|^2 in MATRIX1 (gg only),
                  W <= 350 GeV, p* <= 50 GeV (and E >= -20 GeV, the table)     F_GREEN_FACTOR
  Fuks (19)-(20)  colour-singlet projection: CF1(I,J) = 2/3 -> |M_1|^2         F_DECLARATIONS, F_CF1
  Fuks (21)       rotations and boost to the t-tbar rest frame -> E, p*        F_KINEMATICS, fuks_pstar
  Fuks Sec. 3     GREENDATA module, CALCGREEN interpolation (our table, one
                  uniform grid, no off-by-one index)                           F_MODULE
  JKT (13)-(23)   the table: G/G0 from LSSolver(JKTPotential(alpha_fn=CT25))  table
                  run 7: LSSolver(JKTPotential()) at the run's m_t, Gamma_t     table_mt
  Fuks (12), C_F -> -1/6   octet tables of run 5: JKTPotential.octet_frozen
                  (freeze Q0) via LSSolver; CoulombExact (coulomb_as25)        table_window
"""

import glob
import os
import re
import shutil
import subprocess
import sys
import textwrap
import time

import numpy as np
from scipy.interpolate import CubicSpline

from toponium import Constants, JKTPotential, LSSolver, CoulombExact, Events, M_T


# =====================================================================================
# paths and run plan
# =====================================================================================
HERE = os.path.dirname(os.path.abspath(__file__))
MG_ROOT = os.path.dirname(HERE)
OUTPUT = os.path.join(MG_ROOT, "output")
LHAPDF_DIR = os.path.join(MG_ROOT, "lhapdf_data")
RES = os.path.join(HERE, "results", "6_tfactors")
TABLE = os.path.join(RES, "green_jkt_ct25as.dat")
OCTET_TABLE = os.path.join(RES, "green8_coulomb_ct25as25.dat")   # repulsive octet Coulomb, alpha_s = CT25(25 GeV)
OCTET_Q = 25.0                                                   # GeV: alpha_s scale of the octet (Bohr scale)
ALPHAS_FILE = os.path.join(RES, "alphas_ct25nnlo.dat")
# run 5: tables up to p = 80 GeV; octet variants: Q0 below which alpha_s is frozen (None: tree Coulomb at CT25(25 GeV))
TABLE80 = os.path.join(RES, "green_jkt_ct25as_p80.dat")
OCTETS = {"freeze25": 25.0, "freeze10": 10.0, "freeze50": 50.0, "coulomb_as25": None}
OCTET80 = {v: os.path.join(RES, f"green8_{v}_p80.dat") for v in OCTETS}

CT25_ID = 99999000          # local LHAPDF id for CT25NNLO (its own SetIndex overflows a Fortran integer)
CT25_SOURCES = (os.path.join(HERE, "CT25NNLO"), os.path.expanduser("~/Downloads/CT25NNLO"))
ENERGY, EBEAM = 13, 6500.0

RUNS = {   # Fuks (14), l = e (top side), l' = mu (antitop side); WITHOUT and WITH always the same process
    # runs 1-3: ../output dirs and events deleted 2026-11-04 (results kept in results/6_tfactors/run{1,2,3}_*)
    "run1_pp_fullcolour_13TeV": dict(
        desc="pp -> tt, gg singlet + octet and qqbar", process="p p > t t~ > b e+ ve b~ mu- vm~",
        without="TT6_SM_pp", with_="TT6_Topo_pp", singlet_only=False,
        names=("sm_13TeV", "topo_13TeV_ct25as"), seeds=(1301, 2301), smoke_seeds=(9301, 9302), archived=True,
        octet=False, without_from=None),
    "run2_gg_singlet_13TeV": dict(
        desc="gg -> tt, colour singlet only", process="g g > t t~ > b e+ ve b~ mu- vm~",
        without="TT6_gg1_without", with_="TT6_gg1_with", singlet_only=True,
        names=("without_gg1_13TeV", "with_gg1_13TeV"), seeds=(3301, 4301), smoke_seeds=(9303, 9304), archived=False,
        octet=False, without_from=None),
    "run3_pp_singlet_octet_13TeV": dict(   # WITHOUT = run 1's sample (same process, cards, PDF)
        desc="pp -> tt, gg singlet x |G1/G0|^2, octet (gg octet, qqbar) x |G8/G0|^2",
        process="p p > t t~ > b e+ ve b~ mu- vm~", without="TT6_SM_pp", with_="TT6_Topo18_pp", singlet_only=False,
        names=("sm_13TeV", "with18_13TeV"), seeds=(None, 6301), smoke_seeds=(9301, 9305), archived=False,
        octet=True, without_from="run1_pp_fullcolour_13TeV"),
    "run4_pp_singlet_octet_ll_13TeV": dict(   # as run 3, Fuks's decay-chain syntax, all l+/l- (e, mu) combinations
        desc="pp -> tt -> w+ b w- b~ (w -> l vl, l = e, mu), gg singlet x |G1/G0|^2, octet (gg octet, qqbar) x |G8/G0|^2",
        process="p p > t t~ > w+ b w- b~, w+ > l+ vl, w- > l- vl~", without="TT6_SM_ll", with_="TT6_Topo18_ll",
        singlet_only=False, names=("sm_ll_13TeV", "with18_ll_13TeV"), seeds=(5301, 5401), smoke_seeds=(9306, 9317),
        archived=False, octet=True, without_from=None, seed_step=10, compare_with="run3_pp_singlet_octet_13TeV"),
    "run5_window_kfactors_13TeV": dict(   # the student's design; SM_full = run 4's WITHOUT sample
        desc="K = 1 + Delta/SM_full: SM_win and topo_win runs for gg and qqbar in W <= 350, p* <= 100 GeV, |G/G0|^2 for p* < 80 GeV",
        process="p p > t t~ > w+ b w- b~, w+ > l+ vl, w- > l- vl~",
        process_gg="generate g g > t t~ > w+ b w- b~, w+ > l+ vl, w- > l- vl~",
        process_qq="define q = u c d s\ndefine q~ = u~ c~ d~ s~\ngenerate q q~ > t t~ > w+ b w- b~, w+ > l+ vl, w- > l- vl~",   # q q~ only: no mirror (half of qqbar)
        without="TT6_SM_ll", with_="TT7_T_gg", singlet_only=False, names=("sm_ll_13TeV", None), archived=False,
        octet=True, without_from="run4_pp_singlet_octet_ll_13TeV", window=True, p_rw=80.0,
        # (channel, octet or None for SM_win, seed): his seeds, so the unmodified SM_win runs reproduce his
        window_runs=[("gg", None, 2400), ("gg", "freeze25", 2410), ("qq", None, 2420), ("qq", "freeze25", 2430),
                     ("gg", "freeze10", 2450), ("gg", "freeze50", 2460), ("gg", "coulomb_as25", 2470),
                     ("qq", "freeze10", 2480), ("qq", "freeze50", 2490), ("qq", "coulomb_as25", 2500)]),
    "run6_final_kfactors_13TeV": dict(   # the final run: run 5's design with its own SM_full and new seeds, plus a check
        desc="K = 1 + Delta/SM_full (final): own SM_full (10 x 1M, full phase space); SM_win and topo_win runs for gg and qqbar in W <= 350, "
             "p* <= 100 GeV, |G/G0|^2 for p* < 80 GeV, four octets; check topo_full = 1M full-phase-space WITH (freeze25)",
        process="p p > t t~ > w+ b w- b~, w+ > l+ vl, w- > l- vl~",
        process_gg="generate g g > t t~ > w+ b w- b~, w+ > l+ vl, w- > l- vl~",
        # q q~ AND its beam-swapped mirror q~ q (MIRRORPROCS = .TRUE., as in p p); runs 5 and 6 (first pass,
        # now *_nomirror) had q q~ only, i.e. half of qqbar
        process_qq="define q = u c d s\ndefine q~ = u~ c~ d~ s~\ngenerate q q~ > t t~ > w+ b w- b~, w+ > l+ vl, w- > l- vl~"
                   "\nadd process q~ q > t t~ > w+ b w- b~, w+ > l+ vl, w- > l- vl~",
        without="TT6_SM_ll", with_="TT8_Topo_ll", singlet_only=False, names=("SM_full_ll_13TeV", "topo_full_ll_freeze25"),
        seeds=(7101, 7601), smoke_seeds=(9320, 9321), seed_step=10, archived=False,
        octet=True, without_from=None, window=True, p_rw=80.0, check="freeze25",
        window_runs=[("gg", None, 7400), ("gg", "freeze25", 7410), ("qq", None, 7520), ("qq", "freeze25", 7530),
                     ("gg", "freeze10", 7450), ("gg", "freeze50", 7460), ("gg", "coulomb_as25", 7470),
                     ("qq", "freeze10", 7570), ("qq", "freeze50", 7590), ("qq", "coulomb_as25", 7700)]),
        # qq seeds of the first (no-mirror) pass: 7420, 7430, 7480, 7490, 7500
}


def top_width(m, mW=80.4):
    """LO t -> W b width scaled to Gamma_t(173) = 1.49 GeV: Gamma ~ m^3 (1 - x)^2 (1 + 2x), x = mW^2/m^2."""

    f = lambda m_: m_ ** 3 * (1 - mW ** 2 / m_ ** 2) ** 2 * (1 + 2 * mW ** 2 / m_ ** 2)
    return round(1.49 * f(m) / f(173.0), 4)


MT_SCAN = (170.0, 171.25, 172.5, 173.75, 175.0)
for _i, _m in enumerate(MT_SCAN):   # run 7: run 6's design at five top masses, full JKT singlet, octet coulomb_as25
    _b = 8000 + 200 * _i
    RUNS[f"run7_mt{_m:.2f}".replace(".", "p") + "_13TeV"] = dict(
        RUNS["run6_final_kfactors_13TeV"],
        desc=f"m_t scan, m_t = {_m:g} GeV, Gamma_t = {top_width(_m)} GeV (LO scaling): K = 1 + Delta/SM_full, own SM_full "
             f"(10 x 1M); SM_win and topo_win for gg and qqbar in W <= 2 m_t + 4, p* <= 100 GeV; singlet: full JKT "
             f"(two-loop alpha_s, alpha_s(mZ) = 0.12); octet: coulomb_as25",
        mt=_m, wt=top_width(_m), singlet="jkt", with_="TT7_T_gg", names=("SM_full_ll_13TeV", None),
        seeds=(_b + 1, None), smoke_seeds=(_b + 150, None), check=None,
        # seeds: MadGraph advances them by 3 per survey/refine step; no two of a mass point differ by a multiple of 3
        # below 30 (10 steps, as run 6's SM_full), and the five points are 200 apart
        window_runs=[("gg", None, _b + 104), ("gg", "coulomb_as25", _b + 114),
                     ("qq", None, _b + 124), ("qq", "coulomb_as25", _b + 134)],
        kfiles=[(100.0, "coulomb_as25", ()), (50.0, "coulomb_as25", ())])

RUN, CFG = None, None
SM_DIR = TOPO_DIR = EV_DIR = EV_DIR_WITHOUT = GREEN_FILE = GREEN8_FILE = RES_RUN = None
MT, WT = 173.0, 1.49                  # m_t, Gamma_t of the selected run (run 7: per mass point)
T1_FILE, T8_FILES = None, None        # the p80 singlet table and octet tables the selected run's window runs read


def wlo():
    return 2 * MT - 20          # selection window 2 m_t - 20 <= W <= 2 m_t + 4 (326-350 at m_t = 173)


def whi():
    return 2 * MT + 4


def select(run):
    """Point all run-dependent paths at RUNS[run]: output/<run>/{without, with, events}, results/6_tfactors/<run>."""

    global RUN, CFG, SM_DIR, TOPO_DIR, EV_DIR, EV_DIR_WITHOUT, GREEN_FILE, GREEN8_FILE, RES_RUN, MT, WT, T1_FILE, T8_FILES
    if run not in RUNS:
        sys.exit(f"run must be one of {list(RUNS)}")
    RUN, CFG = run, RUNS[run]
    MT, WT = CFG.get("mt", 173.0), CFG.get("wt", 1.49)
    if "mt" in CFG:                       # run 7: tables for this m_t, Gamma_t (table_mt, written by setup)
        T1_FILE = os.path.join(RES, run, "green1_jkt_p80.dat")
        T8_FILES = {"coulomb_as25": os.path.join(RES, run, "green8_coulomb_as25_p80.dat")}
    else:
        T1_FILE, T8_FILES = TABLE80, OCTET80
    src = CFG["without_from"] or run                                         # where the WITHOUT sample lives
    SM_DIR = os.path.join(OUTPUT, src, CFG["without"])                      # WITHOUT
    TOPO_DIR = os.path.join(OUTPUT, run, CFG["with_"])                      # WITH
    EV_DIR = os.path.join(OUTPUT, run, "events")
    EV_DIR_WITHOUT = os.path.join(OUTPUT, src, "events")
    GREEN_FILE = os.path.join(TOPO_DIR, "SubProcesses", "green_table.dat")  # what MATRIX1 reads
    GREEN8_FILE = os.path.join(TOPO_DIR, "SubProcesses", "green8_table.dat")
    RES_RUN = os.path.join(RES, run)

# table grid: E in [-20, 4] GeV (step 0.025), p in [0, 50] GeV (step 0.125)
E0, DE, NE = -20.0, 0.025, 961
P0, DP, NP = 0.0, 0.125, 401


def window_dir(kind, ch):
    """Runs 5-6: the MadGraph output of SM_win (kind S, unmodified) or topo_win (kind T, with G), channel gg or qq."""

    return os.path.join(OUTPUT, RUN, f"TT7_{kind}_{ch}")


def window_name(ch, octet, seed):
    return f"topo_win_{ch}_{octet}_s{seed}" if octet else f"SM_win_{ch}_s{seed}"


def runs(smoke=False):
    """(name, directory, seed, nevents, octet): 10 WITHOUT + 10 WITH runs of 1M events at 13 TeV, one seed each
    (octet = None); run 5: its 10 window runs (octet = the table the T run reads)."""

    if CFG.get("window"):
        win = [(window_name(ch, o, s) + ("_smoke" if smoke else ""), window_dir("T" if o else "S", ch),
                s + 7000 * smoke, 20000 if smoke else 1000000, o) for ch, o, s in CFG["window_runs"]]
        if CFG["without_from"]:                            # run 5: SM_full is another run's WITHOUT sample
            return win
        # run 6: own SM_full (10 x 1M, as runs 1-4), the window runs, and the check topo_full (1M full phase space, WITH)
        # (run 7: no check sample)
        (n0, n1), (s0, s1), st = CFG["names"], CFG["seeds"], CFG["seed_step"]
        chk = bool(CFG.get("check"))
        if smoke:
            a, b = CFG["smoke_seeds"]
            return ([(f"{n0}_smoke", SM_DIR, a, 20000, None)] + win
                    + [(f"{n1}_smoke", TOPO_DIR, b, 20000, CFG["check"])] * chk)
        return ([(f"{n0}_s{s0 + st * k}", SM_DIR, s0 + st * k, 1000000, None) for k in range(10)] + win
                + [(f"{n1}_s{s1}", TOPO_DIR, s1, 1000000, CFG["check"])] * chk)
    (n0, n1), (s0, s1) = CFG["names"], CFG["seeds"]
    # MadGraph advances the seed by +3 at every survey/refine step of a run (madevent_interface.update_random):
    # seeds 10 apart give runs disjoint seed sequences (s + 3k never equals s + 10 + 3m). Runs 1-3 used step 1.
    st = CFG.get("seed_step", 1)
    own_without = CFG["without_from"] is None          # run 3 reuses run 1's WITHOUT sample
    if smoke:
        a, b = CFG["smoke_seeds"]
        return [(f"{n0}_smoke", SM_DIR, a, 20000, None)] * own_without + [(f"{n1}_smoke", TOPO_DIR, b, 20000, None)]
    without = [(f"{n0}_s{s0 + st * k}", SM_DIR, s0 + st * k, 1000000, None) for k in range(10)] if own_without else []
    return without + [(f"{n1}_s{s1 + st * k}", TOPO_DIR, s1 + st * k, 1000000, None) for k in range(10)]


# =====================================================================================
# 1. CT25NNLO alpha_s(Q) from LHAPDF, and the G/G0 table (JKT (13)-(23) via LSSolver)
# =====================================================================================
AS_DUMP_CC = r"""#include "LHAPDF/LHAPDF.h"
#include <cstdio>
#include <cstdlib>
int main(int argc, char** argv) {
  LHAPDF::setVerbosity(0);
  const LHAPDF::PDF* pdf = LHAPDF::mkPDF(argv[1], 0);
  for (int i = 2; i < argc; ++i) { double Q = atof(argv[i]); printf("%.12e %.14e\n", Q, pdf->alphasQ(Q)); }
  delete pdf;
  return 0;
}
"""


def lhapdf_alphas(Qs):
    """alpha_s(Q) of CT25NNLO member 0 exactly as LHAPDF returns it (compiled once into the results dir)."""

    exe = os.path.join(LHAPDF_DIR, "_as_dump")
    if not os.path.exists(exe):
        src = exe + ".cc"
        open(src, "w").write(AS_DUMP_CC)
        flags = subprocess.run(["lhapdf-config", "--cflags", "--libs"], capture_output=True, text=True).stdout.split()
        subprocess.run([shutil.which("clang++") or "g++", "-std=c++17", "-O2", src, *flags, "-o", exe], check=True)
        os.remove(src)
    out = subprocess.run([exe, "CT25NNLO", *[f"{q:.12e}" for q in Qs]], capture_output=True, text=True, check=True,
                         env=dict(os.environ, LHAPDF_DATA_PATH=LHAPDF_DIR)).stdout
    return np.array([float(l.split()[1]) for l in out.strip().splitlines()])


_AS = {}


def ct25_alphas(Q2):
    """CT25NNLO's alpha_s at Q^2 >= 25 GeV^2: cubic spline in ln Q through the LHAPDF dump (alphas_ct25nnlo.dat).
    The potential needs alpha only above |p| = 5 GeV ((13) and the continuity point of (15))."""

    if "spl" not in _AS:
        d = np.loadtxt(ALPHAS_FILE)
        _AS["spl"] = CubicSpline(np.log(d[:, 0]), d[:, 1])
        _AS["lo"], _AS["hi"] = np.log(d[0, 0]), np.log(d[-1, 0])
    lq = 0.5 * np.log(np.asarray(Q2, dtype=float))
    return _AS["spl"](np.clip(lq, _AS["lo"], _AS["hi"]))


def load_table(path=None):
    path = path or TABLE
    ne, np_ = (int(x) for x in open(path).readline().split()[:2])
    d = np.loadtxt(path, skiprows=1)
    return (d[:, 0] + 1j * d[:, 1]).reshape(ne, np_)


def write_table(path, R):
    with open(path, "w") as f:
        f.write(f"{R.shape[0]} {R.shape[1]} {E0} {DE} {P0} {DP}\n")
        for z in R.ravel():
            f.write(f"{z.real:.10e} {z.imag:.10e}\n")


def octet_ratio(E, p, alpha_s, m_t=M_T, Gamma_t=1.49):
    """G8/G0 for the colour octet: exact Coulomb solution (CoulombExact) for the repulsive potential
    V8 = +alpha_s/(6r), i.e. colour factor -1/6 = C_F x (-1/8). No bound states."""

    return CoulombExact.G(E, p, -alpha_s / 8, m_t, Gamma_t) / Constants.G0(E, p, m_t, Gamma_t)


def interp_table(tab, E, p):
    """Bilinear interpolation exactly as CALCGREEN in matrix1.f (F_MODULE): G/G0 at (E, p)."""

    a, b = (np.asarray(E) - E0) / DE, (np.asarray(p) - P0) / DP
    i = np.clip(np.trunc(a).astype(int), 0, tab.shape[0] - 2)
    j = np.clip(np.trunc(b).astype(int), 0, tab.shape[1] - 2)
    a, b = a - i, b - j
    return ((1 - a) * (1 - b) * tab[i, j] + a * (1 - b) * tab[i + 1, j]
            + (1 - a) * b * tab[i, j + 1] + a * b * tab[i + 1, j + 1])


def table():
    """CT25NNLO alpha_s(Q) from LHAPDF (2000 points, 5 GeV - 1e8 GeV), checked off-grid; then G/G0 on E in [-20, 4]
    (step 0.025), p in [0, 50] (step 0.125) for JKTPotential(alpha_fn=CT25), N = 600 nodes. The ratio is to the
    physical G0(E, p), JKT (6); G includes the energy shift of JKT (20)."""

    os.makedirs(RES, exist_ok=True)
    lines = []
    log = lambda s: (print(s), lines.append(s))
    Q = np.geomspace(5.0, 1e8, 2000)
    np.savetxt(ALPHAS_FILE, np.column_stack([Q, lhapdf_alphas(Q)]), fmt="%.12e",
               header="Q [GeV], alpha_s(Q) of CT25NNLO member 0 from LHAPDF (alphasQ); above QMax = 1e5 GeV LHAPDF "
                      "freezes it")
    _AS.clear()
    Qc = np.geomspace(5.3, 9e4, 97)
    d = np.abs(ct25_alphas(Qc ** 2) / lhapdf_alphas(Qc) - 1)
    log(f"alpha_s(Q) spline vs LHAPDF at 97 off-grid Q in [5.3, 9e4] GeV: max rel. difference {d.max():.1e}")
    Lam = JKTPotential.lambda_msbar_2loop(0.118)
    for q in (5.0, 10.0, 25.0, 91.1876):
        log(f"  alpha_s({q:7.4g} GeV): CT25NNLO {float(ct25_alphas(q * q)):.5f}   JKT (14) two-loop at 0.118 "
            f"{float(JKTPotential.alpha_2loop(np.array(q * q), Lam)):.5f}")

    t0 = time.time()
    pot = JKTPotential(alpha_mz=0.118, alpha_fn=ct25_alphas)
    ref = JKTPotential(alpha_mz=0.118)
    log(f"potential with CT25 alpha_s: C = {pot.C:.5f} (two-loop: {ref.C:.5f}), shift = {pot.shift:.4f} GeV "
        f"(two-loop: {ref.shift:.4f})")
    S, Sref = LSSolver(pot, n=120), LSSolver(ref, n=120)
    Eg, Pg = E0 + DE * np.arange(NE), P0 + DP * np.arange(NP)
    pp = np.maximum(Pg, 1e-3)                           # G is smooth at p -> 0; the spline is in ln p
    R = np.array([S.G(E, pp) / Constants.G0(E, pp) for E in Eg])
    with open(TABLE, "w") as f:
        f.write(f"{NE} {NP} {E0} {DE} {P0} {DP}\n")
        for z in R.ravel():
            f.write(f"{z.real:.10e} {z.imag:.10e}\n")
    rng = np.random.default_rng(1)
    Es, ps = rng.uniform(-20, 4, 300), rng.uniform(0.2, 50, 300)
    direct = np.array([abs(S.G(E, np.array([p]))[0] / Constants.G0(E, p)) ** 2 for E, p in zip(Es, ps)])
    dd = np.abs(np.abs(interp_table(load_table(), Es, ps)) ** 2 / direct - 1)
    log(f"G/G0 table written ({time.time() - t0:.0f} s); |G/G0|^2 table vs solver at 300 random points: median "
        f"{np.median(dd):.1e}, max {dd.max():.1e}")
    for p in (5.0, 15.0, 30.0):
        r1 = abs(S.G(-2.0, np.array([p]))[0] / Constants.G0(-2.0, p)) ** 2
        r2 = abs(Sref.G(-2.0, np.array([p]))[0] / Constants.G0(-2.0, p)) ** 2
        log(f"  |G/G0|^2 at E = -2, p = {p:4.1f}: CT25 alpha_s {r1:.4f}, two-loop at 0.118 {r2:.4f} "
            f"({(r1 / r2 - 1) * 100:+.2f}%)")
    # colour octet: repulsive Coulomb, exact, alpha_s = CT25NNLO's alpha_s at the Bohr scale
    a8 = float(ct25_alphas(OCTET_Q ** 2))
    EE, PP = np.meshgrid(Eg, pp, indexing="ij")
    write_table(OCTET_TABLE, octet_ratio(EE, PP, a8))
    S8 = LSSolver(JKTPotential.coulomb(-a8 / 8), n=48)
    d8 = [abs(CoulombExact.G(np.array(E), np.array(p), -a8 / 8) / S8.G(E, np.array([p]))[0] - 1)
          for E, p in zip(Es[:60], ps[:60])]
    t8 = np.abs(interp_table(load_table(OCTET_TABLE), Es, ps)) ** 2
    x8 = np.abs(octet_ratio(Es, ps, a8)) ** 2
    log(f"octet G8/G0 (repulsive Coulomb, alpha_s = CT25(25 GeV) = {a8:.4f}): exact vs LS solver max {max(d8):.1e}; "
        f"table vs exact max {np.max(np.abs(t8 / x8 - 1)):.1e}")
    for E, p in ((-2.0, 5.0), (0.0, 10.0), (4.0, 30.0), (-2.0, 50.0)):
        log(f"  |G8/G0|^2 at E = {E:+.0f}, p = {p:4.1f}: {abs(octet_ratio(np.array(E), np.array(p), a8)) ** 2:.4f}")
    open(os.path.join(RES, "table_log.txt"), "w").write("\n".join(lines) + "\n")
    for run, cfg in RUNS.items():
        sub = os.path.join(OUTPUT, run, cfg["with_"], "SubProcesses")
        if not cfg["archived"] and os.path.isdir(sub):
            shutil.copy(TABLE, os.path.join(sub, "green_table.dat"))
            if cfg["octet"]:
                shutil.copy(OCTET_TABLE, os.path.join(sub, "green8_table.dat"))


# =====================================================================================
# 2. MadGraph set-up: LHAPDF, outputs, MATRIX1 patch (WITH only), cards
# =====================================================================================
def fortran_string(name, s, width=50):
    """PARAMETER for a long character constant, split over fixed-form continuation lines."""

    parts = [s[i:i + width] for i in range(0, len(s), width)]
    body = "//\n".join(f"     $ '{p}'" for p in parts)
    return f"      CHARACTER*(*) {name}\n      PARAMETER ({name}=\n{body})\n"


def F_MODULE(octet=False):
    """Fuks Sec. 3: GREENDATA module and CALCGREEN. Placed at the top of matrix1.f (Fuks put the
    module in driver.f) so that it is compiled before MATRIX1 uses it. octet=True adds the colour-octet
    table (GRE8/GIM8 from green8_table.dat, same grid) and CALCGREEN8."""

    o8 = octet
    return ("C BEGIN toponium addition: GREENDATA and CALCGREEN (Fuks Sec. 3), our G/G0 table\n"
            "      MODULE GREENDATA\n"
            "      IMPLICIT NONE\n"
            "      INTEGER NEG, NPG\n"
            "      REAL*8 E0G, DEG, P0G, DPG\n"
            "      REAL*8, ALLOCATABLE :: GRE(:,:), GIM(:,:)\n"
            + ("      REAL*8, ALLOCATABLE :: GRE8(:,:), GIM8(:,:)\n" if o8 else "") +
            "      LOGICAL :: LOADED = .FALSE.\n"
            + fortran_string("GFILE", GREEN_FILE)
            + (fortran_string("GFILE8", GREEN8_FILE) if o8 else "") +
            "      CONTAINS\n"
            "      SUBROUTINE LOADGREEN\n"
            "      INTEGER I, J, U\n"
            "      OPEN(NEWUNIT=U, FILE=GFILE, STATUS='OLD')\n"
            "      READ(U,*) NEG, NPG, E0G, DEG, P0G, DPG\n"
            "      ALLOCATE(GRE(NEG,NPG), GIM(NEG,NPG))\n"
            "      DO I = 1, NEG\n"
            "        DO J = 1, NPG\n"
            "          READ(U,*) GRE(I,J), GIM(I,J)\n"
            "        ENDDO\n"
            "      ENDDO\n"
            "      CLOSE(U)\n"
            + ("C     colour octet, same grid\n"
               "      OPEN(NEWUNIT=U, FILE=GFILE8, STATUS='OLD')\n"
               "      READ(U,*)\n"
               "      ALLOCATE(GRE8(NEG,NPG), GIM8(NEG,NPG))\n"
               "      DO I = 1, NEG\n"
               "        DO J = 1, NPG\n"
               "          READ(U,*) GRE8(I,J), GIM8(I,J)\n"
               "        ENDDO\n"
               "      ENDDO\n"
               "      CLOSE(U)\n" if o8 else "") +
            "      LOADED = .TRUE.\n"
            "      END SUBROUTINE LOADGREEN\n"
            "C     G/G0 at (EOMX, MOMX) by bilinear interpolation on the uniform grid\n"
            "      SUBROUTINE CALCGREEN(EOMX, MOMX, GREEN)\n"
            "      REAL*8 EOMX, MOMX, AA, BB\n"
            "      COMPLEX*16 GREEN\n"
            "      INTEGER II, JJ\n"
            "      IF(.NOT.LOADED) CALL LOADGREEN\n"
            "      AA = (EOMX - E0G)/DEG\n"
            "      BB = (MOMX - P0G)/DPG\n"
            "      II = MIN(MAX(INT(AA), 0), NEG-2)\n"
            "      JJ = MIN(MAX(INT(BB), 0), NPG-2)\n"
            "      AA = AA - II\n"
            "      BB = BB - JJ\n"
            "      GREEN = DCMPLX(\n"
            "     $ (1-AA)*(1-BB)*GRE(II+1,JJ+1) + AA*(1-BB)*GRE(II+2,JJ+1)\n"
            "     $ + (1-AA)*BB*GRE(II+1,JJ+2) + AA*BB*GRE(II+2,JJ+2),\n"
            "     $ (1-AA)*(1-BB)*GIM(II+1,JJ+1) + AA*(1-BB)*GIM(II+2,JJ+1)\n"
            "     $ + (1-AA)*BB*GIM(II+1,JJ+2) + AA*BB*GIM(II+2,JJ+2))\n"
            "      END SUBROUTINE CALCGREEN\n"
            + ("C     G8/G0 (colour octet, repulsive) at (EOMX, MOMX), same interpolation\n"
               "      SUBROUTINE CALCGREEN8(EOMX, MOMX, GREEN)\n"
               "      REAL*8 EOMX, MOMX, AA, BB\n"
               "      COMPLEX*16 GREEN\n"
               "      INTEGER II, JJ\n"
               "      IF(.NOT.LOADED) CALL LOADGREEN\n"
               "      AA = (EOMX - E0G)/DEG\n"
               "      BB = (MOMX - P0G)/DPG\n"
               "      II = MIN(MAX(INT(AA), 0), NEG-2)\n"
               "      JJ = MIN(MAX(INT(BB), 0), NPG-2)\n"
               "      AA = AA - II\n"
               "      BB = BB - JJ\n"
               "      GREEN = DCMPLX(\n"
               "     $ (1-AA)*(1-BB)*GRE8(II+1,JJ+1) + AA*(1-BB)*GRE8(II+2,JJ+1)\n"
               "     $ + (1-AA)*BB*GRE8(II+1,JJ+2) + AA*BB*GRE8(II+2,JJ+2),\n"
               "     $ (1-AA)*(1-BB)*GIM8(II+1,JJ+1) + AA*(1-BB)*GIM8(II+2,JJ+1)\n"
               "     $ + (1-AA)*BB*GIM8(II+1,JJ+2) + AA*BB*GIM8(II+2,JJ+2))\n"
               "      END SUBROUTINE CALCGREEN8\n" if o8 else "") +
            "      END MODULE GREENDATA\n"
            "C END toponium addition\n\n")


F_USE = "C BEGIN toponium addition (Fuks Sec. 3)\n      USE GREENDATA\nC END toponium addition\n"

F_DECLARATIONS = ("C BEGIN toponium addition: toponium kinematics (Fuks Sec. 3), singlet colour matrix\n"
                  "      REAL*8 PMOM(0:3,3)\n"
                  "      REAL*8 GXOM, GYOM, GZOM, GEOM\n"
                  "      REAL*8 GMOM, GTHE, GPHI, GBET, GGAM, GMAS\n"
                  "      REAL*8 EXXX, MXXX, MSING\n"
                  "      COMPLEX*16 GREEN, GREEN8\n"
                  "      REAL*8 CF1(NCOLOR,NCOLOR)\n"
                  "C END toponium addition\n")

F_CF1 = ("C BEGIN toponium addition: colour-singlet projection, Fuks (19)-(20); CF above is kept\n"
         "      DATA (CF1(I,  1),I=  1,  2) /6.666666666666666D-01,\n"
         "     $ 6.666666666666666D-01/\n"
         "      DATA (CF1(I,  2),I=  1,  2) /6.666666666666666D-01\n"
         "     $ ,6.666666666666666D-01/\n"
         "C END toponium addition\n")

F_KINEMATICS = """C BEGIN toponium addition: P(.,3..5) = b e+ ve (top), P(.,6..8) = b~ mu- vm~ (antitop)
C     PMOM(.,1) = toponium, PMOM(.,2) = top, PMOM(.,3) = antitop; Fuks Sec. 3
      DO M = 0, 3
        PMOM(M,2) = 0D0
        PMOM(M,3) = 0D0
        DO I = 3, 5
          PMOM(M,2) = PMOM(M,2) + P(M,I)
          PMOM(M,3) = PMOM(M,3) + P(M,I+3)
        ENDDO
        PMOM(M,1) = PMOM(M,2) + PMOM(M,3)
      ENDDO
C     Lorentz transformation factors, Fuks (21)
      GMOM = SQRT(PMOM(1,1)**2 + PMOM(2,1)**2 + PMOM(3,1)**2)
      GMAS = SQRT(PMOM(0,1)**2 - GMOM**2)
      GTHE = DACOS(PMOM(3,1)/(GMOM+1.D-8))
      GPHI = ATAN2(PMOM(2,1), PMOM(1,1))
      GBET = GMOM/PMOM(0,1)
      GGAM = PMOM(0,1)/GMAS
C     Rotating around Z-axis
      DO M = 1, 3
        GXOM = PMOM(1,M)
        GYOM = PMOM(2,M)
        PMOM(1,M) = COS(GPHI)*GXOM + SIN(GPHI)*GYOM
        PMOM(2,M) = COS(GPHI)*GYOM - SIN(GPHI)*GXOM
      ENDDO
C     Rotating around Y-axis
      DO M = 1, 3
        GXOM = PMOM(1,M)
        GZOM = PMOM(3,M)
        PMOM(1,M) = COS(GTHE)*GXOM - SIN(GTHE)*GZOM
        PMOM(3,M) = COS(GTHE)*GZOM + SIN(GTHE)*GXOM
      ENDDO
C     Boost along Z-axis
      DO M = 1, 3
        GEOM = PMOM(0,M)
        GZOM = PMOM(3,M)
        PMOM(0,M) = GGAM*GEOM - GGAM*GBET*GZOM
        PMOM(3,M) = GGAM*GZOM - GGAM*GBET*GEOM
      ENDDO
C END toponium addition
"""

def table_window():
    """Run 5 tables on E in [-20, 4] (step 0.025), p in [0, 80] (step 0.125): the singlet G/G0 as in table() (JKT
    (13)-(23), CT25 alpha_s), and G8/G0 for the four octet potentials: freeze Q0 = JKTPotential.octet_frozen
    (Fuks (12) with C_F -> -1/6, alpha_s(max(q, Q0)), one-loop bracket) via LSSolver; coulomb_as25 = CoulombExact
    at alpha_s = CT25(25 GeV), tree level (octet_ratio, as runs 3-4)."""

    os.makedirs(RES, exist_ok=True)
    lines = []
    log = lambda s: (print(s), lines.append(s))
    Eg, Pg = E0 + DE * np.arange(NE), P0 + DP * np.arange(int(round(80 / DP)) + 1)
    pp = np.maximum(Pg, 1e-3)
    rng = np.random.default_rng(5)
    Es, ps = rng.uniform(-20, 4, 200), rng.uniform(0.2, 80, 200)
    a8 = float(ct25_alphas(OCTET_Q ** 2))

    # code check of the general kernel path: octet_frozen with a constant alpha is exactly Coulomb (bracketed alpha)
    a_eff = a8 * (1 + (31 / 3 - 10 * 5 / 9) * a8 / (4 * np.pi))
    Sc = LSSolver(JKTPotential.octet_frozen(lambda Q2: np.full_like(np.asarray(Q2, dtype=float), a8), 25.0), n=120)
    d = [abs(Sc.G(E, np.array([p]))[0] / CoulombExact.G(np.array(E), np.array(p), -a_eff / 8) - 1)
         for E, p in zip(Es[:40], ps[:40])]
    log(f"octet_frozen with constant alpha = {a8:.5f} vs CoulombExact (alpha x bracket): max rel. diff of G {max(d):.1e}")

    jobs = [("singlet", TABLE80, LSSolver(JKTPotential(alpha_mz=0.118, alpha_fn=ct25_alphas), n=120))]
    jobs += [(v, OCTET80[v], LSSolver(JKTPotential.octet_frozen(ct25_alphas, Q0), n=120))
             for v, Q0 in OCTETS.items() if Q0 is not None]
    for v, path, S in jobs:
        t0 = time.time()
        write_table(path, np.array([S.G(E, pp) / Constants.G0(E, pp) for E in Eg]))
        tab = load_table(path)
        direct = np.array([abs(S.G(E, np.array([p]))[0] / Constants.G0(E, p)) ** 2 for E, p in zip(Es, ps)])
        dt = np.abs(np.abs(interp_table(tab, Es, ps)) ** 2 / direct - 1)
        S2 = LSSolver(S.pot, n=160)
        dn = [abs(abs(S2.G(E, np.array([p]))[0] / Constants.G0(E, p)) ** 2 / x - 1)
              for E, p, x in zip(Es[:40], ps[:40], direct[:40])]
        log(f"{v:12s} {os.path.basename(path)} ({time.time() - t0:.0f} s): |G/G0|^2 table vs solver median "
            f"{np.median(dt):.1e}, max {dt.max():.1e}; n = 120 vs 160 max {max(dn):.1e}")
    old = load_table(TABLE)
    log(f"singlet p80 table vs the runs 1-4 table on p <= 50: max |difference| "
        f"{np.abs(load_table(TABLE80)[:, :old.shape[1]] - old).max():.1e}")
    EE, PP = np.meshgrid(Eg, pp, indexing="ij")
    write_table(OCTET80["coulomb_as25"], octet_ratio(EE, PP, a8))
    log(f"coulomb_as25 {os.path.basename(OCTET80['coulomb_as25'])}: exact repulsive Coulomb, alpha_s = CT25(25 GeV) "
        f"= {a8:.6f}, tree level")
    for E, p in ((-2.0, 5.0), (0.0, 10.0), (4.0, 30.0), (-2.0, 50.0), (2.0, 75.0)):
        r = {v: abs(interp_table(load_table(OCTET80[v]), np.array(E), np.array(p))) ** 2 for v in OCTETS}
        log(f"  |G8/G0|^2 at E = {E:+.0f}, p = {p:4.1f}: " + ", ".join(f"{v} {x:.4f}" for v, x in r.items()))
    out = os.path.join(RES, "run5_window_kfactors_13TeV")
    os.makedirs(out, exist_ok=True)
    open(os.path.join(out, "table_log.txt"), "w").write("\n".join(lines) + "\n")


def table_mt():
    """Run 7 tables for the selected run's m_t, Gamma_t on E in [-20, 4] (step 0.025), p in [0, 80] (step 0.125):
      singlet  G1/G0 for the full JKT potential (13)-(23): JKT's two-loop alpha_s (14) at alpha_s(mZ) = 0.12, NLO
               bracket, IR form (15) below 5 GeV (JKTPotential() defaults), LSSolver n = 120 with this m_t, Gamma_t
      octet    coulomb_as25: exact repulsive Coulomb at alpha_s = CT25(25 GeV), tree level (octet_ratio)
    both divided by G0(E, p; m_t, Gamma_t), JKT (6)."""

    os.makedirs(os.path.dirname(T1_FILE), exist_ok=True)
    lines = [f"{RUN}: m_t = {MT} GeV, Gamma_t = {WT} GeV"]
    log = lambda s: (print(s), lines.append(s))
    Eg, Pg = E0 + DE * np.arange(NE), P0 + DP * np.arange(int(round(80 / DP)) + 1)
    pp = np.maximum(Pg, 1e-3)
    t0 = time.time()
    pot = JKTPotential()
    S = LSSolver(pot, n=120, m_t=MT, Gamma_t=WT)
    write_table(T1_FILE, np.array([S.G(E, pp) / Constants.G0(E, pp, MT, WT) for E in Eg]))
    rng = np.random.default_rng(7)
    Es, ps = rng.uniform(-20, 4, 200), rng.uniform(0.2, 80, 200)
    direct = np.array([abs(S.G(E, np.array([p]))[0] / Constants.G0(E, p, MT, WT)) ** 2 for E, p in zip(Es, ps)])
    dt = np.abs(np.abs(interp_table(load_table(T1_FILE), Es, ps)) ** 2 / direct - 1)
    log(f"singlet full JKT (alpha_s(mZ) = {pot.alpha_mz}, two-loop, NLO bracket; shift {pot.shift:.4f} GeV): "
        f"{os.path.basename(T1_FILE)} ({time.time() - t0:.0f} s); table vs solver median {np.median(dt):.1e}, max {dt.max():.1e}")
    a8 = float(ct25_alphas(OCTET_Q ** 2))
    EE, PP = np.meshgrid(Eg, pp, indexing="ij")
    write_table(T8_FILES["coulomb_as25"], octet_ratio(EE, PP, a8, MT, WT))
    log(f"octet coulomb_as25: exact repulsive Coulomb, alpha_s = CT25(25 GeV) = {a8:.6f}")
    for E, p in ((-2.0, 5.0), (0.0, 10.0), (4.0, 30.0)):
        r1 = abs(interp_table(load_table(T1_FILE), np.array(E), np.array(p))) ** 2
        r8 = abs(interp_table(load_table(T8_FILES["coulomb_as25"]), np.array(E), np.array(p))) ** 2
        log(f"  E = {E:+.0f}, p = {p:4.1f}: |G1/G0|^2 {r1:.4f}, |G8/G0|^2 {r8:.4f}")
    open(os.path.join(os.path.dirname(T1_FILE), "table_log.txt"), "w").write("\n".join(lines) + "\n")


F_GREEN_FACTOR = """C BEGIN toponium addition: Fuks (16) for W <= 350 GeV, p* <= 50 GeV (and E >= -20 GeV)
C     colour-singlet part of |M|^2 with CF1 (Fuks (19)-(20)), same normalisation as MATRIX1
      MSING = 0.D0
      DO I = 1, NCOLOR
        ZTEMP = (0.D0,0.D0)
        DO J = 1, NCOLOR
          ZTEMP = ZTEMP + CF1(J,I)*JAMP(J,1)
        ENDDO
        MSING = MSING + DBLE(ZTEMP*DCONJG(JAMP(I,1)))
      ENDDO
C     |M|^2 -> |M|^2 + (|G/G0|^2 - 1) |M_1|^2 = |M_8|^2 + |G/G0|^2 |M_1|^2
      IF(.NOT.LOADED) CALL LOADGREEN
      EXXX = PMOM(0,1) - 2.0D0*173.D0
      MXXX = SQRT(PMOM(1,2)**2 + PMOM(2,2)**2 + PMOM(3,2)**2)
      IF(GMAS.LE.350.D0 .AND. MXXX.LE.50.D0 .AND. EXXX.GE.E0G) THEN
        CALL CALCGREEN(EXXX, MXXX, GREEN)
        MATRIX1 = MATRIX1 + (DBLE(GREEN*DCONJG(GREEN)) - 1.D0)*MSING
      ENDIF
C END toponium addition

"""


F_CF_SINGLET = ("C BEGIN toponium modification: colour-singlet projection, Fuks (19)-(20) (no octet)\n"
                "      DATA (CF(I,  1),I=  1,  2) /6.666666666666666D-01,\n"
                "     $ 6.666666666666666D-01/\n"
                "C     1 T(1,2,3,6)\n"
                "      DATA (CF(I,  2),I=  1,  2) /6.666666666666666D-01\n"
                "     $ ,6.666666666666666D-01/\n"
                "C     1 T(2,1,3,6)\n"
                "C END toponium modification\n")


def patch_singlet(path):
    """Run 2 (both samples): CF replaced by the singlet projection, so MATRIX1 = |M_1|^2 (idempotent)."""

    src = open(path).read()
    if "toponium modification" in src:
        return "already singlet"
    if not os.path.exists(path + ".orig"):
        shutil.copy(path, path + ".orig")
    new, n = re.subn(r"      DATA \(CF\(I,  1\),I=  1,  2\).*?C     1 T\(2,1,\d+,\d+\)\n", lambda _: F_CF_SINGLET, src,
                     count=1, flags=re.S)
    assert n == 1 and "5.333333333333333D+00" not in new.split("FUNCTION MATRIX1")[1].split("\n      END\n")[0], \
        "colour matrix not replaced"
    open(path, "w").write(new)
    return "singlet"


F_GREEN_FACTOR_OCTET = """C BEGIN toponium addition: singlet and octet Green's functions, W <= 350, p* <= 50 GeV (E >= -20 GeV)
C     colour-singlet part of |M|^2 with CF1 (Fuks (19)-(20)); the octet part is |M|^2 - |M_1|^2
      MSING = 0.D0
      DO I = 1, NCOLOR
        ZTEMP = (0.D0,0.D0)
        DO J = 1, NCOLOR
          ZTEMP = ZTEMP + CF1(J,I)*JAMP(J,1)
        ENDDO
        MSING = MSING + DBLE(ZTEMP*DCONJG(JAMP(I,1)))
      ENDDO
C     |M|^2 -> |G1/G0|^2 |M_1|^2 + |G8/G0|^2 |M_8|^2
      IF(.NOT.LOADED) CALL LOADGREEN
      EXXX = PMOM(0,1) - 2.0D0*173.D0
      MXXX = SQRT(PMOM(1,2)**2 + PMOM(2,2)**2 + PMOM(3,2)**2)
      IF(GMAS.LE.350.D0 .AND. MXXX.LE.50.D0 .AND. EXXX.GE.E0G) THEN
        CALL CALCGREEN(EXXX, MXXX, GREEN)
        CALL CALCGREEN8(EXXX, MXXX, GREEN8)
        MATRIX1 = MATRIX1 + (DBLE(GREEN*DCONJG(GREEN)) - 1.D0)*MSING
     $   + (DBLE(GREEN8*DCONJG(GREEN8)) - 1.D0)*(MATRIX1 - MSING)
      ENDIF
C END toponium addition

"""

F_GREEN_FACTOR_QQ = """C BEGIN toponium addition: qqbar -> ttbar is pure colour octet at LO: |M|^2 -> |G8/G0|^2 |M|^2
C     for W <= 350, p* <= 50 GeV (E >= -20 GeV)
      IF(.NOT.LOADED) CALL LOADGREEN
      EXXX = PMOM(0,1) - 2.0D0*173.D0
      MXXX = SQRT(PMOM(1,2)**2 + PMOM(2,2)**2 + PMOM(3,2)**2)
      IF(GMAS.LE.350.D0 .AND. MXXX.LE.50.D0 .AND. EXXX.GE.E0G) THEN
        CALL CALCGREEN8(EXXX, MXXX, GREEN8)
        MATRIX1 = DBLE(GREEN8*DCONJG(GREEN8))*MATRIX1
      ENDIF
C END toponium addition

"""


def rw_region(text):
    """Run 5: the |G/G0|^2 region of the Fortran above becomes p* < 80 GeV (the student's: R = 1 for p* >= 80);
    runs 1-4 keep p* <= 50 unchanged."""

    if "mt" in CFG:                  # run 7: E = W - 2 m_t and W <= 2 m_t + 4 at this run's m_t
        assert text.count("2.0D0*173.D0") == 1 and text.count("GMAS.LE.350.D0") == 1
        text = (text.replace("2.0D0*173.D0", f"2.0D0*{MT:g}D0").replace("GMAS.LE.350.D0", f"GMAS.LE.{whi():g}D0")
                .replace("W <= 350", f"W <= {whi():g}"))
    if not CFG.get("p_rw"):
        return text
    assert text.count("MXXX.LE.50.D0") == 1
    return (text.replace("MXXX.LE.50.D0", f"MXXX.LT.{CFG['p_rw']:.0f}.D0")
            .replace("p* <= 50 GeV", f"p* < {CFG['p_rw']:.0f} GeV"))


def patch_qq(path):
    """Run 3 WITH: the qqbar subprocess (pure octet) gets |G8/G0|^2 in the window (idempotent)."""

    src = open(path).read()
    if "toponium addition" in src:
        return "already patched"
    if not os.path.exists(path + ".orig"):
        shutil.copy(path, path + ".orig")
    head, sep, body = src.partition("      REAL*8 FUNCTION MATRIX1(P,NHEL,IC, IHEL)\n")
    assert sep, "MATRIX1 not found"
    end = body.index("\n      END\n") + len("\n      END\n")
    m1, rest = body[:end], body[end:]
    for old, new in (("      INTEGER IHEL\n", "      INTEGER IHEL\n" + F_DECLARATIONS),
                     ("      MATRIX1 = 0.D0\n", F_KINEMATICS + "      MATRIX1 = 0.D0\n"),
                     ("      IF(SDE_STRAT.EQ.1)THEN\n", rw_region(F_GREEN_FACTOR_QQ) + "      IF(SDE_STRAT.EQ.1)THEN\n")):
        assert m1.count(old) == 1, f"anchor not unique in qq MATRIX1: {old!r}"
        m1 = m1.replace(old, new)
    open(path, "w").write(F_MODULE(octet=True) + head + sep + F_USE + m1 + rest)
    return "patched (octet)"


def patch_matrix1(path):
    """WITH: Fuks Sec. 3 additions in MATRIX1 of the gg subprocess (idempotent; original kept as matrix1.f.orig).
    With the singlet CF (run 2) MATRIX1 = MSING, so the same lines give |G/G0|^2 |M_1|^2."""

    src = open(path).read()
    if "toponium addition" in src:
        return "already patched"
    if not os.path.exists(path + ".orig"):
        shutil.copy(path, path + ".orig")
    head, sep, body = src.partition("      REAL*8 FUNCTION MATRIX1(P,NHEL,IC, IHEL)\n")
    assert sep, "MATRIX1 not found"

    def once(text, old, new):
        assert text.count(old) == 1, f"anchor not unique in MATRIX1: {old!r}"
        return text.replace(old, new)

    end = body.index("\n      END\n") + len("\n      END\n")
    m1, rest = body[:end], body[end:]
    m1 = F_USE + m1                                                     # right after the FUNCTION line
    m1 = once(m1, "      INTEGER IHEL\n", "      INTEGER IHEL\n" + F_DECLARATIONS)
    m1, n_cf = re.subn(r"(C     1 T\(2,1,\d+,\d+\)\n)", lambda m_: m_.group(1) + F_CF1, m1, count=1)
    assert n_cf == 1, "colour-flow comment T(2,1,..) not found in MATRIX1"   # (3,6) here; (5,8) in decay chains
    m1 = once(m1, "      MATRIX1 = 0.D0\n", F_KINEMATICS + "      MATRIX1 = 0.D0\n")
    factor = F_GREEN_FACTOR_OCTET if CFG["octet"] else F_GREEN_FACTOR
    m1 = once(m1, "      IF(SDE_STRAT.EQ.1)THEN\n", rw_region(factor) + "      IF(SDE_STRAT.EQ.1)THEN\n")
    assert ("5.333333333333333D+00" in m1) != CFG["singlet_only"], "colour matrix not as this run expects"
    open(path, "w").write(F_MODULE(octet=CFG["octet"]) + head + sep + m1 + rest)
    return "patched"


F_WINDOW_DECL = """C BEGIN toponium addition: generation window (run 5)
      REAL*8 PT1(0:3), PT2(0:3), W2, M1S, M2S, XLAM
      INTEGER M
C END toponium addition
"""

F_WINDOW_CUT = """C BEGIN toponium addition: generation window W <= 350 GeV, p* <= 100 GeV
C     t = P(.,3..5) = l+ vl b, tbar = P(.,6..8) = l- vl~ b~ (as F_KINEMATICS)
C     p* = sqrt(lambda(W^2, m_t1^2, m_t2^2))/(2W) <= 100  <=>  lambda <= 4 W^2 100^2
      DO M = 0, 3
        PT1(M) = P(M,3) + P(M,4) + P(M,5)
        PT2(M) = P(M,6) + P(M,7) + P(M,8)
      ENDDO
      M1S = PT1(0)**2 - PT1(1)**2 - PT1(2)**2 - PT1(3)**2
      M2S = PT2(0)**2 - PT2(1)**2 - PT2(2)**2 - PT2(3)**2
      W2 = (PT1(0)+PT2(0))**2 - (PT1(1)+PT2(1))**2
     $   - (PT1(2)+PT2(2))**2 - (PT1(3)+PT2(3))**2
      XLAM = (W2 - M1S - M2S)**2 - 4D0*M1S*M2S
      IF (W2 .GT. 350D0**2 .OR. XLAM .GT. 4D0*W2*100D0**2) THEN
        dummy_cuts = .false.
        RETURN
      ENDIF
C END toponium addition
"""


def patch_window(sub, tau=True, wmax=None):
    """Run 5, all four outputs: events only in W <= 350 GeV (tau window in genps.f GENCMS, as Stage 4; no lower
    edge) and p* <= 100 GeV (dummy_cuts). tau=False (run 6 arbiter): dummy_cuts only, genps.f untouched (W <= 350
    then by the cut alone). wmax (run 6 zoom test): dummy_cuts W <= wmax only, no p* cut. Idempotent; originals kept
    as *.orig."""

    cut = F_WINDOW_CUT if wmax is None else (
        F_WINDOW_CUT.replace("350D0**2", f"{wmax:.0f}D0**2").replace(" .OR. XLAM .GT. 4D0*W2*100D0**2", "")
        .replace("W <= 350 GeV, p* <= 100 GeV", f"W <= {wmax:.0f} GeV (run 6 zoom test, no p* cut)"))
    assert wmax is None or ("XLAM .GT." not in cut and f"{wmax:.0f}D0**2" in cut)
    w = f"{whi():g}D0" if "mt" in CFG else "350D0"         # run 7: W <= 2 m_t + 4 at this run's m_t
    if "mt" in CFG:
        assert wmax is None and cut.count("350D0**2") == 1
        cut = cut.replace("350D0**2", f"{w}**2").replace("W <= 350 GeV", f"W <= {whi():g} GeV")

    done = []
    for f, edits in ((("genps.f", [("      TAUMAX = 1D0\n",
                                     f"      TAUMAX = {w}**2/S ! toponium window W <= {w}\n")]),) * tau +
                     (("dummy_fct.f", [("     . IS_A_ONIUM, do_cuts\n", "     . IS_A_ONIUM, do_cuts\n" + F_WINDOW_DECL),
                                      ("      dummy_cuts=.true.\n", "      dummy_cuts=.true.\n" + cut)]),)):
        path = os.path.join(sub, f)
        src = open(path).read()
        if "toponium" in src:
            continue
        if not os.path.exists(path + ".orig"):
            shutil.copy(path, path + ".orig")
        for old, new in edits:
            assert src.count(old) == 1, f"anchor not unique in {f}: {old!r}"
            src = src.replace(old, new)
        open(path, "w").write(src)
        done.append(f)
    return "window patched: " + (", ".join(done) or "already")


def write_text(path, text, mode=None):
    """Write text to path only if it differs, atomically (running MadGraph jobs read lhapdf_data/ meanwhile)."""

    if os.path.exists(path) and not os.path.islink(path) and open(path).read() == text:
        return
    open(path + ".tmp", "w").write(text)
    if mode:
        os.chmod(path + ".tmp", mode)
    os.replace(path + ".tmp", path)


def put(src, dst):
    """Copy src to dst only if they differ, atomically (concurrent cluster jobs read these tables)."""

    import filecmp
    if os.path.exists(dst) and filecmp.cmp(src, dst, shallow=False):
        return
    shutil.copy(src, dst + ".tmp")
    os.replace(dst + ".tmp", dst)


def set_card(path, values):
    """Set 'value = key' entries of a MadGraph run_card."""

    src = open(path).read()
    for key, val in values.items():
        pat = re.compile(rf"^(\s*)\S+(\s*=\s*{re.escape(key)}\b)", re.M)
        src, n = pat.subn(lambda m: f"{m.group(1)}{val}{m.group(2)}", src, count=1)
        assert n == 1, f"{key} not found in {path}"
    open(path, "w").write(src)


def set_config(path, values):
    """me5_configuration.txt: 'key = value' lines (commented or not)."""

    src = open(path).read()
    for key, val in values.items():
        pat = re.compile(rf"^#?\s*{re.escape(key)}\s*=.*$", re.M)
        src, n = pat.subn(f"{key} = {val}", src, count=1)
        assert n == 1, f"{key} not found in {path}"
    open(path, "w").write(src)


def mg_env():
    return dict(os.environ, LHAPDF_DATA_PATH=LHAPDF_DIR)


def setup():
    # LHAPDF: local data dir with CT25NNLO under a 32-bit id
    os.makedirs(LHAPDF_DIR, exist_ok=True)
    real = shutil.which("lhapdf-config")
    sysdir = subprocess.run([real, "--datadir"], capture_output=True, text=True).stdout.strip()
    if sysdir == LHAPDF_DIR:
        sys.exit("lhapdf-config in PATH is the local wrapper; run with the Homebrew one first in PATH")
    index = open(os.path.join(sysdir, "pdfsets.index")).read().rstrip("\n")
    write_text(os.path.join(LHAPDF_DIR, "pdfsets.index"), index + f"\n{CT25_ID} CT25NNLO 1\n")
    # CT25NNLO: a full copy of the distributed set (so it does not depend on where that is kept); the .info gets
    # the keys LHAPDF 6.5 needs and the file lacks: FlavorScheme/NumFlavors for alpha_s set to its own values;
    # AlphaS_Lambda4/5 are read by MadGraph's LHAGLUE layer into legacy common blocks it never uses (alpha_s
    # comes from the set's own ipol table) and are filled from JKT (14) at the set's alpha_s(mZ), mb = 4.75
    src = next((d for d in CT25_SOURCES if os.path.exists(os.path.join(d, "CT25NNLO.info"))), None)
    if src is None:
        sys.exit(f"CT25NNLO not found in {CT25_SOURCES}")
    ct25 = os.path.join(LHAPDF_DIR, "CT25NNLO")
    if os.path.islink(ct25):
        os.remove(ct25)
    os.makedirs(ct25, exist_ok=True)
    for f in glob.glob(os.path.join(src, "*.dat")):
        dst = os.path.join(ct25, os.path.basename(f))
        if os.path.islink(dst):
            os.remove(dst)
        if not os.path.exists(dst):
            shutil.copy(f, dst)
    info = open(os.path.join(src, "CT25NNLO.info")).read().rstrip("\n")
    amz = float(re.search(r"^AlphaS_MZ:\s*(\S+)", info, re.M).group(1))
    lam5 = JKTPotential.lambda_msbar_2loop(amz, mZ=91.187, nF=5)
    lam4 = JKTPotential.lambda_msbar_2loop(float(JKTPotential.alpha_2loop(np.array(4.75 ** 2), lam5, nF=5)),
                                           mZ=4.75, nF=4)
    write_text(os.path.join(ct25, "CT25NNLO.info"),
        info + "\n# added for LHAPDF 6.5 / MadGraph (not in the distributed file; Lambdas unused by MadGraph):\n"
        f"AlphaS_FlavorScheme: variable\nAlphaS_NumFlavors: 5\nAlphaS_Lambda4: {lam4:.4f}\nAlphaS_Lambda5: {lam5:.4f}\n")
    dst = os.path.join(LHAPDF_DIR, "lhapdf.conf")
    if not os.path.lexists(dst):
        os.symlink(os.path.join(sysdir, "lhapdf.conf"), dst)
    # MadGraph reads the set list from `lhapdf-config --datadir` (not LHAPDF_DATA_PATH): a wrapper that
    # reports this dir and passes everything else to the real lhapdf-config
    wrapper = os.path.join(LHAPDF_DIR, "lhapdf-config")
    # the real lhapdf-config by name, from the PATH of whoever calls it: on the MSU HPCC each node's `module load`
    # gives the build for its own CPU (a fixed path would pin one CPU type: SIGILL on the others)
    write_text(wrapper, f'#!/bin/bash\nif [ "$1" = "--datadir" ]; then echo {LHAPDF_DIR}; '
                        f'else exec lhapdf-config "$@"; fi\n', mode=0o755)
    print(f"LHAPDF data dir {LHAPDF_DIR}: CT25NNLO (alpha_s(mZ) = {amz}) as id {CT25_ID}")
    if "mt" in CFG and not (os.path.exists(T1_FILE) and all(map(os.path.exists, T8_FILES.values()))):
        table_mt()                    # run 7: G tables at this m_t, Gamma_t
    if CFG.get("window"):
        setup_window(wrapper)
        if CFG["without_from"]:
            return
        if not CFG.get("check"):      # run 7: only its own SM_full besides the window outputs
            select(RUN)
            if not os.path.isdir(SM_DIR):
                mg_output(SM_DIR, f"generate {CFG['process']}")
            set_cards(SM_DIR, wrapper)
            print(f"SM_full {SM_DIR}: cards m_t = {MT}, Gamma_t = {WT}")
            return
        select(RUN)   # run 6: also its own SM_full (WITHOUT) and the check topo_full (WITH, full phase space), as runs 1-4 below

    # MadGraph outputs (WITHOUT and WITH: the same process)
    if CFG["archived"]:
        sys.exit(f"{RUN} is archived (its WITH matrix1.f reads the G table from its old path): not set up again")
    own = [TOPO_DIR] if CFG["without_from"] else [SM_DIR, TOPO_DIR]   # never touch another run's directory
    for d in own:
        if os.path.isdir(d):
            print(f"{d} exists")
            continue
        mg_output(d, f"generate {CFG['process']}")

    # matrix elements: run 2 = singlet only in both samples; WITH = Fuks (16) in the gg subprocess (qqbar unchanged)
    # every matrix file of every gg / qq subprocess folder (several lepton combinations may share or split folders)
    mats = lambda d, kind: sorted(f for f in glob.glob(os.path.join(d, "SubProcesses", f"P*_{kind}_*", "matrix*.f")))
    folders = [p for p in glob.glob(os.path.join(TOPO_DIR, "SubProcesses", "P*")) if os.path.isdir(p)]
    stray = [p for p in folders if not re.match(r"P\d+_(gg|qq)_", os.path.basename(p))]
    assert not stray, f"subprocess folders the patch would miss: {stray}"
    for kind in ("gg", "qq"):
        for f in mats(TOPO_DIR, kind):
            assert os.path.basename(f) == "matrix1.f", f"{f}: only one matrix file per folder is supported"
    assert mats(TOPO_DIR, "gg"), "no gg subprocess"
    if CFG["singlet_only"]:
        for f in mats(SM_DIR, "gg") + mats(TOPO_DIR, "gg"):
            print(f"singlet colour {os.path.relpath(f, OUTPUT)}:", patch_singlet(f))
    for f in mats(TOPO_DIR, "gg"):
        print(f"WITH {os.path.relpath(f, OUTPUT)}:", patch_matrix1(f))
    if CFG["octet"]:
        for f in mats(TOPO_DIR, "qq"):
            print(f"WITH {os.path.relpath(f, OUTPUT)}:", patch_qq(f))
    table, table8 = (TABLE80, OCTET80[CFG["check"]]) if CFG.get("check") else (TABLE, OCTET_TABLE)   # run 6: p <= 80
    if os.path.exists(table):
        shutil.copy(table, GREEN_FILE)
    if CFG["octet"] and os.path.exists(table8):
        shutil.copy(table8, GREEN8_FILE)

    # cards (identical in WITHOUT and WITH)
    for d in own:
        set_cards(d, wrapper)
    print("cards: 13 TeV, CT25NNLO, no cuts, m_t = 173, Gamma_t = 1.49, event_norm = average, multicore 10")


def set_cards(d, wrapper):
    """Cards of every run (WITHOUT and WITH identical): 13 TeV, CT25NNLO, no cuts, m_t = MT, Gamma_t = WT (173, 1.49
    except run 7)."""

    set_card(os.path.join(d, "Cards", "run_card.dat"), {
        "pdlabel": "lhapdf", "lhaid": CT25_ID, "use_syst": "False", "event_norm": "average",
        "ebeam1": EBEAM, "ebeam2": EBEAM, "ptl": "0.0", "etal": "-1.0", "drll": "0.0"})   # no cuts
    pc = os.path.join(d, "Cards", "param_card.dat")
    src, n = re.subn(r"^DECAY\s+6\s+\S+\s+# WT", f"DECAY   6 {WT:.6e} # WT", open(pc).read(), flags=re.M)
    src, m = re.subn(r"^(\s+6 )\S+( # MT)", lambda g: f"{g.group(1)}{MT:.6e}{g.group(2)}", src, flags=re.M)
    assert n == 1 and m == 1, "param_card MT/WT"   # m_t, Gamma_t of the run (173, 1.49 except run 7)
    open(pc, "w").write(src)
    set_config(os.path.join(d, "Cards", "me5_configuration.txt"), {
        "run_mode": 2, "nb_core": int(os.environ.get("SLURM_CPUS_PER_TASK", 10)), "lhapdf_py3": wrapper, "automatic_html_opening": "False"})


def mg_output(d, process):
    """MadGraph output d of `process` (mg5 command lines), helicity recycling off (MATRIX1 must stay ours)."""

    os.makedirs(os.path.dirname(d), exist_ok=True)
    script = os.path.join(OUTPUT, f"_{os.path.basename(d)}.mg5")
    open(script, "w").write(f"set automatic_html_opening False\nimport model sm\n{process}\n"
                            f"output {d} --hel_recycling=False\n")
    subprocess.run([sys.executable, "-W", "ignore", os.path.join(MG_ROOT, "bin", "mg5_aMC"), script],
                   cwd=MG_ROOT, check=True, env=mg_env())
    os.remove(script)


def setup_window(wrapper):
    """Runs 5-6: SM_win (TT7_S_*) and topo_win (TT7_T_*) outputs for gg and qqbar, each with the generation window; topo_win gets the Stage 6 MATRIX1 patches
    (gg: singlet x |G1/G0|^2 + octet x |G8/G0|^2; qqbar: x |G8/G0|^2) for p* < 80 GeV; cards as run 4, and as the
    student mc_grouped_subproc = False in the gg runs."""

    global GREEN_FILE, GREEN8_FILE
    for ch in ("gg", "qq"):
        for kind in ("S", "T"):
            d = window_dir(kind, ch)
            if not os.path.isdir(d):
                mg_output(d, CFG[f"process_{ch}"])
            sub = os.path.join(d, "SubProcesses")
            folders = [os.path.basename(p) for p in glob.glob(os.path.join(sub, "P*")) if os.path.isdir(p)]
            assert folders and all(re.match(rf"P\d+_{ch}_", f) for f in folders), f"{d}: subprocesses {folders}"
            print(f"{kind}_{ch}: {patch_window(sub)}")
            if kind == "T":
                GREEN_FILE = os.path.join(sub, "green_table.dat")
                GREEN8_FILE = os.path.join(sub, "green8_table.dat")
                for f in sorted(glob.glob(os.path.join(sub, "P*", "matrix*.f"))):
                    assert os.path.basename(f) == "matrix1.f", f"{f}: only one matrix file per folder is supported"
                    print(f"T_{ch} {os.path.relpath(f, OUTPUT)}:", (patch_matrix1 if ch == "gg" else patch_qq)(f))
                if os.path.exists(T1_FILE):
                    put(T1_FILE, GREEN_FILE)
                if "mt" in CFG:               # run 7: its one octet table, in place before any job starts
                    put(T8_FILES["coulomb_as25"], GREEN8_FILE)
            if ch == "qq" and "add process q~ q" in CFG["process_qq"]:   # the beam-swapped mirror must be there
                ids = set()
                for f in glob.glob(os.path.join(sub, "P*", "leshouche.inc")):
                    ids |= {int(x.split(",")[0]) for x in re.findall(r"IDUP\(I,\s*\d+,\s*\d+\),I=1,\d+\)\s*/\s*([-\d,\s]+)/", open(f).read())}
                assert ids & {1, 2, 3, 4} and ids & {-1, -2, -3, -4}, f"{d}: no q~ q mirror (beam-1 ids {sorted(ids)})"
                print(f"{kind}_qq: beam-1 ids {sorted(ids)} (q q~ and the mirror q~ q)")
            set_cards(d, wrapper)
            rc = os.path.join(d, "Cards", "run_card.dat")
            if ch == "gg" and "mc_grouped_subproc" not in open(rc).read():   # hidden parameter: add the line
                open(rc, "a").write(" False = mc_grouped_subproc ! as the student's gg runs\n")
    GREEN_FILE = GREEN8_FILE = None
    print(f"cards: as run 4 (13 TeV, CT25NNLO, no cuts, m_t = {MT}, Gamma_t = {WT}, event_norm = average)")


# =====================================================================================
# 3. runs: generate, keep per-event t and tbar 4-momenta, delete the LHE file
# =====================================================================================
def ev_path(name, smoke):
    return os.path.join(EV_DIR, "smoke" if smoke else "", f"{name}.npz")


def do_runs(smoke=False, only=None):
    """All samples of runs(smoke) in turn, or only the one named `only` (one sample per cluster job)."""

    if CFG.get("window"):
        for ch in ("gg", "qq"):
            sub = os.path.join(window_dir("T", ch), "SubProcesses")
            put(T1_FILE, os.path.join(sub, "green_table.dat"))
        if CFG.get("check"):                          # run 6's F: singlet table here, octet table per run below
            put(T1_FILE, GREEN_FILE)
    elif not os.path.exists(GREEN_FILE) or (CFG["octet"] and not os.path.exists(GREEN8_FILE)):
        sys.exit("no G table in the WITH directory: run 'table' after 'setup'")
    os.makedirs(os.path.join(EV_DIR, "smoke" if smoke else "", "logs"), exist_ok=True)
    todo = [r for r in runs(smoke) if only in (None, r[0])]
    if not todo:
        sys.exit(f"no sample {only}; samples: {[r[0] for r in runs(smoke)]}")
    for name, d, seed, nev, octet in todo:
        if CFG.get("window") and d in [window_dir(k, c) for k in "ST" for c in ("gg", "qq")]:
            sample = octet or "SM_win"                # window run (run 5's TOPO_DIR is its topo_win_gg dir)
        else:
            sample = "with" if d == TOPO_DIR else "without"
        generate(name, d, seed, nev, octet, sample, ev_path(name, smoke))


def generate(name, d, seed, nev, octet, sample, out):
    """One MadGraph run in directory d -> npz out (w, t, tb per event); skipped if out exists; LHE deleted.
    In a SLURM job the run happens in a private copy d__<name>__<job id> (deleted afterwards), so several jobs can use
    the same output at once, even two jobs of the same sample; MATRIX1 reads its G tables from d itself (absolute paths), which nothing changes during runs."""

    if os.path.exists(out):
        print(f"{name}: done")
        return
    if octet:                                     # runs 5-7: MATRIX1 reads this octet table at run time
        put(T8_FILES[octet], os.path.join(d, "SubProcesses", "green8_table.dat"))
    jid = os.environ.get("SLURM_JOB_ID", "")
    job = bool(jid)
    if job:
        jd = f"{d}__{name}__{jid}"                # own job id: never another running job's copy
        if os.path.isdir(jd):
            shutil.rmtree(jd)                     # this job, requeued after preemption
        shutil.copytree(d, jd, symlinks=True)
        d = jd
    try:
        set_card(os.path.join(d, "Cards", "run_card.dat"), {"iseed": seed, "nevents": nev})
        if job:                                   # as many MadGraph cores as SLURM gave the job
            set_config(os.path.join(d, "Cards", "me5_configuration.txt"),
                       {"nb_core": int(os.environ.get("SLURM_CPUS_PER_TASK", 1))})
        if os.path.isdir(os.path.join(d, "Events", name)):
            shutil.rmtree(os.path.join(d, "Events", name))           # an interrupted earlier attempt
        t0 = time.time()
        log = os.path.join(os.path.dirname(out), "logs", f"{name}{'_' + jid if job else ''}.log")
        print(f"{name}: generating {nev} events (seed {seed}) in {d} ...", flush=True)
        with open(log, "w") as f:
            r = subprocess.run([sys.executable, "-W", "ignore", "bin/generate_events", name, "-f"],
                               cwd=d, stdout=f, stderr=subprocess.STDOUT, env=mg_env())
        lhe = os.path.join(d, "Events", name, "unweighted_events.lhe.gz")
        if r.returncode != 0 or not os.path.exists(lhe):
            sys.exit(f"{name}: MadGraph failed, see {log}")
        w, pt, ptb = Events.read_lhe(lhe)
        tmp = f"{out}.{jid or 'local'}.tmp.npz"
        np.savez(tmp, w=w, t=pt, tb=ptb, seed=seed, sample=sample, run=RUN, mt=MT, wt=WT)
        os.replace(tmp, out)
        for f in glob.glob(os.path.join(d, "Events", name, "*.lhe*")):
            os.remove(f)
        print(f"{name}: {len(w)} events, sigma = {w.mean():.5g} pb, {(time.time() - t0) / 60:.1f} min", flush=True)
    finally:
        if job:
            shutil.rmtree(d, ignore_errors=True)

ARB_RUNS = [("gg", 7701), ("qq", 7711)]   # run 6 arbiter: SM, window by dummy_cuts only (genps.f untouched)
ARB_EVENTS = 500000
ZOOM = dict(wmax=500.0, seed=7801, events=1000000, dir="TT9_SMcut500_ll", name="SM_cut500_ll")   # pp, W <= 500 only


def arbiter_dir(ch):
    return os.path.join(OUTPUT, RUN, f"TT9_SMcut_{ch}")


def arbiter():
    """Run 6 arbiter of SM_win vs SM_full: SM_cut_gg, SM_cut_qq = the SM_win processes and cards with the generation
    window imposed only by dummy_cuts (W <= 350 GeV, p* <= 100 GeV, F_WINDOW_CUT), the tau window of genps.f NOT
    applied, so MadGraph samples the full tau range; SM_cut500 = SM_full's p p process and cards with only W <= 500 GeV
    (dummy_cuts). In the selection window (326 <= W <= 350, p* < 100 or 50)
    sigma(SM_cut) is compared with sigma(SM_full) and sigma(SM_win), in total and per W and p* bin."""

    wrapper = os.path.join(LHAPDF_DIR, "lhapdf-config")
    for ch, seed in ARB_RUNS:
        d = arbiter_dir(ch)
        if not os.path.isdir(d):
            mg_output(d, CFG[f"process_{ch}"])
            print(f"SM_cut_{ch}: {patch_window(os.path.join(d, 'SubProcesses'), tau=False)}")
            set_cards(d, wrapper)
            rc = os.path.join(d, "Cards", "run_card.dat")
            if ch == "gg" and "mc_grouped_subproc" not in open(rc).read():   # as SM_win_gg
                open(rc, "a").write(" False = mc_grouped_subproc ! as the student's gg runs\n")
        os.makedirs(os.path.join(EV_DIR, "logs"), exist_ok=True)
        generate(f"SM_cut_{ch}_s{seed}", d, seed, ARB_EVENTS, None, "SM_cut", ev_path(f"SM_cut_{ch}_s{seed}", False))
    # zoom test: SM_full's own process and cards (p p, grouped), only W <= 500 GeV imposed (dummy_cuts): the
    # selection window becomes a few % of the sample instead of 0.7%
    d = os.path.join(OUTPUT, RUN, ZOOM["dir"])
    if not os.path.isdir(d):
        mg_output(d, f"generate {CFG['process']}")
        print(f"{ZOOM['name']}: {patch_window(os.path.join(d, 'SubProcesses'), tau=False, wmax=ZOOM['wmax'])}")
        set_cards(d, wrapper)
    zname = f"{ZOOM['name']}_s{ZOOM['seed']}"
    generate(zname, d, ZOOM["seed"], ZOOM["events"], None, "SM_cut500", ev_path(zname, False))

    lines = [f"STAGE 6, {RUN}: arbiter, SM_full vs SM_win vs SM_cut vs SM_cut500 in the selection window (sigma in fb,"
             " pulls in sigma)", ""]
    log = lambda s_="": (print(s_), lines.append(s_))
    grp = {"SM_full": [load_group(f"{CFG['names'][0]}_s*.npz")],
           "SM_win": [load_group(window_name(ch, None, s) + ".npz") for ch, o, s in CFG["window_runs"] if o is None],
           "SM_cut": [load_group(f"SM_cut_{ch}_s{s}.npz") if ch == "gg" else reflect_z(load_group(f"SM_cut_{ch}_s{s}.npz"))
                      for ch, s in ARB_RUNS],   # SM_cut_qq: q q~ only, + its z-reflection (the mirror)
           "SM_cut500": [load_group(f"{zname}.npz")]}
    for k, gs in grp.items():
        for g in gs:
            E, ps = E_pstar(g["t"], g["tb"])
            g["W"], g["ps"] = E + 2 * MT, ps
        log(f"{k:8s} " + ", ".join(f"{len(g['files'])} run(s) {g['files'][0]}.. sigma {g['w'].sum():.6g} pb ({len(g['w'])} events)"
                                     for g in gs))
    def sig(m):                                        # {sample: (sigma, err)} [fb] of the events passing m(W, p*)
        out = {}
        for k, gs in grp.items():
            x = [(g["w"][m(g["W"], g["ps"])].sum(), (g["w"][m(g["W"], g["ps"])] ** 2).sum()) for g in gs]
            out[k] = (1e3 * sum(a for a, _ in x), 1e3 * np.sqrt(sum(b for _, b in x)))
        return out
    pull = lambda a, b: (a[0] - b[0]) / np.hypot(a[1], b[1])
    def row(lab, r):
        f, w, c, z = r["SM_full"], r["SM_win"], r["SM_cut"], r["SM_cut500"]
        log(f"  {lab:12s} SM_full {f[0]:7.3f}+-{f[1]:.3f}  SM_win {w[0]:7.3f}+-{w[1]:.3f}  SM_cut {c[0]:7.3f}+-{c[1]:.3f}  "
            f"SM_cut500 {z[0]:7.3f}+-{z[1]:.3f} fb | cut/win {c[0] / w[0]:.3f} ({pull(c, w):+5.1f})  "
            f"cut500/full {z[0] / f[0]:.3f} ({pull(z, f):+5.1f})  cut500/win {z[0] / w[0]:.3f} ({pull(z, w):+5.1f})")
    for pcut in (100.0, 50.0):
        log("")
        log(f"selection window 326 <= W <= 350 GeV, p* < {pcut:.0f} GeV:")
        row("total", sig(lambda W, p: (W >= 326) & (W <= 350) & (p < pcut)))
        for a, b in zip([326, 332, 338, 342, 344, 346, 348], [332, 338, 342, 344, 346, 348, 350]):
            row(f"W {a}-{b}", sig(lambda W, p: (W >= a) & (W < b + (b == 350) * 1e-9) & (p < pcut)))
    log("")
    log("p* bins, 326 <= W <= 350:")
    for a in range(0, 100, 10):
        row(f"p* {a}-{a + 10}", sig(lambda W, p: (W >= 326) & (W <= 350) & (p >= a) & (p < a + 10)))
    os.makedirs(RES_RUN, exist_ok=True)
    open(os.path.join(RES_RUN, "arbiter.txt"), "w").write("\n".join(lines) + "\n")


# =====================================================================================
# 4. observables, T-factors, plots
# =====================================================================================
def pT(p):
    return np.hypot(p[:, 1], p[:, 2])


def rap(p):
    return 0.5 * np.log((p[:, 0] + p[:, 3]) / (p[:, 0] - p[:, 3]))


def observables(t, tb):
    tt = t + tb
    return {"m_tt": np.sqrt(np.maximum(Events.minv2(tt), 0)), "pT_t": pT(t), "absy_t": np.abs(rap(t)),
            "y_tt": rap(tt), "absy_tt": np.abs(rap(tt)), "absyboost_tt": np.abs(rap(t) + rap(tb)) / 2,
            "HT_tt": pT(t) + pT(tb)}


def E_pstar(t, tb):
    """E = W - 2 m_t and p* = sqrt(lambda(W^2, m_t1^2, m_t2^2))/(2W) (Lorentz-invariant form of Fuks (21))."""

    m1s, m2s, W2 = Events.minv2(t), Events.minv2(tb), Events.minv2(t + tb)
    W = np.sqrt(W2)
    lam = (W2 - m1s - m2s) ** 2 - 4 * m1s * m2s
    return W - 2 * MT, np.sqrt(np.maximum(lam, 0)) / (2 * W)


def fuks_pstar(t, tb):
    """Fuks (21) transcribed from F_KINEMATICS (rotations about z and y, boost along z): |p_t| in the
    t-tbar rest frame. Used only to check E_pstar against the code that runs inside MATRIX1."""

    P = t + tb
    gmom = np.sqrt(P[:, 1] ** 2 + P[:, 2] ** 2 + P[:, 3] ** 2)
    gmas = np.sqrt(P[:, 0] ** 2 - gmom ** 2)
    the, phi = np.arccos(P[:, 3] / (gmom + 1e-8)), np.arctan2(P[:, 2], P[:, 1])
    bet, gam = gmom / P[:, 0], P[:, 0] / gmas
    x, y, z, e = t[:, 1], t[:, 2], t[:, 3], t[:, 0]
    x, y = np.cos(phi) * x + np.sin(phi) * y, np.cos(phi) * y - np.sin(phi) * x
    x, z = np.cos(the) * x - np.sin(the) * z, np.cos(the) * z + np.sin(the) * x
    e, z = gam * e - gam * bet * z, gam * z - gam * bet * e
    return np.sqrt(x ** 2 + y ** 2 + z ** 2)


DISTRIBUTIONS = {   # name -> [(observable, edges)] (2 entries = 2D); all binnings as provided, in the order given
    "m_ttbar_345-1600": [("m_tt", [345, 400, 470, 550, 650, 800, 1100, 1600])],
    "pT_top_0-500": [("pT_t", [0, 60, 100, 150, 200, 260, 320, 400, 500])],
    "abs_y_top_0-2.5": [("absy_t", [0, 0.35, 0.85, 1.45, 2.5])],
    "pT_top_0-600": [("pT_t", [0, 80, 150, 250, 600])],
    "abs_y_top_x_pT_top": [("absy_t", [0, 0.35, 0.85, 1.45, 2.5]), ("pT_t", [0, 80, 150, 250, 600])],
    "abs_y_ttbar_0-2.4": [("absy_tt", [0.0, 0.12, 0.24, 0.36, 0.49, 0.62, 0.76, 0.91, 1.06, 1.21, 1.39, 1.59, 2.4])],
    "y_ttbar_-2.6-2.6": [("y_tt", [-2.6, -1.6, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.6, 2.6])],
    "m_ttbar_250-3500": [("m_tt", [250, 400, 480, 560, 640, 720, 800, 900, 1000, 1150, 1300, 1500, 1700, 2000,
                                   2300, 3500])],
    "m_ttbar_325-2000": [("m_tt", [325, 400, 480, 580, 700, 860, 1020, 1250, 1500, 2000])],
    "abs_y_ttbar_0-2.5": [("absy_tt", [0, 0.25, 0.5, 0.8, 1.1, 1.4, 1.8, 2.5])],
    "abs_y_boost_0-2.5": [("absyboost_tt", [0, 0.25, 0.5, 0.8, 1.05, 1.3, 1.6, 1.85, 2.15, 2.5])],
    "HT_ttbar_0-2000": [("HT_tt", [0, 90, 170, 260, 360, 470, 580, 710, 830, 2000])],
}


def load_group(pattern, base=None):
    """Concatenate runs; per-event weight w/N_total (event_norm = average: each LHE weight is that run's sigma)."""

    files = sorted(glob.glob(os.path.join(base or EV_DIR, pattern)))
    if not files:
        return None
    ds = [np.load(f) for f in files]
    n = sum(len(d["w"]) for d in ds)
    return dict(w=np.concatenate([d["w"] for d in ds]) / n, t=np.concatenate([d["t"] for d in ds]),
                tb=np.concatenate([d["tb"] for d in ds]), files=[os.path.basename(f) for f in files],
                seeds=[int(d["seed"]) for d in ds], sigmas=[float(d["w"].mean()) for d in ds])


def cells(spec):
    """[(edges of observable 1), (edges of observable 2) or None] for each bin / 2D cell, in the order of binned()."""

    if len(spec) == 1:
        e = spec[0][1]
        return [((e[k], e[k + 1]), None) for k in range(len(e) - 1)]
    (_, e1), (_, e2) = spec
    return [((e1[a], e1[a + 1]), (e2[b], e2[b + 1])) for a in range(len(e1) - 1) for b in range(len(e2) - 1)]


def binned(obs, w, spec):
    """sum of w and sqrt(sum of w^2) in each bin (2D: each cell, first observable slowest); outside is dropped."""

    idx, shape = [], []
    for name, edges in spec:
        idx.append(np.digitize(obs[name], edges) - 1)
        shape.append(len(edges) - 1)
    ok = np.all([(i >= 0) & (i < n) for i, n in zip(idx, shape)], axis=0)
    flat = np.ravel_multi_index([i[ok] for i in idx], shape)
    size = int(np.prod(shape))
    return (np.bincount(flat, weights=w[ok], minlength=size),
            np.sqrt(np.bincount(flat, weights=w[ok] ** 2, minlength=size)))


LABEL = {"m_tt": (r"$m_{t\bar t}$", "GeV"), "absy_tt": (r"$|y_{t\bar t}|$", ""), "y_tt": (r"$y_{t\bar t}$", ""),
         "absyboost_tt": (r"$|y_t + y_{\bar t}|\,/\,2$", ""), "HT_tt": (r"$p_T^{\,t} + p_T^{\,\bar t}$", "GeV"),
         "pT_t": (r"$p_T^{\,t}$", "GeV"), "absy_t": (r"$|y_t|$", "")}


def rng(a, b):
    return f"{a:g}–{b:g}"


def plot_title(dname, spec):
    note = ""
    if len(spec) == 2:
        (o1, e1), (o2, e2) = spec
        return (f"{LABEL[o1][0]} x {LABEL[o2][0]}: {len(e1) - 1} x {len(e2) - 1} cells" + note)
    (obs, edges), = spec
    sym, unit = LABEL[obs]
    return f"{sym}: {len(edges) - 1} bins, {edges[0]:g} to {edges[-1]:g}" + (f" {unit}" if unit else "") + note


def draw_ratio(ax, dname, spec, T, Te, title=None, ylabel=r"$T = \sigma_{\rm with} \,/\, \sigma_{\rm without}$"):
    """T = sigma_with / sigma_without per bin: one equal-width bar per bin, drawn from T = 1, with its MC error and
    the bin edges written below. 2D: the cells grouped by the first observable, groups labelled on top."""

    n = len(T)
    x = np.arange(n)
    ax.bar(x, T - 1, bottom=1, width=0.8, color=np.where(T >= 1, "#eb6834", "#2a78d6"), alpha=0.85)
    ax.errorbar(x, T, yerr=Te, fmt="none", ecolor="#222222", elinewidth=1.1, capsize=3)
    ax.axhline(1, color="#222222", lw=0.8)
    if len(spec) == 1:
        (obs, edges), = spec
        labels = [rng(a, b) for a, b in zip(edges[:-1], edges[1:])]
        sym, unit = LABEL[obs]
        ax.set_xlabel(f"bin of {sym}" + (f" [{unit}]" if unit else ""))
    else:
        (o1, e1), (o2, e2) = spec
        m = len(e2) - 1
        labels = [rng(e2[k % m], e2[k % m + 1]) for k in range(n)]
        for g in range(len(e1) - 1):
            if g:
                ax.axvline(g * m - 0.5, color="#8a8a85", lw=0.9)
            ax.text((g + 0.5) * m - 0.5, 1.0, f"{LABEL[o1][0]} {rng(e1[g], e1[g + 1])}", transform=ax.get_xaxis_transform(),
                    ha="center", va="bottom", fontsize=8.5)
        sym, unit = LABEL[o2]
        ax.set_xlabel(f"bin of {sym}" + (f" [{unit}]" if unit else "") + f", in slices of {LABEL[o1][0]}")
    ax.set_xticks(x, labels, rotation=45 if n > 8 else 0, ha="right" if n > 8 else "center", fontsize=8.5)
    ax.set_ylabel(ylabel)
    ax.ticklabel_format(axis="y", useOffset=False, style="plain")
    ax.set_xlim(-0.6, n - 0.4)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#e6e6e3", lw=0.6)
    ax.set_axisbelow(True)
    if title:
        ax.set_title(title, fontsize=10, pad=22 if len(spec) == 2 else 6)


def plot_distribution(dname, spec, T, Te, n0, n1, outdir):
    """One figure per binning: T = sigma_with / sigma_without per bin."""

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(outdir, exist_ok=True)
    fig, ax = plt.subplots(figsize=(max(7.0, 0.62 * len(T) + 3), 4.8), layout="constrained")
    draw_ratio(ax, dname, spec, T, Te)
    fig.suptitle(f"T = sigma(with toponium) / sigma(without) per bin\n{plot_title(dname, spec)}\n"
                 f"{CFG['desc']}; 13 TeV, LO, CT25NNLO; {n0 / 1e6:g}M + {n1 / 1e6:g}M events; error bars: MC statistics",
                 fontsize=10)
    fig.savefig(os.path.join(outdir, f"{dname}.png"), dpi=130)
    plt.close(fig)


def plot_all(results, n0, n1, path):
    """All binnings on one page (same content as the single figures)."""

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    nrow = (len(results) + 2) // 2
    fig, axes = plt.subplots(nrow, 2, figsize=(14, 3.9 * nrow), layout="constrained")
    for ax, (dname, (spec, T, Te)) in zip(axes.flat, results.items()):
        draw_ratio(ax, dname, spec, T, Te, title=plot_title(dname, spec))
    for ax in axes.flat[len(results):]:
        ax.axis("off")
    axes.flat[-1].text(0.0, 0.9, "Each bar: T = sigma(with toponium) / sigma(without)\nin that bin (T = 1: no change).\n"
                       "Orange: toponium adds events; blue: fewer.\nError bars: Monte Carlo statistics.\n\n"
                       + "\n".join(textwrap.wrap(CFG['desc'], 60)) + f";\n13 TeV, LO MadGraph, CT25NNLO with its alpha_s(Q);\n{n0 / 1e6:g}M events without + "
                       f"{n1 / 1e6:g}M with toponium.",
                       fontsize=11, va="top", transform=axes.flat[-1].transAxes)
    fig.suptitle(f"T-factors per bin: {CFG['desc']}", fontsize=14)
    fig.savefig(path, dpi=110)
    plt.close(fig)


def analyse(subdir=""):
    out = os.path.join(RES_RUN, subdir)
    os.makedirs(out, exist_ok=True)
    lines = [f"STAGE 6, {RUN}: {CFG['desc']}; toponium T-factors at 13 TeV, T = sigma_with / sigma_without per bin",
             f"process: {CFG['process']}" + ("; colour singlet only (CF = 2/3, Fuks (19)-(20))" if CFG["singlet_only"]
                                              else "; full colour"),
             "WITH: gg |M|^2 -> |M|^2 + (|G/G0|^2 - 1)|M_1|^2 (Fuks (16), singlet (19)-(20)) for W <= 350, p* <= 50 GeV;"
             + (" octet: gg |M_8|^2 and qqbar |M|^2 x |G8/G0|^2 (repulsive Coulomb, alpha_s = CT25(25 GeV));"
                if CFG["octet"] else ""),
             "G from our JKT solver with CT25NNLO's alpha_s(Q); MadGraph LO, CT25NNLO, m_t = 173, Gamma_t = 1.49", ""]
    log = lambda s="": (print(s), lines.append(s))
    tag = "smoke" if subdir else "s*"
    g0 = load_group(os.path.join(subdir, f"{CFG['names'][0]}_{tag}.npz"), base=EV_DIR_WITHOUT)
    g1 = load_group(os.path.join(subdir, f"{CFG['names'][1]}_{tag}.npz"))
    if g0 is None or g1 is None:
        sys.exit("need both WITHOUT and WITH events")
    for lab, g in (("without", g0), ("with", g1)):
        log(f"{lab:8s}: {len(g['w'])} events in {len(g['files'])} runs, sigma = {g['w'].sum():.6g} pb")
        for f, sd, sg in zip(g["files"], g["seeds"], g["sigmas"]):
            log(f"            {f:34s} seed {sd:5d}  sigma {sg:.6g} pb")
    log(f"total sigma_with / sigma_without = {g1['w'].sum() / g0['w'].sum():.5f}")
    log("")
    o0, o1 = observables(g0["t"], g0["tb"]), observables(g1["t"], g1["tb"])
    rows, results = [], {}
    for dname, spec in DISTRIBUTIONS.items():
        s0, s0e = binned(o0, g0["w"], spec)
        s1, s1e = binned(o1, g1["w"], spec)
        T = s1 / s0
        Te = T * np.hypot(s1e / s1, s0e / s0)
        for k, (b1, b2) in enumerate(cells(spec)):
            rows.append([dname, k, *b1, *(b2 if b2 else ("", "")), s0[k], s0e[k], s1[k], s1e[k], T[k], Te[k]])
        plot_distribution(dname, spec, T, Te, len(g0["w"]), len(g1["w"]), os.path.join(out, "plots"))
        results[dname] = (spec, T, Te)
        log(f"{dname}: first bin T = {T[0]:.5f} +- {Te[0]:.5f}; max |T - 1| = {np.max(np.abs(T - 1)):.4f}")
    path = os.path.join(out, "tfactors_13TeV.csv")
    with open(path, "w") as f:
        f.write(f"# {RUN}: {CFG['desc']} ({CFG['process']}).\n"
                "# T = sigma_with / sigma_without per bin (13 TeV, LO MadGraph, CT25NNLO with its alpha_s(Q)).\n"
                "# WITH: gg |M|^2 -> |M|^2 + (|G/G0|^2 - 1)|M_singlet|^2 for W <= 350 GeV, p* <= 50 GeV (Fuks (16), "
                "(19)-(20)); G from our JKT solver with CT25NNLO alpha_s(Q)."
                + (" Octet (gg octet, qqbar) x |G8/G0|^2, repulsive Coulomb at alpha_s = CT25(25 GeV)." if CFG["octet"]
                   else "") +
                " Parton-level tops (t = b l+ nu, tbar = b~ l- nu~). Errors: MC statistics of the two independent "
                "samples.\n"
                "# HT_ttbar = pT(t) + pT(tbar); abs_y_boost = |y_t + y_tbar|/2; abs_y_top_x_pT_top: x1 = |y_t|, x2 = pT_t.\n")
        f.write("distribution,bin,x1_lo,x1_hi,x2_lo,x2_hi,sigma_without_pb,sigma_without_err,sigma_with_pb,"
                "sigma_with_err,T,T_err\n")
        for r in rows:
            f.write(",".join(x if isinstance(x, str) else f"{x:.8g}" for x in r) + "\n")
    plot_all(results, len(g0["w"]), len(g1["w"]), os.path.join(out, "plots", "all_bins.png"))
    log(f"-> {os.path.relpath(path, HERE)}, plots/ ({len(DISTRIBUTIONS)} figures + all_bins.png)")
    if CFG["octet"] and CFG["without_from"] == "run1_pp_fullcolour_13TeV":   # run 3: shares run 1's WITHOUT sample
        octet_effect(g0, g1, o0, o1, subdir, out, log)
    if CFG.get("compare_with") and not subdir:
        compare_runs(CFG["compare_with"], path, log)
    open(os.path.join(out, "summary.txt"), "w").write("\n".join(lines) + "\n")


def read_tfactors(path):
    """{(distribution, bin): (T, T_err)} from a tfactors_13TeV.csv."""

    rows = [r for r in open(path) if not r.startswith("#")]
    head = rows[0].strip().split(",")
    out = {}
    for r in rows[1:]:
        x = dict(zip(head, r.strip().split(",")))
        out[(x["distribution"], int(x["bin"]))] = (float(x["T"]), float(x["T_err"]))
    return out


def compare_runs(other, path, log):
    """Per-bin T of this run vs another run with the same physics (e.g. run 4 vs run 3: only the lepton content
    differs, so T must agree within the MC errors). Pull = (T_this - T_other) / sqrt(err_this^2 + err_other^2)."""

    a, b = read_tfactors(path), read_tfactors(os.path.join(RES, other, "tfactors_13TeV.csv"))
    keys = [k for k in a if k in b and np.hypot(a[k][1], b[k][1]) > 0]
    pull = np.array([(a[k][0] - b[k][0]) / np.hypot(a[k][1], b[k][1]) for k in keys])
    log("")
    log(f"COMPARISON WITH {other} (same physics, different lepton content): T_this - T_other per bin")
    log(f"  {len(keys)} bins with T != 1: pull mean {pull.mean():+.2f}, rms {pull.std():.2f}, max |pull| "
        f"{np.abs(pull).max():.2f}; |pull| > 2 in {np.sum(np.abs(pull) > 2)} bins (expect ~{0.046 * len(keys):.1f})")
    for k in [k for k in keys if k[1] == 0]:
        log(f"  {k[0]:22s} bin 0: T = {a[k][0]:.5f} +- {a[k][1]:.5f}   {other.split('_')[0]}: {b[k][0]:.5f} +- {b[k][1]:.5f}")


def octet_effect(g0, g1, o0, o1, subdir, out, log):
    """Run 3 vs run 1 (same WITHOUT sample): Delta T = T(singlet + octet) - T(singlet only) per bin, the octet change as
    a fraction of the toponium effect, and an alpha_s band (0.12, 0.18) estimated by rescaling with |G8/G0|^2 on the
    WITHOUT events in the window, normalised to the measured octet change (labelled an estimate)."""

    r1 = RUNS["run1_pp_fullcolour_13TeV"]
    tag = "smoke" if subdir else "s*"
    ga = load_group(os.path.join(subdir, f"{r1['names'][1]}_{tag}.npz"),
                    base=os.path.join(OUTPUT, "run1_pp_fullcolour_13TeV", "events"))
    if ga is None:
        log("no run-1 WITH events for the comparison")
        return
    oa = observables(ga["t"], ga["tb"])
    s0tot, s1tot, satot = g0["w"].sum(), g1["w"].sum(), ga["w"].sum()
    # error on the totals: scatter of the per-run MadGraph cross sections (the total is not a count of events)
    serr = lambda g: np.std(g["sigmas"], ddof=1) / np.sqrt(len(g["sigmas"])) if len(g["sigmas"]) > 1 else float("nan")
    d_oct, d_err = s1tot - satot, np.hypot(serr(g1), serr(ga))
    log("")
    log("COLOUR-OCTET EFFECT (run 3 - run 1, same WITHOUT sample)")
    log(f"  sigma_with: singlet only (run 1) {satot:.6g} pb, singlet + octet (run 3) {s1tot:.6g} pb")
    log(f"  toponium excess (singlet only): {satot - s0tot:.5f} pb; octet change: {d_oct:+.5f} +- {d_err:.5f} pb "
        f"= {d_oct / (satot - s0tot) * 100:+.1f}% of the toponium excess; total T {satot / s0tot:.5f} -> {s1tot / s0tot:.5f}")
    # alpha_s band from the WITHOUT events in the window
    E, ps = E_pstar(g0["t"], g0["tb"])
    win = (E + 2 * MT <= 350) & (ps <= 50) & (E >= E0)
    a8 = float(ct25_alphas(OCTET_Q ** 2))
    shape = {}                       # w (|G8/G0|^2 - 1) per event, nonzero only in the window (evaluated there only)
    for a in (0.12, a8, 0.18):
        v = np.zeros(len(E))
        idx = np.flatnonzero(win)
        for c in range(0, len(idx), 20000):                     # chunks: the exact G8 uses a 160-point quadrature
            k = idx[c:c + 20000]
            v[k] = np.abs(octet_ratio(E[k], ps[k], a)) ** 2 - 1
        shape[a] = v * g0["w"]
    f8 = d_oct / shape[a8].sum()
    log(f"  effective octet fraction of the window rate (measured / [sum w (|G8/G0|^2 - 1)]): {f8:.3f}")
    for a in (0.12, 0.18):
        log(f"  estimate at alpha_s = {a}: octet change {f8 * shape[a].sum():+.5f} pb = "
            f"{f8 * shape[a].sum() / (satot - s0tot) * 100:+.1f}% of the toponium excess")
    rows, res = [], {}
    for dname, spec in DISTRIBUTIONS.items():
        s0, s0e = binned(o0, g0["w"], spec)
        s1, s1e = binned(o1, g1["w"], spec)
        sa, sae = binned(oa, ga["w"], spec)
        Ta, T1 = sa / s0, s1 / s0
        Tae, T1e = Ta * np.hypot(sae / sa, s0e / s0), T1 * np.hypot(s1e / s1, s0e / s0)
        dT, dTe = (s1 - sa) / s0, np.hypot(s1e, sae) / s0
        band = {a: f8 * binned(o0, shape[a], spec)[0] / s0 for a in (0.12, 0.18)}
        for k, (b1, b2) in enumerate(cells(spec)):
            frac = dT[k] / (Ta[k] - 1) * 100 if abs(Ta[k] - 1) > 3 * Tae[k] else float("nan")
            rows.append([dname, k, *b1, *(b2 if b2 else ("", "")), Ta[k], Tae[k], T1[k], T1e[k], dT[k], dTe[k], frac,
                         band[0.12][k], band[0.18][k]])
        res[dname] = (spec, Ta, Tae, T1, T1e)
        log(f"  {dname}: first bin T {Ta[0]:.5f} -> {T1[0]:.5f}, Delta T = {dT[0]:+.5f} +- {dTe[0]:.5f}")
    path = os.path.join(out, "octet_effect.csv")
    with open(path, "w") as f:
        f.write("# Colour-octet threshold effect: T_singlet = run 1 (gg singlet x |G1/G0|^2); T_singlet_octet = run 3 "
                "(and gg octet, qqbar x |G8/G0|^2, repulsive Coulomb, alpha_s = CT25(25 GeV)); same WITHOUT sample.\n"
                "# dT = T_singlet_octet - T_singlet; dT_pct_of_effect = dT/(T_singlet - 1) in % (nan where T_singlet - 1 "
                "< 3 sigma); dT_est_a0.12/0.18: alpha_s band, ESTIMATE from rescaling (see summary).\n")
        f.write("distribution,bin,x1_lo,x1_hi,x2_lo,x2_hi,T_singlet,T_singlet_err,T_singlet_octet,T_singlet_octet_err,"
                "dT,dT_err,dT_pct_of_effect,dT_est_a0.12,dT_est_a0.18\n")
        for r in rows:
            f.write(",".join(x if isinstance(x, str) else f"{x:.8g}" for x in r) + "\n")
    plot_compare(res, len(g0["w"]), len(ga["w"]), len(g1["w"]), os.path.join(out, "plots"))
    log(f"  -> {os.path.relpath(path, HERE)}, plots/compare_*.png")


def plot_compare(res, n0, na, n1, outdir):
    """Per binning: T with the singlet only (run 1) and with singlet + octet (run 3), paired bars from T = 1."""

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(outdir, exist_ok=True)

    def draw(ax, dname, spec, Ta, Tae, T1, T1e, title=None):
        n = len(Ta)
        x = np.arange(n)
        for dx, T, Te, c, lab in ((-0.2, Ta, Tae, "#eb6834", "singlet only (run 1)"),
                                  (0.2, T1, T1e, "#7a4fc9", "singlet + octet (run 3)")):
            ax.bar(x + dx, T - 1, bottom=1, width=0.38, color=c, alpha=0.85, label=lab)
            ax.errorbar(x + dx, T, yerr=Te, fmt="none", ecolor="#222222", elinewidth=1, capsize=2)
        ax.axhline(1, color="#222222", lw=0.8)
        if len(spec) == 1:
            (obs, edges), = spec
            labels = [rng(a, b) for a, b in zip(edges[:-1], edges[1:])]
            sym, unit = LABEL[obs]
            ax.set_xlabel(f"bin of {sym}" + (f" [{unit}]" if unit else ""))
        else:
            (o1, e1), (o2, e2) = spec
            m = len(e2) - 1
            labels = [rng(e2[k % m], e2[k % m + 1]) for k in range(n)]
            for g in range(len(e1) - 1):
                if g:
                    ax.axvline(g * m - 0.5, color="#8a8a85", lw=0.9)
                ax.text((g + 0.5) * m - 0.5, 1.0, f"{LABEL[o1][0]} {rng(e1[g], e1[g + 1])}",
                        transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=8.5)
            sym, unit = LABEL[o2]
            ax.set_xlabel(f"bin of {sym}" + (f" [{unit}]" if unit else "") + f", in slices of {LABEL[o1][0]}")
        ax.set_xticks(x, labels, rotation=45 if n > 8 else 0, ha="right" if n > 8 else "center", fontsize=8.5)
        ax.set_ylabel(r"$T = \sigma_{\rm with} \,/\, \sigma_{\rm without}$")
        ax.ticklabel_format(axis="y", useOffset=False, style="plain")
        ax.set_xlim(-0.6, n - 0.4)
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", color="#e6e6e3", lw=0.6)
        ax.set_axisbelow(True)
        ax.legend(frameon=False, fontsize=8.5)
        if title:
            ax.set_title(title, fontsize=10, pad=22 if len(spec) == 2 else 6)

    note = (f"pp -> tt, 13 TeV, LO, CT25NNLO; {n0 / 1e6:g}M without, {na / 1e6:g}M with (singlet), "
            f"{n1 / 1e6:g}M with (singlet + octet); octet: repulsive Coulomb, alpha_s(25 GeV)")
    for dname, (spec, Ta, Tae, T1, T1e) in res.items():
        fig, ax = plt.subplots(figsize=(max(7.5, 0.75 * len(Ta) + 3), 4.8), layout="constrained")
        draw(ax, dname, spec, Ta, Tae, T1, T1e)
        fig.suptitle(f"T per bin: toponium in the colour singlet only vs singlet + octet\n{plot_title(dname, spec)}\n{note}",
                     fontsize=10)
        fig.savefig(os.path.join(outdir, f"compare_{dname}.png"), dpi=130)
        plt.close(fig)
    nrow = (len(res) + 1) // 2
    fig, axes = plt.subplots(nrow, 2, figsize=(15, 3.9 * nrow), layout="constrained")
    for ax, (dname, (spec, Ta, Tae, T1, T1e)) in zip(axes.flat, res.items()):
        draw(ax, dname, spec, Ta, Tae, T1, T1e, title=plot_title(dname, spec))
    for ax in axes.flat[len(res):]:
        ax.axis("off")
    fig.suptitle("T per bin: toponium in the colour singlet only (run 1) vs singlet + octet (run 3)\n" + note, fontsize=12)
    fig.savefig(os.path.join(outdir, "compare_all.png"), dpi=110)
    plt.close(fig)

# =====================================================================================
# 5. runs 5 and 6: K = 1 + Delta/SM_full in the student's CSV layout
# =====================================================================================
KFILES = [   # (selection p* cut, octet in column T, octets of the extra columns): his two files
    (100.0, "freeze25", ("freeze10", "freeze50", "coulomb_as25")),
    (50.0, "coulomb_as25", ("freeze25", "freeze10", "freeze50")),
]


def kfiles():
    """The K-factor files of the selected run: KFILES (runs 5-6), or the run's own list (run 7: one octet)."""

    return CFG.get("kfiles", KFILES)


def kstem(pcut, head):
    return f"kfactors_13TeV_W{wlo():g}-{whi():g}_pstar{pcut:.0f}_octet-{head}"


def kfile(pcut, head):
    """File name of a K-factor CSV: runs 5-6 end in '_ours' (vs the student's files of the same name); run 7 does not."""

    return kstem(pcut, head) + ("" if "mt" in CFG else "_ours") + ".csv"


def read_his(path):
    """His CSV: {(distribution, bin): row dict} and {seed: sigma [pb]} from the RUNS lines of the header."""

    rows, sig = {}, {}
    lines = open(path).read().splitlines()
    for l in lines:
        m = re.search(r"sigma ([0-9.eE+-]+) pb\s+iseed (\d+)", l)
        if l.startswith("#") and m:
            sig[int(m.group(2))] = float(m.group(1))
    data = [l for l in lines if l and not l.startswith("#")]
    head = data[0].split(",")
    for l in data[1:]:
        x = dict(zip(head, l.split(",")))
        rows[(x["distribution"], int(x["bin"]))] = x
    return rows, sig


def _draw_k(ax, dname, spec, rows, title=None):
    """K = 1 + Delta/SM_full per bin (column T) as draw_ratio bars, the other octet versions (columns T_<octet>) as ticks."""

    n = len(cells(spec))
    K = np.array([float(rows[(dname, k)]["T"]) for k in range(n)])
    Ke = np.array([float(rows[(dname, k)]["T_err"]) for k in range(n)])
    draw_ratio(ax, dname, spec, K, Ke, title=title, ylabel=r"$K = 1 + \Delta / \mathrm{SM_{full}}$")
    others = [c for c in rows[(dname, 0)] if c.startswith("T_") and c != "T_err"]
    on = K != 1   # bins with no toponium effect have K = 1 exactly for every octet version
    for i, c in enumerate(others):
        ax.scatter(np.arange(n)[on], np.array([float(rows[(dname, k)][c]) for k in range(n)])[on], marker="_", s=110,
                   lw=1.6, color="#3d3d3a", zorder=4, label="other octet versions" if i == 0 else None)
    return K


def plot_kfactors(path, outdir):
    """Run 5-7 figures from one K-factor CSV (columns as written by write_kfactors), one per binning of DISTRIBUTIONS:
      top     d sigma / dX per bin for SM_full (LO, column sigma_without_pb) and SM_full + Delta (column sigma_with_pb), each
              divided by its bin width (2D binning: by the width in the second observable, per slice of the first)
      bottom  K = 1 + Delta/SM_full (column T) with its MC error, and the other octet versions (columns T_<octet>)
    plus all_bins.png (the K panels of all binnings on one page). Bins are equal-width slots labelled by their edges."""

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows, _ = read_his(path)
    m = re.search(r"pstar(\d+)_octet-(\w+?)(?:_ours)?\.csv$", os.path.basename(path))
    pcut, head = m.group(1), m.group(2)
    info = (f"m_t = {MT:g} GeV, window {wlo():g} <= W <= {whi():g} GeV, p* < {pcut} GeV; octet '{head}'; pp -> tt -> b l nu b l nu (l = e, mu), 13 TeV, LO, "
            f"CT25NNLO\nSM_full: 10M LO events, full phase space; Delta: window runs SM_win, topo_win for gg and qqbar, 1M events each; "
            f"error bars: MC statistics")
    os.makedirs(outdir, exist_ok=True)
    for dname, spec in DISTRIBUTIONS.items():
        n = len(cells(spec))
        get = lambda col: np.array([float(rows[(dname, k)][col]) for k in range(n)])
        width = np.array([(b2 or b1)[1] - (b2 or b1)[0] for b1, b2 in cells(spec)])
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(max(7.0, 0.62 * n + 3), 7.6), sharex=True, layout="constrained",
                                       gridspec_kw=dict(height_ratios=[1.15, 1]))
        groups = [range(n)] if len(spec) == 1 else [range(g, g + len(spec[1][1]) - 1)
                                                     for g in range(0, n, len(spec[1][1]) - 1)]
        for col, ecol, color, ls, lw, lab in (
                ("sigma_with_pb", "sigma_with_err", "#eb6834", "-", 2.4, "with threshold effects (SM_full + Delta)"),
                ("sigma_without_pb", "sigma_without_err", "#3d3d3a", "--", 1.3, "LO (SM_full)")):
            y, ye = get(col) / width, get(ecol) / width
            for i, g in enumerate(groups):
                g = np.array(g)
                ax1.stairs(y[g], np.r_[g - 0.5, g[-1] + 0.5], baseline=None, color=color, ls=ls, lw=lw,
                           label=lab if i == 0 else None)
            ax1.errorbar(np.arange(n), y, yerr=ye, fmt="none", ecolor=color, elinewidth=1.0, capsize=2)
        for g in groups[1:]:
            ax1.axvline(g[0] - 0.5, color="#8a8a85", lw=0.9)
        ax1.set_yscale("log")
        sym, unit = LABEL[spec[-1][0]]
        ax1.set_ylabel(r"$d\sigma/d$" + sym + (f" [pb/{unit}]" if unit else " [pb]"))
        ax1.legend(frameon=False, fontsize=9)
        ax1.spines[["top", "right"]].set_visible(False)
        ax1.grid(axis="y", color="#e6e6e3", lw=0.6)
        ax1.set_axisbelow(True)
        _draw_k(ax2, dname, spec, rows)
        ax2.legend(frameon=False, fontsize=9, loc="upper right")
        fig.suptitle(f"Toponium K-factors per bin: {plot_title(dname, spec)}\n{info}", fontsize=9.5)
        fig.savefig(os.path.join(outdir, f"{dname}.png"), dpi=130)
        plt.close(fig)

    nrow = (len(DISTRIBUTIONS) + 2) // 2
    fig, axes = plt.subplots(nrow, 2, figsize=(14, 3.9 * nrow), layout="constrained")
    for ax, (dname, spec) in zip(axes.flat, DISTRIBUTIONS.items()):
        _draw_k(ax, dname, spec, rows, title=plot_title(dname, spec))
    for ax in axes.flat[len(DISTRIBUTIONS):]:
        ax.axis("off")
    axes.flat[-1].text(0.0, 0.9, "Each bar: K = 1 + Delta/SM_full in that bin (K = 1: no change),\nwith its MC error. "
                       "Dark ticks: the same K with the other\noctet potentials.\n\n" + info.replace("; ", ";\n"),
                       fontsize=10.5, va="top", transform=axes.flat[-1].transAxes)
    fig.suptitle(f"Toponium K-factors per bin ({RUN}): p* < {pcut} GeV, octet '{head}'", fontsize=14)
    fig.savefig(os.path.join(outdir, "all_bins.png"), dpi=110)
    plt.close(fig)
    return outdir


PLAIN = {"m_tt": "m_ttbar_GeV", "pT_t": "pT_top_GeV", "absy_t": "abs_y_top", "y_tt": "y_ttbar", "absy_tt": "abs_y_ttbar",
         "absyboost_tt": "abs_y_boost", "HT_tt": "HT_ttbar_GeV"}


def export(subdir=""):
    """Easy-to-read K-factors, one CSV per binning of DISTRIBUTIONS (+ all_kfactors.txt), from the two kfactors CSVs:
      K          = 1 + Delta/SM_full, p* < 100 GeV, octet freeze25 (column T of the p* < 100 file)
      stat_err   = its MC error (column T_err)
      model_unc  = sqrt(dK_pstar^2 + dK_octet^2): dK_pstar = K(p* < 100) - K(p* < 50), same octet (column T_freeze25
                   of the p* < 50 file); dK_octet = largest |K_other octet - K| (columns T_<octet> of the p* < 100 file)
      K_norm     = K / K_range, K_range = sum(sigma_with)/sum(sigma_without) over the binning's bins: for distributions
                   normalised to their own integral."""

    src = os.path.join(RES_RUN, subdir)
    (p1, h1, _), (p2, h2, _) = kfiles()
    r100, _ = read_his(os.path.join(src, kfile(p1, h1)))
    r50, _ = read_his(os.path.join(src, kfile(p2, h2)))
    col50 = f"T_{h1}" if h1 != h2 else "T"          # K at p* < 50 with the octet of K (run 7: the same file column T)
    out = os.path.join(src, "kfactors_final")
    os.makedirs(out, exist_ok=True)
    nocts = len(kfiles()[0][2]) + 1
    head = (f"# Toponium K-factors (colour singlet + octet), p p -> t t~, 13 TeV, LO MadGraph5_aMC@NLO 2.9.27, CT25NNLO "
            f"(run {RUN}; m_t = {MT:g} GeV, Gamma_t = {WT:g} GeV)\n# K = sigma(with toponium) / sigma(SM) per bin; multiply the SM prediction of the bin by K.\n"
            f"# K_norm = K / K over the whole binning: use it for distributions normalised to their own integral.\n"
            f"# stat_err: MC statistics. model_unc: p* cut (100 vs 50 GeV)"
            + (f" and octet potential ({nocts} versions), in quadrature.\n" if nocts > 1 else ".\n") +
            f"# Bins: lower edge included, upper excluded; tops = parton level (b + l + nu).\n")
    allt = [head.rstrip("\n")]
    for dname, spec in DISTRIBUTIONS.items():
        n = len(cells(spec))
        rows = [r100[(dname, k)] for k in range(n)]
        sw = sum(float(r["sigma_with_pb"]) for r in rows) / sum(float(r["sigma_without_pb"]) for r in rows)
        cols = [f"{PLAIN[o]}_lo,{PLAIN[o]}_hi" for o, _ in spec]
        lines = [f"# {dname}",
                 ",".join(cols) + ",K,stat_err,model_unc,K_norm"]
        for k, ((a1, b1), c2) in enumerate(cells(spec)):
            r = rows[k]
            K, dK = float(r["T"]), float(r["T_err"])
            dp = K - float(r50[(dname, k)][col50])
            do = max([abs(float(r[c]) - K) for c in r if c.startswith("T_") and c != "T_err"], default=0.0)
            edges = [a1, b1] + (list(c2) if c2 else [])
            lines.append(",".join(f"{e:g}" for e in edges) + f",{K:.4f},{dK:.4f},{np.hypot(dp, do):.4f},{K / sw:.4f}")
        open(os.path.join(out, f"{dname}.csv"), "w").write(head + "\n".join(lines) + "\n")
        allt += [""] + lines
    open(os.path.join(out, "all_kfactors.txt"), "w").write("\n".join(allt) + "\n")
    return out


def reflect_z(g):
    """The beam-swapped mirror of a q q~ sample (q~ from beam 1, q from beam 2), which a 'generate q q~ > ...' output
    lacks (MIRRORPROCS = .FALSE.): at LO in symmetric p p collisions it is the same sample with p_z -> -p_z. Returns
    the sample plus its mirror (sigma doubled; each event once as generated and once reflected)."""

    flip = lambda p: p * np.array([1.0, 1.0, 1.0, -1.0])
    return dict(g, w=np.concatenate([g["w"], g["w"]]), t=np.concatenate([g["t"], flip(g["t"])]),
                tb=np.concatenate([g["tb"], flip(g["tb"])]))


def kfactors(subdir="", his=(), reflect_qq=False, tag=""):
    """Run 5, his definitions: in each bin
        sigma_bin(run) = sum of w over the run's events in the bin and in the file's window (326 <= W <= 350 GeV,
                         p* < cut), w = sigma_run / N_run (unweighted events, event_norm = average),
        Delta = (topo_win_gg - SM_win_gg) + (topo_win_qq - SM_win_qq),   K = T = 1 + Delta/C,   sigma_with = C + Delta,
        C = SM_full, the full phase-space LO sample (run 5: run 4's WITHOUT; run 6: its own), no window; errors MC statistics: delta_err from the four runs in quadrature,
        T_err = sqrt((delta_err/C)^2 + (Delta C_err/C^2)^2), sigma_with_err = hypot(C_err, delta_err).
    his: his CSV files (data only) to compare with, matched to ours by file name.
    reflect_qq: preview of the q q~ mirror fix: every q q~ window sample + its z-reflection (reflect_z); output in tag/."""

    out = os.path.join(RES_RUN, subdir, tag)
    os.makedirs(out, exist_ok=True)
    lines = [f"STAGE 6, {RUN}: {CFG['desc']}", ""]
    log = lambda s="": (print(s), lines.append(s))
    own = CFG["without_from"] is None                 # run 6: own SM_full (and the check topo_full) among runs()
    C = load_group(os.path.join(subdir, f"{CFG['names'][0]}_s*.npz") if not subdir else
                   os.path.join(subdir, f"{CFG['names'][0]}_smoke.npz") if own else f"{CFG['names'][0]}_s*.npz",
                   base=EV_DIR_WITHOUT)
    if C is None:
        sys.exit(f"missing the SM_full sample in {EV_DIR_WITHOUT}: run 'run {RUN}' first")
    oC = observables(C["t"], C["tb"])
    log(f"SM_full: {os.path.relpath(EV_DIR_WITHOUT, MG_ROOT)}, {len(C['w'])} events in {len(C['files'])} runs, "
        f"sigma = {C['w'].sum():.6g} pb")
    win, F = {}, None                 # (channel, octet) -> events, observables, W, p*;  F: run 6's check sample
    for name, d, seed, nev, octet in runs(smoke=bool(subdir)):
        if d == SM_DIR:
            continue
        g = load_group(os.path.join(subdir, f"{name}.npz"))
        if g is None:
            sys.exit(f"missing {name}.npz: run 'run {RUN}' first")
        E, ps = E_pstar(g["t"], g["tb"])
        if d == TOPO_DIR and CFG.get("check"):
            F = dict(g=g, o=observables(g["t"], g["tb"]), W=E + 2 * MT, ps=ps, seed=seed, name=name, octet=octet)
            log(f"{name:28s} seed {seed}  sigma {g['w'].sum():.6g} pb  {len(g['w'])} events (check topo_full, full phase space)")
            continue
        ch = os.path.basename(d).split("_")[-1]
        if reflect_qq and ch == "qq":
            g = reflect_z(g)
            E, ps = E_pstar(g["t"], g["tb"])
        win[(ch, octet)] = dict(g=g, o=observables(g["t"], g["tb"]), W=E + 2 * MT, ps=ps, seed=seed, name=name)
        log(f"{name:28s} seed {seed}  sigma {g['w'].sum():.6g} pb  {len(g['w'])} events; max W {(E + 2 * MT).max():.2f}"
            f", max p* {ps.max():.2f} GeV")
    his_by_name = {os.path.basename(h).replace(".csv", ""): h for h in his}
    for (ch, o), v in win.items():               # q q~ + mirror q~ q: y_tt symmetric (without the mirror <y_tt> ~ +0.2)
        if ch == "qq":
            y, w_ = v["o"]["y_tt"], v["g"]["w"]
            m, e = (w_ * y).sum() / w_.sum(), np.sqrt((w_ ** 2 * (y - (w_ * y).sum() / w_.sum()) ** 2).sum()) / w_.sum()
            log(f"{v['name']:28s} mirror check: <y_tt> = {m:+.4f} +- {e:.4f}" + ("" if abs(m) < 5 * e + 0.01 else
                "   ** ASYMMETRIC: q~ q mirror missing? **"))

    for pcut, head, extra in kfiles():
        stem = kstem(pcut, head)
        sel = {k: (v["W"] >= wlo()) & (v["W"] <= whi()) & (v["ps"] < pcut) for k, v in win.items()}
        rows, delta = [], {}
        for dname, spec in DISTRIBUTIONS.items():
            c, ce = binned(oC, C["w"], spec)
            K, D = {}, {}
            for o in (head,) + extra:
                parts = [binned(win[(ch, oo)]["o"], win[(ch, oo)]["g"]["w"] * sel[(ch, oo)], spec)
                         for ch in ("gg", "qq") for oo in (o, None)]
                d_ = parts[0][0] - parts[1][0] + parts[2][0] - parts[3][0]
                de = np.sqrt(sum(p[1] ** 2 for p in parts))
                K[o], D[o] = 1 + d_ / c, (d_, de)
            d_, de = D[head]
            Te = np.sqrt((de / c) ** 2 + (d_ * ce / c ** 2) ** 2)
            delta[dname] = d_.sum()
            for k, (b1, b2) in enumerate(cells(spec)):
                rows.append([dname, k, *b1, *(b2 if b2 else ("", "")), c[k], ce[k], c[k] + d_[k], np.hypot(ce[k], de[k]),
                             K[head][k], Te[k], d_[k], de[k], *[K[o][k] for o in extra]])
        path = os.path.join(out, kfile(pcut, head))
        write_kfactors(path, pcut, head, extra, C, win, sel, rows)
        plot_kfactors(path, os.path.join(out, "plots", f"pstar{pcut:.0f}"))
        log("")
        log(f"{os.path.relpath(path, HERE)}")
        log(f"  window {wlo():g} <= W <= {whi():g}, p* < {pcut:.0f}: total Delta = {delta['m_ttbar_250-3500']:.5f} pb "
            f"(K_total = {1 + delta['m_ttbar_250-3500'] / C['w'].sum():.5f}); Delta in m_ttbar_345-1600: "
            f"{delta['m_ttbar_345-1600']:.5f} pb")
        if stem in his_by_name:
            compare_his(his_by_name[stem], rows, extra, win, log)
        if F is not None:
            window_check(pcut, C, oC, win, F, log)
    open(os.path.join(out, "summary.txt"), "w").write("\n".join(lines) + "\n")


def window_check(pcut, C, oC, win, F, log):
    """Run 6: the window design against full phase-space samples, in this file's window (326 <= W <= 350 GeV,
    p* < pcut), with sigma_X = sum of w over the events of X in the window:
      (a) SM      sigma(SM_full) vs sigma(SM_win_gg) + sigma(SM_win_qq): the generation window (genps tau, dummy_cuts) loses nothing
      (b) topo    sigma(topo_full) vs sigma(topo_win_gg) + sigma(topo_win_qq), same octet: the same with the matrix-element patch
      (c) Delta   Delta_full = topo_full - SM_full per bin (both selected in the window) vs Delta of the K-factor file, and
                   K_full = 1 + Delta_full/C in the threshold bins (independent samples: errors in quadrature)."""

    if "W" not in C:
        E, ps = E_pstar(C["t"], C["tb"])
        C["W"], C["ps"] = E + 2 * MT, ps
    inw = lambda v: (v["W"] >= wlo()) & (v["W"] <= whi()) & (v["ps"] < pcut)
    sig = lambda w, m: (w[m].sum(), np.sqrt((w[m] ** 2).sum()))
    pull = lambda a, b: f"ratio {a[0] / b[0]:.4f}, pull {(a[0] - b[0]) / np.hypot(a[1], b[1]):+.2f}"
    o = F["octet"]
    log(f"  window check (p* < {pcut:.0f} GeV; topo_full: {F['name']}, octet '{o}'):")
    sC, sF = sig(C["w"], inw(C)), sig(F["g"]["w"], inw(F))
    sS = [sig(win[(ch, None)]["g"]["w"], inw(win[(ch, None)])) for ch in ("gg", "qq")]
    sT = [sig(win[(ch, o)]["g"]["w"], inw(win[(ch, o)])) for ch in ("gg", "qq")]
    add = lambda xs: (sum(x[0] for x in xs), np.sqrt(sum(x[1] ** 2 for x in xs)))
    sS, sT = add(sS), add(sT)
    log(f"    (a) SM_full in window {sC[0]:.5f} +- {sC[1]:.5f} pb  vs  SM_win_gg + SM_win_qq {sS[0]:.5f} +- {sS[1]:.5f} pb: {pull(sC, sS)}")
    log(f"    (b) topo_full in window {sF[0]:.5f} +- {sF[1]:.5f} pb  vs  topo_win_gg + topo_win_qq {sT[0]:.5f} +- {sT[1]:.5f} pb: {pull(sF, sT)}")
    dF, dW = (sF[0] - sC[0], np.hypot(sF[1], sC[1])), (sT[0] - sS[0], np.hypot(sT[1], sS[1]))
    log(f"    (c) total Delta: full samples {dF[0]:.5f} +- {dF[1]:.5f} pb  vs  window runs {dW[0]:.5f} +- {dW[1]:.5f} pb: "
        f"{pull(dF, dW)}")
    oF = F["o"]
    for dname, k in (("m_ttbar_250-3500", 0), ("m_ttbar_345-1600", 0), ("HT_ttbar_0-2000", 0), ("pT_top_0-500", 0)):
        spec = DISTRIBUTIONS[dname]
        cb, cbe = binned(oC, C["w"] * inw(C), spec)
        fb, fbe = binned(oF, F["g"]["w"] * inw(F), spec)
        call, _ = binned(oC, C["w"], spec)
        parts = [binned(win[(ch, oo)]["o"], win[(ch, oo)]["g"]["w"] * inw(win[(ch, oo)]), spec)
                 for ch in ("gg", "qq") for oo in (o, None)]
        dw = parts[0][0][k] - parts[1][0][k] + parts[2][0][k] - parts[3][0][k]
        dwe = np.sqrt(sum(p_[1][k] ** 2 for p_ in parts))
        Kw, Kwe = 1 + dw / call[k], dwe / call[k]                 # window K for F's octet (C errors common: left out)
        Kf, Kfe = 1 + (fb[k] - cb[k]) / call[k], np.hypot(fbe[k], cbe[k]) / call[k]
        log(f"        {dname:18s} bin {k}: K_full {Kf:.5f} +- {Kfe:.5f}  vs  K_window {Kw:.5f} +- {Kwe:.5f} "
            f"(pull {(Kf - Kw) / np.hypot(Kfe, Kwe):+.2f})")


def compare_his(path, rows, extra, win, log):
    """Ours vs his, same layout: run cross sections by seed (his seeds; the SM_win runs are unmodified MadGraph), total
    Delta, per-bin K: pull (K_ours - K_his)/hypot(errors) for column T, largest |difference| for the other columns."""

    h, hsig = read_his(path)
    log(f"  vs {os.path.basename(path)}:")
    for v in win.values():
        if v["seed"] in hsig:
            s = v["g"]["w"].sum()
            log(f"    {v['name']:28s} sigma ours {s:.6g} pb, his {hsig[v['seed']]:.6g} pb ({(s / hsig[v['seed']] - 1) * 100:+.2f}%)")
    dh = sum(float(x["delta_pb"]) for (dn, _), x in h.items() if dn == "m_ttbar_250-3500")
    do = sum(r[12] for r in rows if r[0] == "m_ttbar_250-3500")
    log(f"    total Delta: ours {do:.5f} pb, his {dh:.5f} pb ({(do / dh - 1) * 100:+.1f}%)")
    pulls, worst = [], []
    for r in rows:
        x = h.get((r[0], r[1]))
        if x is None or float(x["T_err"]) == 0 and r[11] == 0:
            continue
        pulls.append((r[10] - float(x["T"])) / np.hypot(r[11], float(x["T_err"])))
        worst.append((abs(pulls[-1]), r[0], r[1], r[10], r[11], float(x["T"]), float(x["T_err"])))
    pulls = np.array(pulls)
    log(f"    column T, {len(pulls)} bins with a toponium effect: pull mean {pulls.mean():+.2f}, rms "
        f"{np.sqrt((pulls ** 2).mean()):.2f}, max |pull| {np.abs(pulls).max():.2f}")
    for _, dn, b, k, ke, kh, khe in sorted(worst, reverse=True)[:6]:
        log(f"      {dn:22s} bin {b:2d}: ours {k:.5f} +- {ke:.5f}, his {kh:.5f} +- {khe:.5f}")
    for i, o in enumerate(extra):
        d = [abs(r[14 + i] - float(h[(r[0], r[1])][f"T_{o}"])) for r in rows if (r[0], r[1]) in h]
        log(f"    column T_{o}: max |K_ours - K_his| = {max(d):.5f}")


def write_kfactors(path, pcut, head, extra, C, win, sel, rows):
    """The CSV in his layout (same columns, distributions and rows), header describing our set-up."""

    a8 = float(ct25_alphas(OCTET_Q ** 2))
    sing = ("full JKT potential (13)-(15): JKT's two-loop alpha_s (14) at alpha_s(mZ) = 0.12, NLO bracket, IR form (15) "
            "below 5 GeV" if CFG.get("singlet") == "jkt" else "JKT potential (13)-(15) with CT25 alpha_s")
    desc = {o: (f"V~8 = +(1/6) 4 pi a_s(Q)/q^2 [1 + a1 a_s(Q)/(4 pi)], Q = max(q, {Q0:g} GeV), for ALL q "
                f"(JKTPotential.octet_frozen, LSSolver n = 120)" if Q0 else
                f"V~8 = +(1/6) 4 pi a_s/q^2, a_s = a_s(25 GeV) = {a8:.6f} FIXED, tree level (CoulombExact, exact)")
            for o, Q0 in OCTETS.items()}
    run_lines = []
    for (ch, o), v in win.items():
        what = (f"{'g g ' if ch == 'gg' else 'q q~'} SM |M|^2, window only" if o is None else
                f"{'g g ' if ch == 'gg' else 'q q~'} with G, octet '{o}'")
        run_lines.append(f"#   {v['name']:28s} {what:38s} sigma {v['g']['w'].sum():.6g} pb  iseed {v['seed']}  "
                         f"events {int(sel[(ch, o)].sum())}")
    own = CFG["without_from"] is None
    origin = ("#   our Green's functions (toponium.py: JKTPotential, LSSolver, CoulombExact) and our MATRIX1 patch; the student's design\n"
              "#   and definitions, with an independent SM_full sample and new seeds." if own else
              "#   to compare with the student's file of the same name without '_ours': same design and definitions as his;\n"
              "#   our Green's functions (toponium.py: JKTPotential, LSSolver, CoulombExact) and our MATRIX1 patch.")
    csrc = "this run's own LO sample" if own else "run 4's WITHOUT sample"
    seeds = "new seeds (listed below)" if own else "his seeds"
    qqmir = "; add process q~ q > (same), the beam-swapped mirror" if "add process q~ q" in CFG["process_qq"] else \
            " (no mirror q~ q: half of qqbar)"
    hdr = f"""# TOPONIUM K-FACTORS{"" if "mt" in CFG else " (OURS)"}  p p -> t t~ -> b l+ nu b~ l- nu~ (l = e, mu), 13 TeV, LO, MadGraph5_aMC@NLO 2.9.27, PDF CT25NNLO with its own alpha_s(Q).
# Produced by toponium_study/6_tfactors.py (run {RUN}):
{origin}
#
# DEFINITIONS (as his)
#   W = m(t t~) of the six decay products, t = b + l+ + nu, t~ = b~ + l- + nu~ (parton level, off-shell tops); p* = |p_t| in the t t~ rest frame;
#   E = W - 2 m_t, m_t = {MT:g} GeV, Gamma_t = {WT:g} GeV.
#   WINDOW (this file): {wlo():g} GeV <= W <= {whi():g} GeV and p* < {pcut:.0f} GeV, applied as an event SELECTION on the window runs.
#   GENERATION WINDOW (all SM_win and topo_win runs): W <= {whi():g} GeV (TAUMAX = {whi():g}^2/S in genps.f) and p* <= 100 GeV (dummy_cuts), no lower W edge.
#   C = SM_full = sigma_bin of the p p plain-SM sample (no window, no re-weighting) = column sigma_without: {csrc},
#       {len(C['w'])} events in {len(C['files'])} runs (seeds {C['seeds'][0]}..{C['seeds'][-1]}), sigma {C['w'].sum():.6g} pb.
#   SM_win_gg, topo_win_gg, SM_win_qq, topo_win_qq: the window runs below, restricted to this file's window.
#   Delta = (topo_win_gg - SM_win_gg) + (topo_win_qq - SM_win_qq).
#   T (column) = K = 1 + Delta/C;  sigma_with = C + Delta;  sigma_bin(run) = sigma_run x (events in bin and window)/(all events of the run).
#
# GENERATION: MadGraph 2.9.27, model sm. g g runs: generate g g > t t~ > w+ b w- b~, w+ > l+ vl, w- > l- vl~;
#   q q~ runs: define q = u c d s; define q~ = u~ c~ d~ s~; generate q q~ > (same){qqmir}. run_card: nevents = 1000000, iseed = {seeds},
#   pdlabel = lhapdf (CT25NNLO member 0), use_syst = False, event_norm = average, lepton cuts off, hel_recycling = False (output flag),
#   mc_grouped_subproc = False in the g g runs; param_card: m_t = {MT:g}, Gamma_t = {WT:g} GeV.
#
# MATRIX-ELEMENT MODIFICATION (topo_win and topo_full runs; Fuks et al. arXiv:2411.18962 (16), (19)-(20))
#   g g : |M|^2 -> |M1|^2 |G1/G0|^2 + |M8|^2 |G8/G0|^2,  |M1|^2 from MadGraph's colour sum with CF1(I,J) = 2/3, |M8|^2 = |M|^2 - |M1|^2.
#   q q~: |M|^2 -> |M|^2 |G8/G0|^2.  Applied for W <= {whi():g} GeV, p* < {CFG['p_rw']:.0f} GeV, E >= -20 GeV (R = 1 elsewhere, as his).
#   G/G0 by bilinear interpolation (our CALCGREEN) on E = -20..4 GeV step 0.025, p = 0..80 GeV step 0.125.
#
# GREEN'S FUNCTIONS: S-wave Lippmann-Schwinger equation solved as in JKT Sec. 3 (LSSolver, n = 120 per subinterval);
#   G0(E,p) = 1/(E + i Gamma_t - p^2/m_t); alpha_s(Q) = CT25NNLO's (LHAPDF).
#   SINGLET (R1): {sing}, cut (19) at q_cut = 50 MeV, d3 term (20) as an energy shift; R1 with the unshifted G0.
#   OCTET (R8), colour factor +1/6:
#     [used for T]  {head}: {desc[head]}
""" + "".join(f"#     [column T_{o}] {o}: {desc[o]}\n" for o in extra) + """#
# RUNS (1,000,000 unweighted events each; 'events' = events passing THIS file's window)
""" + "\n".join(run_lines) + """
#
# OBSERVABLES AND BINS (as his): lower edge included, upper excluded. HT_ttbar = pT(t) + pT(tbar); abs_y_boost = |y_t + y_tbar|/2;
#   abs_y_top_x_pT_top: x1 = |y_t|, x2 = pT_t.
# ERRORS: MC statistics only; delta_err from the four window runs; T_err = sqrt((delta_err/C)^2 + (Delta*C_err/C^2)^2).
"""
    with open(path, "w") as f:
        f.write(hdr)
        f.write("distribution,bin,x1_lo,x1_hi,x2_lo,x2_hi,sigma_without_pb,sigma_without_err,sigma_with_pb,sigma_with_err,"
                "T,T_err,delta_pb,delta_err," + ",".join(f"T_{o}" for o in extra) + "\n")
        for r in rows:
            f.write(",".join(x if isinstance(x, str) else f"{x:.8g}" for x in r) + "\n")


OCTETCHECK = dict(name="topo_win_gg_none_s7850", seed=7850, table=os.path.join(RES, "green8_none_p80.dat"),
                  # LO gg -> tt -> 6 bodies over full phase space (denominators of the gg rows): full colour, and the
                  # colour-singlet projection of Fuks (19)-(20) (Fuks et al.'s own sample)
                  gg_lo=[("SM_gg_full_s7860", "TT10_gg_full", 7860, False), ("SM_gg1_full_s7870", "TT10_gg_singlet", 7870, True)])


def octetcheck():
    """Run 6: the octet effect measured directly with MadGraph.
        Fuks (16) singlet only:  |M|^2 -> |M|^2 + (|G1/G0|^2 - 1)|M_1|^2   (gg; qqbar unchanged)
        run 6 (singlet + octet): |M|^2 -> |G1/G0|^2 |M_1|^2 + |G8/G0|^2 |M_8|^2   (gg),  |G8/G0|^2 |M|^2   (qqbar)
    The run-6 topo_win_gg output (F_GREEN_FACTOR_OCTET) with an octet table G8/G0 = 1 IS the singlet-only matrix
    element, so one extra gg window run (1M events, seed 7850) gives it; singlet-only qqbar = SM_win_qq. In 326 <= W <= 350:
        Delta_1   = topo_win_gg_none - SM_win_gg
        Delta_1+8 = (topo_win_gg_freeze25 - SM_win_gg) + (topo_win_qq_freeze25 - SM_win_qq)     (the run-6 K-factors)
        octet effect = Delta_1+8 - Delta_1;  K = 1 + Delta/SM_full in the first bin of m_tt (CMS, ATLAS) and H_T.
    -> results/6_tfactors/run6_final_kfactors_13TeV/octet_check.txt"""

    assert RUN == "run6_final_kfactors_13TeV", "octetcheck uses run 6's outputs and samples"
    if not os.path.exists(OCTETCHECK["table"]):
        shape = load_table(T8_FILES["freeze25"]).shape
        write_table(OCTETCHECK["table"], np.ones(shape, dtype=complex))
    T8_FILES["none"] = OCTETCHECK["table"]
    d = window_dir("T", "gg")
    put(T1_FILE, os.path.join(d, "SubProcesses", "green_table.dat"))
    os.makedirs(os.path.join(EV_DIR, "logs"), exist_ok=True)
    try:
        generate(OCTETCHECK["name"], d, OCTETCHECK["seed"], 1000000, "none", "singlet_only",
                 ev_path(OCTETCHECK["name"], False))
    finally:                                  # leave the octet table run 6 used in place
        put(T8_FILES["freeze25"], os.path.join(d, "SubProcesses", "green8_table.dat"))
    wrapper = os.path.join(LHAPDF_DIR, "lhapdf-config")
    for name, dname, seed, singlet in OCTETCHECK["gg_lo"]:
        dd = os.path.join(OUTPUT, RUN, dname)
        if not os.path.isdir(dd):
            mg_output(dd, CFG["process_gg"])
            for f in sorted(glob.glob(os.path.join(dd, "SubProcesses", "P*_gg_*", "matrix*.f"))):
                if singlet:
                    print(f"singlet colour {os.path.relpath(f, OUTPUT)}:", patch_singlet(f))
            set_cards(dd, wrapper)
        generate(name, dd, seed, 1000000, None, name.rsplit("_s", 1)[0], ev_path(name, False))

    seeds = {(ch, o): s for ch, o, s in CFG["window_runs"]}
    grp = {"SM_gg": window_name("gg", None, seeds[("gg", None)]),
           "T_gg_18": window_name("gg", "freeze25", seeds[("gg", "freeze25")]),
           "SM_qq": window_name("qq", None, seeds[("qq", None)]),
           "T_qq_18": window_name("qq", "freeze25", seeds[("qq", "freeze25")]),
           "T_gg_1": OCTETCHECK["name"]}
    C = load_group(f"{CFG['names'][0]}_s*.npz")
    oC = observables(C["t"], C["tb"])
    ev = {}
    for k, n in grp.items():
        g = load_group(f"{n}.npz")
        E, ps = E_pstar(g["t"], g["tb"])
        ev[k] = dict(g=g, o=observables(g["t"], g["tb"]), W=E + 2 * MT, ps=ps)
    lines = [f"STAGE 6, {RUN}: octet check -- singlet only (Fuks (16)) vs singlet + octet (run 6), MadGraph window runs",
             f"samples: " + ", ".join(f"{k} = {n}" for k, n in grp.items()),
             f"SM_full: {len(C['w'])} events, sigma = {C['w'].sum():.6g} pb", ""]
    log = lambda s="": (print(s), lines.append(s))
    firsts = [("m_ttbar_250-3500", "first CMS m_tt bin [250, 400]"), ("m_ttbar_325-2000", "first ATLAS m_tt bin [325, 400]"),
              ("HT_ttbar_0-2000", "first H_T bin [0, 90]")]
    for pcut in (100.0, 50.0):
        sel = {k: (v["W"] >= wlo()) & (v["W"] <= whi()) & (v["ps"] < pcut) for k, v in ev.items()}
        s = {k: (v["g"]["w"][sel[k]].sum(), np.sqrt((v["g"]["w"][sel[k]] ** 2).sum())) for k, v in ev.items()}
        diff = lambda a, b: (s[a][0] - s[b][0], np.hypot(s[a][1], s[b][1]))
        d1 = diff("T_gg_1", "SM_gg")
        dgg, dqq = diff("T_gg_18", "SM_gg"), diff("T_qq_18", "SM_qq")
        d18 = (dgg[0] + dqq[0], np.hypot(dgg[1], dqq[1]))
        # octet effect = (T_gg_18 - T_gg_1) + (T_qq_18 - SM_qq): SM_gg cancels; four independent samples
        oc = (d18[0] - d1[0], np.sqrt(s["T_gg_18"][1] ** 2 + s["T_gg_1"][1] ** 2 + s["T_qq_18"][1] ** 2 + s["SM_qq"][1] ** 2))
        log(f"window {wlo():g} <= W <= {whi():g}, p* < {pcut:.0f} GeV  [pb]")
        log(f"  SM_win gg {s['SM_gg'][0]:.5f}, qq {s['SM_qq'][0]:.5f}; topo gg singlet only {s['T_gg_1'][0]:.5f}, "
            f"gg singlet+octet {s['T_gg_18'][0]:.5f}, qq octet {s['T_qq_18'][0]:.5f}")
        log(f"  Delta singlet only      {d1[0]:+.5f} +- {d1[1]:.5f}")
        log(f"  Delta singlet + octet   {d18[0]:+.5f} +- {d18[1]:.5f}   (gg {dgg[0]:+.5f}, qq {dqq[0]:+.5f})")
        log(f"  octet effect            {oc[0]:+.5f} +- {oc[1]:.5f}  = {oc[0] / d1[0]:+.1%} of the singlet-only excess")
        for dname, label in firsts:
            spec = DISTRIBUTIONS[dname]
            c = binned(oC, C["w"], spec)[0][0]
            part = lambda k: binned(ev[k]["o"], ev[k]["g"]["w"] * sel[k], spec)[0][0]
            k1 = 1 + (part("T_gg_1") - part("SM_gg")) / c
            k18 = 1 + (part("T_gg_18") - part("SM_gg") + part("T_qq_18") - part("SM_qq")) / c
            log(f"  {label:32s} K singlet only {k1:.4f} -> singlet + octet {k18:.4f}   (Delta K {k18 - k1:+.4f})")
        log()

    # the staircase: the same singlet-only / singlet+octet changes over ever larger denominators
    lo = {k: load_group(f"{n}.npz") for k, n in (("gg", OCTETCHECK["gg_lo"][0][0]), ("gg1", OCTETCHECK["gg_lo"][1][0]))}
    olo = {k: observables(g["t"], g["tb"]) for k, g in lo.items()}
    log(f"LO gg (full colour) {lo['gg']['w'].sum():.5g} pb, gg colour singlet {lo['gg1']['w'].sum():.5g} pb, "
        f"pp (SM_full) {C['w'].sum():.5g} pb; singlet share of gg {lo['gg1']['w'].sum() / lo['gg']['w'].sum():.3f}")
    for pcut in (100.0, 50.0):
        sel = {k: (v["W"] >= wlo()) & (v["W"] <= whi()) & (v["ps"] < pcut) for k, v in ev.items()}
        log(f"staircase, p* < {pcut:.0f} GeV: K = 1 + Delta/LO in the first bin (+- MC)")
        for dname, label in firsts:
            spec = DISTRIBUTIONS[dname]
            part = lambda k: binned(ev[k]["o"], ev[k]["g"]["w"] * sel[k], spec)
            den = {"gg1": binned(olo["gg1"], lo["gg1"]["w"], spec), "gg": binned(olo["gg"], lo["gg"]["w"], spec),
                   "pp": binned(oC, C["w"], spec)}
            P = {k: (part(k)[0][0], part(k)[1][0]) for k in ("T_gg_1", "T_gg_18", "SM_gg", "T_qq_18", "SM_qq")}
            err = lambda *ks: float(np.sqrt(sum(P[k][1] ** 2 for k in ks)))
            d1 = (P["T_gg_1"][0] - P["SM_gg"][0], err("T_gg_1", "SM_gg"))
            d18g = (P["T_gg_18"][0] - P["SM_gg"][0], err("T_gg_18", "SM_gg"))
            d18 = (d18g[0] + P["T_qq_18"][0] - P["SM_qq"][0], err("T_gg_18", "SM_gg", "T_qq_18", "SM_qq"))
            rows = [("gg -> tt, colour-singlet events only (Fuks et al.)", d1, "gg1"),
                    ("gg -> tt, all colours, singlet modified", d1, "gg"),
                    ("gg -> tt, singlet + octet", d18g, "gg"),
                    ("pp -> tt, singlet modified", d1, "pp"),
                    ("pp -> tt, singlet + octet (final)", d18, "pp")]
            log(f"  {label}")
            for txt, (dl, de), dk in rows:
                c, ce = den[dk][0][0], den[dk][1][0]
                ke = np.hypot(de / c, dl * ce / c ** 2)          # as kfactors: T_err
                log(f"    {txt:52s} K = {1 + dl / c:.4f} +- {ke:.4f}   (LO {c:.5f} pb)")
        log()
    path = os.path.join(RES_RUN, "octet_check.txt")
    open(path, "w").write("\n".join(lines) + "\n")
    print("->", path)
    return path


def mt_summary():
    """Run 7: K (p* < 100) per bin against m_t, from each mass point's kfactors_final/, for the bins where K != 1 at
    some mass -> results/6_tfactors/run7_mtscan_summary.txt."""

    names = [n for n in RUNS if n.startswith("run7_mt")]
    tabs = {}
    for n in names:
        d = os.path.join(RES, n, "kfactors_final")
        if os.path.isdir(d):
            tabs[RUNS[n]["mt"]] = {f: [l for l in open(os.path.join(d, f)).read().splitlines() if l and l[0] not in "#"
                                       and not l[0].isalpha()] for f in sorted(os.listdir(d)) if f.endswith(".csv")}
    if not tabs:
        sys.exit("no run7 kfactors_final/ yet: kfactors + export each run7_mt*_13TeV first")
    ms = sorted(tabs)
    out = ["K = 1 + Delta/SM_full per bin (p* < 100 GeV, octet coulomb_as25, full JKT singlet) against m_t; "
           "bins with K != 1 at some mass; +- = MC error", "m_t [GeV]:      " + "".join(f"{m:>17g}" for m in ms),
           "Gamma_t [GeV]:  " + "".join(f"{top_width(m):>17g}" for m in ms), ""]
    for f in tabs[ms[0]]:
        for i, l in enumerate(tabs[ms[0]][f]):
            rows = [tabs[m][f][i].split(",") for m in ms if f in tabs[m]]
            if all(r[-4] == "1.0000" for r in rows):
                continue
            edges = ",".join(rows[0][:-4])
            out.append(f"{f[:-4]:28s} [{edges}]".ljust(48) + "".join(f"{r[-4]:>9s} +-{r[-3]:>6s}" for r in rows))
    path = os.path.join(RES, "run7_mtscan_summary.txt")
    open(path, "w").write("\n".join(out) + "\n")
    print("\n".join(out))
    return path


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    if cmd in ("setup", "smoke", "run", "one", "list", "analyse", "kfactors", "plotk", "arbiter", "export", "octetcheck"):
        if len(sys.argv) < 3:
            sys.exit(f"usage: 6_tfactors.py {cmd} RUN   with RUN one of {list(RUNS)}")
        select(sys.argv[2])
    if cmd == "mtscan":                            # run 7: K vs m_t over the five mass points
        mt_summary()
    elif cmd == "mtscan_names":                    # run 7 run names, one per line (for submit_mtscan.sh)
        print("\n".join(n for n in RUNS if n.startswith("run7_mt")))
    elif cmd == "table":
        table_window() if sys.argv[2:3] == ["run5"] else table()
    elif cmd == "kfactors":
        if "--reflect-qq" in sys.argv:            # preview of the q q~ mirror fix -> RUN/preview_qq_mirror/
            kfactors(his=[a for a in sys.argv[3:] if a != "--reflect-qq"], reflect_qq=True, tag="preview_qq_mirror")
        else:
            kfactors(his=sys.argv[3:])
    elif cmd == "export":                          # easy-to-read K-factors -> RUN/[SUBDIR/]kfactors_final/
        print(export(sys.argv[3] if len(sys.argv) > 3 else ""))
    elif cmd == "arbiter":
        arbiter()
    elif cmd == "octetcheck":                      # run 6: singlet only vs singlet + octet, one extra gg window run
        octetcheck()
    elif cmd == "plotk":
        for f in sorted(glob.glob(os.path.join(RES_RUN, "kfactors_*.csv"))):
            print(plot_kfactors(f, os.path.join(RES_RUN, "plots", re.search(r"pstar\d+", f).group(0))))
    elif cmd == "smoke" and CFG.get("window"):
        do_runs(smoke=True)
        kfactors(subdir="smoke")
    elif cmd == "setup":
        setup()
    elif cmd == "smoke":
        do_runs(smoke=True)
        analyse(subdir="smoke")
    elif cmd == "run":
        do_runs()
    elif cmd == "one":                             # one sample by name (cluster jobs): one RUN NAME
        do_runs(smoke=sys.argv[3].endswith("_smoke"), only=sys.argv[3])
    elif cmd == "list":                            # the sample names of RUN, one per line
        print("\n".join(r[0] for r in runs()))
    elif cmd == "analyse":
        analyse()
    else:
        sys.exit(__doc__)
