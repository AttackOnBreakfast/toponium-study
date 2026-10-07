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
alpha_s everywhere = CT25NNLO's alpha_s(Q), as CTEQ use it: MadGraph takes alpha_s from the PDF (LHAPDF), and the
JKT potential behind G uses the same LHAPDF alpha_s(Q) in (13) instead of JKT's two-loop (14).

Run from toponium_study/ with MadGraph's extra dependency (RUN = one of the names above):
    uv run --with six python 6_tfactors.py table           # CT25 alpha_s(Q) from LHAPDF, G/G0 table (all runs)
    uv run --with six python 6_tfactors.py setup RUN       # LHAPDF dir, MadGraph outputs, MATRIX1 patches, cards
    uv run --with six python 6_tfactors.py smoke RUN       # 20k WITHOUT + 20k WITH, sanity check
    uv run --with six python 6_tfactors.py run RUN         # 20 runs of 1M events (resumable)       (~4-5 h)
    uv run --with six python 6_tfactors.py analyse RUN     # T-factors per bin, CSV, plots
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
"""

import glob
import os
import re
import shutil
import subprocess
import sys
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

CT25_ID = 99999000          # local LHAPDF id for CT25NNLO (its own SetIndex overflows a Fortran integer)
CT25_SOURCES = (os.path.join(HERE, "CT25NNLO"), os.path.expanduser("~/Downloads/CT25NNLO"))
ENERGY, EBEAM = 13, 6500.0

RUNS = {   # Fuks (14), l = e (top side), l' = mu (antitop side); WITHOUT and WITH always the same process
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
}
RUN, CFG = None, None
SM_DIR = TOPO_DIR = EV_DIR = EV_DIR_WITHOUT = GREEN_FILE = GREEN8_FILE = RES_RUN = None


def select(run):
    """Point all run-dependent paths at RUNS[run]: output/<run>/{without, with, events}, results/6_tfactors/<run>."""

    global RUN, CFG, SM_DIR, TOPO_DIR, EV_DIR, EV_DIR_WITHOUT, GREEN_FILE, GREEN8_FILE, RES_RUN
    if run not in RUNS:
        sys.exit(f"run must be one of {list(RUNS)}")
    RUN, CFG = run, RUNS[run]
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


def runs(smoke=False):
    """(name, directory, seed, nevents): 10 WITHOUT + 10 WITH runs of 1M events at 13 TeV, one seed each."""

    (n0, n1), (s0, s1) = CFG["names"], CFG["seeds"]
    own_without = CFG["without_from"] is None          # run 3 reuses run 1's WITHOUT sample
    if smoke:
        a, b = CFG["smoke_seeds"]
        return [(f"{n0}_smoke", SM_DIR, a, 20000)] * own_without + [(f"{n1}_smoke", TOPO_DIR, b, 20000)]
    without = [(f"{n0}_s{s0 + k}", SM_DIR, s0 + k, 1000000) for k in range(10)] if own_without else []
    return without + [(f"{n1}_s{s1 + k}", TOPO_DIR, s1 + k, 1000000) for k in range(10)]


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
        subprocess.run(["clang++", "-std=c++17", "-O2", src, *flags, "-o", exe], check=True)
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
    d = np.loadtxt(path or TABLE, skiprows=1)
    return (d[:, 0] + 1j * d[:, 1]).reshape(NE, NP)


def write_table(path, R):
    with open(path, "w") as f:
        f.write(f"{NE} {NP} {E0} {DE} {P0} {DP}\n")
        for z in R.ravel():
            f.write(f"{z.real:.10e} {z.imag:.10e}\n")


def octet_ratio(E, p, alpha_s):
    """G8/G0 for the colour octet: exact Coulomb solution (CoulombExact) for the repulsive potential
    V8 = +alpha_s/(6r), i.e. colour factor -1/6 = C_F x (-1/8). No bound states."""

    return CoulombExact.G(E, p, -alpha_s / 8) / Constants.G0(E, p)


def interp_table(tab, E, p):
    """Bilinear interpolation exactly as CALCGREEN in matrix1.f (F_MODULE): G/G0 at (E, p)."""

    a, b = (np.asarray(E) - E0) / DE, (np.asarray(p) - P0) / DP
    i = np.clip(np.trunc(a).astype(int), 0, NE - 2)
    j = np.clip(np.trunc(b).astype(int), 0, NP - 2)
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
    new, n = re.subn(r"      DATA \(CF\(I,  1\),I=  1,  2\).*?C     1 T\(2,1,3,6\)\n", lambda _: F_CF_SINGLET, src,
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
                     ("      IF(SDE_STRAT.EQ.1)THEN\n", F_GREEN_FACTOR_QQ + "      IF(SDE_STRAT.EQ.1)THEN\n")):
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
    m1 = once(m1, "C     1 T(2,1,3,6)\n", "C     1 T(2,1,3,6)\n" + F_CF1)
    m1 = once(m1, "      MATRIX1 = 0.D0\n", F_KINEMATICS + "      MATRIX1 = 0.D0\n")
    factor = F_GREEN_FACTOR_OCTET if CFG["octet"] else F_GREEN_FACTOR
    m1 = once(m1, "      IF(SDE_STRAT.EQ.1)THEN\n", factor + "      IF(SDE_STRAT.EQ.1)THEN\n")
    assert ("5.333333333333333D+00" in m1) != CFG["singlet_only"], "colour matrix not as this run expects"
    open(path, "w").write(F_MODULE(octet=CFG["octet"]) + head + sep + m1 + rest)
    return "patched"


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
    open(os.path.join(LHAPDF_DIR, "pdfsets.index"), "w").write(index + f"\n{CT25_ID} CT25NNLO 1\n")
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
    open(os.path.join(ct25, "CT25NNLO.info"), "w").write(
        info + "\n# added for LHAPDF 6.5 / MadGraph (not in the distributed file; Lambdas unused by MadGraph):\n"
        f"AlphaS_FlavorScheme: variable\nAlphaS_NumFlavors: 5\nAlphaS_Lambda4: {lam4:.4f}\nAlphaS_Lambda5: {lam5:.4f}\n")
    dst = os.path.join(LHAPDF_DIR, "lhapdf.conf")
    if not os.path.lexists(dst):
        os.symlink(os.path.join(sysdir, "lhapdf.conf"), dst)
    # MadGraph reads the set list from `lhapdf-config --datadir` (not LHAPDF_DATA_PATH): a wrapper that
    # reports this dir and passes everything else to the real lhapdf-config
    wrapper = os.path.join(LHAPDF_DIR, "lhapdf-config")
    open(wrapper, "w").write(f'#!/bin/bash\nif [ "$1" = "--datadir" ]; then echo {LHAPDF_DIR}; '
                             f'else exec {real} "$@"; fi\n')
    os.chmod(wrapper, 0o755)
    print(f"LHAPDF data dir {LHAPDF_DIR}: CT25NNLO (alpha_s(mZ) = {amz}) as id {CT25_ID}")

    # MadGraph outputs (WITHOUT and WITH: the same process)
    if CFG["archived"]:
        sys.exit(f"{RUN} is archived (its WITH matrix1.f reads the G table from its old path): not set up again")
    own = [TOPO_DIR] if CFG["without_from"] else [SM_DIR, TOPO_DIR]   # never touch another run's directory
    for d in own:
        if os.path.isdir(d):
            print(f"{d} exists")
            continue
        os.makedirs(os.path.dirname(d), exist_ok=True)
        script = os.path.join(OUTPUT, f"_{os.path.basename(d)}.mg5")
        open(script, "w").write("set automatic_html_opening False\nimport model sm\n"
                                f"generate {CFG['process']}\noutput {d} --hel_recycling=False\n")
        subprocess.run([sys.executable, "-W", "ignore", os.path.join(MG_ROOT, "bin", "mg5_aMC"), script],
                       cwd=MG_ROOT, check=True, env=mg_env())
        os.remove(script)

    # matrix elements: run 2 = singlet only in both samples; WITH = Fuks (16) in the gg subprocess (qqbar unchanged)
    gg = lambda d: os.path.join(d, "SubProcesses", "P1_gg_blvlbxlvl", "matrix1.f")
    if CFG["singlet_only"]:
        print("WITHOUT matrix1.f (gg):", patch_singlet(gg(SM_DIR)))
        print("WITH    matrix1.f (gg):", patch_singlet(gg(TOPO_DIR)))
    print("WITH    matrix1.f (gg):", patch_matrix1(gg(TOPO_DIR)))
    if CFG["octet"]:
        print("WITH    matrix1.f (qq):", patch_qq(os.path.join(TOPO_DIR, "SubProcesses", "P1_qq_blvlbxlvl", "matrix1.f")))
    if os.path.exists(TABLE):
        shutil.copy(TABLE, GREEN_FILE)
    if CFG["octet"] and os.path.exists(OCTET_TABLE):
        shutil.copy(OCTET_TABLE, GREEN8_FILE)

    # cards (identical in WITHOUT and WITH)
    for d in own:
        set_card(os.path.join(d, "Cards", "run_card.dat"), {
            "pdlabel": "lhapdf", "lhaid": CT25_ID, "use_syst": "False", "event_norm": "average",
            "ebeam1": EBEAM, "ebeam2": EBEAM, "ptl": "0.0", "etal": "-1.0", "drll": "0.0"})   # no cuts
        pc = os.path.join(d, "Cards", "param_card.dat")
        src, n = re.subn(r"^DECAY\s+6\s+\S+\s+# WT", "DECAY   6 1.490000e+00 # WT", open(pc).read(), flags=re.M)
        assert n == 1 and re.search(r"^\s+6 1.730000e\+02 # MT", src, re.M), "param_card MT/WT"
        open(pc, "w").write(src)
        set_config(os.path.join(d, "Cards", "me5_configuration.txt"), {
            "run_mode": 2, "nb_core": 10, "lhapdf_py3": wrapper, "automatic_html_opening": "False"})
    print("cards: 13 TeV, CT25NNLO, no cuts, m_t = 173, Gamma_t = 1.49, event_norm = average, multicore 10")


# =====================================================================================
# 3. runs: generate, keep per-event t and tbar 4-momenta, delete the LHE file
# =====================================================================================
def ev_path(name, smoke):
    return os.path.join(EV_DIR, "smoke" if smoke else "", f"{name}.npz")


def do_runs(smoke=False):
    if not os.path.exists(GREEN_FILE) or (CFG["octet"] and not os.path.exists(GREEN8_FILE)):
        sys.exit("no G table in the WITH directory: run 'table' after 'setup'")
    os.makedirs(os.path.join(EV_DIR, "smoke" if smoke else "", "logs"), exist_ok=True)
    for name, d, seed, nev in runs(smoke):
        out = ev_path(name, smoke)
        if os.path.exists(out):
            print(f"{name}: done")
            continue
        set_card(os.path.join(d, "Cards", "run_card.dat"), {"iseed": seed, "nevents": nev})
        if os.path.isdir(os.path.join(d, "Events", name)):
            shutil.rmtree(os.path.join(d, "Events", name))           # an interrupted earlier attempt
        t0 = time.time()
        log = os.path.join(os.path.dirname(out), "logs", f"{name}.log")
        print(f"{name}: generating {nev} events (seed {seed}) ...", flush=True)
        with open(log, "w") as f:
            r = subprocess.run([sys.executable, "-W", "ignore", "bin/generate_events", name, "-f"],
                               cwd=d, stdout=f, stderr=subprocess.STDOUT, env=mg_env())
        lhe = os.path.join(d, "Events", name, "unweighted_events.lhe.gz")
        if r.returncode != 0 or not os.path.exists(lhe):
            sys.exit(f"{name}: MadGraph failed, see {log}")
        w, pt, ptb = Events.read_lhe(lhe)
        np.savez(out, w=w, t=pt, tb=ptb, seed=seed, sample="with" if d == TOPO_DIR else "without", run=RUN)
        for f in glob.glob(os.path.join(d, "Events", name, "*.lhe*")):
            os.remove(f)
        print(f"{name}: {len(w)} events, sigma = {w.mean():.5g} pb, {(time.time() - t0) / 60:.1f} min", flush=True)


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
    return W - 2 * M_T, np.sqrt(np.maximum(lam, 0)) / (2 * W)


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


def draw_ratio(ax, dname, spec, T, Te, title=None):
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
    ax.set_ylabel(r"$T = \sigma_{\rm with} \,/\, \sigma_{\rm without}$")
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
                       f"{CFG['desc']};\n13 TeV, LO MadGraph, CT25NNLO with its alpha_s(Q);\n{n0 / 1e6:g}M events without + "
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
             "WITH: gg |M|^2 -> |M|^2 + (|G/G0|^2 - 1)|M_1|^2 (Fuks (16), singlet (19)-(20)) for W <= 350, p* <= 50 GeV;",
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
                "(19)-(20)); G from our JKT solver with CT25NNLO alpha_s(Q). Parton-level tops (t = b e+ ve, "
                "tbar = b~ mu- vm~). Errors: MC statistics of the two independent samples.\n"
                "# HT_ttbar = pT(t) + pT(tbar); abs_y_boost = |y_t + y_tbar|/2; abs_y_top_x_pT_top: x1 = |y_t|, x2 = pT_t.\n")
        f.write("distribution,bin,x1_lo,x1_hi,x2_lo,x2_hi,sigma_without_pb,sigma_without_err,sigma_with_pb,"
                "sigma_with_err,T,T_err\n")
        for r in rows:
            f.write(",".join(x if isinstance(x, str) else f"{x:.8g}" for x in r) + "\n")
    plot_all(results, len(g0["w"]), len(g1["w"]), os.path.join(out, "plots", "all_bins.png"))
    log(f"-> {os.path.relpath(path, HERE)}, plots/ ({len(DISTRIBUTIONS)} figures + all_bins.png)")
    if CFG["octet"]:
        octet_effect(g0, g1, o0, o1, subdir, out, log)
    open(os.path.join(out, "summary.txt"), "w").write("\n".join(lines) + "\n")


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
    win = (E + 2 * M_T <= 350) & (ps <= 50) & (E >= E0)
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

if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    if cmd in ("setup", "smoke", "run", "analyse"):
        if len(sys.argv) < 3:
            sys.exit(f"usage: 6_tfactors.py {cmd} RUN   with RUN one of {list(RUNS)}")
        select(sys.argv[2])
    if cmd == "table":
        table()
    elif cmd == "setup":
        setup()
    elif cmd == "smoke":
        do_runs(smoke=True)
        analyse(subdir="smoke")
    elif cmd == "run":
        do_runs()
    elif cmd == "analyse":
        analyse()
    else:
        sys.exit(__doc__)
