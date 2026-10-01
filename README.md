# Toponium threshold study

A self-contained physics study layered on MadGraph (not part of the framework): toponium /
threshold effects in off-shell `t t~`, following Fuks, Hagiwara, Ma, Zheng, arXiv:2411.18962
(EPJC 85 (2025) 157), with the Green's function computed as in Jezabek, Kuhn, Teubner,
Z. Phys. C 56 (1992) 653 (JKT).

**Read `writeup.tex` / `writeup.pdf` first.** It explains what Fuks et al. do, what we do,
how they are compared, the solver step by step with JKT equation numbers, and all results.

## Files

Seven Python files: one library and six scripts you run.

`toponium.py` is the only library. It holds six classes, top to bottom:

| Class | Side | What it is |
|---|---|---|
| `Constants` | both | m_t = 173, Gamma_t = 1.49, C_F, n_F, and the free Green's function `G0` |
| `JKTPotential` | ours | JKT's potential (Eqs. 13-15, cut by 19-20); `JKTPotential.coulomb(alpha)` gives fixed-alpha Coulomb |
| `LSSolver` | ours | Lippmann-Schwinger equation solved as in JKT Sec. 3 (Eqs. 16-23): `LSSolver(pot).G(E, p)` |
| `CoulombExact` | ours | exact analytic G for fixed-alpha Coulomb: the ground truth used to test `LSSolver` |
| `FuksTables` | Fuks's | their published tables (swData) and their `CALCGREEN` routine, used as-is |
| `Events` | both | MadGraph events (from the caches) and the per-event weight \|G/G0\|^2 (Fuks 15) |

Every class opens with an EQUATION MAP (paper equation -> method). `reference/` holds both
PDFs and swData.

## Run order (from inside this folder: `cd toponium_study`)

Python 3.11 environment managed by uv (`pyproject.toml`, `uv.lock`, `.python-version`; the
environment itself is `.venv/`). First time only: `uv sync`. Then run each stage as
`uv run python 0_solver_checks.py`, and so on.

| Stage | Script | What it does | Time |
|---|---|---|---|
| 0 | `0_solver_checks.py` | our solver vs exact Coulomb, convergence, q_cut independence | ~3 min |
| 1 | `1_validate_jkt.py` | our solver reproduces JKT's own Figs. 2-5 | ~30 s |
| 2 | `2_identify_swdata.py` | which potential is in swData (answer: fixed-alpha Coulomb, alpha = 0.15) | ~3 min |
| 3 | `3_fuks_fig1.py` | Fuks Fig. 1 from swData, ours alongside | ~10 s |
| 4 | `4_fuks_figs2to5.py` | Fuks Figs. 2-5: same events re-weighted with Fuks's G and ours | ~40 s |
| 5 | `5_comparison_tables.py` | CSV tables (parameters, V vs Q^2, G on E-p grids; colour singlet) for comparing with an independent solver | ~1 min |

Each stage writes `results/<stage>/summary.txt` (+ figures). 

## Conventions

`m_t = 173 GeV`, `Gamma_t = 1.49 GeV`, `C_F = 4/3`, `n_F = 5`; running potentials anchored at
`alpha_s(m_Z) = 0.12` as the paper states. swData itself corresponds to fixed alpha = 0.15.
G everywhere has the sign of Fuks (11) / JKT (6); swData stores -G.

## Stage 4 event samples (regeneration recipe)

`../output/` was cleared on 2026-09-24. Stage 4 runs from its per-event kinematics caches
(`results/4_fuks_figs2to5/_cache_singlet330.npz`, `_cache_full.npz`: W, E, p*, m_tL, m_tH,
weight). **Keep these files**: they are the only copy. To regenerate the samples (~10 min each):

1. `./bin/mg5_aMC` from the MG5 root:
   `import model sm` / `generate g g > t t~ > b e+ ve b~ mu- vm~` /
   `output output/GGtoTTbar_Singlet_Fuks13TeV_W330 --hel_recycling=False`
   (and again as `output/GGtoTTbar_Full_Fuks13TeV` for the full-colour sample).
2. Singlet only: in `SubProcesses/P1_gg_blvlbxlvl/matrix1.f` set every colour-matrix entry
   `CF(I,J)` to `6.666666666666666D-01` (Fuks (19)-(20): the colour-singlet projection).
3. `SubProcesses/genps.f`, subroutine `GENCMS`: `TAUMIN = 330D0**2/S`, `TAUMAX = 350D0**2/S`
   (340 for the full-colour sample), so only 330 <= W = sqrt(shat) <= 350 GeV is sampled.
4. `SubProcesses/dummy_fct.f`, function `dummy_cuts`: with t = P(3)+P(4)+P(5) and
   tbar = P(6)+P(7)+P(8), keep events with 330 <= W <= 350 and
   p* = sqrt(lambda(W^2, m_t^2, m_tbar^2))/(2W) < 100 GeV.
5. `Cards/run_card.dat`: `ebeam1 = ebeam2 = 6500`, `pdlabel = lhapdf`, `lhaid = 14400`
   (CT18NLO), `nevents = 500000`, `use_syst = False`. If `Cards/me5_configuration.txt` has a
   stale `mg5_path`, point it at this MG5 checkout.
6. `./bin/generate_events run_01 -f`; delete the old `_cache_*.npz` so stage 4 re-reads the events.
   Check: sigma(singlet)/sigma(full) in 340-350 GeV should be about 2/7 (0.281 obtained).
