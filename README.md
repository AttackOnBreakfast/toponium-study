# Toponium threshold study

A self-contained physics study layered on MadGraph (not part of the framework): toponium /
threshold effects in off-shell `t t~`, following Fuks, Hagiwara, Ma, Zheng, arXiv:2411.18962
(EPJC 85 (2025) 157), with the Green's function computed as in Jezabek, Kuhn, Teubner,
Z. Phys. C 56 (1992) 653 (JKT).

**Read `writeup.tex` / `writeup.pdf` first.** It explains what Fuks et al. do, what we do,
how they are compared, the solver step by step with JKT equation numbers, and all results.

## Files

Eight Python files: one library and seven scripts you run.

`toponium.py` is the only library. It holds six classes, top to bottom:

| Class | Side | What it is |
|---|---|---|
| `Constants` | both | m_t = 173, Gamma_t = 1.49, C_F, n_F, and the free Green's function `G0` |
| `JKTPotential` | ours | JKT's potential (Eqs. 13-15, cut by 19-20); `JKTPotential.coulomb(alpha)` gives fixed-alpha Coulomb |
| `LSSolver` | ours | Lippmann-Schwinger equation solved as in JKT Sec. 3 (Eqs. 16-23): `LSSolver(pot).G(E, p)` |
| `CoulombExact` | ours | exact analytic G for fixed-alpha Coulomb: the ground truth used to test `LSSolver` |
| `FuksTables` | Fuks's | their published tables (swData) and their `CALCGREEN` routine, used as-is |
| `Events` | both | MadGraph events (from the caches) and the per-event weight \|G/G0\|^2 (Fuks 15) |

Every class opens with an EQUATION MAP (paper equation -> method). `reference/` holds the papers
and swData:

| File | Paper |
|---|---|
| `BF01474740.pdf` | Jezabek, Kuhn, Teubner, Z. Phys. C 56 (1992) 653 (JKT): potential and solver |
| `2411.18962v2.pdf` | Fuks, Hagiwara, Ma, Zheng, arXiv:2411.18962: matrix-element reweighting; `fuks2024/` = their swData tables |
| `0804.1014v2.pdf` | Hagiwara, Sumino, Yokoya, arXiv:0804.1014 (HSY): colour-octet potential, factor -1/6 |
| `2503.22382v2.pdf` | CMS, arXiv:2503.22382: threshold excess in dilepton tt |
| `ATLAS-CONF-2025-008.pdf` | ATLAS-CONF-2025-008: threshold excess in dilepton tt |

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
| 4 | `4_fuks_figs2to5.py` | Fuks Figs. 2-5: same events re-weighted with Fuks's G and ours. `--sample fuksproc_ct25`: Fuks's exact process (all e/mu combinations, CT25NNLO; `--generate` makes the samples with MadGraph); `--octet`: the figures with the octet in the signal (`fuksproc_ct25_octet/`); the octet study with the Green's functions in MATRIX1 writes `fuks_octet_mg/` (samples LO, A1/A18 = Fuks's tables singlet / + octet, B1/B18 = ours) | ~40 s (MadGraph: ~40 min) |
| 5 | `5_comparison_tables.py` | CSV tables (parameters, V vs Q^2, G on E-p grids; colour singlet) for comparing with an independent solver | ~1 min |
| 6 | `6_tfactors.py` | toponium K-factors at 13 TeV in the ATLAS/CMS binnings: LO MadGraph with \|G1/G0\|^2 on the gg colour singlet and \|G8/G0\|^2 on the octet (gg and qqbar) in MATRIX1, our G with CT25NNLO's alpha_s(Q). `uv run --with six python 6_tfactors.py <cmd> RUN`; runs below | ~4-5 h per run |

Each stage writes `results/<stage>/summary.txt` (+ figures).

### Stage 6 runs (`results/6_tfactors/RUN/`)

| Run | What | Commands |
|---|---|---|
| `run1_pp_fullcolour_13TeV`, `run2_gg_singlet_13TeV`, `run3_pp_singlet_octet_13TeV`, `run4_pp_singlet_octet_ll_13TeV` | earlier full-sample runs (two independent 10M samples, p* <= 50): pp singlet only; gg singlet only; pp singlet + octet; same with all e/mu combinations. Runs 1-3: results only | `setup`/`smoke`/`run`/`analyse` |
| `run5_window_kfactors_13TeV` | the window design (K = 1 + Delta/SM_full, Delta from samples generated only in W <= 350, p* <= 100; \|G/G0\|^2 for p* < 80; four octet potentials), with the seeds of M. Yanez Reyes's independent calculation, for the comparison | `table run5`, `setup`/`smoke`/`run`, `kfactors run5_window_kfactors_13TeV [CSV ...]` |
| **`run6_final_kfactors_13TeV`** | **final K-factors**: own SM_full (10M), q q~ and q~ q window samples, full-phase-space check. Plain CSVs per binning in `kfactors_final/` | `setup`/`smoke`/`run`/`kfactors`/`export`/`arbiter` |
| run 6, `octetcheck` | singlet only vs singlet + octet, one extra gg window sample + gg LO samples (`octet_check.txt`). First CMS m_tt bin [250, 400]: gg singlet events only (Fuks's sample) K = 1.229; gg all colours 1.061 (+ octet 1.055); pp 1.049; pp + octet 1.042 (final). The octet repulsion removes 13.5% of the singlet excess | `octetcheck run6_final_kfactors_13TeV` |
| `run7_mt{170p00,171p25,172p50,173p75,175p00}_13TeV` | m_t scan: run 6's design at each m_t, Gamma_t scaled at LO (`top_width`), its own G0(E, p; m_t, Gamma_t) and window 2 m_t - 20 <= W <= 2 m_t + 4; singlet full JKT (two-loop alpha_s, alpha_s(mZ) = 0.12), octet coulomb_as25. Run on the MSU HPCC (`submit_mtscan.sh`, `slurm_*.sh`). K vs m_t: `run7_mtscan_summary.txt`; per-mass CSVs with a README: `toponium_kfactors_mtscan_13TeV/` (+ `.zip`) | per run: `setup`/`smoke`/`one RUN NAME`/`kfactors`/`export`; then `mtscan` |

## Reproduction

In the repo: all code, the summaries, K-factor CSVs and plots of every run, the stage-4 event caches of
Fuks's figures (`results/4_fuks_figs2to5/_cache_*.npz`: per-event kinematics from MadGraph, the only copy of
those events), and the Green's-function tables of the early runs.

Not in the repo, regenerated by the code:
- The event caches of the two further stage-4 studies (MadGraph, ~40 min each, needs `uv run --with six`):
  `4_fuks_figs2to5.py --generate` for `fuksproc_ct25/` (seeds 7301, 7302) and `--octet-mg-generate` for
  `fuks_octet_mg/` (seeds 8401-8405). Their figures and summaries are in the repo.
- Smoke tests (`results/6_tfactors/*/smoke/`): rerun with `smoke RUN`.
- Green's-function tables on the p* <= 80 grid (`results/6_tfactors/*_p80.dat`, ~20 MB each, ~1 min each):
  `6_tfactors.py table run5` (singlet + four octets), `octetcheck` (all-ones octet table), and `setup` of each
  `run7_*` run (that mass's tables). Each run's `table_log.txt` records the checks against the solver.
- The per-event K-factor samples (`../output/<run>/events/*.npz`, up to ~2 GB per run): `setup`, `smoke`, `run` (or
  `one RUN NAME` per sample on the cluster) regenerate them with the seeds fixed in `RUNS` in `6_tfactors.py`;
  `kfactors` and `export` then rebuild the CSVs. Needs MadGraph (this checkout), `six`, and the local LHAPDF
  set `../lhapdf_data/` (CT25NNLO as id 99999000).

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
