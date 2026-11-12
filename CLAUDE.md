# Toponium study

A self-contained physics-analysis project layered on top of MadGraph, not part of the
framework: a from-scratch NRQCD Green's-function solver for toponium/threshold effects in
off-shell `t t~` (following Fuks et al., arXiv:2411.18962, solved as in JKT, Z. Phys. C 56 (1992) 653).
`README.md` has the layout and run order; `writeup.tex` is the methodology note for the group.
One library, `toponium.py` (six classes). The user wants
few files and easy navigation: don't split it into modules or packages again. Stages `0_solver_checks.py` .. `6_tfactors.py` run from this
folder and write `results/<stage>/`. Keep the equation-mapped style: EQUATION MAP headers,
docstrings starting with the equation implemented.
Python: the uv project in this folder. Run stages with
`uv run python <stage>.py`. ~/miniconda3 no longer exists (CONDA_PREFIX in the shell is stale).

The MadGraph event samples in `../output/` were deleted to free disk space. Stage 4
(`4_fuks_figs2to5.py`) now runs only from `results/4_fuks_figs2to5/_cache_*.npz`, which are the
ONLY remaining copy of those events: never delete or overwrite them. Regenerating the samples
follows the recipe in `README.md`.

Stage 6 (`6_tfactors.py`) drives MadGraph: run it as `uv run --with six python 6_tfactors.py <cmd>` (MadGraph needs `six`;
Homebrew python3 lacks it). Each run lives in `../output/<run>/` (its MadGraph dirs, `events/`, README) and `results/6_tfactors/<run>/`.
The final K-factors are run 6 (`run6_final_kfactors_13TeV`; plain per-binning CSVs in `results/6_tfactors/run6_final_kfactors_13TeV/kfactors_final/`).
Runs 1-3 exist only as results (`results/6_tfactors/run{1,2,3}_*`): their `../output/` dirs and events were deleted on 2026-11-04.
The MadGraph integration scratch (`SubProcesses/P*/G*/`) was deleted from runs 4-6; a rerun recreates it. It patches MATRIX1
in the WITH/topo dirs and uses a local LHAPDF dir `../lhapdf_data/` (CT25NNLO as id 99999000, `lhapdf-config` wrapper).
Per-event npz files are in `../output/<run>/events/` (outside the repo).
Run 7 (`run7_mt{170p00,171p25,172p50,173p75,175p00}_13TeV`) is the m_t scan: run 6's design per mass, Gamma_t scaled at LO
(`top_width`), full-JKT singlet (alpha_s(mZ) = 0.12, two-loop) and octet coulomb_as25 only; tables per run in
`results/6_tfactors/<run>/`. It runs on the MSU HPCC (`submit_mtscan.sh`, `slurm_*.sh`): only `/mnt/home/filippo3`
and `/mnt/scratch/filippo3` may be touched there; the user runs every cluster command.

Git: this folder is the public repo `AttackOnBreakfast/toponium-study`. `talk/` (PIKIMO 21 slides) is local only and
gitignored, as are the regenerable `results/6_tfactors/*_p80.dat` tables (README "Reproduction" says how to rebuild them).
