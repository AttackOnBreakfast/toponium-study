# Toponium study

A self-contained physics-analysis project layered on top of MadGraph, not part of the
framework: a from-scratch NRQCD Green's-function solver for toponium/threshold effects in
off-shell `t t~` (following Fuks et al., arXiv:2411.18962, solved as in JKT, Z. Phys. C 56 (1992) 653).
`README.md` has the layout and run order; `writeup.tex` is the methodology note for the group.
One library, `toponium.py`, with six classes: `Constants`, `JKTPotential` (JKT Sec. 2 +
cut (19)-(20)), `LSSolver` (JKT Sec. 3), `CoulombExact` (exact fixed-alpha Coulomb), `FuksTables`
(Fuks's swData + CALCGREEN, as-is), `Events` (events, |G/G0|^2 re-weighting). The user wants
few files and easy navigation: don't split it into modules or packages again. Stages `0_solver_checks.py` .. `6_tfactors.py` run from this
folder and write `results/<stage>/`. Keep the equation-mapped style: EQUATION MAP headers,
docstrings starting with the equation implemented.
Python: the uv project in this folder (Python 3.11, numpy/scipy/matplotlib; `.venv/`). Run stages with
`uv run python <stage>.py`. ~/miniconda3 no longer exists (CONDA_PREFIX in the shell is stale).

The MadGraph event samples in `../output/` were deleted to free disk space. Stage 4
(`4_fuks_figs2to5.py`) now runs only from `results/4_fuks_figs2to5/_cache_*.npz`, which are the
ONLY remaining copy of those events: never delete or overwrite them. Regenerating the samples
follows the recipe in `README.md`.

Stage 6 (`6_tfactors.py`) drives MadGraph: run it as `uv run --with six python 6_tfactors.py <cmd>` (MadGraph needs `six`;
Homebrew python3 lacks it). Each run lives in `../output/<run>/` (its two MadGraph dirs, `events/`, README) and `results/6_tfactors/<run>/`;
`run1_pp_fullcolour_13TeV` is archived (do not re-run it in place). It patches the gg MATRIX1 of each WITH dir and uses a
local LHAPDF dir `../lhapdf_data/` (CT25NNLO as id 99999000, `lhapdf-config` wrapper). Per-event npz files are in
`../output/<run>/events/` (outside the repo).
