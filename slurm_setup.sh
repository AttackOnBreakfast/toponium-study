#!/bin/bash
# Set-up of one Stage 6 run on the MSU HPCC (G tables, MadGraph outputs, patches, cards):  sbatch slurm_setup.sh RUN
#SBATCH --job-name=topo-setup
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/mnt/scratch/filippo3/output/slurm-%x-%j.log

module load LHAPDF/6.5.4-GCC-12.3.0
source ~/.local/bin/env
export UV_PYTHON_PREFERENCE=only-managed   # uv's portable Python, not the module one (built for one CPU type)
cd /mnt/home/filippo3/MG5_aMC_v2_9_27/toponium_study
uv run --with six python -W ignore 6_tfactors.py setup "$1"
