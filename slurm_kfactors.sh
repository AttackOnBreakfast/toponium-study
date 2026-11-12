#!/bin/bash
# K-factors of one Stage 6 run on the MSU HPCC (kfactors CSVs, plots, kfactors_final/):  sbatch slurm_kfactors.sh RUN
# or, after all five run 7 mass points:  sbatch slurm_kfactors.sh mtscan   (K vs m_t summary)
#SBATCH --job-name=topo-kfac
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/mnt/scratch/filippo3/output/slurm-%x-%j.log

module load LHAPDF/6.5.4-GCC-12.3.0
source ~/.local/bin/env
export UV_PYTHON_PREFERENCE=only-managed   # uv's portable Python, not the module one (built for one CPU type)
cd /mnt/home/filippo3/MG5_aMC_v2_9_27/toponium_study
if [ "$1" = "mtscan" ]; then
  uv run --with six python -W ignore 6_tfactors.py mtscan
else
  uv run --with six python -W ignore 6_tfactors.py kfactors "$1" && uv run --with six python -W ignore 6_tfactors.py export "$1"
fi
