#!/bin/bash
# One MadGraph sample of a Stage 6 run on the MSU HPCC:  sbatch slurm_one.sh RUN SAMPLE
# (sample names: uv run --with six python 6_tfactors.py list RUN). The job runs MadGraph in a private copy of the
# output directory, so jobs of the same run can run at the same time. Resubmitting a finished sample does nothing.
#SBATCH --job-name=topo-one
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=/mnt/scratch/filippo3/output/slurm-%x-%j.log

module load LHAPDF/6.5.4-GCC-12.3.0
source ~/.local/bin/env
export UV_PYTHON_PREFERENCE=only-managed   # uv's portable Python, not the module one (built for one CPU type)
cd /mnt/home/filippo3/MG5_aMC_v2_9_27/toponium_study
uv run --with six python -W ignore 6_tfactors.py one "$1" "$2"
