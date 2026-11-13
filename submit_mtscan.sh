#!/bin/bash
# Run 7 (m_t scan, five masses) on the MSU HPCC, all jobs at once with dependencies:
#   per mass: setup (the five set-ups one after another: they share lhapdf_data/) -> its 14 samples -> kfactors + export
#   at the end: the K vs m_t summary (results/6_tfactors/run7_mtscan_summary.txt)
# Usage: bash submit_mtscan.sh                   (from toponium_study/, on a dev node)
#        QOS=scavenger bash submit_mtscan.sh     the 70 sample jobs in the scavenger QOS (preemptible: a preempted
#        sample is requeued and restarts from scratch; finished samples are kept). Setup and kfactors jobs stay normal.
#        CONSTRAINT=amd24 ...                    all Python jobs only on nodes with that feature (the uv Python was built on
#                                                dev-amd24: other CPU types die with SIGILL, exit code 132)
#        SKIP_SETUP=1 ...                        no setup jobs (set-ups already done): samples start at once
# A failed sample leaves its later jobs pending (DependencyNeverSatisfied): resubmit that sample with slurm_one.sh,
# then the kfactors job by hand (finished samples are skipped, so resubmitting is cheap).
set -e
cd /mnt/home/filippo3/MG5_aMC_v2_9_27/toponium_study
module load LHAPDF/6.5.4-GCC-12.3.0
source ~/.local/bin/env
export UV_PYTHON_PREFERENCE=only-managed
PY="uv run --with six python -W ignore 6_tfactors.py"
SAMPLE_OPTS=""
[ -n "$QOS" ] && SAMPLE_OPTS="--account=$QOS --qos=$QOS --requeue"   # MSU: scavenger is an account with its own QOS
CONS=""
[ -n "$CONSTRAINT" ] && CONS="--constraint=$CONSTRAINT"
SAMPLE_OPTS="$SAMPLE_OPTS $CONS"
prev=""
finals=""
for RUN in $($PY mtscan_names); do
  if [ -n "$SKIP_SETUP" ]; then
    s=""; DEP=""
  else
    s=$(sbatch --parsable ${prev:+--dependency=afterok:$prev} $CONS slurm_setup.sh "$RUN"); DEP="--dependency=afterok:$s"
  fi
  ids=""
  for NAME in $($PY list "$RUN"); do
    j=$(sbatch --parsable $SAMPLE_OPTS $DEP slurm_one.sh "$RUN" "$NAME")
    ids="$ids:$j"
  done
  k=$(sbatch --parsable $CONS --dependency=afterok$ids slurm_kfactors.sh "$RUN")
  finals="$finals:$k"
  prev=$s
  echo "$RUN: setup ${s:-skipped}, 14 samples ${ids#:}, kfactors $k"
done
m=$(sbatch --parsable $CONS --dependency=afterok$finals slurm_kfactors.sh mtscan)
echo "summary $m"
