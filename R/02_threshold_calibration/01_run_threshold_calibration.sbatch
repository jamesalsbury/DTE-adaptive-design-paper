#!/bin/bash
#
# PP threshold calibration (Table 3, kappa*): 15,000 replicates per scenario.
# Paper scale: 50 array tasks x NSIMS=300 per task.
# Submit from this folder:  sbatch 01_run_threshold_calibration.sbatch
# For a quick local test before submitting at full scale, reduce NSIMS to ~5 and array size to 1-2 tasks.
#
#SBATCH --job-name=threshold_calib
#SBATCH --array=1-50
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

module load rjags/4-13-foss-2022a-R-4.2.1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export TMPDIR=${TMPDIR:-/tmp/$USER/$SLURM_JOB_ID}
mkdir -p "$TMPDIR"

NSIMS=300
SEED=$SLURM_ARRAY_TASK_ID

Rscript 01_run_threshold_calibration.R $NSIMS $SEED
