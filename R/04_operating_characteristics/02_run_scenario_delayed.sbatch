#!/bin/bash
#
# Operating characteristics, scenario S2 (delayed effect): 100,000 replicates.
# Paper scale: 50 array tasks x NSIMS=2000 per task.
# Submit from this folder:  sbatch 02_run_scenario_delayed.sbatch
# For a quick local test before submitting at full scale, reduce NSIMS to ~5 and array size to 1-2 tasks.
#
#SBATCH --job-name=table4_S2
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

NSIMS=2000
SEED=$SLURM_ARRAY_TASK_ID

Rscript 02_run_scenario_delayed.R $NSIMS $SEED
