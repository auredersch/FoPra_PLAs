#!/bin/bash
#SBATCH --job-name=compare_qc_predictions
#SBATCH --output=/nfs/home/students/f.mathis/FoPra_PLAs/logs/compare_qc_predictions_%A_%a.out
#SBATCH --error=/nfs/home/students/f.mathis/FoPra_PLAs/logs/compare_qc_predictions_%A_%a.err
#SBATCH --chdir=/nfs/home/students/f.mathis/FoPra_PLAs
#SBATCH --time=00:30:00
#SBATCH --mem=128G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-40%8

set -euo pipefail
source ~/.bashrc
conda activate FoPra

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# One task per row of combined/run_config.csv (currently 40 runs).
# Optional first argument: a different QC grid directory.
INPUT_DIR="${1:-results/benchmarking/qc_grid_no_na_full}"
Rscript src/benchmarking/search_qc_params/compare_qc_predictions.R "$INPUT_DIR" "$SLURM_ARRAY_TASK_ID"
