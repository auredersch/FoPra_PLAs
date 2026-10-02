#!/bin/bash
#SBATCH --job-name=frequency_qc_grid
#SBATCH --output=/nfs/home/students/f.mathis/FoPra_PLAs/logs/frequency_qc_grid_%A_%a.out
#SBATCH --error=/nfs/home/students/f.mathis/FoPra_PLAs/logs/frequency_qc_grid_%A_%a.err
#SBATCH --chdir=/nfs/home/students/f.mathis/FoPra_PLAs
#SBATCH --time=02:00:00
#SBATCH --mem=128G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=0-4

set -euo pipefail
source ~/.bashrc
conda activate FoPra

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# Usage: sbatch run_frequency_qc_grid.sh [exclude_na_metadata] [adtnorm_input_dir] [output_root] [grid_csv]
# One task per dataset: 0=heart, 1=vaccine, 2=immune_aging, 3=sepsis, 4=skin.
# Defaults: retain NA/Unassigned, full tuned ADTnorm inputs, standard frequency grid.
if (( $# > 4 )); then
  echo "Expected: [TRUE|FALSE] [adtnorm_input_dir] [output_root] [grid_csv]" >&2
  exit 2
fi

EXCLUDE_NA="${1:-FALSE}"
EXCLUDE_NA="${EXCLUDE_NA^^}"
case "${EXCLUDE_NA}" in
  TRUE) DEFAULT_OUTPUT_ROOT="results/benchmarking/frequency_qc_full_no_na" ;;
  FALSE) DEFAULT_OUTPUT_ROOT="results/benchmarking/frequency_qc_full" ;;
  *) echo "exclude_na_metadata must be TRUE or FALSE." >&2; exit 2 ;;
esac

case "${SLURM_ARRAY_TASK_ID:-}" in
  0|1|2|3|4) ;;
  *) echo "Submit as an array job with task IDs 0-4." >&2; exit 2 ;;
esac

DATASETS=(heart vaccine immune_aging sepsis skin)
DATASET="${DATASETS[$SLURM_ARRAY_TASK_ID]}"
INPUT_DIR="${2:-}"
OUTPUT_ROOT="${3:-${DEFAULT_OUTPUT_ROOT}}"
OUTPUT_DIR="${OUTPUT_ROOT}/by_dataset/${DATASET}"
GRID_CSV="${4:-}"

export PLA_PROJECT_ROOT="${PWD}"
echo "Dataset: ${DATASET}; exclude_na_metadata: ${EXCLUDE_NA}"
echo "Output: ${OUTPUT_DIR}"
Rscript src/benchmarking/search_qc_params/Benchmarking_Frequency_QC_Grid.R \
  "${INPUT_DIR}" "${OUTPUT_DIR}" "${GRID_CSV}" "${DATASET}" "${EXCLUDE_NA}"
