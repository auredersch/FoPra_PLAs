#!/bin/bash

#SBATCH --job-name=analyze_liana
#SBATCH --array=0-4
#SBATCH --output=slurm_logs/analyze_liana_%A_%a.out
#SBATCH --error=slurm_logs/analyze_liana_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2

set -euo pipefail

cd "$SLURM_SUBMIT_DIR"

INPUT_DIR="/nfs/home/students/i.kaciran/FoPra_PLAs/results/cell_signaling"
OUTPUT_DIR="/nfs/home/students/i.kaciran/FoPra_PLAs/results/cell_signaling/analysis"
SCRIPT="/nfs/home/students/i.kaciran/FoPra_PLAs/src/cell_signaling/analyze_liana.R"
RSCRIPT="/nfs/home/students/i.kaciran/.conda/envs/liana_r/bin/Rscript"

FILES=(
  "immune_aging_liana_results_pla_status.rds"
  "heart_liana_results_pla_status.rds"
  "sepsis_liana_results_pla_status.rds"
  "vaccine_liana_results_pla_status.rds"
  "our_data_liana_results_pla_status.rds"
)

TASK_ID="$SLURM_ARRAY_TASK_ID"
INPUT_FILE="${INPUT_DIR}/${FILES[$TASK_ID]}"

mkdir -p "$OUTPUT_DIR"

echo "Job ID: $SLURM_JOB_ID"
echo "Array job ID: $SLURM_ARRAY_JOB_ID"
echo "Array task ID: $TASK_ID"
echo "Rscript: $RSCRIPT"
echo "Input file: $INPUT_FILE"
echo "Output directory: $OUTPUT_DIR"
echo "Node: $(hostname)"
echo "Start time: $(date)"

if [[ ! -x "$RSCRIPT" ]]; then
  echo "ERROR: Rscript does not exist or is not executable: $RSCRIPT" >&2
  exit 1
fi

if [[ ! -f "$SCRIPT" ]]; then
  echo "ERROR: R script does not exist: $SCRIPT" >&2
  exit 1
fi

if [[ ! -f "$INPUT_FILE" ]]; then
  echo "ERROR: Input file does not exist: $INPUT_FILE" >&2
  exit 1
fi

"$RSCRIPT" "$SCRIPT" \
  "$INPUT_FILE" \
  "$OUTPUT_DIR"

echo "Finished analysis."
echo "End time: $(date)"