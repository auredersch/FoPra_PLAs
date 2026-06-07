#!/bin/bash

#SBATCH --job-name=liana
#SBATCH --array=0-4
#SBATCH --output=slurm_logs/liana_%A_%a.out
#SBATCH --error=slurm_logs/liana_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=4

set -euo pipefail

cd "$SLURM_SUBMIT_DIR"

DATA_DIR="/nfs/home/students/i.kaciran/FoPra_PLAs/data/datasets"
OUTPUT_DIR="/nfs/home/students/i.kaciran/FoPra_PLAs/results/cell_signaling"
SCRIPT="/nfs/home/students/i.kaciran/FoPra_PLAs/src/cell_signaling/run_liana.R"

FILES=(
  "gated_ImmuneAging.rds"
  "gated_heart_processed.rds"
  "gated_sepsis_processed.rds"
  "gated_vaccine_processed.rds"
  "gated_our_dataset_processed.rds"
)

DATASET_NAMES=(
  "immune_aging"
  "heart"
  "sepsis"
  "vaccine"
  "our_data"
)

TASK_ID="$SLURM_ARRAY_TASK_ID"

FILE_NAME="${FILES[$TASK_ID]}"
DATASET_NAME="${DATASET_NAMES[$TASK_ID]}"
INPUT_FILE="${DATA_DIR}/${FILE_NAME}"

mkdir -p "$OUTPUT_DIR"

echo "Job ID: $SLURM_JOB_ID"
echo "Array job ID: $SLURM_ARRAY_JOB_ID"
echo "Array task ID: $TASK_ID"
echo "Dataset: $DATASET_NAME"
echo "Input file: $INPUT_FILE"
echo "Node: $(hostname)"
echo "Start time: $(date)"

if [[ ! -f "$INPUT_FILE" ]]; then
  echo "ERROR: Input file does not exist: $INPUT_FILE" >&2
  exit 1
fi

conda run --no-capture-output -n liana_r \
  Rscript "$SCRIPT" \
  "$INPUT_FILE" \
  "$DATASET_NAME" \
  "$OUTPUT_DIR"

echo "Finished dataset: $DATASET_NAME"
echo "End time: $(date)"