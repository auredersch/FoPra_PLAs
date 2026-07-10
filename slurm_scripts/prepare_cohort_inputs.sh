#!/bin/bash

#SBATCH --job-name=prep_dccs_inputs
#SBATCH --array=0-9
#SBATCH --output=/nfs/home/students/i.kaciran/FoPra_PLAs/slurm_logs/prep_%A_%a.out
#SBATCH --error=/nfs/home/students/i.kaciran/FoPra_PLAs/slurm_logs/prep_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=2

set -euo pipefail

PROJECT_DIR="/nfs/home/students/i.kaciran/FoPra_PLAs"
SCRIPT_DIR="${PROJECT_DIR}/src/cell_signaling/dccs"
DATA_DIR="${PROJECT_DIR}/data/datasets"

RESULTS_DIR="${PROJECT_DIR}/results/differential_ccs"
PREP_DIR="${RESULTS_DIR}/prepared_inputs"
LOG_DIR="${RESULTS_DIR}/slurm_logs"

R_BIN="/nfs/home/students/i.kaciran/.conda/envs/liana_r/bin/Rscript"
PREP_SCRIPT="${SCRIPT_DIR}/prepare_cohort_rds.R"

mkdir -p "${RESULTS_DIR}" "${PREP_DIR}" "${LOG_DIR}"

DATASETS=(
  "${DATA_DIR}/gated_heart_processed.rds|heart"
  "${DATA_DIR}/gated_ImmuneAging.rds|immune_aging"
  "${DATA_DIR}/gated_sepsis_processed.rds|sepsis"
  "${DATA_DIR}/gated_vaccine_processed.rds|vaccine"
  "${DATA_DIR}/gated_our_dataset_processed.rds|our_data"
)

MODES=(
  "withHealthy"
  "noHealthy"
)

N_DATASETS=${#DATASETS[@]}
N_MODES=${#MODES[@]}
TOTAL_TASKS=$((N_DATASETS * N_MODES))

TASK_ID=${SLURM_ARRAY_TASK_ID}

if (( TASK_ID >= TOTAL_TASKS )); then
  echo "Task ID ${TASK_ID} exceeds total tasks ${TOTAL_TASKS}. Exiting."
  exit 0
fi

MODE_INDEX=$((TASK_ID % N_MODES))
DATASET_INDEX=$((TASK_ID / N_MODES))

DATASET_ENTRY="${DATASETS[$DATASET_INDEX]}"
MODE="${MODES[$MODE_INDEX]}"

INPUT_RDS="${DATASET_ENTRY%%|*}"
DATASET_TYPE="${DATASET_ENTRY##*|}"
DATASET_BASE="$(basename "${INPUT_RDS}" .rds)"

OUTPUT_RDS="${PREP_DIR}/${DATASET_BASE}_${MODE}.rds"

echo "SLURM job ID: ${SLURM_JOB_ID}"
echo "SLURM array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Dataset: ${DATASET_BASE}"
echo "Dataset type: ${DATASET_TYPE}"
echo "Mode: ${MODE}"
echo "Input RDS: ${INPUT_RDS}"
echo "Output RDS: ${OUTPUT_RDS}"

if [[ ! -f "${INPUT_RDS}" ]]; then
  echo "Missing dataset: ${INPUT_RDS}"
  exit 1
fi

if [[ ! -f "${PREP_SCRIPT}" ]]; then
  echo "Missing preparation script: ${PREP_SCRIPT}"
  exit 1
fi

"${R_BIN}" "${PREP_SCRIPT}" \
  "${INPUT_RDS}" \
  "${OUTPUT_RDS}" \
  "${DATASET_TYPE}" \
  "${MODE}"

echo "Finished preparing: ${OUTPUT_RDS}"