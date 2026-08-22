#!/bin/bash

#SBATCH --job-name=multiniche_preflight
#SBATCH --array=0-17%4
#SBATCH --output=/nfs/home/students/i.kaciran/FoPra_PLAs/slurm_logs/multiniche_preflight_%A_%a.out
#SBATCH --error=/nfs/home/students/i.kaciran/FoPra_PLAs/slurm_logs/multiniche_preflight_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2

set -euo pipefail

PROJECT_DIR="/nfs/home/students/i.kaciran/FoPra_PLAs"

SCRIPT_DIR="${PROJECT_DIR}/src/cell_signaling/dccs"
DATA_DIR="${PROJECT_DIR}/data/datasets"
PREP_DIR="${DATA_DIR}/prepared_inputs"
RESULTS_DIR="${PROJECT_DIR}/results/differential_ccs"
PREFLIGHT_OUT_DIR="${RESULTS_DIR}/multinichetr_preflight"
LOG_DIR="${PROJECT_DIR}/slurm_logs"

R_BIN="/nfs/home/students/i.kaciran/.conda/envs/liana_r/bin/Rscript"
SCRIPT="${SCRIPT_DIR}/multinichenet_preflight.R"

mkdir -p \
  "${PREFLIGHT_OUT_DIR}" \
  "${LOG_DIR}"

DATASETS=(
  "gated_heart_processed"
  "gated_ImmuneAging"
  "gated_sepsis_processed"
  "gated_vaccine_processed"
  "gated_our_dataset_processed"
  "gated_skin_processed"
)

MODES=(
  "all"
  "diseasedOnly"
  "healthyOnly"
)

N_MODES=${#MODES[@]}
TASK_ID=${SLURM_ARRAY_TASK_ID}

MODE_INDEX=$((TASK_ID % N_MODES))
DATASET_INDEX=$((TASK_ID / N_MODES))

DATASET_BASE="${DATASETS[$DATASET_INDEX]}"
MODE="${MODES[$MODE_INDEX]}"

if [[ "${MODE}" == "all" ]]; then
  DATASET_NAME="${DATASET_BASE}"
  INPUT_RDS="${DATA_DIR}/${DATASET_BASE}.rds"
else
  DATASET_NAME="${DATASET_BASE}_${MODE}"
  INPUT_RDS="${PREP_DIR}/${DATASET_NAME}.rds"
fi

DATASET_OUT_DIR="${PREFLIGHT_OUT_DIR}/${DATASET_NAME}"

mkdir -p "${DATASET_OUT_DIR}"

echo "============================================================"
echo "MULTINICHENET PREFLIGHT ARRAY TASK"
echo "============================================================"
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "SLURM job ID: ${SLURM_JOB_ID}"
echo "SLURM array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Dataset: ${DATASET_BASE}"
echo "Mode: ${MODE}"
echo "Input RDS: ${INPUT_RDS}"
echo "Output directory: ${DATASET_OUT_DIR}"
echo "============================================================"

if [[ ! -f "${INPUT_RDS}" ]]; then
  echo "ERROR: Missing RDS: ${INPUT_RDS}"
  exit 1
fi

if [[ ! -f "${SCRIPT}" ]]; then
  echo "ERROR: Missing preflight script: ${SCRIPT}"
  exit 1
fi

if [[ ! -x "${R_BIN}" ]]; then
  echo "ERROR: R binary is missing or not executable: ${R_BIN}"
  exit 1
fi

"${R_BIN}" "${SCRIPT}" \
  "${INPUT_RDS}" \
  "${DATASET_OUT_DIR}"

echo "Finished preflight for ${DATASET_NAME}"