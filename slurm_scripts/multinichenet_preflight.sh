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
RESULTS_DIR="${PROJECT_DIR}/results/differential_ccs"
PREP_DIR="${RESULTS_DIR}/prepared_inputs"
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

N_DATASETS=${#DATASETS[@]}
N_MODES=${#MODES[@]}

TOTAL_TASKS=$((N_DATASETS * N_MODES))
TASK_ID=${SLURM_ARRAY_TASK_ID}

MODE_INDEX=$((TASK_ID % N_MODES))
DATASET_INDEX=$((TASK_ID / N_MODES))

DATASET_BASE="${DATASETS[$DATASET_INDEX]}"
MODE="${MODES[$MODE_INDEX]}"
DATASET_NAME="${DATASET_BASE}_${MODE}"

PREPARED_RDS="${PREP_DIR}/${DATASET_NAME}.rds"
DATASET_OUT_DIR="${PREFLIGHT_OUT_DIR}/${DATASET_NAME}"

mkdir -p "${DATASET_OUT_DIR}"

echo "============================================================"
echo "MULTINICHENET PREFLIGHT ARRAY TASK"
echo "============================================================"
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "Working directory: $(pwd)"
echo "SLURM job ID: ${SLURM_JOB_ID}"
echo "SLURM array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Dataset index: ${DATASET_INDEX}"
echo "Mode index: ${MODE_INDEX}"
echo "Dataset: ${DATASET_BASE}"
echo "Mode: ${MODE}"
echo "Dataset name: ${DATASET_NAME}"
echo "Prepared RDS: ${PREPARED_RDS}"
echo "Script: ${SCRIPT}"
echo "Output directory: ${DATASET_OUT_DIR}"
echo "============================================================"

if [[ ! -f "${PREPARED_RDS}" ]]; then
  echo "ERROR: Missing prepared RDS: ${PREPARED_RDS}"
  echo "Run prepare_cohort_inputs_slurm.sh first."
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

echo "Checking R binary:"
ls -lh "${R_BIN}"
"${R_BIN}" --version

echo "Checking preflight script:"
ls -lh "${SCRIPT}"

echo "Checking prepared RDS:"
ls -lh "${PREPARED_RDS}"

echo "Starting preflight command:"
echo "${R_BIN} ${SCRIPT} ${PREPARED_RDS} ${DATASET_OUT_DIR}"
echo "============================================================"

set +e

"${R_BIN}" "${SCRIPT}" \
  "${PREPARED_RDS}" \
  "${DATASET_OUT_DIR}"

EXIT_CODE=$?

set -e

echo "============================================================"
echo "R preflight exit code: ${EXIT_CODE}"
echo "Finished at: $(date)"
echo "============================================================"

if [[ "${EXIT_CODE}" -ne 0 ]]; then
  echo "ERROR: Preflight failed for ${DATASET_NAME}"
  exit "${EXIT_CODE}"
fi

echo "Finished preflight for ${DATASET_NAME}"