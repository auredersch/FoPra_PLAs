#!/bin/bash

#SBATCH --job-name=dccs_batch
#SBATCH --array=0-35%4
#SBATCH --output=/nfs/home/students/i.kaciran/FoPra_PLAs/slurm_logs/dccs_%A_%a.out
#SBATCH --error=/nfs/home/students/i.kaciran/FoPra_PLAs/slurm_logs/dccs_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=4

set -euo pipefail

PROJECT_DIR="/nfs/home/students/i.kaciran/FoPra_PLAs"
SCRIPT_DIR="${PROJECT_DIR}/src/cell_signaling/dccs"

DATA_DIR="${PROJECT_DIR}/data/datasets"
PREP_DIR="${DATA_DIR}/prepared_inputs"

RESULTS_DIR="${PROJECT_DIR}/results/differential_ccs"
LOG_DIR="${PROJECT_DIR}/slurm_logs"

R_BIN="/nfs/home/students/i.kaciran/.conda/envs/liana_r/bin/Rscript"

mkdir -p "${RESULTS_DIR}" "${LOG_DIR}"

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

METHODS=(
  "multinichetr|${SCRIPT_DIR}/multinichetr.R"
  "scDiffCom|${SCRIPT_DIR}/scDiffCom.R"
)

N_DATASETS=${#DATASETS[@]}
N_MODES=${#MODES[@]}
N_METHODS=${#METHODS[@]}

TASK_ID=${SLURM_ARRAY_TASK_ID}

METHOD_INDEX=$((TASK_ID % N_METHODS))
MODE_INDEX=$(((TASK_ID / N_METHODS) % N_MODES))
DATASET_INDEX=$((TASK_ID / (N_METHODS * N_MODES)))

DATASET_BASE="${DATASETS[$DATASET_INDEX]}"
MODE="${MODES[$MODE_INDEX]}"
METHOD_ENTRY="${METHODS[$METHOD_INDEX]}"

METHOD="${METHOD_ENTRY%%|*}"
SCRIPT="${METHOD_ENTRY##*|}"

# Full datasets:
#   data/datasets/gated_heart_processed.rds
#
# Subset datasets:
#   data/datasets/prepared_inputs/gated_heart_processed_diseasedOnly.rds
#   data/datasets/prepared_inputs/gated_heart_processed_healthyOnly.rds
if [[ "${MODE}" == "all" ]]; then
  INPUT_RDS="${DATA_DIR}/${DATASET_BASE}.rds"
else
  INPUT_RDS="${PREP_DIR}/${DATASET_BASE}_${MODE}.rds"
fi

METHOD_OUT_DIR="${RESULTS_DIR}/${METHOD}"

mkdir -p "${METHOD_OUT_DIR}"

echo "============================================================"
echo "DCCS ARRAY TASK"
echo "============================================================"
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "Working directory: $(pwd)"
echo "SLURM job ID: ${SLURM_JOB_ID}"
echo "SLURM array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Dataset: ${DATASET_BASE}"
echo "Mode: ${MODE}"
echo "Method: ${METHOD}"
echo "Input RDS: ${INPUT_RDS}"
echo "Script: ${SCRIPT}"
echo "Output dir: ${METHOD_OUT_DIR}"
echo "Log dir: ${LOG_DIR}"
echo "============================================================"

if [[ ! -f "${INPUT_RDS}" ]]; then
  echo "ERROR: Missing input RDS: ${INPUT_RDS}"

  if [[ "${MODE}" == "all" ]]; then
    echo "Expected the full dataset under: ${DATA_DIR}"
  else
    echo "Expected the prepared cohort dataset under: ${PREP_DIR}"
    echo "Run prepare_cohort_inputs_slurm.sh if it has not been created."
  fi

  exit 1
fi

if [[ ! -f "${SCRIPT}" ]]; then
  echo "ERROR: Missing analysis script: ${SCRIPT}"
  exit 1
fi

if [[ ! -x "${R_BIN}" ]]; then
  echo "ERROR: R binary is missing or not executable: ${R_BIN}"
  exit 1
fi

echo "Checking R binary:"
ls -lh "${R_BIN}"
"${R_BIN}" --version

echo "Checking R script:"
ls -lh "${SCRIPT}"

# Useful for confirming that SLURM used the current script version
echo "Script checksum:"
md5sum "${SCRIPT}"

echo "Checking input RDS:"
ls -lh "${INPUT_RDS}"

echo "Checking output directory:"
ls -ld "${METHOD_OUT_DIR}"

echo "Starting R analysis command:"
echo "${R_BIN} ${SCRIPT} ${INPUT_RDS} ${METHOD_OUT_DIR}"
echo "============================================================"

set +e

"${R_BIN}" "${SCRIPT}" \
  "${INPUT_RDS}" \
  "${METHOD_OUT_DIR}"

EXIT_CODE=$?

set -e

echo "============================================================"
echo "R analysis exit code: ${EXIT_CODE}"
echo "Finished at: $(date)"
echo "============================================================"

if [[ "${EXIT_CODE}" -ne 0 ]]; then
  echo "ERROR: R analysis failed for ${METHOD} on ${DATASET_BASE} (${MODE})"
  exit "${EXIT_CODE}"
fi

echo "Finished ${METHOD} on ${DATASET_BASE} (${MODE})"