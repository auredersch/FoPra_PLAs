#!/bin/bash

#SBATCH --job-name=fix_pla_status
#SBATCH --output=/nfs/home/students/i.kaciran/FoPra_PLAs/slurm_logs/fix_pla_status_%j.out
#SBATCH --error=/nfs/home/students/i.kaciran/FoPra_PLAs/slurm_logs/fix_pla_status_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1

set -euo pipefail

PROJECT_DIR="/nfs/home/students/i.kaciran/FoPra_PLAs"
SCRIPT="${PROJECT_DIR}/src/dataset/adapt_dataset.R"
LOG_DIR="${PROJECT_DIR}/slurm_logs"

R_BIN="/nfs/home/students/i.kaciran/.conda/envs/liana_r/bin/Rscript"

echo "============================================================"
echo "PLA STATUS UPDATE"
echo "============================================================"
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "Working directory: $(pwd)"
echo "SLURM job ID: ${SLURM_JOB_ID:-unknown}"
echo "R binary: ${R_BIN}"
echo "R script: ${SCRIPT}"
echo "Log directory: ${LOG_DIR}"
echo "============================================================"

if [[ ! -x "${R_BIN}" ]]; then
  echo "ERROR: Rscript is missing or not executable:"
  echo "${R_BIN}"
  exit 1
fi

if [[ ! -f "${SCRIPT}" ]]; then
  echo "ERROR: R script does not exist:"
  echo "${SCRIPT}"
  exit 1
fi

echo "R version:"
"${R_BIN}" --version

echo "R script:"
ls -lh "${SCRIPT}"

echo "Starting pla_status update..."
echo "============================================================"

"${R_BIN}" "${SCRIPT}"

EXIT_CODE=$?

echo "============================================================"
echo "Exit code: ${EXIT_CODE}"
echo "Finished at: $(date)"
echo "============================================================"

if (( EXIT_CODE != 0 )); then
  echo "ERROR: pla_status update failed."
  exit "${EXIT_CODE}"
fi

echo "pla_status update completed successfully."