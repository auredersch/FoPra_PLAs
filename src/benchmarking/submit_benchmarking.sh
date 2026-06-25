#!/bin/bash
#SBATCH --job-name=PLA_Benchmarking
#SBATCH --output=slurm_logs/benchmark_%A_%a.out
#SBATCH --error=slurm_logs/benchmark_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=80G
#SBATCH --time=06:00:00
#SBATCH --array=1-5

mkdir -p slurm_logs

# --- DATASET CONFIGURATION ---
DATASETS=(
  "gated_heart_processed.rds"
  "gated_sepsis_processed.rds"
  "gated_vaccine_processed.rds"
  "gated_ImmuneAging.rds"
  "gated_our_dataset_processed.rds"
)

INDEX=$((SLURM_ARRAY_TASK_ID - 1))
CURRENT_FILE=${DATASETS[$INDEX]}

# Zwinge Slurm, das R 4.4.2 vom Cluster zu laden (das mit dem funktionierenden Seurat!)
export PATH="/nfs/data/cluster/software/R/4.4.2/lib/R/bin:$PATH"

# Wir zeigen R deinen neuen 4.4-Ordner. Den System-Ordner findet R 4.4.2 von alleine!
export R_LIBS_USER="/cmnfs/home/students/a.dersch/R/x86_64-pc-linux-gnu-library/4.4"

echo "======================================================"
echo "Starte Task $SLURM_ARRAY_TASK_ID für $CURRENT_FILE mit R 4.4.2"
echo "======================================================"

# --- PIPELINE STARTEN ---
Rscript Benchmarking_Master_v4.R \
  "AUCell" \
  "MANNE_DN" \
  "MANNE_COVID19_COMBINED_COHORT_VS_HEALTHY_DONOR_PLATELETS_DN.v2025.1.Hs" \
  TRUE \
  "gmm_dist_dual" \
  "$CURRENT_FILE"