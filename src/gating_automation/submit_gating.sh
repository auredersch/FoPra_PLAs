#!/bin/bash
#SBATCH --job-name=PLA_Auto_Gating
#SBATCH --output=slurm_logs/gating_%A_%a.out
#SBATCH --error=slurm_logs/gating_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10         
#SBATCH --mem=60G                  
#SBATCH --time=03:00:00            
#SBATCH --array=1-15               # 5 Datensätze * 3 Modi = 15 parallele Tasks
#SBATCH --chdir=/nfs/home/students/a.dersch/FoPra_PLAs # Fixiert den Startordner

mkdir -p slurm_logs

# --- DATENSÄTZE DEFINIEREN ---
DATASETS=(
  "gated_heart_processed.rds"
  "gated_sepsis_processed.rds"
  "gated_vaccine_processed.rds"
  "gated_ImmuneAging.rds"
  "gated_our_dataset_processed.rds"
)

ZERO_INDEX=$((SLURM_ARRAY_TASK_ID - 1))

DATASET_INDEX=$((ZERO_INDEX / 3))
CURRENT_FILE=${DATASETS[$DATASET_INDEX]}

MODUS_INDEX=$((ZERO_INDEX % 3))

if [ $MODUS_INDEX -eq 0 ]; then
  CURRENT_MODE="raw"
elif [ $MODUS_INDEX -eq 1 ]; then
  CURRENT_MODE="qc_tolerant"
else
  CURRENT_MODE="qc_strict"
fi

export PATH="/nfs/data/cluster/software/R/4.4.2/lib/R/bin:$PATH"
export R_LIBS_USER="/cmnfs/home/students/a.dersch/R/x86_64-pc-linux-gnu-library/4.4"

echo "======================================================"
echo "Slurm Task ID:   $SLURM_ARRAY_TASK_ID"
echo "Verarbeite Datei: $CURRENT_FILE"
echo "Verwende Modus:  $CURRENT_MODE"
echo "PROJEKT-PFAD:    $(pwd)"
echo "======================================================"

Rscript src/gating_automation/03_PLA_Automated_Gating.R "$CURRENT_FILE" "$CURRENT_MODE"

echo "Task $SLURM_ARRAY_TASK_ID erfolgreich beendet."