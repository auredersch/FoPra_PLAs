#!/bin/bash
#SBATCH --job-name=PLA_Benchmarking
#SBATCH --output=slurm_logs/benchmark_%A_%a.out
#SBATCH --error=slurm_logs/benchmark_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20        
#SBATCH --mem=120G             
#SBATCH --time=02:00:00            
#SBATCH --array=1-5              
#SBATCH --chdir=/nfs/home/students/a.dersch/FoPra_PLAs

mkdir -p slurm_logs

COHORTEN=("heart" "sepsis" "vaccine" "immune_aging" "impact")

ZERO_INDEX=$((SLURM_ARRAY_TASK_ID - 1))
COHORTE=${COHORTEN[$ZERO_INDEX]}

MODUS="qc_tolerant"
GT_SOURCE="gmm_dual"  

CURRENT_FILE="${COHORTE}_${MODUS}_automated_gating.rds"

# Environment Setup
export PATH="/nfs/data/cluster/software/R/4.4.2/lib/R/bin:$PATH"
export R_LIBS_USER="/cmnfs/home/students/a.dersch/R/x86_64-pc-linux-gnu-library/4.4"

echo "======================================================"
echo "TEST LAUF - Slurm Task ID: $SLURM_ARRAY_TASK_ID"
echo "Verarbeite Kohorte:        $COHORTE"
echo "Nutze Datei:               $CURRENT_FILE"
echo "Modus (Filter):            $MODUS"
echo "Ground Truth Quelle:       $GT_SOURCE"
echo "======================================================"

# --- PIPELINE AUSFÜHRUNG ---
Rscript src/benchmarking/Benchmarking_Master_v4.R \
  "AUCell" \
  "MANNE_DN" \
  "MANNE_COVID19_COMBINED_COHORT_VS_HEALTHY_DONOR_PLATELETS_DN.v2025.1.Hs" \
  FALSE \
  "gmm_dist_dual" \
  "$CURRENT_FILE" \
  "$GT_SOURCE"

echo "Task $SLURM_ARRAY_TASK_ID (Kohorte: $COHORTE) erfolgreich beendet."