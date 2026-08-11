#!/bin/bash
#SBATCH --job-name=Deep_analysis
#SBATCH --output=slurm_logs/analysis_%A_%a.out
#SBATCH --error=slurm_logs/analysis_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10        
#SBATCH --mem=100G             
#SBATCH --time=02:00:00            
#SBATCH --array=1-6              
#SBATCH --chdir=/nfs/home/students/a.dersch/FoPra_PLAs

mkdir -p slurm_logs

COHORTEN=("heart" "sepsis" "vaccine" "immune_aging" "impact" "skin")
#COHORTEN=("vaccine")

ZERO_INDEX=$((SLURM_ARRAY_TASK_ID - 1))
COHORTE=${COHORTEN[$ZERO_INDEX]}

MODUS="raw"
GT_SOURCE="gmm_dual"  

CURRENT_FILE="pbmc_benchmarked_${COHORTE}_AUCell_gmm_dist_dual_${MODUS}_GT_${GT_SOURCE}.rds"

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


Rscript src/benchmarking/deep_analysis.R "$CURRENT_FILE" "$GT_SOURCE" "MANNE_DN"

echo "Task $SLURM_ARRAY_TASK_ID (Kohorte: $COHORTE) erfolgreich beendet."