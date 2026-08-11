#!/bin/bash
#SBATCH --job-name=PLA_Benchmarking
#SBATCH --output=slurm_logs/benchmark_%A_%a.out
#SBATCH --error=slurm_logs/benchmark_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=120G
#SBATCH --time=03:00:00
#SBATCH --array=1-66            # 5 Kohorten x 12 Signaturen = 60 Jobs
#SBATCH --chdir=/nfs/home/students/a.dersch/FoPra_PLAs

mkdir -p slurm_logs

COHORTEN=("heart" "sepsis" "vaccine" "immune_aging" "impact" "skin")
#COHORTEN=("skin")
N_COHORTEN=${#COHORTEN[@]}   # = 6

SIGNATUREN=(
  "MANNE_DN:MANNE_COVID19_COMBINED_COHORT_VS_HEALTHY_DONOR_PLATELETS_DN.v2025.1.Hs"
  "MANNE_UP:MANNE_COVID19_COMBINED_COHORT_VS_HEALTHY_DONOR_PLATELETS_UP.v2025.1.Hs"
  "GOBP:GOBP_REGULATION_OF_PLATELET_ACTIVATION.v2025.1.Hs"
  "REACTOME:REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION.v2025.1.Hs"
  "GNATENKO:GNATENKO_PLATELET_SIGNATURE.v2025.1.Hs"
  "WP:WP_PLATELETMEDIATED_INTERACTIONS_WITH_VASCULAR_AND_CIRCULATING_CELLS.v2025.1.Hs"
  "HP:HP_ABNORMAL_PLATELET_MEMBRANE_PROTEIN_EXPRESSION.v2025.1.Hs"
  "OVERLAP:OVERLAP_GOBP_MANNEDN_REACTOME"
  "UNION_ALL:UNION_ALL"
  "OVERLAP_MIN2:OVERLAP_MIN2"
  "HUMAN:human_immune_markers_unique"
)
#SIGNATUREN=("HUMAN:human_immune_markers_unique" "LEUKOCYTE:leukocyte_activation")
N_SIGS=${#SIGNATUREN[@]}     # = 11

ZERO_INDEX=$((SLURM_ARRAY_TASK_ID - 1))
COHORTE_IDX=$((ZERO_INDEX % N_COHORTEN))
SIG_IDX=$((ZERO_INDEX / N_COHORTEN))

COHORTE=${COHORTEN[$COHORTE_IDX]}
SIG_ENTRY=${SIGNATUREN[$SIG_IDX]}

SIG_NAME="${SIG_ENTRY%%:*}"
SIG_FILE="${SIG_ENTRY##*:}"

MODUS="raw"
GT_SOURCE="biologist"
METHOD="AUCell"
USE_EXT="FALSE"
THRESH="gmm_dist_dual"

CURRENT_FILE="${COHORTE}_${MODUS}_automated_gating.rds"

export PATH="/nfs/data/cluster/software/R/4.4.2/lib/R/bin:$PATH"
export R_LIBS_USER="/cmnfs/home/students/a.dersch/R/x86_64-pc-linux-gnu-library/4.4"

echo "======================================================"
echo "Array Task ID:   $SLURM_ARRAY_TASK_ID"
echo "Kohorte:         $COHORTE  (Index $COHORTE_IDX)"
echo "Signatur:        $SIG_NAME  (Index $SIG_IDX)"
echo "Datei:           $CURRENT_FILE"
echo "GT-Quelle:       $GT_SOURCE"
echo "Modus:           $MODUS"
echo "======================================================"

Rscript src/benchmarking/Benchmarking_Master_v4.R \
  "$METHOD"       \
  "$SIG_NAME"     \
  "$SIG_FILE"     \
  "$USE_EXT"      \
  "$THRESH"       \
  "$CURRENT_FILE" \
  "$GT_SOURCE"

EXIT_CODE=$?
if [ $EXIT_CODE -ne 0 ]; then
  echo "FEHLER: Task $SLURM_ARRAY_TASK_ID (${COHORTE} / ${SIG_NAME}) ist mit Exit-Code $EXIT_CODE fehlgeschlagen!"
  exit $EXIT_CODE
fi

echo "Task $SLURM_ARRAY_TASK_ID (${COHORTE} / ${SIG_NAME}) erfolgreich beendet."