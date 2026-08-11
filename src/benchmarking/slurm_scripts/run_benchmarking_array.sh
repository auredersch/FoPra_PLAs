#!/bin/bash
#SBATCH --job-name=PLA_Benchmarking
#SBATCH --output=slurm_logs/benchmark_%A_%a.out
#SBATCH --error=slurm_logs/benchmark_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=120G
#SBATCH --time=03:00:00
#SBATCH --array=1-66             # Fallback-Groesse fuer Direktaufruf (siehe unten) - 6 Kohorten x 11 Signaturen x 1 Filter_Mode x 1 GT_Source
#SBATCH --chdir=/nfs/home/students/a.dersch/FoPra_PLAs

# -------------------------------------------------------------------
# WICHTIG: Filter_Mode(s), GT_Source(s) und Signatur(en) werden ueber die
# Umgebungsvariablen FILTER_MODES_LIST, GT_SOURCES_LIST und SIGNATURES_LIST
# gesteuert (kommagetrennt bzw. "ALL"), die vom Wrapper-Skript
# submit_benchmarking.sh per --export gesetzt werden. Die --array-Direktive oben
# ist dabei nur ein sicherer Fallback fuer einen (nicht empfohlenen) direkten
# "sbatch run_benchmarking_array.sh"-Aufruf ohne Wrapper - bitte IMMER ueber
# submit_benchmarking.sh submitten, da nur dort die Array-Groesse zur gewaehlten
# Kombination passend berechnet wird.
# -------------------------------------------------------------------

mkdir -p slurm_logs

source "/nfs/home/students/a.dersch/FoPra_PLAs/src/benchmarking/slurm_scripts/benchmarking_lib.sh"

# Kohorten filtern (Default: alle 6, siehe benchmarking_lib.sh)
DATASETS_LIST="${DATASETS_LIST:-ALL}"
filter_cohorts "$DATASETS_LIST"
COHORTEN=("${FILTERED_COHORTEN[@]}")
N_COHORTEN=${#COHORTEN[@]}

# Signaturen filtern (Default: alle 11, siehe benchmarking_lib.sh)
SIGNATURES_LIST="${SIGNATURES_LIST:-ALL}"
filter_signatures "$SIGNATURES_LIST"
SIGNATUREN=("${FILTERED_SIGNATUREN[@]}")
N_SIGS=${#SIGNATUREN[@]}

# Fallback, falls das Skript versehentlich ohne Wrapper direkt submittet wird:
# bewusst NUR eine Kombination (raw + biologist), nicht alle 6 - damit ein
# "vergessener" --export nicht aus Versehen wieder 396 Jobs auf einmal auslöst.
FILTER_MODES_LIST="${FILTER_MODES_LIST:-raw}"
GT_SOURCES_LIST="${GT_SOURCES_LIST:-biologist}"

# Leerzeichen entfernen, dann in Bash-Arrays splitten
IFS=',' read -ra FILTER_MODES <<< "$(echo "$FILTER_MODES_LIST" | tr -d ' ')"
IFS=',' read -ra GT_SOURCES   <<< "$(echo "$GT_SOURCES_LIST"   | tr -d ' ')"

N_MODUS=${#FILTER_MODES[@]}
N_GT=${#GT_SOURCES[@]}

# Validierung gegen bekannte, gueltige Werte - faengt Tippfehler beim Submit ab,
# statt dass Benchmarking_Master_v4.R spaeter mit "Fehler: Ungueltige GT_SOURCE
# uebergeben!" abbricht und man erst in den Slurm-Logs suchen muss, warum.
VALID_FILTER_MODES=("raw" "qc_tolerant" "qc_strict")
VALID_GT_SOURCES=("biologist" "gmm_dual" "gmm_single")

for fm in "${FILTER_MODES[@]}"; do
  if [[ ! " ${VALID_FILTER_MODES[*]} " =~ " ${fm} " ]]; then
    echo "FEHLER: Ungueltiger Filter_Mode '$fm'. Erlaubt: ${VALID_FILTER_MODES[*]}"
    exit 1
  fi
done
for gt in "${GT_SOURCES[@]}"; do
  if [[ ! " ${VALID_GT_SOURCES[*]} " =~ " ${gt} " ]]; then
    echo "FEHLER: Ungueltige GT_Source '$gt'. Erlaubt: ${VALID_GT_SOURCES[*]}"
    exit 1
  fi
done

TOTAL_JOBS=$((N_COHORTEN * N_SIGS * N_MODUS * N_GT))

# Konsistenzcheck: stimmt die tatsaechlich gebrauchte Array-Groesse mit der beim
# Submit gesetzten ueberein? SLURM_ARRAY_TASK_MAX ist die groesste Task-ID des
# aktuell laufenden Arrays.
if [ -n "${SLURM_ARRAY_TASK_MAX:-}" ] && [ "$SLURM_ARRAY_TASK_MAX" -ne "$TOTAL_JOBS" ]; then
  echo "WARNUNG: SLURM_ARRAY_TASK_MAX ($SLURM_ARRAY_TASK_MAX) weicht von der aus" \
       "FILTER_MODES_LIST/GT_SOURCES_LIST/SIGNATURES_LIST/DATASETS_LIST berechneten" \
       "Job-Anzahl ($TOTAL_JOBS) ab. Wurde ueber submit_benchmarking.sh submittet?"
fi

ZERO_INDEX=$((SLURM_ARRAY_TASK_ID - 1))

COHORTE_IDX=$((ZERO_INDEX % N_COHORTEN))
REST1=$((ZERO_INDEX / N_COHORTEN))

SIG_IDX=$((REST1 % N_SIGS))
REST2=$((REST1 / N_SIGS))

MODUS_IDX=$((REST2 % N_MODUS))
GT_IDX=$((REST2 / N_MODUS))

COHORTE=${COHORTEN[$COHORTE_IDX]}
SIG_ENTRY=${SIGNATUREN[$SIG_IDX]}
MODUS=${FILTER_MODES[$MODUS_IDX]}
GT_SOURCE=${GT_SOURCES[$GT_IDX]}

SIG_NAME="${SIG_ENTRY%%:*}"
SIG_FILE="${SIG_ENTRY##*:}"

METHOD="AUCell"
USE_EXT="FALSE"
THRESH="gmm_dist_dual"

CURRENT_FILE="${COHORTE}_${MODUS}_automated_gating.rds"

export PATH="/nfs/data/cluster/software/R/4.4.2/lib/R/bin:$PATH"
export R_LIBS_USER="/cmnfs/home/students/a.dersch/R/x86_64-pc-linux-gnu-library/4.4"

# Standardmaessig FALSE - RDS-Objekte nur speichern, wenn explizit fuer diese Tranche
# gewuenscht (siehe submit_benchmarking.sh). Verhindert, dass die NFS-Quota durch
# RDS-Dateien fuer viele Kombinationen gefuellt wird, obwohl nur eine Handvoll davon
# tatsaechlich per deep_analysis.R untersucht wird.
SAVE_RDS="${SAVE_RDS:-FALSE}"

echo "======================================================"
echo "Array Task ID:   $SLURM_ARRAY_TASK_ID  (von $TOTAL_JOBS)"
echo "Kohorte:         $COHORTE  (Index $COHORTE_IDX)"
echo "Signatur:        $SIG_NAME  (Index $SIG_IDX, aus: $(for e in "${SIGNATUREN[@]}"; do echo -n "${e%%:*} "; done))"
echo "Filter_Mode:     $MODUS  (Index $MODUS_IDX, aus: ${FILTER_MODES[*]})"
echo "GT_Source:       $GT_SOURCE  (Index $GT_IDX, aus: ${GT_SOURCES[*]})"
echo "SAVE_RDS:        $SAVE_RDS"
echo "Datei:           $CURRENT_FILE"
echo "======================================================"

Rscript src/benchmarking/Benchmarking_Master_v4.R \
  "$METHOD"       \
  "$SIG_NAME"     \
  "$SIG_FILE"     \
  "$USE_EXT"      \
  "$THRESH"       \
  "$CURRENT_FILE" \
  "$GT_SOURCE"     \
  "$SAVE_RDS"

EXIT_CODE=$?
if [ $EXIT_CODE -ne 0 ]; then
  echo "FEHLER: Task $SLURM_ARRAY_TASK_ID (${COHORTE} / ${SIG_NAME} / ${MODUS} / ${GT_SOURCE}) ist mit Exit-Code $EXIT_CODE fehlgeschlagen!"
  exit $EXIT_CODE
fi

echo "Task $SLURM_ARRAY_TASK_ID (${COHORTE} / ${SIG_NAME} / ${MODUS} / ${GT_SOURCE}) erfolgreich beendet."
