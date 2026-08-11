#!/bin/bash
# -------------------------------------------------------------------
# Wrapper zum tranchenweisen Submitten des Benchmarking-Arrays.
# Berechnet die noetige Array-Groesse fuer die gewaehlte(n) Filter_Mode/GT_Source/
# Signatur/Dataset-Kombination(en) und uebergibt sie per --array an sbatch (kann
# NICHT im Job-Skript selbst stehen, da #SBATCH-Direktiven beim Einreichen
# statisch gelesen werden).
#
# Nutzung (benannte Flags, Reihenfolge beliebig - bleibt stabil, auch wenn spaeter
# weitere Filter-Dimensionen dazukommen):
#   ./submit_benchmarking.sh --filter <Filter_Modes> --gt <GT_Sources> \
#       [--sig <Signaturen>] [--datasets <Kohorten>] [--save-rds TRUE|FALSE] [-y]
#
# Alle Werte kommagetrennt fuer mehrere; "ALL" (Default) = alle verfuegbaren Werte.
#
# Beispiele:
#   ./run_benchmarking.sh --filter raw --gt biologist
#       -> alle 6 Kohorten x alle 11 Signaturen = 66 Jobs
#
#   ./run_benchmarking.sh --filter raw --gt biologist --sig HP
#       -> alle 6 Kohorten x nur HP = 6 Jobs
#
#   ./run_benchmarking.sh --filter raw --gt gmm_dual --datasets sepsis,skin --sig HP
#       -> nur sepsis+skin x nur HP = 2 Jobs
#
#   ./run_benchmarking.sh --filter raw --gt gmm_dual --datasets sepsis --save-rds TRUE -y
#       -> nur sepsis, alle Signaturen, inkl. RDS-Speicherung, ohne Rueckfrage
#
# Verfuegbare Kuerzel: siehe benchmarking_lib.sh (SIGNATUREN, COHORTEN).
# -------------------------------------------------------------------
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "/nfs/home/students/a.dersch/FoPra_PLAs/src/benchmarking/slurm_scripts/benchmarking_lib.sh"

FILTER_MODES_LIST=""
GT_SOURCES_LIST=""
SIGNATURES_ARG="ALL"
DATASETS_ARG="ALL"
SAVE_RDS_ARG="FALSE"
AUTO_YES=""

usage() {
  echo "Usage: $0 --filter <Filter_Modes> --gt <GT_Sources> [--sig <Signaturen>] [--datasets <Kohorten>] [--save-rds TRUE|FALSE] [-y]"
  echo "Beispiel: $0 --filter raw --gt gmm_dual --datasets sepsis --sig HP"
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --filter)    FILTER_MODES_LIST="$2"; shift 2 ;;
    --gt)        GT_SOURCES_LIST="$2"; shift 2 ;;
    --sig)       SIGNATURES_ARG="$2"; shift 2 ;;
    --datasets)  DATASETS_ARG="$2"; shift 2 ;;
    --save-rds)  SAVE_RDS_ARG="$2"; shift 2 ;;
    -y|--yes)    AUTO_YES="-y"; shift ;;
    -h|--help)   usage ;;
    *)
      echo "FEHLER: Unbekanntes Argument '$1'."
      usage
      ;;
  esac
done

if [[ -z "$FILTER_MODES_LIST" || -z "$GT_SOURCES_LIST" ]]; then
  echo "FEHLER: --filter und --gt sind Pflichtargumente."
  usage
fi

if [[ "$SAVE_RDS_ARG" != "TRUE" && "$SAVE_RDS_ARG" != "FALSE" ]]; then
  echo "FEHLER: --save-rds muss TRUE oder FALSE sein, nicht '$SAVE_RDS_ARG'."
  exit 1
fi

IFS=',' read -ra FM_ARR <<< "$(echo "$FILTER_MODES_LIST" | tr -d ' ')"
IFS=',' read -ra GT_ARR <<< "$(echo "$GT_SOURCES_LIST"   | tr -d ' ')"

VALID_FILTER_MODES=("raw" "qc_tolerant" "qc_strict")
VALID_GT_SOURCES=("biologist" "gmm_dual" "gmm_single")

for fm in "${FM_ARR[@]}"; do
  if [[ ! " ${VALID_FILTER_MODES[*]} " =~ " ${fm} " ]]; then
    echo "FEHLER: Ungueltiger Filter_Mode '$fm'. Erlaubt: ${VALID_FILTER_MODES[*]}"
    exit 1
  fi
done
for gt in "${GT_ARR[@]}"; do
  if [[ ! " ${VALID_GT_SOURCES[*]} " =~ " ${gt} " ]]; then
    echo "FEHLER: Ungueltige GT_Source '$gt'. Erlaubt: ${VALID_GT_SOURCES[*]}"
    exit 1
  fi
done

# Signaturen und Kohorten validieren/filtern ueber dieselben Funktionen, die auch
# run_benchmarking_array.sh nutzt (aus benchmarking_lib.sh) - garantiert, dass
# Wrapper (Array-Groessen-Berechnung) und Array-Job (tatsaechliche Ausfuehrung) immer
# dieselbe Auswahl zugrunde legen.
filter_signatures "$SIGNATURES_ARG"
N_SIGS=${#FILTERED_SIGNATUREN[@]}
SIG_SHORT_NAMES="$(for e in "${FILTERED_SIGNATUREN[@]}"; do echo -n "${e%%:*} "; done)"

filter_cohorts "$DATASETS_ARG"
N_COHORTEN=${#FILTERED_COHORTEN[@]}

N_MODUS=${#FM_ARR[@]}
N_GT=${#GT_ARR[@]}
TOTAL_JOBS=$((N_COHORTEN * N_SIGS * N_MODUS * N_GT))

echo "======================================================"
echo "Naechste Tranche:"
echo "  Filter_Mode(s): $FILTER_MODES_LIST"
echo "  GT_Source(s):   $GT_SOURCES_LIST"
echo "  Dataset(s):     $DATASETS_ARG  ->  ${FILTERED_COHORTEN[*]}"
echo "  Signatur(en):   $SIGNATURES_ARG  ->  $SIG_SHORT_NAMES"
echo "  SAVE_RDS:       $SAVE_RDS_ARG"
echo "  Jobs gesamt:    $TOTAL_JOBS  ($N_COHORTEN Kohorte(n) x $N_SIGS Signatur(en) x $N_MODUS Filter_Mode(s) x $N_GT GT_Source(s))"
if [[ "$SAVE_RDS_ARG" == "TRUE" ]]; then
  echo "  WARNUNG: SAVE_RDS=TRUE speichert bis zu $TOTAL_JOBS einzelne RDS-Dateien -"
  echo "  bei grossen Datasets kann das die NFS-Quota fuellen."
fi
echo "======================================================"

if [[ "$AUTO_YES" != "-y" ]]; then
  read -r -p "Fortfahren und $TOTAL_JOBS Jobs submitten? [y/N] " CONFIRM
  if [[ "$CONFIRM" != "y" && "$CONFIRM" != "Y" ]]; then
    echo "Abgebrochen."
    exit 0
  fi
fi

mkdir -p slurm_logs

JOBID=$(sbatch --parsable \
  --array=1-"$TOTAL_JOBS" \
  --export=ALL,FILTER_MODES_LIST="$FILTER_MODES_LIST",GT_SOURCES_LIST="$GT_SOURCES_LIST",SIGNATURES_LIST="$SIGNATURES_ARG",DATASETS_LIST="$DATASETS_ARG",SAVE_RDS="$SAVE_RDS_ARG" \
  "$SCRIPT_DIR/run_benchmarking_array.sh")

echo "-> Submittet als Job-ID $JOBID ($TOTAL_JOBS Tasks)."
echo "-> Status pruefen mit: squeue -j $JOBID"