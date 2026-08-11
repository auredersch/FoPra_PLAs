#!/bin/bash
# -------------------------------------------------------------------
# Geteilte Konfiguration fuer submit_benchmarking.sh und run_benchmarking_array.sh.
# Zentral an EINER Stelle gepflegt, damit der Wrapper (der die Array-Groesse fuer
# sbatch berechnet) und der Array-Job (der die Signaturen tatsaechlich verwendet)
# nie auseinanderlaufen koennen - vorher haette man die Liste an zwei Stellen
# synchron halten muessen, was genau die Art von stillem Bug ist, die schwer
# aufzuspueren ist (falsche Job-Anzahl vs. tatsaechlich verwendete Signaturen).
# -------------------------------------------------------------------

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

COHORTEN=("heart" "sepsis" "vaccine" "immune_aging" "impact" "skin")

# filter_signatures <comma-getrennte Kuerzel oder "ALL"/leer>
# Setzt das globale Array FILTERED_SIGNATUREN. Bricht mit klarer Fehlermeldung ab,
# wenn ein angegebenes Kuerzel nicht existiert (faengt Tippfehler ab, z.B. "Hp"
# statt "HP", statt dass eine leere Signaturliste durchrutscht).
filter_signatures() {
  local requested="$1"
  FILTERED_SIGNATUREN=()

  if [[ -z "$requested" || "$requested" == "ALL" ]]; then
    FILTERED_SIGNATUREN=("${SIGNATUREN[@]}")
    return 0
  fi

  local req_clean
  req_clean="$(echo "$requested" | tr -d ' ')"
  local REQUESTED_SIGS
  IFS=',' read -ra REQUESTED_SIGS <<< "$req_clean"

  local all_short_names=""
  for entry in "${SIGNATUREN[@]}"; do
    all_short_names+="${entry%%:*} "
  done

  for req in "${REQUESTED_SIGS[@]}"; do
    local found=false
    for entry in "${SIGNATUREN[@]}"; do
      if [[ "${entry%%:*}" == "$req" ]]; then
        FILTERED_SIGNATUREN+=("$entry")
        found=true
        break
      fi
    done
    if [[ "$found" == false ]]; then
      echo "FEHLER: Signatur-Kuerzel '$req' nicht gefunden. Verfuegbar: $all_short_names"
      exit 1
    fi
  done
}

# filter_cohorts <comma-getrennte Kohorten oder "ALL"/leer>
# Setzt das globale Array FILTERED_COHORTEN. Gleiche Logik wie filter_signatures,
# nur dass COHORTEN einfache Strings sind (kein "Kuerzel:Datei"-Paar).
filter_cohorts() {
  local requested="$1"
  FILTERED_COHORTEN=()

  if [[ -z "$requested" || "$requested" == "ALL" ]]; then
    FILTERED_COHORTEN=("${COHORTEN[@]}")
    return 0
  fi

  local req_clean
  req_clean="$(echo "$requested" | tr -d ' ')"
  local REQUESTED_COHORTS
  IFS=',' read -ra REQUESTED_COHORTS <<< "$req_clean"

  for req in "${REQUESTED_COHORTS[@]}"; do
    local found=false
    for c in "${COHORTEN[@]}"; do
      if [[ "$c" == "$req" ]]; then
        FILTERED_COHORTEN+=("$c")
        found=true
        break
      fi
    done
    if [[ "$found" == false ]]; then
      echo "FEHLER: Kohorte '$req' nicht gefunden. Verfuegbar: ${COHORTEN[*]}"
      exit 1
    fi
  done
}