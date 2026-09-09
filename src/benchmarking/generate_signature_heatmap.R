library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

METRICS_DIR <- "/nfs/home/students/a.dersch/FoPra_PLAs/results/benchmarking/metrics"
OUTPUT_DIR  <- "/nfs/home/students/a.dersch/FoPra_PLAs/results/benchmarking"
HEATMAP_DIR <- file.path(OUTPUT_DIR, "signature_heatmaps")
dir.create(HEATMAP_DIR, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)

FILTER_MODES_ARG      <- if (length(args) >= 1 && nzchar(args[1])) args[1] else "raw,qc_tolerant"
GT_SOURCES_ARG         <- if (length(args) >= 2 && nzchar(args[2])) args[2] else "biologist,gmm_dual,gmm_single"
METRICS_ARG            <- if (length(args) >= 3 && nzchar(args[3])) args[3] else "F1,Prec,Rec"
FIXED_THRESHOLD_MODE   <- if (length(args) >= 4 && nzchar(args[4])) args[4] else "gmm_dist_dual"

FILTER_MODES <- trimws(strsplit(FILTER_MODES_ARG, ",")[[1]])
GT_SOURCES   <- trimws(strsplit(GT_SOURCES_ARG, ",")[[1]])
METRICS      <- trimws(strsplit(METRICS_ARG, ",")[[1]])

VALID_FILTER_MODES <- c("raw", "qc_tolerant", "qc_strict")
VALID_GT_SOURCES   <- c("biologist", "gmm_dual", "gmm_single")
VALID_METRICS      <- c("F1", "Prec", "Rec", "SRI")

invalid_fm <- setdiff(FILTER_MODES, VALID_FILTER_MODES)
invalid_gt <- setdiff(GT_SOURCES, VALID_GT_SOURCES)
invalid_metric <- setdiff(METRICS, VALID_METRICS)

if (length(invalid_fm) > 0) {
  stop("Ungueltige(r) Filter_Mode: ", paste(invalid_fm, collapse = ", "),
       ". Erlaubt: ", paste(VALID_FILTER_MODES, collapse = ", "))
}
if (length(invalid_gt) > 0) {
  stop("Ungueltige(r) GT_Source: ", paste(invalid_gt, collapse = ", "),
       ". Erlaubt: ", paste(VALID_GT_SOURCES, collapse = ", "))
}
if (length(invalid_metric) > 0) {
  stop("Ungueltige Metrik: ", paste(invalid_metric, collapse = ", "),
       ". Erlaubt: ", paste(VALID_METRICS, collapse = ", "))
}

print("======================================================")
print(paste("Filter_Modes:    ", paste(FILTER_MODES, collapse = ", ")))
print(paste("GT_Sources:      ", paste(GT_SOURCES, collapse = ", ")))
print(paste("Metriken:        ", paste(METRICS, collapse = ", ")))
print(paste("Threshold_Mode:  ", FIXED_THRESHOLD_MODE))
print("======================================================")

METRIC_LABELS <- c(F1 = "F1-Score", Prec = "Precision", Rec = "Recall", SRI = "Signal-weighted Recall")

all_files <- list.files(METRICS_DIR, pattern = "\\.csv$", full.names = TRUE)

if (length(all_files) == 0) {
  stop("Fehler: Keine Metrik-CSVs in ", METRICS_DIR, " gefunden!")
}

raw_list <- lapply(all_files, function(f) {
  df <- tryCatch(read_csv(f, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(df)) return(NULL)

  required <- c("Dataset", "Signature", "Filter_Mode", "GT_Source", "Threshold_Mode",
                "F1", "Prec", "Rec")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    warning("Datei ", basename(f), " uebersprungen - fehlende Spalten: ",
            paste(missing, collapse = ", "))
    return(NULL)
  }

  if (!"SRI" %in% colnames(df)) {
    warning("Datei ", basename(f), " enthaelt noch keine SRI-Spalte (vermutlich vor ",
            "Einfuehrung der Metrik gelaufen) - SRI wird als NA gesetzt.")
    df$SRI <- NA_real_
  }

  if (length(unique(df$Filter_Mode)) > 1 || length(unique(df$GT_Source)) > 1 ||
      length(unique(df$Threshold_Mode)) > 1) {
    warning("Datei ", basename(f), " enthaelt gemischte Filter_Mode/GT_Source/",
            "Threshold_Mode-Werte - wird uebersprungen, um stille Vermischung zu vermeiden.")
    return(NULL)
  }

  df %>% select(Dataset, Signature, Filter_Mode, GT_Source, Threshold_Mode, F1, Prec, Rec, SRI)
})

combined_all <- bind_rows(raw_list)

if (nrow(combined_all) == 0) {
  stop("Fehler: Nach dem Einlesen sind keine gueltigen Metrik-Zeilen uebrig!")
}

combined_all <- combined_all %>% filter(Threshold_Mode == FIXED_THRESHOLD_MODE)

if (nrow(combined_all) == 0) {
  stop("Fehler: Keine Zeilen mit Threshold_Mode == '", FIXED_THRESHOLD_MODE, "' gefunden!")
}

generate_single_heatmap <- function(data, filter_mode, gt_source, metric) {
  subset_df <- data %>%
    filter(Filter_Mode == filter_mode, GT_Source == gt_source) %>%
    group_by(Dataset, Signature) %>%
    summarise(Value = mean(.data[[metric]], na.rm = TRUE), .groups = "drop")

  if (nrow(subset_df) == 0) {
    warning("Keine Daten fuer Filter_Mode=", filter_mode, ", GT_Source=", gt_source,
            ", Metric=", metric, " - Heatmap wird uebersprungen.")
    return(NULL)
  }

  metric_label <- METRIC_LABELS[[metric]]

  plot_obj <- ggplot(subset_df, aes(x = Signature, y = Dataset, fill = Value)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = sprintf("%.2f", Value)), color = "black", size = 4, fontface = "bold") +
    scale_fill_gradient2(low = "#4B8BFF", mid = "#FFF9A6", high = "#FF4B4B",
                         midpoint = 0.5, limits = c(0, 1)) +
    theme_minimal() +
    labs(
      title = paste("Signature Benchmarking Performance -", metric_label),
      subtitle = paste0("Filter_Mode: ", filter_mode, " | GT_Source: ", gt_source,
                        " | Threshold_Mode: ", FIXED_THRESHOLD_MODE),
      x = "Platelet Signatures",
      y = "Datasets / Cohorts",
      fill = metric_label
    ) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 11, face = "bold"),
      axis.text.y = element_text(size = 11, face = "bold"),
      panel.grid = element_blank()
    )

  out_file <- file.path(HEATMAP_DIR,
    paste0("heatmap_", tolower(metric), "_", filter_mode, "_", gt_source, ".png"))
  ggsave(filename = out_file, plot = plot_obj, width = 10, height = 6, dpi = 300)
  print(paste("-> Heatmap gespeichert:", out_file))

  subset_df %>% mutate(Filter_Mode = filter_mode, GT_Source = gt_source, Metric = metric)
}

n_combos_expected <- length(FILTER_MODES) * length(GT_SOURCES) * length(METRICS)
print(paste("Erzeuge bis zu", n_combos_expected, "Einzel-Heatmaps (",
            length(FILTER_MODES), "Filter_Mode(s) x", length(GT_SOURCES), "GT_Source(s) x",
            length(METRICS), "Metrik(en) )..."))

all_combo_data <- list()

for (metric in METRICS) {
  for (filter_mode in FILTER_MODES) {
    for (gt_source in GT_SOURCES) {
      key <- paste(metric, filter_mode, gt_source, sep = "_")
      all_combo_data[[key]] <- generate_single_heatmap(combined_all, filter_mode, gt_source, metric)
    }
  }
}

combo_df <- bind_rows(all_combo_data)

if (nrow(combo_df) == 0) {
  stop("Fehler: Es konnten keine Heatmaps erzeugt werden - bitte Metrik-Daten pruefen.")
}

print("Gefundene Datenkombinationen (Dataset x Signature x Filter_Mode x GT_Source x Metric):")
print(table(combo_df$Dataset, combo_df$Metric))

generate_faceted_overview <- function(data, metric) {
  subset_df <- data %>% filter(Metric == metric)
  if (nrow(subset_df) == 0) return(NULL)

  metric_label <- METRIC_LABELS[[metric]]

  plot_obj <- ggplot(subset_df, aes(x = Signature, y = Dataset, fill = Value)) +
    geom_tile(color = "white", size = 0.4) +
    geom_text(aes(label = sprintf("%.2f", Value)), color = "black", size = 3, fontface = "bold") +
    scale_fill_gradient2(low = "#4B8BFF", mid = "#FFF9A6", high = "#FF4B4B",
                         midpoint = 0.5, limits = c(0, 1)) +
    facet_grid(Filter_Mode ~ GT_Source) +
    theme_minimal() +
    labs(
      title = paste("Signature Benchmarking Performance -", metric_label, "- Uebersicht"),
      subtitle = paste0("Alle Filter_Mode x GT_Source Kombinationen | Threshold_Mode: ", FIXED_THRESHOLD_MODE),
      x = "Platelet Signatures",
      y = "Datasets / Cohorts",
      fill = metric_label
    ) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9, face = "bold"),
      axis.text.y = element_text(size = 9, face = "bold"),
      strip.text = element_text(size = 10, face = "bold"),
      panel.grid = element_blank()
    )

  out_file <- file.path(HEATMAP_DIR, paste0("heatmap_overview_", tolower(metric), "_all_modes.png"))
  ggsave(filename = out_file, plot = plot_obj, width = 16, height = 9, dpi = 300)
  print(paste("-> Uebersichts-Heatmap gespeichert:", out_file))
}

for (metric in METRICS) {
  generate_faceted_overview(combo_df, metric)
}

message("Fertig! ", nrow(distinct(combo_df, Filter_Mode, GT_Source, Metric)),
        " Einzel- + ", length(METRICS), " Uebersichts-PNGs in: ", HEATMAP_DIR)