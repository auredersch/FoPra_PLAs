# Detaillierte FP/FN-Analyse pro Datensatz
# Input: benchmarked RDS aus Benchmarking_Master_v4

library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(pheatmap)

project_root <- "/nfs/home/students/a.dersch/FoPra_PLAs"

args <- commandArgs(trailingOnly = TRUE)
RDS_FILE   <- if(length(args) >= 1) args[1] else "pbmc_benchmarked_vaccine_AUCell_gmm_dist_dual_raw_GT_biologist.rds"
GT_SOURCE  <- if(length(args) >= 2) args[2] else "gmm_dual"
SIG_NAME   <- if(length(args) >= 3) args[3] else "MANNE_DN"
# Optional (opt-in, da AUCell-Rankings dafuer neu aufgebaut werden muessen - kostet
# spuerbar Zeit): vergleicht die aktuell verwendete Immune_Score-Signatur mit einer
# alternativen Leukozyten-Signatur, um zu testen, ob eine spezifischere Signatur die
# myeloide FP-Haeufung reduzieren wuerde (Hypothese aus der Impact-Analyse).
RUN_LEUKOCYTE_COMPARISON <- if(length(args) >= 4) as.logical(args[4]) else FALSE

RDS_DIR    <- file.path(project_root, "data/datasets/benchmarked_objects")

# ADT-Platelet-Marker mit Alias-Liste: robust gegenueber unterschiedlicher Feature-
# Benennung zwischen Datasets (z.B. Gen-Symbol vs. CITE-seq-Panel-Name vs. Seurat-
# Suffix ".1" bei Namenskollision zwischen RNA- und ADT-Assay). Reihenfolge pro
# Marker ist Prioritaet: spezifischere/eindeutigere Aliase zuerst, generische zuletzt.
# Neue Marker koennen einfach als weiterer Listen-Eintrag ergaenzt werden, der Rest
# des Skripts (Suche, Plot, Benennung) passt sich automatisch an.
marker_aliases <- list(
  CD41  = c("CD41", "ITGA2B.1", "ITGA2B", "GPIIB"),
  CD61  = c("CD61", "ITGB3.1", "ITGB3", "GPIIIA"),
  CD62P = c("CD62P", "SELP.1", "SELP", "PSELECTIN", "CD62")
)

get_dataset_short <- function(filename) {
  if (grepl("heart",        filename, ignore.case = TRUE)) return("heart")
  if (grepl("sepsis",       filename, ignore.case = TRUE)) return("sepsis")
  if (grepl("vaccine",      filename, ignore.case = TRUE)) return("vaccine")
  if (grepl("immune_aging", filename, ignore.case = TRUE)) return("immune_aging")
  if (grepl("impact",       filename, ignore.case = TRUE)) return("impact")
  if (grepl("skin",         filename, ignore.case = TRUE)) return("skin")
  return("unknown")
}
DATASET_SHORT <- get_dataset_short(RDS_FILE)

OUT_DIR <- file.path(project_root, "results/benchmarking/fp_fn_analysis",
                     DATASET_SHORT, paste0(SIG_NAME, "_GT_", GT_SOURCE))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

print(paste("Lade:", RDS_FILE))
pbmc <- readRDS(file.path(RDS_DIR, RDS_FILE))

if (GT_SOURCE == "biologist") {
  GT_COLUMN <- if ("pla_status" %in% colnames(pbmc@meta.data)) "pla_status" else "pla.status"
} else if (GT_SOURCE == "gmm_dual") {
  GT_COLUMN <- "automative_gating_double"
} else {
  GT_COLUMN <- "automative_gating_single"
}
POSITIVE_VAL <- "PLA"

required_cols <- c("Error_Type", "Z_Score", "Raw_Score", "celltype_clean", GT_COLUMN)
missing <- setdiff(required_cols, colnames(pbmc@meta.data))
if (length(missing) > 0) stop(paste("Fehlende Spalten:", paste(missing, collapse = ", ")))

meta <- pbmc@meta.data %>% filter(!is.na(Error_Type))
print(paste("Zellen für Analyse:", nrow(meta)))
print(table(meta$Error_Type))

png(file.path(OUT_DIR, "1_Score_Distribution_GT.png"), 1200, 600)
meta$GT_Label <- ifelse(meta[[GT_COLUMN]] == POSITIVE_VAL, "GT Positiv (PLA)", "GT Negativ (platelet-free)")
print(
  ggplot(meta, aes(x = Z_Score, fill = GT_Label)) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(values = c("GT Positiv (PLA)" = "#FF4B4B",
                                 "GT Negativ (platelet-free)" = "#4B8BFF")) +
    facet_wrap(~ celltype_clean, scales = "free_y") +
    theme_minimal(base_size = 11) +
    labs(title = paste("Score-Verteilung GT+ vs GT- |", DATASET_SHORT, "|", SIG_NAME),
         subtitle = "Trennbarkeit pro Zelltyp — gut trennbar = zwei klare Peaks",
         x = "Z-Score", fill = "Ground Truth") +
    theme(legend.position = "bottom")
)
dev.off()

fn_meta <- meta %>% filter(Error_Type == "FN")
tp_meta <- meta %>% filter(Error_Type == "TP")

png(file.path(OUT_DIR, "2_FN_Score_Distribution.png"), 1000, 600)
if (nrow(fn_meta) > 0) {
  bind_rows(
    fn_meta %>% mutate(Gruppe = "FN (missed PLAs)"),
    tp_meta %>% mutate(Gruppe = "TP (detected PLAs)")
  ) %>%
  ggplot(aes(x = Z_Score, fill = Gruppe)) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(values = c("FN (missed PLAs)" = "#1E90FF",
                                 "TP (detected PLAs)" = "#228B22")) +
    theme_minimal(base_size = 12) +
    labs(title = paste("Z-Score der FN vs. TP Zellen |", DATASET_SHORT),
         subtitle = "FN mit hohem Score = Threshold-Problem | FN mit niedrigem Score = biologischer Grenzfall",
         x = "Z-Score", fill = NULL) +
    theme(legend.position = "bottom")
} else {
  plot.new(); text(0.5, 0.5, "Keine FN-Zellen vorhanden", cex = 1.5)
}
dev.off()
print(paste("Median Z-Score FN:", round(median(fn_meta$Z_Score, na.rm=TRUE), 3)))
print(paste("Median Z-Score TP:", round(median(tp_meta$Z_Score, na.rm=TRUE), 3)))

error_summary <- meta %>%
  group_by(celltype_clean, Error_Type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(celltype_clean) %>%
  mutate(frac = n / sum(n),
         total = sum(n))

png(file.path(OUT_DIR, "3a_ErrorType_Barplot_Relative.png"), 1200, 600)
print(
  ggplot(error_summary, aes(x = reorder(celltype_clean, -total), y = frac, fill = Error_Type)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("TP"="#228B22","FP"="#FF4500","FN"="#1E90FF","TN"="#D3D3D3")) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste("Error-Typ Verteilung pro Zelltyp (relativ) |", DATASET_SHORT),
         x = "Zelltyp", y = "Anteil", fill = "Error Type")
)
dev.off()

png(file.path(OUT_DIR, "3b_ErrorType_Barplot_Absolute.png"), 1200, 600)
print(
  ggplot(error_summary, aes(x = reorder(celltype_clean, -total), y = n, fill = Error_Type)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("TP"="#228B22","FP"="#FF4500","FN"="#1E90FF","TN"="#D3D3D3")) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste("Error-Typ Verteilung pro Zelltyp (absolut) |", DATASET_SHORT),
         x = "Zelltyp", y = "Anzahl Zellen", fill = "Error Type")
)
dev.off()

cm_wide <- error_summary %>%
  select(celltype_clean, Error_Type, n) %>%
  pivot_wider(names_from = Error_Type, values_from = n, values_fill = 0)

# Sicherstellen dass alle 4 Spalten vorhanden
for (col in c("TP","FP","FN","TN")) {
  if (!col %in% colnames(cm_wide)) cm_wide[[col]] <- 0
}

cm_mat <- as.matrix(cm_wide[, c("TP","FP","FN","TN")])
rownames(cm_mat) <- cm_wide$celltype_clean

annotation_row <- data.frame(
  FPR = round(cm_mat[,"FP"] / pmax(cm_mat[,"FP"] + cm_mat[,"TN"], 1), 2),
  FNR = round(cm_mat[,"FN"] / pmax(cm_mat[,"FN"] + cm_mat[,"TP"], 1), 2),
  row.names = rownames(cm_mat)
)

png(file.path(OUT_DIR, "4_Confusion_Matrix_Heatmap.png"), 1000, 700)
pheatmap(cm_mat,
         cluster_rows = FALSE, cluster_cols = FALSE,
         display_numbers = TRUE, number_format = "%d",
         color = colorRampPalette(c("white", "#FF6B35"))(50),
         annotation_row = annotation_row,
         main = paste("Confusion Matrix pro Zelltyp |", DATASET_SHORT, "|", SIG_NAME),
         fontsize = 11)
dev.off()

adt_assay <- if ("ADT_corrected" %in% Assays(pbmc)) "ADT_corrected" else "ADT"

if (adt_assay %in% Assays(pbmc)) {
  adt_features <- rownames(pbmc[[adt_assay]])
  
  # Sucht pro Marker den ersten passenden Feature-Namen entlang der Alias-Prioritaet.
  # fixed = TRUE, da Aliase Literalstrings sind (kein Regex-Interpretationsrisiko
  # bei Punkten in Suffixen wie "ITGA2B.1"). Gibt zusaetzlich zurueck, ueber welchen
  # Alias-Rang der Treffer erfolgte, um generische/mehrdeutige Treffer zu markieren.
  find_adt <- function(features, patterns) {
    for (i in seq_along(patterns)) {
      hit <- grep(patterns[i], features, value = TRUE, ignore.case = TRUE, fixed = TRUE)
      if (length(hit) > 0) return(list(feature = hit[1], alias_rank = i, alias_used = patterns[i]))
    }
    return(NULL)
  }

  matched_markers <- lapply(names(marker_aliases), function(marker_name) {
    res <- find_adt(adt_features, marker_aliases[[marker_name]])
    if (is.null(res)) {
      print(paste0("-> Marker '", marker_name, "': kein passendes Feature gefunden (Aliase: ",
                   paste(marker_aliases[[marker_name]], collapse = ", "), ")."))
      return(NULL)
    }
    # Warnung, falls nur der letzte (generischste) Alias in der Liste gegriffen hat -
    # z.B. koennte "CD62" statt "CD62P" auch "CD62L" (L-Selectin, biologisch andere
    # Funktion) treffen. Das Skript nutzt den Treffer trotzdem, aber macht auf das
    # Risiko einer Fehlzuordnung aufmerksam, statt es stillschweigend zu uebernehmen.
    if (res$alias_rank == length(marker_aliases[[marker_name]]) && length(marker_aliases[[marker_name]]) > 1) {
      warning(paste0("Marker '", marker_name, "' wurde nur ueber den generischsten Alias ('",
                     res$alias_used, "') als Feature '", res$feature, "' gefunden - ",
                     "bitte manuell pruefen, ob das tatsaechlich der gewuenschte Marker ist ",
                     "(z.B. 'CD62' koennte auch CD62L statt CD62P treffen)."))
    }
    res$feature
  })
  names(matched_markers) <- names(marker_aliases)

  features_to_plot <- unlist(matched_markers, use.names = TRUE)
  
  if (length(features_to_plot) > 0) {
    old_assay <- DefaultAssay(pbmc)
    DefaultAssay(pbmc) <- adt_assay
    
    # Nur FP, TN, TP, FN Zellen
    pbmc_sub <- subset(pbmc, cells = rownames(meta))
    
    png(file.path(OUT_DIR, "5_ADT_Expression_FP_vs_TN.png"), 1400, 500 * length(features_to_plot))
    print(
      VlnPlot(pbmc_sub,
              features = features_to_plot,
              group.by = "Error_Type",
              pt.size = 0,
              ncol = length(features_to_plot)) &
        scale_fill_manual(values = c("TP"="#228B22","FP"="#FF4500",
                                     "FN"="#1E90FF","TN"="#D3D3D3")) &
        labs(subtitle = paste(DATASET_SHORT, "— ADT Expression der Error-Typen")) &
        theme(legend.position = "none")
    )
    dev.off()
    
    DefaultAssay(pbmc) <- old_assay
  } else {
    print("-> Keine ADT-Platelet-Marker gefunden. Plot 5 übersprungen.")
  }
} else {
  print("-> Kein ADT-Assay vorhanden. Plot 5 übersprungen.")
}

threshold_z <- min(meta$Z_Score[meta$Error_Type %in% c("TP","FP")], na.rm = TRUE)

png(file.path(OUT_DIR, "6_Threshold_vs_Distribution.png"), 1200, 500)
print(
  ggplot(meta, aes(x = Z_Score, fill = Error_Type)) +
    geom_histogram(bins = 80, alpha = 0.7, position = "stack") +
    scale_fill_manual(values = c("TP"="#228B22","FP"="#FF4500",
                                 "FN"="#1E90FF","TN"="#D3D3D3")) +
    geom_vline(xintercept = threshold_z, color = "red", linetype = "dashed", size = 1) +
    annotate("text", x = threshold_z + 0.05, y = Inf, vjust = 2,
             label = paste("Threshold =", round(threshold_z, 3)), color = "red", size = 4) +
    theme_minimal(base_size = 12) +
    labs(title = paste("Threshold-Position in der Score-Verteilung |", DATASET_SHORT),
         subtitle = "FN links vom Threshold = zu hoch gesetzt | FP rechts = zu sensitiv",
         x = "Z-Score", y = "Anzahl Zellen", fill = "Error Type")
)
dev.off()

# --- 7. 2D-THRESHOLD-SCATTER (bivariate Entscheidungsgrenze bei gmm_dist_dual) ---
# Wichtig: Plot 6 (oben) zeichnet nur eine 1D-Linie auf Basis von Z_Score, obwohl die
# tatsaechliche Klassifikation bei THRESH_MODE == "gmm_dist_dual" bivariat ist
# (Platelet_High & Immune_High, siehe Benchmarking_Master_v4.R). Eine Zelle mit hohem
# Z_Score kann trotzdem FN sein, wenn Immune_High fehlschlaegt - das wuerde in Plot 6
# faelschlich als "Threshold zu hoch gesetzt" interpretiert werden. Dieser Plot zeigt
# die reale 2D-Grenze und macht sichtbar, an welcher Dimension eine Zelle "scheitert".
has_bivariate_cols <- all(c("Platelet_High", "Immune_High", "Immune_Z") %in% colnames(pbmc@meta.data))

if (has_bivariate_cols) {
  meta$Immune_Z      <- pbmc@meta.data[rownames(meta), "Immune_Z"]
  meta$Platelet_High <- pbmc@meta.data[rownames(meta), "Platelet_High"]
  meta$Immune_High   <- pbmc@meta.data[rownames(meta), "Immune_High"]

  # Rekonstruktion identisch zur Logik in Benchmarking_Master_v4.R (THRESH_MODE == "gmm_dist_dual"):
  # THRESHOLD_Z = kleinster Z_Score innerhalb der Platelet_High-Gruppe
  # THRESHOLD_I = kleinster Immune_Z innerhalb der Immune_High-Gruppe
  thresh_z_2d <- suppressWarnings(min(meta$Z_Score[meta$Platelet_High], na.rm = TRUE))
  thresh_i_2d <- suppressWarnings(min(meta$Immune_Z[meta$Immune_High], na.rm = TRUE))

  if (is.finite(thresh_z_2d) && is.finite(thresh_i_2d)) {
    png(file.path(OUT_DIR, "7_Threshold_2D_Scatter.png"), 1100, 900)
    print(
      ggplot(meta, aes(x = Z_Score, y = Immune_Z, color = Error_Type)) +
        geom_point(alpha = 0.35, size = 0.8) +
        geom_vline(xintercept = thresh_z_2d, color = "black", linetype = "dashed", linewidth = 0.8) +
        geom_hline(yintercept = thresh_i_2d, color = "black", linetype = "dashed", linewidth = 0.8) +
        scale_color_manual(values = c("TP" = "#228B22", "FP" = "#FF4500",
                                     "FN" = "#1E90FF", "TN" = "#D3D3D3")) +
        annotate("text", x = thresh_z_2d, y = max(meta$Immune_Z, na.rm = TRUE),
                 label = paste("Platelet-Cutoff =", round(thresh_z_2d, 2)),
                 vjust = -0.5, size = 3.5) +
        annotate("text", x = max(meta$Z_Score, na.rm = TRUE), y = thresh_i_2d,
                 label = paste("Immune-Cutoff =", round(thresh_i_2d, 2)),
                 hjust = 1, vjust = -0.5, size = 3.5) +
        theme_minimal(base_size = 12) +
        labs(title = paste("Bivariate Entscheidungsgrenze (gmm_dist_dual) |", DATASET_SHORT, "|", SIG_NAME),
             subtitle = paste("FN oben-rechts vom Kreuzpunkt = Immune-Dimension verantwortlich |",
                              "FN unten-links = Platelet-Dimension verantwortlich"),
             x = "Z_Score (Platelet)", y = "Immune_Z", color = "Error Type")
    )
    dev.off()
  } else {
    print("-> Konnte 2D-Thresholds nicht rekonstruieren (keine gueltigen Platelet_High/Immune_High-Werte). Plot 7 uebersprungen.")
  }
} else {
  print("-> Spalten Platelet_High/Immune_High/Immune_Z fehlen (vermutlich anderer Threshold_Mode als gmm_dist_dual). Plot 7 uebersprungen, Plot 6 (1D) bleibt als Fallback bestehen.")
}

# --- 8. LEUKOCYTE-SIGNATUR-VERGLEICH (optional, --> RUN_LEUKOCYTE_COMPARISON) ---
# Testet Hypothese 2 aus der Impact-Analyse: liefert eine spezifischere Leukozyten-
# Aktivierungssignatur (leukocyte_activation.csv) eine bessere Trennschaerfe fuer die
# Immune-Dimension als die aktuell in Benchmarking_Master_v4.R verwendete generische
# GOBP_LEUKOCYTE_ACTIVATION_INVOLVED_IN_INFLAMMATORY_RESPONSE-Signatur? Bewusst nur
# als Vergleichsplot umgesetzt (keine volle Re-Klassifikation mit neuem GMM-Fit) -
# haelt den Umfang proportional zur explorativen Fragestellung, statt eine zweite
# vollstaendige Benchmarking-Pipeline zu duplizieren.
if (RUN_LEUKOCYTE_COMPARISON) {
  if (!"Immune_Z" %in% colnames(pbmc@meta.data)) {
    print("-> Immune_Z fehlt im Objekt (vermutlich anderer Threshold_Mode) - Leukocyte-Vergleich uebersprungen.")
  } else {
    library(AUCell)
    source(file.path(project_root, "src", "benchmarking", "read_and_extend_gene_list.R"))

    CURRENT_IMMUNE_SIG <- "GOBP_LEUKOCYTE_ACTIVATION_INVOLVED_IN_INFLAMMATORY_RESPONSE.v2025.1.Hs"
    PATH_CURRENT_IMMUNE_SIG <- file.path(project_root, "data", "signatures",
                                          paste0(CURRENT_IMMUNE_SIG, ".csv"))
    PATH_ALT_IMMUNE_SIG <- file.path(project_root, "data", "signatures", "leukocyte_activation.csv")

    if (!file.exists(PATH_ALT_IMMUNE_SIG)) {
      print(paste("-> Alternative Signatur nicht gefunden unter", PATH_ALT_IMMUNE_SIG,
                  "- Leukocyte-Vergleich uebersprungen."))
    } else {
      current_genes <- read_gene_list(PATH_CURRENT_IMMUNE_SIG)
      alt_genes     <- read_gene_list(PATH_ALT_IMMUNE_SIG)

      # Sanity-Check VOR der teuren AUCell-Neuberechnung: falls beide Listen
      # (fast) identisch sind, waere der Vergleich wenig aussagekraeftig - lieber
      # jetzt abbrechen/warnen, als erst nach mehreren Minuten Rechenzeit zu merken,
      # dass beide Signaturen praktisch dasselbe Gen-Set sind.
      overlap_frac <- length(intersect(current_genes, alt_genes)) /
                        length(union(current_genes, alt_genes))
      print(paste0("-> Gen-Set-Ueberlappung (Jaccard) zwischen aktueller und alternativer ",
                   "Immune-Signatur: ", round(100 * overlap_frac, 1), "%",
                   " (", length(current_genes), " vs. ", length(alt_genes), " Gene)"))
      if (overlap_frac > 0.9) {
        print("-> WARNUNG: Ueberlappung > 90% - die beiden Signaturen sind fast identisch, ")
        print("   ein Vergleich wird vermutlich kaum unterscheidbare Ergebnisse liefern.")
      }

      alt_genes_present <- intersect(alt_genes, rownames(pbmc[["RNA"]]))
      if (length(alt_genes_present) < 5) {
        print(paste("-> Nur", length(alt_genes_present), "Gene der alternativen Signatur im",
                    "Datensatz vorhanden - Leukocyte-Vergleich uebersprungen (zu wenige Gene)."))
      } else {
        print("-> Baue AUCell-Rankings neu auf, um alternative Immune-Signatur zu scoren (kann dauern)...")
        expression_matrix <- GetAssayData(pbmc, assay = "RNA", layer = "data")
        rankings_alt <- AUCell_buildRankings(expression_matrix, plotStats = FALSE)
        auc_alt <- AUCell_calcAUC(list(Immune_Score_Alt = alt_genes_present), rankings_alt)
        pbmc$Immune_Score_Alt <- as.numeric(getAUC(auc_alt)[1, ])
        pbmc$Immune_Z_Alt     <- as.vector(scale(pbmc$Immune_Score_Alt))

        meta$Immune_Z_Alt <- pbmc@meta.data[rownames(meta), "Immune_Z_Alt"]

        png(file.path(OUT_DIR, "8_Leukocyte_Signature_Comparison.png"), 1300, 1000)
        print(
          ggplot(meta, aes(x = Immune_Z, y = Immune_Z_Alt, color = Error_Type)) +
            geom_point(alpha = 0.3, size = 0.6) +
            geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
            facet_wrap(~ celltype_clean, scales = "free") +
            scale_color_manual(values = c("TP" = "#228B22", "FP" = "#FF4500",
                                         "FN" = "#1E90FF", "TN" = "#D3D3D3")) +
            theme_minimal(base_size = 11) +
            labs(
              title = paste("Aktuelle vs. alternative Immune-Signatur |", DATASET_SHORT, "|", SIG_NAME),
              subtitle = paste0("x: ", CURRENT_IMMUNE_SIG, " (aktuell)  |  y: leukocyte_activation (alternativ) | ",
                                "Gen-Set-Overlap: ", round(100 * overlap_frac, 1), "% | ",
                                "Punkte UNTER der Diagonale = alternative Signatur bewertet Zelle niedriger ",
                                "(koennte myeloide FPs unter den Immune-Cutoff druecken)"),
              x = "Immune_Z (aktuell)", y = "Immune_Z (alternativ)", color = "Error Type"
            ) +
            theme(strip.text = element_text(size = 9, face = "bold"))
        )
        dev.off()
        print("-> Plot 8 (Leukocyte_Signature_Comparison) gespeichert.")

        # Kurze quantitative Zusammenfassung fuer genau die Zellen, die uns
        # interessieren: myeloide FPs. Wie viele wuerden mit der alternativen
        # Signatur einen NIEDRIGEREN Immune_Z bekommen (= Hinweis, dass die
        # alternative Signatur sie eher aussortieren wuerde)?
        fp_cells <- meta %>% filter(Error_Type == "FP")
        if (nrow(fp_cells) > 0) {
          frac_lower <- mean(fp_cells$Immune_Z_Alt < fp_cells$Immune_Z, na.rm = TRUE)
          print(paste0("-> Bei ", round(100 * frac_lower, 1), "% der FP-Zellen (n=", nrow(fp_cells),
                       ") liegt Immune_Z_Alt niedriger als Immune_Z - moegliches Signal, ",
                       "dass die alternative Signatur einen Teil der FPs aussortieren wuerde. ",
                       "Das ist eine Tendenzaussage, KEIN Beweis - dafuer waere eine echte ",
                       "Re-Klassifikation mit neu gefittetem GMM-Threshold noetig."))
        }
      }
    }
  }
}

summary_table <- cm_wide %>%
  mutate(
    Recall    = round(TP / pmax(TP + FN, 1), 3),
    Precision = round(TP / pmax(TP + FP, 1), 3),
    FPR       = round(FP / pmax(FP + TN, 1), 3),
    FNR       = round(FN / pmax(FN + TP, 1), 3),
    F2        = round((1 + 4) * (Precision * Recall) / pmax((4 * Precision) + Recall, 1e-6), 3),
    N_total   = TP + FP + FN + TN,
    N_PLA_GT  = TP + FN
  ) %>%
  arrange(desc(FNR))

write.csv(summary_table,
          file.path(OUT_DIR, "summary_per_celltype.csv"),
          row.names = FALSE)

print("=== Zusammenfassung pro Zelltyp (sortiert nach FNR) ===")
print(summary_table %>% select(celltype_clean, N_PLA_GT, Recall, FPR, FNR, F2))

message("Fertig! Alle Plots in: ", OUT_DIR)