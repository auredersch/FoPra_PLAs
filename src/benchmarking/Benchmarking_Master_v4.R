# -------------------------------------------------------------------
# Benchmarking_Master_v5_MultiCohort.R - FINAL GROUND TRUTH SELECTOR WITH QC
# -------------------------------------------------------------------
start_time <- Sys.time() 

library(Seurat)
library(AUCell)
library(pROC)
#library(PRROC)     
library(dplyr)
library(ggplot2)
library(pheatmap)
library(tidyr)
library(mclust)

project_root <- "/nfs/home/students/a.dersch/FoPra_PLAs"

args <- commandArgs(trailingOnly = TRUE)
METHOD_NAME    <- if(length(args) >= 1) args[1] else "AUCell"
SIG_NAME       <- if(length(args) >= 2) args[2] else "MANNE_DN"
SIG_FILE_BASE  <- if(length(args) >= 3) args[3] else "MANNE_COVID19_COMBINED_COHORT_VS_HEALTHY_DONOR_PLATELETS_DN.v2025.1.Hs"
USE_EXTENSION  <- if(length(args) >= 4) as.logical(args[4]) else FALSE
THRESH_MODE    <- if(length(args) >= 5) args[5] else "gmm_dist_dual" 
CURRENT_FILE   <- if(length(args) >= 6) args[6] else "sepsis_qc_strict_automated_gating.rds" 
GT_SOURCE      <- if(length(args) >= 7) args[7] else "biologist"    # other: "gmm_dual" "gmm_single"
# Standardmaessig FALSE: RDS-Objekte werden nur fuer Kombinationen gebraucht, die per
# deep_analysis.R untersucht werden sollen (siehe Kommentar bei der Speicherung weiter
# unten) - fuer die grosse Heatmap-Matrix reicht die Metrics-CSV. Explizit auf TRUE
# setzen fuer die paar Kombinationen, die ihr wirklich deep-analysen wollt.
SAVE_RDS       <- if(length(args) >= 8) as.logical(args[8]) else FALSE

get_dataset_short <- function(filename) {
  if (grepl("heart", filename, ignore.case = TRUE)) return("heart")
  if (grepl("sepsis", filename, ignore.case = TRUE)) return("sepsis")
  if (grepl("vaccine", filename, ignore.case = TRUE)) return("vaccine")
  if (grepl("immune_aging", filename, ignore.case = TRUE)) return("immune_aging")
  if (grepl("impact", filename, ignore.case = TRUE)) return("impact")
  if (grepl("skin", filename, ignore.case = TRUE)) return("skin")
  return("unknown")
}

EXTRACTED_MODE <- case_when(
  grepl("raw", CURRENT_FILE) ~ "raw",
  grepl("qc_tolerant", CURRENT_FILE) ~ "qc_tolerant",
  grepl("qc_strict", CURRENT_FILE) ~ "qc_strict",
  TRUE ~ "raw"
)

DATASET_SHORT <- get_dataset_short(CURRENT_FILE)
if(DATASET_SHORT == "unknown") stop("Fehler: Datensatz-Kürzel konnte nicht ermittelt werden!")

INPUT_DIR   <- "/nfs/home/students/a.dersch/FoPra_PLAs/data/datasets_automated/"
PATH_DATA   <- file.path(INPUT_DIR, CURRENT_FILE)

RDS_OUT_DIR <- "/nfs/home/students/a.dersch/FoPra_PLAs/data/datasets/benchmarked_objects"
PLOT_BASE   <- "/nfs/home/students/a.dersch/FoPra_PLAs/results/benchmarking/"
OUT_DIR <- file.path(
  PLOT_BASE, 
  DATASET_SHORT, 
  EXTRACTED_MODE, 
  paste0(METHOD_NAME, "_Ext", USE_EXTENSION, "_", THRESH_MODE, "_GT_", GT_SOURCE),
  SIG_NAME, 
  "/"
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RDS_OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(project_root, "results/benchmarking/metrics"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(project_root, "results/benchmarking/celltype_data"), recursive = TRUE, showWarnings = FALSE)

print(paste("Lade Datensatz:", CURRENT_FILE, "als", DATASET_SHORT))
pbmc <- readRDS(PATH_DATA)

print("Erste 10 Gennamen im RNA-Assay:")
print(head(rownames(pbmc[["RNA"]]), 10))

if(!"celltype_clean" %in% colnames(pbmc@meta.data) && "lineage" %in% colnames(pbmc@meta.data)) {
    pbmc$celltype_clean <- pbmc$lineage
}
if(!"celltype.l3" %in% colnames(pbmc@meta.data)) {
    pbmc$celltype.l3 <- pbmc$celltype_clean
}

# --- QC-FILTERUNG ---
QC_BASE <- "/nfs/home/students/a.dersch/FoPra_PLAs/results/sample_qc/"

if (EXTRACTED_MODE != "raw") {
  qc_table_file <- file.path(QC_BASE, DATASET_SHORT, "19_final_sample_lineage_QC_table.csv")
  
  if (file.exists(qc_table_file)) {
    print(paste("-> Filtere Zellen für Benchmarking basierend auf Modus:", EXTRACTED_MODE))
    kollege_qc <- read.csv(qc_table_file)
    
    if (EXTRACTED_MODE == "qc_strict") {
      trusted_pairs <- kollege_qc %>% 
        filter(final_pair_status %in% c("trusted", "trusted_but_extreme_PLA"))
    } else if (EXTRACTED_MODE == "qc_tolerant") {
      trusted_pairs <- kollege_qc %>% 
        filter(final_pair_status %in% c("trusted", "trusted_but_extreme_PLA", 
                                        "sample_ADT_suspicious", "lineage_ADT_suspicious"))
    }
    
    actual_sample_col <- case_when(
      "donor_id" %in% colnames(pbmc@meta.data) ~ "donor_id",
      "sample_id" %in% colnames(pbmc@meta.data) ~ "sample_id",
      TRUE ~ "sample"
    )
    
    pbmc$csv_celltype_name <- case_when(
      pbmc$celltype_clean == "CD4 T" ~ "CD4 T cells",
      pbmc$celltype_clean == "CD8 T cells" ~ "CD8 T cells",
      TRUE ~ pbmc$celltype_clean
    )
    
    pbmc$match_key <- paste0(pbmc@meta.data[[actual_sample_col]], "_", pbmc$csv_celltype_name)
    trusted_keys <- paste0(trusted_pairs$sample_id, "_", trusted_pairs$celltype_id)
    
    keep_cells <- colnames(pbmc)[pbmc$match_key %in% trusted_keys]
    
    if(length(keep_cells) < 20) {
       stop("Kritischer Fehler: Zu wenige Zellen nach QC-Filter übrig!")
    }
    
    pbmc <- subset(pbmc, cells = keep_cells)
    print(paste("-> Filterung abgeschlossen. Verbleibende Zellen:", ncol(pbmc)))
  } else {
    warning("-> Keine QC-Tabelle gefunden. Benchmarking läuft auf ungefiltertem Objekt!")
  }
}

# --- GROUND TRUTH SELEKTION ---
if (GT_SOURCE == "biologist") {
    print("-> Nutze die manuelle Annotation der Biologin (pla_status) als Ground Truth.")
    if("pla_status" %in% colnames(pbmc@meta.data)) {
        GT_COLUMN <- "pla_status"
    } else if("pla.status" %in% colnames(pbmc@meta.data)) {
        GT_COLUMN <- "pla.status"
    } else {
        stop("Fehler: Spalte 'pla_status' wurde im Objekt nicht gefunden!")
    }
    POSITIVE_VAL <- "PLA"
    
} else if (GT_SOURCE == "gmm_dual") {
    print("-> Nutze automatisches 2D GMM Gating (automative_gating_double) als Ground Truth.")
    if("automative_gating_double" %in% colnames(pbmc@meta.data)) {
        GT_COLUMN <- "automative_gating_double"
    } else {
        stop("Fehler: Spalte 'automative_gating_double' fehlt!")
    }
    POSITIVE_VAL <- "PLA"

} else if (GT_SOURCE == "gmm_single") {
    print("-> Nutze automatisches 1D GMM Gating (automative_gating_single) als Ground Truth.")
    if("automative_gating_single" %in% colnames(pbmc@meta.data)) {
        GT_COLUMN <- "automative_gating_single"
    } else {
        stop("Fehler: Spalte 'automative_gating_single' fehlt!")
    }
    POSITIVE_VAL <- "PLA"

} else {
    stop("Fehler: Ungültige GT_SOURCE übergeben!")
}

# --- GENLISTE LADEN ---
source(file.path(project_root, "src", "benchmarking", "read_and_extend_gene_list.R"))
PATH_SIG <- file.path(project_root, "data", "signatures", paste0(SIG_FILE_BASE, ".csv"))
genes <- read_gene_list(PATH_SIG)

# Overlap-Check: bricht ab wenn zu wenige Gene matchen
gene_overlap <- intersect(genes, rownames(pbmc[["RNA"]]))
overlap_frac <- length(gene_overlap) / length(genes)
print(paste("Gen-Overlap Signatur <-> RNA-Assay:", length(gene_overlap), "von", length(genes), "Genen",
            "(", round(100 * overlap_frac, 1), "% )"))

# FIX: Der alte Check (< 10 absolute Gene) konnte zwei ganz unterschiedliche Faelle
# nicht unterscheiden: (a) echter Gen-ID-Mismatch (Symbol vs. Ensembl) - erkennbar an
# einer NIEDRIGEN Ueberlappungsquote, selbst bei grosser Signatur, und (b) eine von
# Natur aus KLEINE, aber korrekt gematchte Signatur (z.B. HP_ABNORMAL_PLATELET_
# MEMBRANE_PROTEIN_EXPRESSION mit nur 6 Genen total, davon 6/6 = 100% gematcht). Fall
# (b) wurde bisher faelschlich als "kritischer ID-Mismatch-Fehler" abgebrochen, obwohl
# der Gen-Abgleich technisch perfekt funktioniert hat - nur die Signatur ist klein.
MIN_ABSOLUTE_GENES   <- 3    # AUCell-Ranking wird bei extrem wenigen Genen wenig aussagekraeftig
MIN_OVERLAP_FRACTION <- 0.5  # deutlich niedrigere Quote deutet auf ID-Mismatch hin

if (length(gene_overlap) < MIN_ABSOLUTE_GENES) {
  stop(paste("Kritischer Fehler: Nur", length(gene_overlap), "Gene matchen - unter der",
             "absoluten Mindestanzahl von", MIN_ABSOLUTE_GENES, "fuer ein sinnvolles AUCell-Ranking.",
             "\nBeispiel Objekt:", paste(head(rownames(pbmc[["RNA"]]), 3), collapse=", "),
             "\nBeispiel Signatur:", paste(head(genes, 3), collapse=", ")))
} else if (overlap_frac < MIN_OVERLAP_FRACTION) {
  stop(paste("Kritischer Fehler: Nur", round(100 * overlap_frac, 1), "% der Signatur-Gene matchen (",
             length(gene_overlap), "von", length(genes), ") - das deutet auf einen Gen-ID-Mismatch hin",
             "(Symbol vs. Ensembl?), nicht auf eine kleine Signatur.",
             "\nBeispiel Objekt:", paste(head(rownames(pbmc[["RNA"]]), 3), collapse=", "),
             "\nBeispiel Signatur:", paste(head(genes, 3), collapse=", ")))
} else if (length(gene_overlap) < 10) {
  warning(paste0("Nur ", length(gene_overlap), " Gene matchen, ABER das sind ", round(100 * overlap_frac, 1),
                 "% der gesamten Signatur (", length(genes), " Gene total) - kein ID-Mismatch, sondern eine ",
                 "von Natur aus kleine Signatur. AUCell-Scores auf so wenigen Genen koennen verrauschter/",
                 "instabiler sein als bei groesseren Signaturen (z.B. UNION_ALL) - das beim Vergleich in ",
                 "der Heatmap im Hinterkopf behalten, nicht direkt 1:1 mit grossen Signaturen vergleichen."))
}

# --- IMMUNE CONFIG ---
IMMUNE_SIG <- "GOBP_LEUKOCYTE_ACTIVATION_INVOLVED_IN_INFLAMMATORY_RESPONSE.v2025.1.Hs"
PATH_IMMUNE_SIG <- file.path(project_root, "data", "signatures", paste0(IMMUNE_SIG, ".csv"))
immune_genes <- read_gene_list(PATH_IMMUNE_SIG)
immune_genes <- intersect(immune_genes, rownames(pbmc[["RNA"]]))

# --- RNA ASSAY SETZEN ---
if ("RNA" %in% Assays(pbmc)) {
    DefaultAssay(pbmc) <- "RNA"
    print("-> DefaultAssay auf 'RNA' umgestellt.")
} else {
    warning("-> WARNUNG: Kein RNA-Assay im Objekt gefunden!")
}

# --- SCORING ---
print(paste("--- Calculating Scores using", METHOD_NAME, "---"))

if (METHOD_NAME == "AUCell" || METHOD_NAME == "WeightedAUCell") {
    expression_matrix <- GetAssayData(pbmc, assay = "RNA", layer = "data")
    rankings <- AUCell_buildRankings(expression_matrix, plotStats = FALSE)
    
    auc_orig <- AUCell_calcAUC(list(Platelet_Orig = genes), rankings)
    pbmc$Raw_Score_Original <- as.numeric(getAUC(auc_orig)[1, ])

    auc_imm <- AUCell_calcAUC(list(Immune_Score = immune_genes), rankings)
    pbmc$Immune_Score <- as.numeric(getAUC(auc_imm)[1, ])

} else if (METHOD_NAME == "UCell") {
    pbmc <- AddModuleScore_UCell(pbmc, features = list(Platelet_Orig = genes), name = NULL)
    pbmc$Raw_Score_Original <- pbmc$Platelet_Orig
    pbmc <- AddModuleScore_UCell(pbmc, features = list(Immune_Score = immune_genes), name = NULL)
    pbmc$Immune_Score <- pbmc$Immune_Score_UCell  # Fix: korrekter Spaltenname

} else if (METHOD_NAME == "AddModuleScore") {
    pbmc <- AddModuleScore(pbmc, features = list(genes), name = "AMS_Orig")
    pbmc$Raw_Score_Original <- pbmc$AMS_Orig1
    pbmc <- AddModuleScore(pbmc, features = list(immune_genes), name = "AMS_Immune")
    pbmc$Immune_Score <- pbmc$AMS_Immune1
}

# --- GENLISTEN-ERWEITERUNG ---
if (USE_EXTENSION) {
    EXT_DIR <- file.path(project_root, "results/extended_lists")
    dir.create(EXT_DIR, recursive = TRUE, showWarnings = FALSE)
    
    EXT_FILE <- file.path(EXT_DIR, paste0("ext_", DATASET_SHORT, "_", SIG_NAME, "_", METHOD_NAME, ".csv"))
    
    if (file.exists(EXT_FILE)) {
        message("Lade existierende erweiterte Liste...")
        extended_genes <- read.csv(EXT_FILE)$geneName
    } else {
        message("Berechne neue Extension...")
        res_ext <- extend_gene_set(pbmc, base_genes = genes, score_name = "Raw_Score_Original")
        extended_genes <- res_ext$extended_genes
        write.csv(data.frame(geneName = extended_genes), EXT_FILE, row.names = FALSE)
    }
    final_genes <- extended_genes
} else {
    final_genes <- genes
}

# --- FINALES SCORING ---
if (METHOD_NAME == "AUCell" || METHOD_NAME == "WeightedAUCell") {
    pbmc$Raw_Score <- as.numeric(getAUC(AUCell_calcAUC(list(Platelet_Score = final_genes), rankings))[1, ])
} else if (METHOD_NAME == "UCell") {
    pbmc <- AddModuleScore_UCell(pbmc, features = list(Platelet_Score = final_genes), name = NULL)
    pbmc$Raw_Score <- pbmc$Platelet_Score
} else {
    pbmc <- AddModuleScore(pbmc, features = list(final_genes), name = "AMS")
    pbmc$Raw_Score <- pbmc$AMS1
}

pbmc$Z_Score   <- as.vector(scale(pbmc$Raw_Score))
pbmc$Immune_Z  <- as.vector(scale(pbmc$Immune_Score))

n_na_z <- sum(is.na(pbmc$Z_Score))
if (n_na_z > 0) {
  warning(paste("-> WARNUNG:", n_na_z, "Zellen haben NA als Z_Score und werden ignoriert."))
}

# --- GROUND TRUTH RESPONSE ---
pbmc$GT_Response <- ifelse(pbmc[[GT_COLUMN]] == POSITIVE_VAL, 1, 0)

# --- ROC & AU-PR ---
roc_obj <- roc(response = pbmc$GT_Response, predictor = pbmc$Z_Score,
               direction = "<", quiet = TRUE)

pos_scores <- pbmc$Z_Score[pbmc$GT_Response == 1]
neg_scores <- pbmc$Z_Score[pbmc$GT_Response == 0]

pos_scores <- pos_scores[!is.na(pos_scores)]
neg_scores <- neg_scores[!is.na(neg_scores)]

AUPR <- 0
pr_curve_data <- data.frame(recall = c(0, 1), precision = c(1, 0))

if (length(pos_scores) > 0 && length(neg_scores) > 0) {
  all_thresholds <- sort(unique(pbmc$Z_Score[!is.na(pbmc$Z_Score)]), decreasing = TRUE)
  pr_df <- data.frame(threshold = all_thresholds, precision = NA, recall = NA)
  
  for(i in seq_along(all_thresholds)) {
    t_val <- all_thresholds[i]
    cur_tp <- sum(pos_scores > t_val)
    cur_fp <- sum(neg_scores > t_val)
    
    pr_df$recall[i]    <- cur_tp / length(pos_scores)
    pr_df$precision[i] <- if((cur_tp + cur_fp) > 0) cur_tp / (cur_tp + cur_fp) else 1
  }
  
  pr_df <- pr_df[order(pr_df$recall), ]
  pr_curve_data <- pr_df
  
  for(i in 2:nrow(pr_df)) {
    delta_recall <- pr_df$recall[i] - pr_df$recall[i-1]
    mean_precision <- (pr_df$precision[i] + pr_df$precision[i-1]) / 2
    AUPR <- AUPR + (delta_recall * mean_precision)
  }
}
prevalence <- mean(pbmc$GT_Response, na.rm = TRUE)

print(paste("AU-ROC:", round(as.numeric(auc(roc_obj)), 4),
            "| AU-PR (Native):", round(AUPR, 4),
            "| Baseline (Prevalence):", round(prevalence, 4)))

# --- THRESHOLDING ---
THRESHOLD_I <- -Inf

if (THRESH_MODE == "youden") {
    THRESHOLD_Z <- as.numeric(coords(roc_obj, x = "best", best.method = "youden")$threshold)

} else if (THRESH_MODE == "percentile") {
    prob_cutoff <- 0.90
    THRESHOLD_Z <- as.numeric(quantile(pbmc$Z_Score, probs = prob_cutoff))

} else if (THRESH_MODE == "manual") {
    med <- median(pbmc$Z_Score)
    mad_val <- mad(pbmc$Z_Score)
    THRESHOLD_Z <- med + (1.5 * mad_val)

} else if (THRESH_MODE == "null_dist_platelet") {
    z_grid <- unique(quantile(pbmc$Z_Score, probs = seq(0.05, 0.95, 0.05)))
    best_f1 <- -1
    THRESHOLD_Z <- z_grid[1]
    for(tz in z_grid) {
        pred <- pbmc$Z_Score > tz
        tp <- sum(pred & pbmc$GT_Response == 1); fp <- sum(pred & pbmc$GT_Response == 0)
        fn <- sum(!pred & pbmc$GT_Response == 1)
        prec_t <- if((tp+fp) > 0) tp/(tp+fp) else 0
        rec_t  <- if((tp+fn) > 0) tp/(tp+fn) else 0
        f1_t   <- if((prec_t+rec_t) > 0) 2*(prec_t*rec_t)/(prec_t+rec_t) else 0
        if(!is.na(f1_t) && f1_t > best_f1) { best_f1 <- f1_t; THRESHOLD_Z <- tz }
    }

} else if (THRESH_MODE == "null_dist_immune_dual") {
    z_grid <- unique(quantile(pbmc$Z_Score, probs = seq(0.1, 0.9, 0.1)))
    i_grid <- unique(quantile(pbmc$Immune_Z, probs = seq(0.1, 0.9, 0.1)))
    best_f1 <- -1
    THRESHOLD_Z <- z_grid[1]; THRESHOLD_I <- i_grid[1]
    for(tz in z_grid) {
        for(ti in i_grid) {
            pred <- (pbmc$Z_Score > tz) & (pbmc$Immune_Z > ti)
            tp <- sum(pred & pbmc$GT_Response == 1); fp <- sum(pred & pbmc$GT_Response == 0)
            fn <- sum(!pred & pbmc$GT_Response == 1)
            prec_t <- if((tp+fp) > 0) tp/(tp+fp) else 0
            rec_t  <- if((tp+fn) > 0) tp/(tp+fn) else 0
            f1_t   <- if((prec_t+rec_t) > 0) 2*(prec_t*rec_t)/(prec_t+rec_t) else 0
            if(!is.na(f1_t) && f1_t > best_f1) {
                best_f1 <- f1_t; THRESHOLD_Z <- tz; THRESHOLD_I <- ti
            }
        }
    }

} else if (THRESH_MODE %in% c("gmm_dist_platelet", "gmm_dist_dual")) {
    auc_obs <- pbmc$Raw_Score
    auc_imm <- pbmc$Immune_Score

    fit_plat <- Mclust(auc_obs, G = 2)
    plat_high <- which.max(fit_plat$parameters$mean)
    pbmc$Platelet_High <- fit_plat$classification == plat_high
    pbmc$Immune_High <- TRUE

    if(THRESH_MODE == "gmm_dist_dual") {
        idx <- which(pbmc$Platelet_High)
        fit_imm <- Mclust(auc_imm[idx], G = 2)
        imm_high <- which.max(fit_imm$parameters$mean)
        pbmc$Immune_High <- FALSE
        pbmc$Immune_High[idx] <- fit_imm$classification == imm_high
        THRESHOLD_I <- min(pbmc$Immune_Z[pbmc$Immune_High], na.rm = TRUE)
    }
    THRESHOLD_Z <- min(pbmc$Z_Score[pbmc$Platelet_High], na.rm = TRUE)

} else if (THRESH_MODE == "kmeans") {
    set.seed(42)
    km_data <- FetchData(pbmc, vars = c("Z_Score", "Immune_Z")) %>% drop_na()
    km_fit <- kmeans(km_data, centers = 4, nstart = 50, iter.max = 100)
    pbmc$KMeans_Cluster <- NA
    pbmc$KMeans_Cluster[rownames(km_data)] <- km_fit$cluster

    cluster_stats <- km_data %>%
        mutate(Cluster = km_fit$cluster) %>%
        group_by(Cluster) %>%
        summarise(
            Mean_Z = mean(Z_Score), Mean_Immune_Z = mean(Immune_Z),
            Score_Sum = Mean_Z + Mean_Immune_Z,
            Min_Z = min(Z_Score), Min_Immune = min(Immune_Z)
        )
    positive_cluster <- cluster_stats$Cluster[which.max(cluster_stats$Score_Sum)]
    pbmc$Platelet_High <- pbmc$KMeans_Cluster == positive_cluster
    pbmc$Immune_High   <- pbmc$KMeans_Cluster == positive_cluster
    THRESHOLD_Z <- cluster_stats$Min_Z[cluster_stats$Cluster == positive_cluster]
    THRESHOLD_I <- cluster_stats$Min_Immune[cluster_stats$Cluster == positive_cluster]
}

# --- KLASSIFIKATION ---
positive_condition <- if (THRESH_MODE == "kmeans") {
    pbmc@meta.data$KMeans_Cluster == positive_cluster
} else if (THRESH_MODE == "gmm_dist_platelet") {
    pbmc@meta.data$Platelet_High
} else if (THRESH_MODE == "gmm_dist_dual") {
    pbmc@meta.data$Platelet_High & pbmc@meta.data$Immune_High
} else {
    (pbmc@meta.data$Z_Score > THRESHOLD_Z) & (pbmc@meta.data$Immune_Z > THRESHOLD_I)
}

positive_condition <- as.logical(as.vector(positive_condition))
pbmc$Prediction <- factor(ifelse(positive_condition, "Positive", "Negative"),
                          levels = c("Negative", "Positive"))

gt_values    <- as.character(pbmc@meta.data[[GT_COLUMN]])
pos_val_char <- as.character(POSITIVE_VAL)

pbmc$Error_Type <- case_when(
    pbmc$Prediction == "Positive" & gt_values == pos_val_char ~ "TP",
    pbmc$Prediction == "Positive" & gt_values != pos_val_char ~ "FP",
    pbmc$Prediction == "Negative" & gt_values == pos_val_char ~ "FN",
    pbmc$Prediction == "Negative" & gt_values != pos_val_char ~ "TN",
    TRUE ~ NA_character_
)

# --- METRIKEN SCHREIBEN (inkl. AU-PR) ---
# WICHTIG (nach Fehleranalyse): Dieser Block steht jetzt VOR der RDS-Speicherung,
# nicht mehr danach. Grund: saveRDS() ist der mit Abstand fehleranfaelligste und
# teuerste Schritt (grosse Datei, NFS-Schreiblast bei vielen parallelen Array-Jobs).
# Vorher fuehrte ein Schreibfehler bei saveRDS() (z.B. "error writing to connection")
# zu einem sofortigen Skriptabbruch ("Execution halted") - und damit gingen auch die
# laengst berechneten Metriken (F1/Prec/Rec/SRI/AUROC/AUPR) verloren, obwohl die gar
# nichts mit dem RDS-Schreibfehler zu tun hatten. Jetzt werden die Metriken zuerst
# gesichert, damit ein RDS-Problem hoechstens den RDS-Export betrifft, nicht den
# gesamten Benchmarking-Lauf.
runtime_min <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
tp <- sum(pbmc$Error_Type == "TP", na.rm = TRUE)
fp <- sum(pbmc$Error_Type == "FP", na.rm = TRUE)
fn <- sum(pbmc$Error_Type == "FN", na.rm = TRUE)
tn <- sum(pbmc$Error_Type == "TN", na.rm = TRUE)

prec <- if((tp + fp) > 0) tp / (tp + fp) else 0
rec  <- if((tp + fn) > 0) tp / (tp + fn) else 0
f1   <- if((prec + rec) > 0) 2 * prec * rec / (prec + rec) else 0

# --- NEUE METRIK: Signal-weighted Recall (Spaltenname bleibt "SRI") ---
# Gewichteter Recall: statt jede echte PLA-Zelle binaer (gefunden ja/nein) gleich zu
# zaehlen, wird sie mit ihrem eigenen Raw_Score gewichtet. Eine stark positive,
# knapp verpasste PLA "kostet" damit mehr als eine schwach positive, knapp verpasste.
# WICHTIG: Raw_Score (immer >= 0), NICHT Z_Score - Z_Score kann negativ sein
# (Zellen unter dem Mittelwert), was Zaehler/Nenner unsinnig machen wuerde.
gt_positive_mask <- pbmc$Error_Type %in% c("TP", "FN")
detected_mask     <- pbmc$Error_Type == "TP"

sri_denominator <- sum(pbmc$Raw_Score[gt_positive_mask], na.rm = TRUE)
sri <- if (sri_denominator > 0) {
  sum(pbmc$Raw_Score[detected_mask], na.rm = TRUE) / sri_denominator
} else {
  NA_real_
}
if (is.na(sri)) {
  warning("SRI konnte nicht berechnet werden (keine GT-positiven Zellen mit Raw_Score > 0 vorhanden).")
}

# Prevalence als Baseline für AU-PR Interpretation
prevalence <- mean(pbmc$GT_Response, na.rm = TRUE)

METRICS_FILE <- file.path(project_root, "results/benchmarking/metrics",
  paste0("metrics_", DATASET_SHORT, "_", METHOD_NAME, "_", SIG_NAME, "_",
         THRESH_MODE, "_", EXTRACTED_MODE, "_GT_", GT_SOURCE, ".csv"))

write.csv(data.frame(
    Dataset        = DATASET_SHORT,
    Method         = METHOD_NAME,
    Signature      = SIG_NAME,
    Ext            = USE_EXTENSION,
    Filter_Mode    = EXTRACTED_MODE,
    GT_Source      = GT_SOURCE,
    Threshold_Mode = THRESH_MODE,
    # Ranking-Metriken (schwellenwertfrei)
    AUROC          = as.numeric(auc(roc_obj)),
    AUPR           = AUPR,                       
    Prevalence     = prevalence,                 
    F1             = f1,
    Prec           = prec,
    Rec            = rec,
    SRI            = sri,
    TP             = tp,
    FP             = fp,
    FN             = fn,
    TN             = tn,
    Threshold_Value = THRESHOLD_Z,
    Runtime_Min    = runtime_min
), METRICS_FILE, row.names = FALSE)

print(paste("-> Metriken gespeichert:", METRICS_FILE))

# --- SAVE RDS (optional, gesteuert ueber SAVE_RDS) ---
# FIX 1: SIG_NAME wird jetzt mit in den Dateinamen aufgenommen. Vorher fehlte das -
# bei parallelen Array-Jobs (mehrere Signaturen fuer denselben Dataset/Filter_Mode/
# GT_Source) haben sich die Jobs gegenseitig dieselbe RDS-Datei ueberschrieben, und
# je nach Fertigstellungsreihenfolge landete eine beliebige Signatur im File.
# FIX 2: Das ist aber gleichzeitig die Ursache fuer die "error writing to connection"-
# Fehler - vorher wurde pro Dataset/Filter_Mode/GT_Source EINE RDS-Datei ueberschrieben
# (wenig Netto-Speicherbedarf), jetzt wird JEDE Signatur einzeln gespeichert (~11x mehr
# RDS-Dateien als vorher). Das fuellt die NFS-Quota deutlich schneller. Da fuer die
# Heatmaps nur die Metrics-CSV gebraucht wird (siehe oben) und die RDS-Objekte nur fuer
# deep_analysis.R noetig sind, ist SAVE_RDS jetzt standardmaessig FALSE - explizit auf
# TRUE setzen (8. CLI-Argument) nur fuer die Kombinationen, die wirklich per Deep
# Analysis untersucht werden sollen.
if (SAVE_RDS) {
  RDS_FILE_NAME <- file.path(RDS_OUT_DIR,
    paste0("pbmc_benchmarked_", DATASET_SHORT, "_", METHOD_NAME, "_", SIG_NAME, "_",
           THRESH_MODE, "_", EXTRACTED_MODE, "_GT_", GT_SOURCE, ".rds"))
  message("Speichere Seurat-Objekt unter: ", RDS_FILE_NAME)

  save_result <- tryCatch({
    saveRDS(pbmc, file = RDS_FILE_NAME)
    TRUE
  }, error = function(e) {
    message("WARNUNG: RDS-Speicherung fehlgeschlagen (", conditionMessage(e), ").")
    message("Metriken wurden bereits erfolgreich gespeichert (siehe oben) - nur der ",
            "RDS-Export fuer diese Kombination fehlt. Moegliche Ursache: NFS-Quota voll ",
            "(pruefen mit 'quota -s' bzw. 'df -h ", RDS_OUT_DIR, "') oder kurzzeitige ",
            "NFS-Schreiblast bei vielen parallelen Jobs.")
    FALSE
  })

  if (!save_result) {
    # Bewusst KEIN stop()/quit(status=1) hier - der Lauf war inhaltlich erfolgreich
    # (Metriken vorhanden), nur der optionale RDS-Export ist fehlgeschlagen. Ein
    # Nicht-Null-Exit-Code wuerde im Slurm-Array-Log als "FEHLER" markiert werden,
    # obwohl die fuer die Heatmap relevanten Daten vollstaendig sind.
    message("-> Lauf wird als erfolgreich gewertet (Metriken vollstaendig), RDS-Export übersprungen.")
  }
} else {
  message("-> SAVE_RDS = FALSE, RDS-Objekt wird fuer diese Kombination nicht gespeichert.")
}

# --- CELLTYPE-DATEN ---
ct_data <- pbmc@meta.data %>%
    group_by(celltype_clean, celltype.l3) %>%
    summarise(
        TP      = sum(Error_Type == "TP", na.rm = TRUE),
        FP      = sum(Error_Type == "FP", na.rm = TRUE),
        FN      = sum(Error_Type == "FN", na.rm = TRUE),
        TN      = sum(Error_Type == "TN", na.rm = TRUE),
        Mean_Z  = mean(Z_Score, na.rm = TRUE),
        n       = n(),
        .groups = "drop"
    ) %>%
    mutate(
        Accuracy         = (TP + TN) / n,
        Balanced_Accuracy = 0.5 * ((TP / (TP + FN + 1e-6)) + (TN / (TN + FP + 1e-6))),
        Dataset    = DATASET_SHORT,
        Method     = METHOD_NAME,
        Signature  = SIG_NAME,
        Mode       = THRESH_MODE,
        Filter_Mode = EXTRACTED_MODE,
        GT_Source  = GT_SOURCE
    )

CT_FILE <- file.path(project_root, "results/benchmarking/celltype_data",
  paste0("ct_", DATASET_SHORT, "_", METHOD_NAME, "_", SIG_NAME, "_",
         THRESH_MODE, "_", EXTRACTED_MODE, "_GT_", GT_SOURCE, ".csv"))

write.csv(ct_data, CT_FILE, row.names = FALSE)

# --- PLOTS ---
available_reductions <- names(pbmc@reductions)
print(paste("Verfügbare Reduktionen:", paste(available_reductions, collapse = ", ")))

CHOSEN_REDUCTION <- if (DATASET_SHORT == "immune_aging") "GEX_umap_mrvi" else "umap_totalVI"
if (!CHOSEN_REDUCTION %in% names(pbmc@reductions)) {
  warning(paste("Gewünschte Reduktion", CHOSEN_REDUCTION, "fehlt. Weiche auf erste verfügbare aus."))
  CHOSEN_REDUCTION <- names(pbmc@reductions)[1]
}
print(paste("-> Visualisierung mit Reduktion:", CHOSEN_REDUCTION))

png(paste0(OUT_DIR, "1_Reference_UMAP.png"), 1200, 800)
print(DimPlot(pbmc, reduction = CHOSEN_REDUCTION, group.by = "celltype_clean", label = TRUE, repel = TRUE) +
      labs(title = paste("Reference Gating (", CHOSEN_REDUCTION, ") -", DATASET_SHORT)))
dev.off()

png(paste0(OUT_DIR, "2_Score_UMAP.png"), 900, 700)
print(FeaturePlot(pbmc, reduction = CHOSEN_REDUCTION, features = "Raw_Score") +
      scale_colour_viridis_c(option = "magma") +
      labs(title = paste("Raw Score -", SIG_NAME)))
dev.off()

png(paste0(OUT_DIR, "3_Error_Mapping_UMAP.png"), 1000, 800)
print(DimPlot(pbmc, reduction = CHOSEN_REDUCTION, group.by = "Error_Type") +
      scale_color_manual(values = c("TP" = "#228B22", "FP" = "#FF4500", "FN" = "#1E90FF", "TN" = "#D3D3D3")) +
      labs(title = "Error Mapping", subtitle = paste("Threshold Z =", round(THRESHOLD_Z, 2))))
dev.off()

# --- 3b. CONFUSION MATRIX AUF DATASET-EBENE (aggregiert ueber alle Zelltypen) ---
# Ergaenzung zur Error-UMAP: dieselbe Konfiguration, gleicher OUT_DIR, damit beide
# Plots direkt nebeneinander interpretiert werden koennen. Die zelltyp-aufgeloeste
# Confusion Matrix bleibt Aufgabe von deep_analysis.R (andere Detailebene).
cm_dataset <- data.frame(
  Prediction  = factor(c("Positive", "Positive", "Negative", "Negative"),
                        levels = c("Positive", "Negative")),
  GroundTruth = factor(c("Positive", "Negative", "Positive", "Negative"),
                        levels = c("Positive", "Negative")),
  Label = c("TP", "FP", "FN", "TN"),
  Count = c(tp, fp, fn, tn)
)
cm_dataset$Frac <- cm_dataset$Count / sum(cm_dataset$Count)

png(paste0(OUT_DIR, "3b_Confusion_Matrix_Dataset.png"), 800, 700)
print(
  ggplot(cm_dataset, aes(x = GroundTruth, y = Prediction, fill = Label)) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(aes(label = paste0(Label, "\n", Count, "\n(", sprintf("%.1f%%", 100 * Frac), ")")),
              color = "black", size = 6, fontface = "bold") +
    scale_fill_manual(values = c("TP" = "#228B22", "FP" = "#FF4500",
                                 "FN" = "#1E90FF", "TN" = "#D3D3D3")) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "none", panel.grid = element_blank()) +
    labs(title = paste("Confusion Matrix (Dataset-Ebene) -", DATASET_SHORT),
         subtitle = paste0(SIG_NAME, " | ", METHOD_NAME, " | GT: ", GT_SOURCE,
                           " | Filter: ", EXTRACTED_MODE, " | N = ", sum(cm_dataset$Count)),
         x = "Ground Truth", y = "Prediction")
)
dev.off()

png(paste0(OUT_DIR, "4a_Density_ZScore_PLA_Status.png"), 1200, 800)
print(ggplot(pbmc@meta.data, aes(x = Z_Score, fill = !!sym(GT_COLUMN))) +
      geom_density(alpha = 0.5) + theme_minimal() +
      scale_fill_manual(values = c("PLA" = "#FF4B4B", "platelet-free" = "#4B8BFF")) +
      geom_vline(xintercept = THRESHOLD_Z, linetype = "dashed", color = "red", size = 1) +
      labs(title = "Global Z-Score Distribution",
           subtitle = paste("Threshold Z:", round(THRESHOLD_Z, 2),
                            "| F1:", round(f1, 3),
                            "| AUROC:", round(as.numeric(auc(roc_obj)), 3),
                            "| AUPR:", round(AUPR, 3)),
           x = "Z-Score", fill = "PLA Status"))
dev.off()

png(paste0(OUT_DIR, "4b_Density_RawScore_PLA_Status.png"), 1200, 800)
print(ggplot(pbmc@meta.data, aes(x = Raw_Score, fill = !!sym(GT_COLUMN))) +
      geom_density(alpha = 0.5) + theme_minimal() +
      scale_fill_manual(values = c("PLA" = "#FF4B4B", "platelet-free" = "#4B8BFF")) +
      labs(title = "Global Raw Score Distribution"))
dev.off()

png(paste0(OUT_DIR, "6a_Violin_RawScore_Celltypes.png"), 1200, 600)
print(VlnPlot(pbmc, features = "Raw_Score", group.by = "celltype_clean", pt.size = 0))
dev.off()

png(paste0(OUT_DIR, "6b_Split_Violin_PLA_Status.png"), 1400, 700)
print(VlnPlot(pbmc, features = "Raw_Score", group.by = "celltype_clean",
              split.by = GT_COLUMN, pt.size = 0) +
      scale_fill_manual(values = c("PLA" = "#FF4B4B", "platelet-free" = "#4B8BFF")))
dev.off()

# ROC + PR Kurve nebeneinander
png(paste0(OUT_DIR, "9_ROC_PR_Curves.png"), 1400, 700)
par(mfrow = c(1, 2))

plot(roc_obj, col = "#E41A1C", lwd = 3,
     main = paste("ROC  |  AUC:", round(as.numeric(auc(roc_obj)), 3)))

plot(pr_curve_data$recall, pr_curve_data$precision, type = "l", col = "#1E90FF", lwd = 3,
     xlim = c(0, 1), ylim = c(0, 1), xlab = "Recall", ylab = "Precision",
     main = paste("PR Curve  |  AU-PR:", round(AUPR, 3), " | Baseline:", round(prevalence, 3)))
abline(h = prevalence, lty = 2, col = "grey50")

par(mfrow = c(1, 1))
dev.off()

message("Done! Alle Ergebnisse für ", DATASET_SHORT, " (", SIG_NAME, ", ", EXTRACTED_MODE, ") in: ", OUT_DIR)
