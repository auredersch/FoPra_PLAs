# -------------------------------------------------------------------
# Benchmarking_Master_v5_MultiCohort.R - FINAL GROUND TRUTH SELECTOR WITH QC
# -------------------------------------------------------------------
start_time <- Sys.time() 

library(Seurat)
library(AUCell)
library(pROC)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(tidyr)
library(mclust)

project_root <- "/nfs/home/students/a.dersch/FoPra_PLAs"

# --- CLI ARGUMENTE / DEFAULTS ---
args <- commandArgs(trailingOnly = TRUE)
METHOD_NAME    <- if(length(args) >= 1) args[1] else "AUCell"
SIG_NAME       <- if(length(args) >= 2) args[2] else "MANNE_DN"
SIG_FILE_BASE  <- if(length(args) >= 3) args[3] else "MANNE_COVID19_COMBINED_COHORT_VS_HEALTHY_DONOR_PLATELETS_DN.v2025.1.Hs"
USE_EXTENSION  <- if(length(args) >= 4) as.logical(args[4]) else FALSE
THRESH_MODE    <- if(length(args) >= 5) args[5] else "gmm_dist_dual" 
CURRENT_FILE   <- if(length(args) >= 6) args[6] else "sepsis_qc_strict_automated_gating.rds" 
GT_SOURCE      <- if(length(args) >= 7) args[7] else "biologist"    # other: "gmm_dual" "gmm_single"

# --- DATASET MAP (FIXED WITH grepl & impact) ---
get_dataset_short <- function(filename) {
  if (grepl("heart", filename, ignore.case = TRUE)) return("heart")
  if (grepl("sepsis", filename, ignore.case = TRUE)) return("sepsis")
  if (grepl("vaccine", filename, ignore.case = TRUE)) return("vaccine")
  if (grepl("immune_aging", filename, ignore.case = TRUE)) return("immune_aging")
  if (grepl("impact", filename, ignore.case = TRUE)) return("impact")
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

# --- CONFIGURATION ---
INPUT_DIR   <- "/nfs/home/students/a.dersch/FoPra_PLAs/data/datasets_automated/"
PATH_DATA   <- file.path(INPUT_DIR, CURRENT_FILE)

# --- OUTPUT STRUKTUR ---
RDS_OUT_DIR <- "/nfs/home/students/a.dersch/FoPra_PLAs/data/datasets/benchmarked_objects/"
PLOT_BASE   <- "/nfs/home/students/a.dersch/FoPra_PLAs/results/benchmarking/"
OUT_DIR <- file.path(PLOT_BASE, DATASET_SHORT, EXTRACTED_MODE, paste0(METHOD_NAME, "_Ext", USE_EXTENSION, "_", THRESH_MODE, "_GT_", GT_SOURCE), "/")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RDS_OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(project_root, "results/benchmarking/metrics"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(project_root, "results/benchmarking/celltype_data"), recursive = TRUE, showWarnings = FALSE)

# --- DATEN LADEN ---
print(paste("Lade Datensatz:", CURRENT_FILE, "als", DATASET_SHORT))
pbmc <- readRDS(PATH_DATA)

# Test
head(rownames(pbmc), 10)

if(!"celltype_clean" %in% colnames(pbmc@meta.data) && "lineage" %in% colnames(pbmc@meta.data)) {
    pbmc$celltype_clean <- pbmc$lineage
}
if(!"celltype.l3" %in% colnames(pbmc@meta.data)) {
    pbmc$celltype.l3 <- pbmc$celltype_clean
}

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
    print(paste("-> Filterung abgeschlossen. Verbleibende Zellen für Evaluierung:", ncol(pbmc)))
  } else {
    warning("-> Keine QC-Tabelle gefunden. Benchmarking läuft auf ungefiltertem Objekt!")
  }
}

# --- DYNAMISCHE SELEKTION DER GROUND TRUTH ---
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
    print("-> Nutze dein automatisches 2D GMM Gating (automative_gating_double) als Ground Truth.")
    if("automative_gating_double" %in% colnames(pbmc@meta.data)) {
        GT_COLUMN <- "automative_gating_double"
    } else {
        stop("Fehler: Spalte 'automative_gating_double' fehlt!")
    }
    POSITIVE_VAL <- "PLA"

} else if (GT_SOURCE == "gmm_single") {
    print("-> [Automatic 1D GMM Gating (automative_gating_single) als Ground Truth.")
    if("automative_gating_single" %in% colnames(pbmc@meta.data)) {
        GT_COLUMN <- "automative_gating_single"
    } else {
        stop("Fehler: Spalte 'automative_gating_single' fehlt!")
    }
    POSITIVE_VAL <- "PLA"

}else {
    stop("Fehler: Ungültige GT_SOURCE übergeben!")
}

# --- GENLISTE LADEN ---
base_dir <- getwd()
source(file.path(project_root, "src", "benchmarking", "read_and_extend_gene_list.R"))
PATH_SIG <- file.path(project_root, "data", "signatures", paste0(SIG_FILE_BASE, ".csv"))
genes <- read_gene_list(PATH_SIG)

# --- Immune Config --- 
IMMUNE_SIG <- "GOBP_LEUKOCYTE_ACTIVATION_INVOLVED_IN_INFLAMMATORY_RESPONSE.v2025.1.Hs"
PATH_IMMUNE_SIG <- file.path(project_root, "data", "signatures", paste0(IMMUNE_SIG, ".csv"))
immune_genes <- read_gene_list(PATH_IMMUNE_SIG)
#immune_genes <- intersect(immune_genes, rownames(pbmc))
immune_genes <- intersect(immune_genes, rownames(pbmc[["RNA"]]))

# Scoring

if ("RNA" %in% Assays(pbmc)) {
    DefaultAssay(pbmc) <- "RNA"
    print("-> DefaultAssay global für das Scoring auf 'RNA' umgestellt.")
} else {
    warning("-> WARNUNG: Kein RNA-Assay im Objekt gefunden!")
}

print(paste("--- Calculating Scores using", METHOD_NAME, "---"))

if (METHOD_NAME == "AUCell" || METHOD_NAME == "WeightedAUCell") {
    expression_matrix <- GetAssayData(pbmc, assay = "RNA", layer = "data")
    rankings <- AUCell_buildRankings(expression_matrix, plotStats=FALSE)
    auc_orig <- AUCell_calcAUC(list(Platelet_Orig = genes), rankings)
    pbmc$Raw_Score_Original <- as.numeric(getAUC(auc_orig)[1, ])

    auc_imm <- AUCell_calcAUC(list(Immune_Score = immune_genes), rankings)
    pbmc$Immune_Score <- as.numeric(getAUC(auc_imm)[1, ])

} else if (METHOD_NAME == "UCell") {
    pbmc <- AddModuleScore_UCell(pbmc, features = list(Platelet_Orig = genes), name = NULL)
    pbmc$Raw_Score_Original <- pbmc$Platelet_Orig
    pbmc <- AddModuleScore_UCell(pbmc, features = list(Immune_Score = immune_genes), name = NULL)
    pbmc$Immune_Score <- pbmc$Immune_Score 

} else if (METHOD_NAME == "AddModuleScore") {
    pbmc <- AddModuleScore(pbmc, features = list(genes), name = "AMS_Orig")
    pbmc$Raw_Score_Original <- pbmc$AMS_Orig1
    pbmc <- AddModuleScore(pbmc, features = list(immune_genes), name = "AMS_Immune")
    pbmc$Immune_Score <- pbmc$AMS_Immune1
}

# Gensequenz-Erweiterung
if (USE_EXTENSION) {
    EXT_DIR <- file.path(project_root, "results/extended_lists")
    dir.create(EXT_DIR, recursive = TRUE, showWarnings = FALSE) # Sicherheitsnetz!
    
    EXT_FILE <- file.path(EXT_DIR, paste0("ext_", DATASET_SHORT, "_", SIG_NAME, "_", METHOD_NAME, ".csv"))
    
    if (file.exists(EXT_FILE)) {
        message("Lade existierende Liste...")
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

# Finales Scoring
if (METHOD_NAME == "AUCell" || METHOD_NAME == "WeightedAUCell") {
    pbmc$Raw_Score <- as.numeric(getAUC(AUCell_calcAUC(list(Platelet_Score = final_genes), rankings))[1, ])
} else if (METHOD_NAME == "UCell") {
    pbmc <- AddModuleScore_UCell(pbmc, features = list(Platelet_Score = final_genes), name = NULL)
    pbmc$Raw_Score <- pbmc$Platelet_Score
} else {
    pbmc <- AddModuleScore(pbmc, features = list(final_genes), name = "AMS")
    pbmc$Raw_Score <- pbmc$AMS1
}

pbmc$Z_Score <- as.vector(scale(pbmc$Raw_Score))
pbmc$Immune_Z <- as.vector(scale(pbmc$Immune_Score))

# --- THRESHOLDING & EVALUIERUNG ---
pbmc$GT_Response <- ifelse(pbmc[[GT_COLUMN]] == POSITIVE_VAL, 1, 0)
roc_obj <- roc(response = pbmc$GT_Response, predictor = pbmc$Z_Score, direction = "<", quiet = TRUE)

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
    for(tz in z_grid) {
        pred <- pbmc$Z_Score > tz
        tp <- sum(pred & pbmc$GT_Response == 1); fp <- sum(pred & pbmc$GT_Response == 0)
        fn <- sum(!pred & pbmc$GT_Response == 1); prec <- tp/(tp+fp); rec <- tp/(tp+fn)
        f1 <- 2*(prec*rec)/(prec+rec)
        if(!is.na(f1) && f1 > best_f1) { best_f1 <- f1; THRESHOLD_Z <- tz }
    }
} else if (THRESH_MODE == "null_dist_immune_dual") {
    z_grid <- unique(quantile(pbmc$Z_Score, probs = seq(0.1, 0.9, 0.1)))
    i_grid <- unique(quantile(pbmc$Immune_Z, probs = seq(0.1, 0.9, 0.1)))
    best_f1 <- -1
    for(tz in z_grid) {
        for(ti in i_grid) {
            pred <- (pbmc$Z_Score > tz) & (pbmc$Immune_Z > ti)
            tp <- sum(pred & pbmc$GT_Response == 1); fp <- sum(pred & pbmc$GT_Response == 0)
            fn <- sum(!pred & pbmc$GT_Response == 1); prec <- tp/(tp+fp); rec <- tp/(tp+fn)
            f1 <- 2*(prec*rec)/(prec+rec)
            if(!is.na(f1) && f1 > best_f1) { 
                best_f1 <- f1; THRESHOLD_Z <- tz; THRESHOLD_I <- ti 
            }
        }
    }
} else if (THRESH_MODE == "gmm_dist_platelet" || THRESH_MODE == "gmm_dist_dual"){
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
} else if (THRESH_MODE == "kmeans"){
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
    pbmc$Immune_High <- pbmc$KMeans_Cluster == positive_cluster
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
pbmc$Prediction <- factor(ifelse(positive_condition, "Positive", "Negative"), levels = c("Negative","Positive"))

gt_values <- as.character(pbmc@meta.data[[GT_COLUMN]])
pos_val_char <- as.character(POSITIVE_VAL)

pbmc$Error_Type <- case_when(
    pbmc$Prediction == "Positive" & gt_values == pos_val_char ~ "TP",
    pbmc$Prediction == "Positive" & gt_values != pos_val_char ~ "FP",
    pbmc$Prediction == "Negative" & gt_values == pos_val_char ~ "FN",
    pbmc$Prediction == "Negative" & gt_values != pos_val_char ~ "TN",
    TRUE ~ NA_character_
)

# --- SAVE RDS ---
RDS_FILE_NAME <- paste0(RDS_OUT_DIR, "pbmc_benchmarked_", DATASET_SHORT, "_", METHOD_NAME, "_", THRESH_MODE, "_", EXTRACTED_MODE, "_GT_", GT_SOURCE, ".rds")
message("Speichere Seurat-Objekt unter: ", RDS_FILE_NAME)
saveRDS(pbmc, file = RDS_FILE_NAME)

# Metriken schreiben
runtime_min <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
tp <- sum(pbmc$Error_Type == "TP", na.rm = TRUE); fp <- sum(pbmc$Error_Type == "FP", na.rm = TRUE)
fn <- sum(pbmc$Error_Type == "FN", na.rm = TRUE); tn <- sum(pbmc$Error_Type == "TN", na.rm = TRUE)
prec <- if((tp + fp) > 0) tp / (tp + fp) else 0
rec  <- if((tp + fn) > 0) tp / (tp + fn) else 0
f1   <- if((prec + rec) > 0) 2 * prec * rec / (prec + rec) else 0

write.csv(data.frame(
    Dataset = DATASET_SHORT, Method = METHOD_NAME, Signature = SIG_NAME, Ext = USE_EXTENSION, 
    Filter_Mode = EXTRACTED_MODE, GT_Source = GT_SOURCE, Threshold_Mode = THRESH_MODE,
    AUC = as.numeric(auc(roc_obj)), F1 = f1, Prec = prec, Rec = rec, TP = tp, FP = fp, FN = fn, TN = tn, 
    Threshold_Value = THRESHOLD_Z, Runtime_Min = runtime_min
),
paste0("results/metrics/metrics_", DATASET_SHORT, "_", METHOD_NAME, "_", THRESH_MODE, "_", EXTRACTED_MODE, "_GT_", GT_SOURCE, ".csv"), row.names = FALSE)

ct_data <- pbmc@meta.data %>% 
    group_by(celltype_clean, celltype.l3) %>% 
    summarise(
        TP = sum(Error_Type == "TP", na.rm = TRUE), FP = sum(Error_Type == "FP", na.rm = TRUE),
        FN = sum(Error_Type == "FN", na.rm = TRUE), TN = sum(Error_Type == "TN", na.rm = TRUE),
        Mean_Z = mean(Z_Score, na.rm = TRUE), n = n(), .groups = "drop" 
    ) %>%
    mutate(
        Accuracy = (TP + TN) / n,
        Balanced_Accuracy = 0.5 * ((TP / (TP + FN + 1e-6)) + (TN / (TN + FP + 1e-6))),
        Dataset = DATASET_SHORT, Method = METHOD_NAME, Mode = THRESH_MODE, Filter_Mode = EXTRACTED_MODE, GT_Source = GT_SOURCE
    )

write.csv(ct_data, 
          paste0("results/celltype_data/ct_", DATASET_SHORT, "_", METHOD_NAME, "_", THRESH_MODE, "_", EXTRACTED_MODE, "_GT_", GT_SOURCE, ".csv"), 
          row.names = FALSE)

#CHOSEN_REDUCTION <- "GEX_umap" 
available_reductions <- names(pbmc@reductions)
print(paste("Verfügbare UMAPs/PCAs in dieser Kohorte:", paste(available_reductions, collapse=", ")))

CHOSEN_REDUCTION <- if (DATASET_SHORT == "immune_aging") "GEX_umap_mrvi" else "umap_totalVI"

if (!CHOSEN_REDUCTION %in% names(pbmc@reductions)) {
  warning(paste("Gewünschte Reduktion", CHOSEN_REDUCTION, "fehlt. Weiche auf erste verfügbare aus."))
  CHOSEN_REDUCTION <- names(pbmc@reductions)[1]
}

print(paste("-> Für Visualisierung gewählte RNA-Reduktion:", CHOSEN_REDUCTION))

png(paste0(OUT_DIR, "1_Reference_UMAP.png"), 1200, 800)
print(DimPlot(pbmc, reduction = CHOSEN_REDUCTION, group.by="celltype_clean", label=T, repel=T) + 
      labs(title=paste("Reference Gating (", CHOSEN_REDUCTION, ") -", DATASET_SHORT)))
dev.off()

png(paste0(OUT_DIR, "2_Score_UMAP.png"), 900, 700)
print(FeaturePlot(pbmc, reduction = CHOSEN_REDUCTION, features="Raw_Score") + 
      scale_colour_viridis_c(option="magma") + 
      labs(title="Raw Score Intensity"))
dev.off()

png(paste0(OUT_DIR, "3_Error_Mapping_UMAP.png"), 1000, 800)
print(DimPlot(pbmc, reduction = CHOSEN_REDUCTION, group.by="Error_Type") + 
      scale_color_manual(values=c("TP"="#228B22","FP"="#FF4500","FN"="#1E90FF","TN"="#D3D3D3")) + 
      labs(title="Error Mapping", subtitle=paste("Threshold Z =", round(THRESHOLD_Z, 2))))
dev.off()

png(paste0(OUT_DIR, "4a_Density_ZScore_PLA_Status.png"), 1200, 800)
print(ggplot(pbmc@meta.data, aes(x=Z_Score, fill=!!sym(GT_COLUMN))) + 
      geom_density(alpha=0.5) + theme_minimal() + 
      scale_fill_manual(values=c("PLA"="#FF4B4B", "platelet-free"="#4B8BFF")) +
      geom_vline(xintercept=THRESHOLD_Z, linetype="dashed", color="red", size=1) +
      labs(title="Global Z-Score Distribution", subtitle=paste("Threshold Z:", round(THRESHOLD_Z, 2)), x="Z-Score", fill="PLA Status"))
dev.off()

png(paste0(OUT_DIR, "4b_Density_RawScore_PLA_Status.png"), 1200, 800)
print(ggplot(pbmc@meta.data, aes(x=Raw_Score, fill=!!sym(GT_COLUMN))) + 
      geom_density(alpha=0.5) + theme_minimal() + 
      scale_fill_manual(values=c("PLA"="#FF4B4B", "platelet-free"="#4B8BFF")) + labs(title="Global Raw Score Distribution"))
dev.off()

png(paste0(OUT_DIR, "6a_Violin_RawScore_Celltypes.png"), 1200, 600)
print(VlnPlot(pbmc, features="Raw_Score", group.by="celltype_clean", pt.size=0))
dev.off()

png(paste0(OUT_DIR, "6b_Split_Violin_PLA_Status.png"), 1400, 700)
print(VlnPlot(pbmc, features="Raw_Score", group.by="celltype_clean", split.by=GT_COLUMN, pt.size=0) +
      scale_fill_manual(values=c("PLA"="#FF4B4B", "platelet-free"="#4B8BFF")))
dev.off()

png(paste0(OUT_DIR, "9_ROC_Curve.png"), 700, 700)
plot(roc_obj, col="#E41A1C", lwd=3, main=paste("ROC AUC:", round(auc(roc_obj),3)))
dev.off()

message("Done! Alle Plots und Ergebnisse für ", DATASET_SHORT, " (", EXTRACTED_MODE, ") in: ", OUT_DIR)