start_time <- Sys.time() 

library(Seurat)
library(AUCell)
library(dplyr)
library(tidyr)
library(mclust)

project_root <- "/nfs/home/students/f.mathis/FoPra_PLAs"

args <- commandArgs(trailingOnly = TRUE)
METHOD_NAME    <- if(length(args) >= 1) args[1] else "AUCell"
SIG_NAME       <- if(length(args) >= 2) args[2] else "GOBP_REG"
SIG_FILE_BASE  <- if(length(args) >= 3) args[3] else "GOBP_REGULATION_OF_PLATELET_ACTIVATION.v2025.1.Hs.csv"
USE_EXTENSION  <- if(length(args) >= 4) as.logical(args[4]) else FALSE
THRESH_MODE    <- if(length(args) >= 5) args[5] else "gmm_dist_dual" 
CURRENT_FILE   <- if(length(args) >= 6) args[6] else "/nfs/home/students/f.mathis/FoPra_PLAs/data/stemi_raw" 
SAVE_RDS       <- if(length(args) >= 7) as.logical(args[7]) else TRUE
THRESH_SCOPE   <- if(length(args) >= 8) args[8] else "global"       # "global" oder "per_celltype"

set.seed(42)

# Argument 6: Pfad zum Eingabe-RDS; 7: SAVE_RDS; 8: THRESH_SCOPE.
PATH_DATA <- CURRENT_FILE
DATASET_SHORT <- tools::file_path_sans_ext(basename(CURRENT_FILE))
if (!file.exists(PATH_DATA)) stop("Eingabedatei fehlt: ", PATH_DATA)
if (is.na(USE_EXTENSION) || is.na(SAVE_RDS)) stop("Extension und SAVE_RDS müssen TRUE/FALSE sein.")
if (!METHOD_NAME %in% c("AUCell", "UCell", "AddModuleScore")) stop("Unbekannte oder nicht implementierte Scoring-Methode.")
if (!THRESH_MODE %in% c("percentile", "manual", "gmm_dist_platelet", "gmm_dist_dual", "kmeans")) stop("Ungültiger Threshold-Modus ohne GT.")
if (!THRESH_SCOPE %in% c("global", "per_celltype")) stop("Ungültiger THRESH_SCOPE.")
if (THRESH_SCOPE == "per_celltype" && !grepl("^gmm_dist_", THRESH_MODE)) stop("per_celltype ist nur für GMM implementiert.")

PLOT_BASE   <- "/nfs/home/students/f.mathis/FoPra_PLAs/results/benchmarking/"
OUT_DIR <- file.path(
  PLOT_BASE, 
  DATASET_SHORT, 
  paste0(METHOD_NAME, "_Ext", USE_EXTENSION, "_", THRESH_MODE, "_Scope", THRESH_SCOPE),
  SIG_NAME, 
  "/"
)

RDS_OUT_DIR <- file.path(
  project_root,
  "data/datasets/benchmarked_objects",
  DATASET_SHORT,
  paste0(METHOD_NAME, "_Ext", USE_EXTENSION, "_", THRESH_MODE, "_Scope", THRESH_SCOPE),
  SIG_NAME,
  "/"
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RDS_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

print(paste("Lade Datensatz:", CURRENT_FILE, "als", DATASET_SHORT))

counts <- Read10X(
  data.dir = CURRENT_FILE,
  gene.column = 1
)

metadata <- read.delim(
  file.path(CURRENT_FILE, "metadata.tsv.gz"),
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

metadata <- metadata[colnames(counts), , drop = FALSE]

pbmc <- CreateSeuratObject(
  counts = counts,
  project = "stemi",
  meta.data = metadata
)

pbmc$celltype_clean <- "cell_type_lowerres"

# --- GENLISTE LADEN ---
source(file.path(project_root, "src", "benchmarking", "read_and_extend_gene_list.R"))
PATH_SIG <- file.path(project_root, "data", "signatures", paste0(SIG_FILE_BASE, ".csv"))
genes <- read_gene_list(PATH_SIG)

# Overlap-Check: bricht ab wenn zu wenige Gene matchen
gene_overlap <- intersect(genes, rownames(pbmc[["RNA"]]))
overlap_frac <- length(gene_overlap) / length(genes)
print(paste("Gen-Overlap Signatur <-> RNA-Assay:", length(gene_overlap), "von", length(genes), "Genen",
            "(", round(100 * overlap_frac, 1), "% )"))


MIN_ABSOLUTE_GENES   <- 3
#MIN_ABSOLUTE_GENES   <- 1    
MIN_OVERLAP_FRACTION <- 0.5  
#MIN_OVERLAP_FRACTION <- 0.25 

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
                 "von Natur aus kleine Signatur."))
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

if (METHOD_NAME == "AUCell") {
    expression_matrix <- GetAssayData(pbmc, assay = "RNA", layer = "counts")
    rankings <- AUCell_buildRankings(expression_matrix, plotStats = FALSE)
    
    auc_orig <- AUCell_calcAUC(list(Platelet_Orig = genes), rankings)
    pbmc$Raw_Score_Original <- as.numeric(getAUC(auc_orig)[1, ])

    auc_imm <- AUCell_calcAUC(list(Immune_Score = immune_genes), rankings)
    pbmc$Immune_Score <- as.numeric(getAUC(auc_imm)[1, ])

} else if (METHOD_NAME == "UCell") {
    pbmc <- UCell::AddModuleScore_UCell(pbmc, features = list(Platelet_Orig = genes), name = NULL)
    pbmc$Raw_Score_Original <- pbmc$Platelet_Orig
    pbmc <- UCell::AddModuleScore_UCell(pbmc, features = list(Immune_Score = immune_genes), name = NULL)

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
if (METHOD_NAME == "AUCell") {
    pbmc$Raw_Score <- as.numeric(getAUC(AUCell_calcAUC(list(Platelet_Score = final_genes), rankings))[1, ])
} else if (METHOD_NAME == "UCell") {
    pbmc <- UCell::AddModuleScore_UCell(pbmc, features = list(Platelet_Score = final_genes), name = NULL)
    pbmc$Raw_Score <- pbmc$Platelet_Score
} else {
    pbmc <- AddModuleScore(pbmc, features = list(final_genes), name = "AMS")
    pbmc$Raw_Score <- pbmc$AMS1
}

pbmc$Z_Score   <- as.vector(scale(pbmc$Raw_Score))
pbmc$Immune_Z  <- as.vector(scale(pbmc$Immune_Score))

n_na_z <- sum(is.na(pbmc$Z_Score))
if (n_na_z > 0) {
  warning(paste("-> WARNUNG:", n_na_z, "Zellen haben NA als Z_Score."))
}

# --- THRESHOLDING ---
THRESHOLD_I <- -Inf

if (THRESH_MODE == "percentile") {
    prob_cutoff <- 0.90
    THRESHOLD_Z <- as.numeric(quantile(pbmc$Z_Score, probs = prob_cutoff, na.rm = TRUE))

} else if (THRESH_MODE == "manual") {
    med <- median(pbmc$Z_Score, na.rm = TRUE)
    mad_val <- mad(pbmc$Z_Score, na.rm = TRUE)
    THRESHOLD_Z <- med + (1.5 * mad_val)

} else if (THRESH_MODE %in% c("gmm_dist_platelet", "gmm_dist_dual")) {

  pbmc$Platelet_High <- NA
  pbmc$Immune_High   <- NA

  fit_gmm_pair <- function(values) {
    if (length(values) < 3 || any(!is.finite(values)) || length(unique(values)) < 2) return(NULL)
    tryCatch(Mclust(values, G = 2), error = function(e) NULL)
  }

  if (THRESH_SCOPE == "global") {
    fit_plat <- fit_gmm_pair(pbmc$Raw_Score)
    if (is.null(fit_plat)) stop("Mclust (Platelet, global) fehlgeschlagen!")
    plat_high <- which.max(fit_plat$parameters$mean)
    pbmc$Platelet_High <- fit_plat$classification == plat_high

    if (THRESH_MODE == "gmm_dist_dual") {
      idx <- which(pbmc$Platelet_High)
      fit_imm <- fit_gmm_pair(pbmc$Immune_Score[idx])
      if (!is.null(fit_imm)) {
        imm_high <- which.max(fit_imm$parameters$mean)
        pbmc$Immune_High[idx] <- fit_imm$classification == imm_high
      }
    } else {
      pbmc$Immune_High <- TRUE
    }

    THRESHOLD_Z <- min(pbmc$Z_Score[which(pbmc$Platelet_High)], na.rm = TRUE)
    THRESHOLD_I <- if (THRESH_MODE == "gmm_dist_dual") {
      min(pbmc$Immune_Z[which(pbmc$Immune_High)], na.rm = TRUE)
    } else { -Inf }

    # Zusammenfassung der globalen Klassifikation
    threshold_summary <- data.frame(celltype_clean = "GLOBAL",
                                     Threshold_Z = THRESHOLD_Z,
                                     Threshold_I = THRESHOLD_I,
                                     n_cells = ncol(pbmc))

  } else if (THRESH_SCOPE == "per_celltype") {
    MIN_N_PER_CELLTYPE <- 100   
    threshold_rows <- list()

    for (ct in unique(na.omit(pbmc$celltype_clean))) {
      idx_ct <- which(pbmc$celltype_clean == ct)
      if (length(idx_ct) < MIN_N_PER_CELLTYPE) {
        warning(paste0("Zelltyp '", ct, "' hat nur ", length(idx_ct),
                        " Zellen (< ", MIN_N_PER_CELLTYPE, ") - Threshold-Fit übersprungen, ",
                        "Zellen bleiben unklassifiziert (NA)."))
        next
      }

      fit_plat <- fit_gmm_pair(pbmc$Raw_Score[idx_ct])
      if (is.null(fit_plat)) {
        warning(paste0("Mclust (Platelet) fehlgeschlagen für Zelltyp '", ct, "' - übersprungen."))
        next
      }
      plat_high <- which.max(fit_plat$parameters$mean)
      plat_high_idx <- idx_ct[fit_plat$classification == plat_high]
      pbmc$Platelet_High[idx_ct] <- fit_plat$classification == plat_high
      thresh_z_ct <- min(pbmc$Z_Score[plat_high_idx], na.rm = TRUE)

      thresh_i_ct <- -Inf
      if (THRESH_MODE == "gmm_dist_dual" && length(plat_high_idx) >= MIN_N_PER_CELLTYPE) {
        fit_imm <- fit_gmm_pair(pbmc$Immune_Score[plat_high_idx])
        if (!is.null(fit_imm)) {
          imm_high <- which.max(fit_imm$parameters$mean)
          imm_high_idx <- plat_high_idx[fit_imm$classification == imm_high]
          pbmc$Immune_High[plat_high_idx] <- fit_imm$classification == imm_high
          thresh_i_ct <- min(pbmc$Immune_Z[imm_high_idx], na.rm = TRUE)
        }
      } else if (THRESH_MODE == "gmm_dist_platelet") {
        pbmc$Immune_High[idx_ct] <- TRUE
      }

      threshold_rows[[ct]] <- data.frame(celltype_clean = ct,
                                          Threshold_Z = thresh_z_ct,
                                          Threshold_I = thresh_i_ct,
                                          n_cells = length(idx_ct))
    }

    threshold_summary <- if (length(threshold_rows)) bind_rows(threshold_rows) else
      data.frame(celltype_clean = character(), Threshold_Z = numeric(), Threshold_I = numeric(), n_cells = integer())

    THRESHOLD_Z <- median(threshold_summary$Threshold_Z, na.rm = TRUE)
    THRESHOLD_I <- median(threshold_summary$Threshold_I[is.finite(threshold_summary$Threshold_I)], na.rm = TRUE)

  } else {
    stop(paste("Ungueltiger THRESH_SCOPE:", THRESH_SCOPE, "- erlaubt: 'global', 'per_celltype'"))
  }

  THRESH_OUT_DIR <- OUT_DIR
  dir.create(THRESH_OUT_DIR, recursive = TRUE, showWarnings = FALSE)
  write.csv(threshold_summary,
            file.path(THRESH_OUT_DIR, paste0("thresholds_", DATASET_SHORT, "_", SIG_NAME, "_",
                                              THRESH_MODE, "_", THRESH_SCOPE, ".csv")),
            row.names = FALSE)
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
    pbmc@meta.data$Z_Score > THRESHOLD_Z
}

positive_condition <- as.logical(as.vector(positive_condition))
pbmc$Prediction <- factor(ifelse(positive_condition, "PLA", "platelet-free"),
                          levels = c("platelet-free", "PLA"))


if (anyNA(pbmc$Prediction)) warning(sum(is.na(pbmc$Prediction)), " Zellen nicht klassifiziert (NA).")
write.csv(data.frame(cell_barcode = colnames(pbmc), Raw_Score = pbmc$Raw_Score,
                     Prediction = pbmc$Prediction), file.path(OUT_DIR, "predictions.csv"), row.names = FALSE)

# --- SAVE RDS (optional, 7. CLI-Argument) ---
if (SAVE_RDS) {
  RDS_FILE_NAME <- file.path(RDS_OUT_DIR,
    paste0("pbmc_classified_", DATASET_SHORT, "_", METHOD_NAME, "_", SIG_NAME, "_",
           THRESH_MODE, "_Scope", THRESH_SCOPE, ".rds"))
  saveRDS(pbmc, file = RDS_FILE_NAME)
  message("Seurat-Objekt gespeichert: ", RDS_FILE_NAME)
}
