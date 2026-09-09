# ===================================================================
# 03_PLA_Automated_Gating.R
# ===================================================================

library(Seurat)
library(tidyverse)
library(mclust)
library(ggplot2)
library(patchwork)

# --- CLI ARGUMENTE ---
args <- commandArgs(trailingOnly = TRUE)
CURRENT_FILE <- if(length(args) >= 1) args[1] else "gated_sepsis_processed.rds"
FILTER_MODE  <- if(length(args) >= 2) args[2] else "raw" # Optionen: "raw", "qc_tolerant", "qc_strict"

INPUT_DIR   <- "/nfs/home/students/f.mathis/Dataset_PostQC/"
#INPUT_DIR   <-  "/nfs/home/students/a.dersch/data"
OUTPUT_DIR  <- "/nfs/home/students/a.dersch/FoPra_PLAs/data/datasets_automated_postQC/"
BASE_PLOT   <- "/nfs/home/students/a.dersch/FoPra_PLAs/results/gating_automation_postQC/"
QC_BASE     <- "/nfs/home/students/f.mathis/Dataset_PostQC/"

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# --- DATASET MAPPING ---
#datasets_map = list(
#  "gated_heart_processed.rds"        = "heart",
#  "gated_sepsis_processed.rds"       = "sepsis",
#  "gated_vaccine_processed.rds"      = "vaccine",
#  "gated_ImmuneAging.rds"            = "immune_aging",
#  "gated_our_dataset_processed.rds"  = "impact",
#  "gated_skin_processed.rds"          = "skin"
#)

datasets_map = list(
  "heart.rds"        = "heart",
  "sepsis.rds"       = "sepsis",
  "vaccine.rds"      = "vaccine",
  "immune_aging.rds" = "immune_aging",
  "impact.rds"       = "impact",
  "skin.rds"         = "skin"
)

dataset_type <- datasets_map[[CURRENT_FILE]]
if(is.null(dataset_type)) stop("Fehler: Datensatz-Dateiname nicht in datasets_map gefunden!")

#all_lineages <- c("B cells", "DCs", "NK cells", "Unassigned", "CD4 T", "Macrophages", "Neutrophils", "CD8 T cells")

lineage_aliases <- list(
  "B cells"     = c("B cells"),
  "DCs"         = c("DCs"),
  "NK cells"    = c("NK cells"),
  "Unassigned"  = c("Unassigned"),
  "CD4 T cells" = c("CD4 T cells", "CD4 T"),
  "CD8 T cells" = c("CD8 T cells"),
  "Monocytes"   = c("Monocytes"),
  "Macrophages" = c("Macrophages"),
  "Neutrophils" = c("Neutrophils", "neutrophils")
)

marker_aliases <- list(
  CD41  = c("CD41", "ITGA2B.1", "ITGA2B", "GPIIB"),
  CD61  = c("CD61", "ITGB3.1", "ITGB3", "GPIIIA"),
  CD62P = c("CD62P", "SELP.1", "SELP", "PSELECTIN", "CD62")
)

find_marker_name <- function(seurat_features, alias_list) {
  for (alias in alias_list) {
    query <- grep(paste0("(^|[-_])", alias, "(\\.|$)"), seurat_features, value = TRUE, ignore.case = TRUE)
    if (length(query) > 0) return(query[1])
  }
  return(NULL)
}

print(paste("====================================================="))
print(paste("Datensatz:", dataset_type, "(", CURRENT_FILE, ")"))
print(paste("Filter-Modus:", FILTER_MODE))
print(paste("====================================================="))

seurat_obj <- readRDS(file.path(INPUT_DIR, CURRENT_FILE))

handled_lineages <- unlist(lineage_aliases)
unhandled <- setdiff(unique(na.omit(seurat_obj$lineage)), handled_lineages)
if (length(unhandled) > 0) {
  warning(paste("Nicht gegatete Lineage-Level gefunden:", paste(unhandled, collapse = ", ")))
}

adtnorm_assay_name <- paste0("ADT_ADTnorm_", dataset_type)

if (adtnorm_assay_name %in% Assays(seurat_obj)) {
  DefaultAssay(seurat_obj) <- adtnorm_assay_name
  print(paste("-> Nutze ADTNorm-korrigierten Assay:", adtnorm_assay_name))
} else if ("ADT_corrected" %in% Assays(seurat_obj)) {
  DefaultAssay(seurat_obj) <- "ADT_corrected"
  print("-> Nutze 'ADT_corrected' (Fallback, alte immune_aging-Konvention).")
} else {
  warning(paste0("-> WARNUNG: Kein ADTNorm-korrigierter Assay gefunden fuer '", dataset_type,
                  "' (erwartet: '", adtnorm_assay_name, "' oder 'ADT_corrected'). ",
                  "Verfuegbare Assays: ", paste(Assays(seurat_obj), collapse = ", "),
                  ". Falle zurueck auf UNKORRIGIERTEN 'ADT'-Assay - Ergebnisse ggf. ",
                  "nicht vergleichbar mit ADTNorm-korrigierten Kohorten!"))
  DefaultAssay(seurat_obj) <- "ADT"
}

all_features <- rownames(seurat_obj)

seurat_obj$automative_gating_single <- "platelet-free"
seurat_obj$automative_gating_double <- "platelet-free"

cd41_name <- find_marker_name(all_features, marker_aliases$CD41)
cd61_name <- find_marker_name(all_features, marker_aliases$CD61)
cd62p_name <- find_marker_name(all_features, marker_aliases$CD62P)

if (is.null(cd41_name)) {
  stop(paste("Kritischer Fehler: Kein CD41-Marker in", dataset_type, "gefunden!"))
}

second_marker_name <- NULL
second_marker_label <- "Missing"

if (!is.null(cd61_name)) {
  second_marker_name <- cd61_name
  second_marker_label <- "CD61"
  print("-> Nutze CD41 + CD61 für das Dual Gating.")
} else if (!is.null(cd62p_name)) {
  second_marker_name <- cd62p_name
  second_marker_label <- "CD62P"
  print("-> CD61 fehlt. Nutze stattdessen CD41 + CD62P für das Dual Gating.")
} else {
  warning("-> Weder CD61 noch CD62P gefunden! Dual Gating wird auf Single Gating beschränkt.")
}

trusted_pairs <- NULL

if (FILTER_MODE != "raw") {
  #qc_table_file <- file.path(QC_BASE, dataset_type, "19_final_sample_lineage_QC_table.csv")
  qc_table_file <- file.path(QC_BASE, paste0(dataset_type, ".csv"))
  
  if (file.exists(qc_table_file)) {
    print(paste("-> Integriere Proben-QC-Tabelle im Modus:", FILTER_MODE))
    qc <- read.csv(qc_table_file)
    
    if (FILTER_MODE == "qc_strict") {
      trusted_pairs <- qc %>% 
        filter(final_pair_status %in% c("trusted", "trusted_but_extreme_PLA")) %>%
        select(sample_id, celltype_id)
        
    } else if (FILTER_MODE == "qc_tolerant") {
    
      trusted_pairs <- qc %>% 
        filter(final_pair_status %in% c("trusted", "trusted_but_extreme_PLA", 
                                        "sample_ADT_suspicious", "lineage_ADT_suspicious")) %>%
        select(sample_id, celltype_id)
    }
  } else {
    stop(paste("-> Fehler: Keine QC-Tabelle für", FILTER_MODE, "gefunden! Pipeline bricht ab."))
  }
}

OUT_DIR <- file.path(BASE_PLOT, dataset_type, FILTER_MODE, "/")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

get_gmm_threshold <- function(gmm_fit) {
  means <- gmm_fit$parameters$mean
  vars  <- gmm_fit$parameters$variance$sigmasq
  pros  <- gmm_fit$parameters$pro
  
  neg_idx <- which.min(means)
  pos_idx <- which.max(means)
  
  # Wenn means identisch oder extrem nah aneinander ist
  if (abs(means[neg_idx] - means[pos_idx]) < 1e-4) {
    return(means[neg_idx])
  }
  
  #Kreuzungspunkt zwischen means
  search_grid <- seq(means[neg_idx], means[pos_idx], length.out = 1000)
  
  var_neg <- if(length(vars) >= neg_idx) vars[neg_idx] else vars[1]
  var_pos <- if(length(vars) >= pos_idx) vars[pos_idx] else vars[1]
  
  d_neg <- pros[neg_idx] * dnorm(search_grid, means[neg_idx], sqrt(var_neg))
  d_pos <- pros[pos_idx] * dnorm(search_grid, means[pos_idx], sqrt(var_pos))
  
  cross_idx <- which.min(abs(d_neg - d_pos))
  thresh <- search_grid[cross_idx]
  
  # Falls thresh Länge 0 oder NA ist
  if (length(thresh) == 0 || is.na(thresh)) {
    # arithmetischer Mittelpunkt 
    thresh <- (means[neg_idx] + means[pos_idx]) / 2
  }
  
  return(thresh)
}

actual_sample_col <- NULL
if (!is.null(trusted_pairs)) {
  candidate_cols <- grep("sample|donor", colnames(seurat_obj@meta.data), value = TRUE, ignore.case = TRUE)

  if (length(candidate_cols) == 0) {
    stop(paste0("Kritischer Fehler: Keine Spalte mit 'sample' oder 'donor' im Namen ",
                "in '", dataset_type, "' gefunden. Verfuegbare Spalten: ",
                paste(colnames(seurat_obj@meta.data), collapse = ", ")))
  }

  match_counts <- vapply(candidate_cols, function(col) {
    length(intersect(unique(as.character(seurat_obj@meta.data[[col]])), trusted_pairs$sample_id))
  }, FUN.VALUE = integer(1))

  actual_sample_col <- candidate_cols[which.max(match_counts)]
  print(paste("-> Sample-ID-Spalte automatisch gewaehlt:", actual_sample_col,
              "(", max(match_counts), "uebereinstimmende Sample-IDs von",
              length(unique(trusted_pairs$sample_id)), "in der QC-Tabelle, geprueft gegen:",
              paste(candidate_cols, collapse=", "), ")"))
  if (max(match_counts) == 0) {
    warning("Keine der gefundenen Sample/Donor-Spalten hat Overlap mit trusted_pairs$sample_id!")
  }
}
# Gating
for (canonical_name in names(lineage_aliases)) {

  aliases <- lineage_aliases[[canonical_name]]
  csv_lineage_name <- canonical_name 

  match_idx <- which(seurat_obj$lineage %in% aliases)
  if (length(match_idx) < 50) next

  Zellen_im_Zelltyp <- match_idx
  
  if (!is.null(trusted_pairs)) {
    
    cell_samples <- seurat_obj@meta.data[[actual_sample_col]][Zellen_im_Zelltyp]
    valid_samples <- trusted_pairs %>% 
      filter(celltype_id == csv_lineage_name) %>% 
      pull(sample_id)
    
    Zellen_im_Zelltyp <- Zellen_im_Zelltyp[cell_samples %in% valid_samples]
    
    if (length(Zellen_im_Zelltyp) < 20) {
      print(paste("   -> [QC Filter] Überspringe Lineage", canonical_name, "- Zu wenige verlässliche Zellen übrig."))
      next
    }
  }
  
  # --- 1D GMM (CD41) ---
  cd41_values <- FetchData(seurat_obj, vars = cd41_name)[Zellen_im_Zelltyp, 1]
  gmm_cd41 <- tryCatch(Mclust(cd41_values, G = 2), error = function(e) NULL)
  
  thresh_cd41 <- NULL
  if (!is.null(gmm_cd41)) {
    #neg_cluster_cd41 <- which.min(gmm_cd41$parameters$mean)
    #thresh_cd41 <- max(cd41_values[gmm_cd41$classification == neg_cluster_cd41])
    thresh_cd41 <- get_gmm_threshold(gmm_cd41)
    seurat_obj$automative_gating_single[Zellen_im_Zelltyp] <- ifelse(cd41_values > thresh_cd41, "PLA", "platelet-free")
  }
  
  # --- 2D GMM (CD41 + Zweitmarker) ---
  if (!is.null(second_marker_name) && !is.null(thresh_cd41)) {
    second_values <- FetchData(seurat_obj, vars = second_marker_name)[Zellen_im_Zelltyp, 1]
    gmm_second <- tryCatch(Mclust(second_values, G = 2), error = function(e) NULL)
    
    if (!is.null(gmm_second)) {
      #neg_cluster_second <- which.min(gmm_second$parameters$mean)
      #thresh_second <- max(second_values[gmm_second$classification == neg_cluster_second])
      thresh_second <- get_gmm_threshold(gmm_second)

      seurat_obj$automative_gating_double[Zellen_im_Zelltyp] <- ifelse(
        (cd41_values > thresh_cd41) & (second_values > thresh_second), "PLA", "platelet-free"
      )
    }
  } else {
    seurat_obj$automative_gating_double[Zellen_im_Zelltyp] <- seurat_obj$automative_gating_single[Zellen_im_Zelltyp]
  }
}

neuer_name <- paste0(dataset_type, "_", FILTER_MODE, "_automated_gating.rds")
saveRDS(seurat_obj, file = file.path(OUTPUT_DIR, neuer_name))
print(paste("-> RDS erfolgreich gespeichert unter:", file.path(OUTPUT_DIR, neuer_name)))

# Plots
print("-> Generiere Validierungs-Plots...")
plot_data <- seurat_obj@meta.data %>% filter(lineage %in% unlist(lineage_aliases))

p1 <- ggplot(plot_data, aes(x = lineage, fill = pla_status)) +
  geom_bar(position = "fill") + theme_minimal() +
  labs(title = "1. Biologist (pla_status)", x = "Lineage", y = "Fraction", fill = "PLA") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p2 <- ggplot(plot_data, aes(x = lineage, fill = automative_gating_single)) +
  geom_bar(position = "fill") + theme_minimal() +
  labs(title = "2. GMM (Only CD41)", x = "Lineage", y = "Fraction", fill = "PLA") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p3 <- ggplot(plot_data, aes(x = lineage, fill = automative_gating_double)) +
  geom_bar(position = "fill") + theme_minimal() +
  labs(title = paste("3. GMM Improved (CD41 +", second_marker_label, ")"), x = "Lineage", y = "Fraction", fill = "PLA") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(OUT_DIR, "barplot_fraction_biologist.png"), plot = p1, width = 6, height = 5)
ggsave(paste0(OUT_DIR, "barplot_fraction_gmm_1d.png"), plot = p2, width = 6, height = 5)
ggsave(paste0(OUT_DIR, "barplot_fraction_gmm_2d.png"), plot = p3, width = 6, height = 5)

comp_p1 <- p1 + theme(axis.title.x = element_blank()) + scale_fill_manual(values = c("PLA" = "#FF4B4B", "platelet-free" = "#4B8BFF"), guide = "none")
comp_p2 <- p2 + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + scale_fill_manual(values = c("PLA" = "#FF4B4B", "platelet-free" = "#4B8BFF"), guide = "none")
comp_p3 <- p3 + theme(axis.title.y = element_blank()) + scale_fill_manual(values = c("PLA" = "#FF4B4B", "platelet-free" = "#4B8BFF"), guide = "none")

composition_plot <- (comp_p1 | comp_p2 | comp_p3) + 
  plot_annotation(
    title = paste("PLA Gating Composition Comparison -", toupper(dataset_type)),
    subtitle = paste("Filter Mode:", toupper(FILTER_MODE)),
    theme = theme(plot.title = element_text(size = 16, face = "bold"), plot.subtitle = element_text(size = 12))
  )
ggsave(
  filename = paste0(OUT_DIR, "combined_composition_plot.png"), 
  plot = composition_plot, 
  width = 16, 
  height = 6, 
  dpi = 300
)  

reduction <- if (dataset_type == "immune_aging") "GEX_umap_mrvi" else "umap_totalVI"

if (reduction %in% names(seurat_obj@reductions)) {
  old_assay <- DefaultAssay(seurat_obj)
  if ("RNA" %in% Assays(seurat_obj)) DefaultAssay(seurat_obj) <- "RNA"
  
  png(paste0(OUT_DIR, "1_Reference_Lineage_UMAP.png"), 1000, 700)
  print(DimPlot(seurat_obj, group.by="lineage", reduction = reduction, label=TRUE, repel=TRUE) + labs(title=paste("Reference Lineage (", reduction, ") -", dataset_type)))
  dev.off()
  
  png(paste0(OUT_DIR, "2_Biologist_Gating_UMAP.png"), 1000, 700)
  print(DimPlot(seurat_obj, group.by="pla_status", reduction = reduction) + labs(title=paste("Biologist Gating (", reduction, ") -", dataset_type)))
  dev.off()
  
  png(paste0(OUT_DIR, "3_GMM_Single_Gating_UMAP.png"), 1000, 700)
  print(DimPlot(seurat_obj, group.by="automative_gating_single", reduction = reduction) + labs(title=paste("1D GMM CD41 (", reduction, ") -", dataset_type)))
  dev.off()
  
  png(paste0(OUT_DIR, "4_GMM_Double_Gating_UMAP.png"), 1000, 700)
  print(DimPlot(seurat_obj, group.by="automative_gating_double", reduction = reduction) + labs(title=paste("2D GMM CD41+", second_marker_label, " (", reduction, ") -", dataset_type)))
  dev.off()
  
  DefaultAssay(seurat_obj) <- old_assay
} else {
  warning(paste("Warnung: Gewünschte Reduktion", reduction, "nicht im Objekt gefunden! Überspringe UMAP-Plots."))
}

print(paste("Fertig mit Datensatz:", dataset_type))