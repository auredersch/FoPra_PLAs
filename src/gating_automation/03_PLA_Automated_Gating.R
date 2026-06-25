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

INPUT_DIR   <- "/nfs/home/students/f.mathis/SysBioMed-PLAs/data/datasets/"
OUTPUT_DIR  <- "/nfs/home/students/a.dersch/FoPra_PLAs/data/datasets_automated/"
BASE_PLOT   <- "/nfs/home/students/a.dersch/FoPra_PLAs/results/gating_automation/"
QC_BASE     <- "/nfs/home/students/a.dersch/FoPra_PLAs/results/sample_qc/"

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# --- DATASET MAPPING ---
datasets_map = list(
  "gated_heart_processed.rds"        = "heart",
  "gated_sepsis_processed.rds"       = "sepsis",
  "gated_vaccine_processed.rds"      = "vaccine",
  "gated_ImmuneAging.rds"            = "immune_aging",
  "gated_our_dataset_processed.rds"  = "impact"
)

dataset_type <- datasets_map[[CURRENT_FILE]]
if(is.null(dataset_type)) stop("Fehler: Datensatz-Dateiname nicht in datasets_map gefunden!")

all_lineages <- c("B cells", "DCs", "NK cells", "Unassigned", "CD4 T", "Macrophages", "Neutrophils", "CD8 T cells")

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
DefaultAssay(seurat_obj) <- "ADT"
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
  qc_table_file <- file.path(QC_BASE, dataset_type, "19_final_sample_lineage_QC_table.csv")
  
  if (file.exists(qc_table_file)) {
    print(paste("-> Integriere Proben-QC-Tabelle im Modus:", FILTER_MODE))
    qc <- read.csv(qc_table_file)
    
    if (FILTER_MODE == "qc_strict") {
      # Nur die absoluten Top-Paare zulassen (Alles andere fliegt raus)
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

# --- GATING SCHLEIFE ---
for (lineage_name in all_lineages) {

  csv_lineage_name <- case_when(
    lineage_name == "CD4 T" ~ "CD4 T cells",
    lineage_name == "CD8 T cells" ~ "CD8 T cells",
    TRUE ~ lineage_name
  )
  
  if (!(lineage_name %in% seurat_obj$lineage)) next
  
  Zellen_im_Zelltyp <- which(seurat_obj$lineage == lineage_name)
  if (length(Zellen_im_Zelltyp) < 50) next
  
  if (!is.null(trusted_pairs)) {
    actual_sample_col <- case_when(
      "donor_id" %in% colnames(seurat_obj@meta.data) ~ "donor_id",
      "sample_id" %in% colnames(seurat_obj@meta.data) ~ "sample_id",
      TRUE ~ "sample"
    )
    
    cell_samples <- seurat_obj@meta.data[[actual_sample_col]][Zellen_im_Zelltyp] 
    
    valid_samples <- trusted_pairs %>% 
      filter(celltype_id == csv_lineage_name) %>% 
      pull(sample_id)
    
    Zellen_im_Zelltyp <- Zellen_im_Zelltyp[cell_samples %in% valid_samples]
    
    if (length(Zellen_im_Zelltyp) < 20) {
      print(paste("   -> [QC Filter] Überspringe Lineage", lineage_name, "- Zu wenige verlässliche Zellen übrig."))
      next
    }
  }
  
  # --- 1D GMM (CD41) ---
  cd41_values <- FetchData(seurat_obj, vars = cd41_name)[Zellen_im_Zelltyp, 1]
  gmm_cd41 <- tryCatch(Mclust(cd41_values, G = 2), error = function(e) NULL)
  
  thresh_cd41 <- NULL
  if (!is.null(gmm_cd41)) {
    neg_cluster_cd41 <- which.min(gmm_cd41$parameters$mean)
    thresh_cd41 <- max(cd41_values[gmm_cd41$classification == neg_cluster_cd41])
    seurat_obj$automative_gating_single[Zellen_im_Zelltyp] <- ifelse(cd41_values > thresh_cd41, "PLA", "platelet-free")
  }
  
  # --- 2D GMM (CD41 + Zweitmarker) ---
  if (!is.null(second_marker_name) && !is.null(thresh_cd41)) {
    second_values <- FetchData(seurat_obj, vars = second_marker_name)[Zellen_im_Zelltyp, 1]
    gmm_second <- tryCatch(Mclust(second_values, G = 2), error = function(e) NULL)
    
    if (!is.null(gmm_second)) {
      neg_cluster_second <- which.min(gmm_second$parameters$mean)
      thresh_second <- max(second_values[gmm_second$classification == neg_cluster_second])
      
      seurat_obj$automative_gating_double[Zellen_im_Zelltyp] <- ifelse(
        (cd41_values > thresh_cd41) & (second_values > thresh_second), "PLA", "platelet-free"
      )
    }
  } else {
    seurat_obj$automative_gating_double[Zellen_im_Zelltyp] <- seurat_obj$automative_gating_single[Zellen_im_Zelltyp]
  }
}

# --- RDS SPEICHERN ---
neuer_name <- paste0(dataset_type, "_", FILTER_MODE, "_automated_gating.rds")
saveRDS(seurat_obj, file = file.path(OUTPUT_DIR, neuer_name))
print(paste("-> RDS erfolgreich gespeichert unter:", file.path(OUTPUT_DIR, neuer_name)))

# --- PLOTS SCHREIBEN ---
print("-> Generiere Validierungs-Plots...")
plot_data <- seurat_obj@meta.data %>% filter(lineage %in% all_lineages)

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

reduction <- "ADT_umap"
if (reduction %in% names(seurat_obj@reductions)) {
  Key(seurat_obj[[reduction]]) <- "adtumap_"
  old_assay <- DefaultAssay(seurat_obj)
  if ("RNA" %in% Assays(seurat_obj)) DefaultAssay(seurat_obj) <- "RNA"
  
  png(paste0(OUT_DIR, "1_Reference_Lineage_UMAP.png"), 1000, 700)
  print(DimPlot(seurat_obj, group.by="lineage", reduction = reduction, label=TRUE, repel=TRUE) + labs(title=paste("Reference Lineage -", dataset_type)))
  dev.off()
  
  png(paste0(OUT_DIR, "2_Biologist_Gating_UMAP.png"), 1000, 700)
  print(DimPlot(seurat_obj, group.by="pla_status", reduction = reduction) + labs(title=paste("Biologist Gating -", dataset_type)))
  dev.off()
  
  png(paste0(OUT_DIR, "3_GMM_Single_Gating_UMAP.png"), 1000, 700)
  print(DimPlot(seurat_obj, group.by="automative_gating_single", reduction = reduction) + labs(title=paste("1D GMM (CD41) -", dataset_type)))
  dev.off()
  
  png(paste0(OUT_DIR, "4_GMM_Double_Gating_UMAP.png"), 1000, 700)
  print(DimPlot(seurat_obj, group.by="automative_gating_double", reduction = reduction) + labs(title=paste("2D GMM (CD41 +", second_marker_label, ") -", dataset_type)))
  dev.off()
  
  DefaultAssay(seurat_obj) <- old_assay
}

print(paste("Fertig mit Datensatz:", dataset_type))