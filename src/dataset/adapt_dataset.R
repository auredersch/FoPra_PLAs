library(Seurat)

dirs_to_update <- c(
  #"/nfs/home/students/i.kaciran/FoPra_PLAs/data/datasets",
  "/nfs/home/students/f.mathis/SysBioMed-PLAs/data/datasets"
)

files <- c(
  "gated_ImmuneAging.rds",
  "gated_heart_processed.rds",
  "gated_sepsis_processed.rds",
  "gated_vaccine_processed.rds",
  "gated_our_dataset_processed.rds"
)

for (data_dir in dirs_to_update) {
  message("\nUpdating directory: ", data_dir)
  
  for (file_name in files) {
    file_path <- file.path(data_dir, file_name)
    
    if (!file.exists(file_path)) {
      warning("File does not exist: ", file_path)
      next
    }
    
    message("Processing: ", file_path)
    
    # Backup first
    backup_path <- paste0(file_path, ".backup_before_pla_status_fix")
    
    if (!file.exists(backup_path)) {
      file.copy(file_path, backup_path)
      message("Backup created: ", backup_path)
    } else {
      message("Backup already exists: ", backup_path)
    }
    
    seurat_obj <- readRDS(file_path)
    
    if (!"pla_status" %in% colnames(seurat_obj@meta.data)) {
      warning("No pla_status column in: ", file_path)
      next
    }
    
    n_na_before <- sum(is.na(seurat_obj$pla_status))
    message("NAs before: ", n_na_before)
    
    seurat_obj$pla_status[is.na(seurat_obj$pla_status)] <- "platelet-free"
    
    n_na_after <- sum(is.na(seurat_obj$pla_status))
    message("NAs after: ", n_na_after)
    
    saveRDS(seurat_obj, file_path)
    
    message("Updated: ", file_path)
  }
}

for (data_dir in dirs_to_update) {
  message("\nChecking directory: ", data_dir)
  
  for (file_name in files) {
    file_path <- file.path(data_dir, file_name)
    
    if (!file.exists(file_path)) {
      next
    }
    
    seurat_obj <- readRDS(file_path)
    
    cat("\n", file_name, "\n")
    print(table(seurat_obj$pla_status, useNA = "ifany"))
  }
}