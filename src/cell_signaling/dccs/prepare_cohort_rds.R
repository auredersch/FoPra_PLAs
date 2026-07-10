#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("Usage: Rscript prepare_cohort_rds.R <input_rds> <output_rds> <dataset_type> <mode>")
}

input_rds <- args[[1]]
output_rds <- args[[2]]
dataset_type <- args[[3]]
mode <- args[[4]]

if (!mode %in% c("withHealthy", "noHealthy")) {
  stop("mode must be either 'withHealthy' or 'noHealthy'")
}

message("Input: ", input_rds)
message("Output: ", output_rds)
message("Dataset type: ", dataset_type)
message("Mode: ", mode)

seurat_obj <- readRDS(input_rds)
metadata <- seurat_obj@meta.data

if (dataset_type == "immune_aging" && "disease" %in% colnames(metadata)) {
  metadata <- metadata %>%
    mutate(
      cohort = if_else(
        str_detect(tolower(disease), "healthy|normal|control"),
        "healthy",
        "diseased"
      )
    )

} else if (dataset_type == "heart" && "HF.etiology" %in% colnames(metadata)) {
  metadata <- metadata %>%
    mutate(
      cohort = if_else(
        tolower(HF.etiology) == "donor",
        "healthy",
        "diseased"
      )
    )

} else if (dataset_type == "sepsis" && "Group" %in% colnames(metadata)) {
  metadata <- metadata %>%
    mutate(
      cohort = if_else(
        Group == "Adult Healthy Control",
        "healthy",
        "diseased"
      )
    )

} else if (dataset_type == "vaccine" && "time" %in% colnames(metadata)) {
  metadata <- metadata %>%
    mutate(
      cohort = if_else(
        time == 0,
        "healthy",
        "diseased"
      )
    )

} else if (dataset_type == "our_data" && "condition" %in% colnames(metadata)) {
  metadata <- metadata %>%
    mutate(
      cohort = if_else(
        condition == "Healthy",
        "healthy",
        "diseased"
      )
    )

} else {
  stop(
    "Could not assign cohort for dataset_type = '", dataset_type,
    "'. Available metadata columns: ",
    paste(colnames(metadata), collapse = ", ")
  )
}

seurat_obj@meta.data$cohort <- metadata$cohort

if (mode == "noHealthy") {
  cells_keep <- rownames(seurat_obj@meta.data)[seurat_obj$cohort != "healthy"]

  message("Keeping cells after healthy exclusion: ", length(cells_keep))

  if (length(cells_keep) == 0) {
    message("No diseased cells found after healthy exclusion.")
    message("Keeping original object so withHealthy and noHealthy outputs are identical.")
    seurat_obj$cohort_filter_note <- "noHealthy_requested_but_no_diseased_cells_found_original_object_kept"
  } else {
    seurat_obj <- subset(
      seurat_obj,
      cells = cells_keep
    )

    seurat_obj$cohort_filter_note <- "healthy_cells_excluded"
  }

  rm(cells_keep)
  gc()
}

if (mode == "withHealthy") {
  seurat_obj$cohort_filter_note <- "original_object_with_healthy_cells_kept"
}

dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)

saveRDS(
  seurat_obj,
  output_rds,
  compress = FALSE
)

message("Saved: ", output_rds)