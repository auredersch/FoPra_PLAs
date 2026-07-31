#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop(
    paste(
      "Usage:",
      "Rscript prepare_cohort_rds.R",
      "<input_rds> <output_rds> <dataset_type>",
      "<healthyOnly|diseasedOnly|all>"
    )
  )
}

input_rds <- args[[1]]
output_rds <- args[[2]]
dataset_type <- args[[3]]
mode <- args[[4]]

valid_modes <- c(
  "healthyOnly",
  "diseasedOnly",
  "all"
)

if (!mode %in% valid_modes) {
  stop(
    "mode must be one of: ",
    paste(valid_modes, collapse = ", ")
  )
}

if (!file.exists(input_rds)) {
  stop("Input RDS does not exist: ", input_rds)
}

message("Input: ", input_rds)
message("Output: ", output_rds)
message("Dataset type: ", dataset_type)
message("Mode: ", mode)

seurat_obj <- readRDS(input_rds)

if (!inherits(seurat_obj, "Seurat")) {
  stop("Input file does not contain a Seurat object: ", input_rds)
}

metadata <- seurat_obj@meta.data

# ------------------------------------------------------------
# Assign healthy/diseased cohort according to dataset
# ------------------------------------------------------------

if (
  dataset_type == "immune_aging" &&
  "disease" %in% colnames(metadata)
) {

  disease_value <- tolower(trimws(as.character(metadata$disease)))

  metadata <- metadata %>%
    mutate(
      cohort = case_when(
        str_detect(disease_value, "healthy|normal|control") ~ "healthy",
        !is.na(disease_value) & disease_value != "" ~ "diseased",
        TRUE ~ NA_character_
      )
    )

} else if (
  dataset_type == "heart" &&
  "HF.etiology" %in% colnames(metadata)
) {

  etiology_value <- tolower(trimws(as.character(metadata$HF.etiology)))

  metadata <- metadata %>%
    mutate(
      cohort = case_when(
        etiology_value == "donor" ~ "healthy",
        !is.na(etiology_value) & etiology_value != "" ~ "diseased",
        TRUE ~ NA_character_
      )
    )

} else if (
  dataset_type == "sepsis" &&
  "Group" %in% colnames(metadata)
) {

  group_value <- trimws(as.character(metadata$Group))

  metadata <- metadata %>%
    mutate(
      cohort = case_when(
        group_value == "Adult Healthy Control" ~ "healthy",
        !is.na(group_value) & group_value != "" ~ "diseased",
        TRUE ~ NA_character_
      )
    )

} else if (
  dataset_type == "vaccine" &&
  "time" %in% colnames(metadata)
) {

  time_value <- trimws(as.character(metadata$time))

  metadata <- metadata %>%
    mutate(
      cohort = case_when(
        time_value == "0" ~ "healthy",
        !is.na(time_value) & time_value != "" ~ "diseased",
        TRUE ~ NA_character_
      )
    )

} else if (
  dataset_type == "our_data" &&
  "condition" %in% colnames(metadata)
) {

  condition_value <- tolower(trimws(as.character(metadata$condition)))

  metadata <- metadata %>%
    mutate(
      cohort = case_when(
        condition_value == "healthy" ~ "healthy",
        !is.na(condition_value) & condition_value != "" ~ "diseased",
        TRUE ~ NA_character_
      )
    )

} else if (
  dataset_type == "skin" &&
  "Status" %in% colnames(metadata)
) {

  status_value <- trimws(as.character(metadata$Status))

  metadata <- metadata %>%
    mutate(
      cohort = case_when(
        status_value == "Healthy" ~ "healthy",
        status_value %in% c(
          "AXI",
          "PSA",
          "PSO",
          "PSX"
        ) ~ "diseased",
        TRUE ~ NA_character_
      )
    )

  message("Skin status-to-cohort mapping:")
  print(
    table(
      Status = metadata$Status,
      cohort = metadata$cohort,
      useNA = "ifany"
    )
  )

} else {
  stop(
    "Could not assign cohort for dataset_type = '",
    dataset_type,
    "'. Available metadata columns: ",
    paste(colnames(metadata), collapse = ", ")
  )
}

# ------------------------------------------------------------
# Validate cohort assignment
# ------------------------------------------------------------

if (!"cohort" %in% colnames(metadata)) {
  stop("Internal error: cohort column was not created.")
}

if (anyNA(metadata$cohort)) {
  unmapped_rows <- is.na(metadata$cohort)

  stop(
    "Cohort assignment produced ",
    sum(unmapped_rows),
    " missing values for dataset_type = '",
    dataset_type,
    "'."
  )
}

unexpected_cohorts <- setdiff(
  unique(as.character(metadata$cohort)),
  c("healthy", "diseased")
)

if (length(unexpected_cohorts) > 0) {
  stop(
    "Unexpected cohort values: ",
    paste(unexpected_cohorts, collapse = ", ")
  )
}

seurat_obj$cohort <- metadata$cohort

message("Cohort counts before filtering:")
print(
  table(
    seurat_obj$cohort,
    useNA = "ifany"
  )
)

# ------------------------------------------------------------
# Filter by requested mode
# ------------------------------------------------------------

if (mode == "healthyOnly") {

  cells_keep <- rownames(seurat_obj@meta.data)[
    seurat_obj$cohort == "healthy"
  ]

  if (length(cells_keep) == 0) {
    stop(
      "No healthy cells found. ",
      "Cannot create healthyOnly output."
    )
  }

  seurat_obj <- subset(
    seurat_obj,
    cells = cells_keep
  )

  seurat_obj$cohort_filter_note <- "healthy_cells_only"

} else if (mode == "diseasedOnly") {

  cells_keep <- rownames(seurat_obj@meta.data)[
    seurat_obj$cohort == "diseased"
  ]

  if (length(cells_keep) == 0) {
    stop(
      "No diseased cells found. ",
      "Cannot create diseasedOnly output."
    )
  }

  seurat_obj <- subset(
    seurat_obj,
    cells = cells_keep
  )

  seurat_obj$cohort_filter_note <- "diseased_cells_only"

} else if (mode == "all") {

  seurat_obj$cohort_filter_note <- "all_cells_kept"
}

rm(metadata)

if (exists("cells_keep")) {
  rm(cells_keep)
}

message("Cohort counts after filtering:")
print(
  table(
    seurat_obj$cohort,
    useNA = "ifany"
  )
)

message("Cells in final object: ", ncol(seurat_obj))

# ------------------------------------------------------------
# Save output
# ------------------------------------------------------------

dir.create(
  dirname(output_rds),
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  seurat_obj,
  output_rds,
  compress = FALSE
)

message("Saved: ", output_rds)