#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(liana)
  library(Seurat)
  library(dplyr)
  library(tibble)
})

# ------------------------------------------------------------
# Command-line arguments
# ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop(
    paste(
      "Usage:",
      "Rscript run_liana.R <input_rds> <dataset_name> <output_dir>"
    )
  )
}

input_file  <- args[[1]]
dataset_name <- args[[2]]
output_dir  <- args[[3]]

# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

if (!file.exists(input_file)) {
  stop("Input file does not exist: ", input_file)
}

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

message("Input file: ", input_file)
message("Dataset: ", dataset_name)
message("Output directory: ", output_dir)
message("R version: ", R.version.string)
message("LIANA version: ", packageVersion("liana"))
message("Seurat version: ", packageVersion("Seurat"))

# ------------------------------------------------------------
# Load Seurat object
# ------------------------------------------------------------

message("Loading Seurat object...")

seurat_obj <- readRDS(input_file)

message("Number of cells: ", ncol(seurat_obj))
message("Number of features: ", nrow(seurat_obj))
message("Default assay: ", DefaultAssay(seurat_obj))
message(
  "Active identities: ",
  paste(levels(Idents(seurat_obj)), collapse = ", ")
)

# ------------------------------------------------------------
# Run LIANA
# ------------------------------------------------------------

set.seed(123)

message("Running LIANA...")

liana_raw <- liana_wrap(seurat_obj)

message(
  "Completed methods: ",
  paste(names(liana_raw), collapse = ", ")
)

# ------------------------------------------------------------
# Consensus aggregation
# ------------------------------------------------------------

message("Aggregating LIANA results...")

liana_consensus <- liana_raw %>%
  liana_aggregate() %>%
  arrange(aggregate_rank)

# ------------------------------------------------------------
# Count cells per identity
# ------------------------------------------------------------

cell_type_counts_table <- table(Idents(seurat_obj))

cell_type_counts <- tibble(
  cell_type = names(cell_type_counts_table),
  n_cells = as.integer(cell_type_counts_table)
) %>%
  arrange(desc(n_cells))

# ------------------------------------------------------------
# Bundle results
# ------------------------------------------------------------

liana_result <- list(
  metadata = list(
    dataset_name = dataset_name,
    input_file = normalizePath(input_file),
    created_at = Sys.time(),
    n_cells = ncol(seurat_obj),
    n_features = nrow(seurat_obj),
    default_assay = DefaultAssay(seurat_obj),
    identity_levels = levels(Idents(seurat_obj)),
    methods = names(liana_raw),
    R_version = R.version.string,
    liana_version = as.character(packageVersion("liana")),
    Seurat_version = as.character(packageVersion("Seurat"))
  ),

  cell_type_counts = cell_type_counts,

  raw = liana_raw,

  consensus = liana_consensus
)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

output_file <- file.path(
  output_dir,
  paste0(dataset_name, "_liana_results.rds")
)

saveRDS(
  liana_result,
  file = output_file,
  compress = "xz"
)

message("Successfully saved LIANA results:")
message(output_file)
message("Finished.")