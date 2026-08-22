#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)

input_file <- args[[1]]
output_dir <- args[[2]]

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(input_file)
dataset_name <- tools::file_path_sans_ext(basename(input_file))

sample_col <- if ("sample_ID" %in% colnames(obj@meta.data)) {
  "sample_ID"
} else {
  "sample"
}

pair_col <- if (grepl("ImmuneAging", dataset_name)) {
  "donor_id"
} else if (grepl("our_dataset|impact", dataset_name, ignore.case = TRUE)) {
  "patient"
} else {
  "sample"
}

min_cells_values <- c(3, 5, 8, 10, 15, 20)

meta <- obj@meta.data %>%
  mutate(
    sample_id = as.character(.data[[sample_col]]),
    pair_id = as.character(.data[[pair_col]]),
    group = as.character(pla_status),
    celltype = as.character(lineage)
  )

# 1. each sample must belong to exactly one group and pair
sample_mapping <- meta %>%
  distinct(sample_id, pair_id, group) %>%
  group_by(sample_id) %>%
  summarise(
    n_groups = n_distinct(group),
    n_pairs = n_distinct(pair_id),
    valid = n_groups == 1 & n_pairs == 1,
    .groups = "drop"
  )

write.csv(
  sample_mapping,
  file.path(output_dir, "sample_mapping.csv"),
  row.names = FALSE
)

# 2. check raw counts and pseudobulk library sizes
counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
cell_libsize <- Matrix::colSums(counts)

pb <- meta %>%
  mutate(cell_libsize = cell_libsize[rownames(meta)]) %>%
  group_by(sample_id, pair_id, group, celltype) %>%
  summarise(
    n_cells = n(),
    library_size = sum(cell_libsize),
    .groups = "drop"
  )

write.csv(
  pb,
  file.path(output_dir, "sample_celltype_counts.csv"),
  row.names = FALSE
)

# 3. actual MultiNicheNet-style eligibility
preflight <- pb %>%
  crossing(min_cells = min_cells_values) %>%
  mutate(
    usable = n_cells > min_cells &
      library_size > 0
  ) %>%
  filter(usable) %>%
  distinct(min_cells, celltype, pair_id, group) %>%
  mutate(present = TRUE) %>%
  complete(
    min_cells,
    celltype,
    pair_id,
    group = c("PLA", "platelet-free"),
    fill = list(present = FALSE)
  ) %>%
  pivot_wider(
    names_from = group,
    values_from = present,
    values_fill = FALSE
  ) %>%
  group_by(min_cells, celltype) %>%
  summarise(
    n_complete_pairs = sum(PLA & `platelet-free`),
    eligible = n_complete_pairs >= 2,
    .groups = "drop"
  )

write.csv(
  preflight,
  file.path(output_dir, "min_cells_preflight.csv"),
  row.names = FALSE
)

summary <- preflight %>%
  group_by(min_cells) %>%
  summarise(
    eligible_lineages = paste(
      celltype[eligible],
      collapse = "; "
    ),
    .groups = "drop"
  )

write.csv(
  summary,
  file.path(output_dir, "min_cells_summary.csv"),
  row.names = FALSE
)

print(summary)

cat(
  "\nInvalid sample mappings:",
  sum(!sample_mapping$valid),
  "\nZero pseudobulk libraries:",
  sum(pb$library_size == 0),
  "\n"
)