#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)

input_file <- args[[1]]
output_dir <- args[[2]]

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(input_file)

dataset_name <- tools::file_path_sans_ext(basename(input_file))

if (grepl("ImmuneAging", dataset_name)) {
  pair_col <- "donor_id"
} else if (grepl("our_dataset", dataset_name)) {
  pair_col <- "patient"
} else {
  pair_col <- "sample"
}

if (!pair_col %in% colnames(obj@meta.data)) {
  stop(
    "Pairing column not found: ", pair_col,
    "\nAvailable columns: ",
    paste(colnames(obj@meta.data), collapse = ", ")
  )
}

message("Using pairing column: ", pair_col)

candidate_min_cells <- c(3, 5, 8, 10, 15, 20)

counts <- obj@meta.data %>%
  filter(
    !is.na(.data[[pair_col]]),
    !is.na(pla_status),
    !is.na(lineage)
  ) %>%
  count(
    pair_id = .data[[pair_col]],
    pla_status,
    lineage,
    name = "n_cells"
  ) %>%
  complete(
    pair_id,
    pla_status = c("PLA", "platelet-free"),
    lineage,
    fill = list(n_cells = 0)
  )

# ============================================================
# dataset overview
# ============================================================

pairing_overview <- counts %>%
  group_by(pair_id) %>%
  summarise(
    has_PLA = any(pla_status == "PLA" & n_cells > 0),
    has_platelet_free = any(pla_status == "platelet-free" & n_cells > 0),
    .groups = "drop"
  )

dataset_overview <- tibble(
  dataset = dataset_name,
  pairing_column = pair_col,
  total_cells = ncol(obj),
  total_pair_ids = n_distinct(counts$pair_id),
  complete_pairs = sum(
    pairing_overview$has_PLA &
      pairing_overview$has_platelet_free
  ),
  PLA_only_ids = sum(
    pairing_overview$has_PLA &
      !pairing_overview$has_platelet_free
  ),
  platelet_free_only_ids = sum(
    !pairing_overview$has_PLA &
      pairing_overview$has_platelet_free
  ),
  total_lineages = n_distinct(counts$lineage)
)

write.csv(
  dataset_overview,
  file.path(output_dir, "dataset_overview.csv"),
  row.names = FALSE
)

# ============================================================
# pairing overview by lineage
# ============================================================

pairing_by_lineage <- counts %>%
  mutate(present = n_cells > 0) %>%
  select(pair_id, lineage, pla_status, present) %>%
  pivot_wider(
    names_from = pla_status,
    values_from = present,
    values_fill = FALSE
  ) %>%
  group_by(lineage) %>%
  summarise(
    PLA_ids = sum(PLA),
    platelet_free_ids = sum(`platelet-free`),
    complete_pairs = sum(PLA & `platelet-free`),
    .groups = "drop"
  )

write.csv(
  pairing_by_lineage,
  file.path(output_dir, "pairing_by_lineage.csv"),
  row.names = FALSE
)

preflight <- counts %>%
  crossing(min_cells = candidate_min_cells) %>%
  mutate(usable = n_cells > min_cells) %>%
  select(pair_id, lineage, pla_status, min_cells, usable) %>%
  pivot_wider(
    names_from = pla_status,
    values_from = usable,
    values_fill = FALSE
  ) %>%
  mutate(
    complete_pair = PLA & `platelet-free`
  ) %>%
  group_by(min_cells, lineage) %>%
  summarise(
    n_complete_pairs = sum(complete_pair),
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
    n_eligible_lineages = sum(eligible),
    eligible_lineages = paste(lineage[eligible], collapse = "; "),
    .groups = "drop"
  )

write.csv(
  summary,
  file.path(output_dir, "min_cells_summary.csv"),
  row.names = FALSE
)

print(summary)