#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(liana)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(forcats)
  library(circlize)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript split_liana.R <input_rds> <output_dir>")
}

input_file <- args[[1]]
base_output_dir <- args[[2]]

dataset_name <- tools::file_path_sans_ext(basename(input_file))

dataset_mode <- stringr::str_extract(dataset_name, "(withHealthy|noHealthy)$")

dataset_clean <- dataset_name %>%
  stringr::str_remove("_(withHealthy|noHealthy)$")

if (is.na(dataset_mode)) {
  dataset_mode <- "unknownMode"
}

plot_title <- function(title) {
  paste0(title, "\n", dataset_clean, " | ", dataset_mode)
}

out_dir <- file.path(base_output_dir, dataset_name)
plot_dir <- file.path(out_dir, "plots")
table_dir <- file.path(out_dir, "tables")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

message("Analyzing: ", dataset_name)
message("Input: ", input_file)
message("Output directory: ", out_dir)
message("Plot directory: ", plot_dir)

save_plot <- function(plot, filename, width = 8, height = 6, dpi = 300) {
  ggsave(
    filename = file.path(plot_dir, filename),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
}

seurat_obj <- readRDS(input_file)
DefaultAssay(seurat_obj)

sample_col    <- "sample"       # technical / sample-level replicate
donor_col     <- "donor_id"     # biological donor
celltype_col  <- "celltype"     # broad cell type
celltype_full <- "celltype_full" # more detailed cell type
condition_col <- "pla_status"   # PLA vs platelet-free
lineage_col   <- "lineage"

table(seurat_obj$pla_status, useNA = "ifany")
table(seurat_obj$pla_status, seurat_obj$lineage)

seurat_obj <- subset(
  seurat_obj,
  subset = !is.na(lineage)
)

# split into pla & platelet-free subsets
seurat_obj$pla_status[is.na(seurat_obj$pla_status)] <- "platelet-free"

seurat_pla <- subset(seurat_obj, subset = pla_status == "PLA")
seurat_pf  <- subset(seurat_obj, subset = pla_status == "platelet-free")

# LIANA should infer communication between lineages, so identities are lineage
Idents(seurat_pla) <- "lineage"
Idents(seurat_pf)  <- "lineage"

# perform liana analysis
liana_pla_raw <- liana_wrap(seurat_pla)
liana_pla <- liana_pla_raw %>%
  liana_aggregate() %>%
  arrange(aggregate_rank)

liana_pf_raw <- liana_wrap(seurat_pf)
liana_pf <- liana_pf_raw %>%
  liana_aggregate() %>%
  arrange(aggregate_rank)

# -------------------------
# prepare LIANA split tables
# -------------------------

liana_pla_plot <- liana_pla %>%
  as_tibble() %>%
  transmute(
    source,
    target,
    ligand = ligand.complex,
    receptor = receptor.complex,
    aggregate_rank_PLA = aggregate_rank,
    mean_rank_PLA = mean_rank,
    LRscore_PLA = sca.LRscore,
    specificity_PLA = natmi.edge_specificity
  )

liana_pf_plot <- liana_pf %>%
  as_tibble() %>%
  transmute(
    source,
    target,
    ligand = ligand.complex,
    receptor = receptor.complex,
    aggregate_rank_pf = aggregate_rank,
    mean_rank_pf = mean_rank,
    LRscore_pf = sca.LRscore,
    specificity_pf = natmi.edge_specificity
  )

liana_split <- full_join(
  liana_pla_plot,
  liana_pf_plot,
  by = c("source", "target", "ligand", "receptor")
) %>%
  mutate(
    interaction = paste(ligand, receptor, sep = " → "),
    lineage_pair = paste(source, target, sep = " → ")
  )

write.csv(liana_pla, file.path(table_dir, "liana_pla_aggregated.csv"), row.names = FALSE)
write.csv(liana_pf, file.path(table_dir, "liana_platelet_free_aggregated.csv"), row.names = FALSE)
write.csv(liana_split, file.path(table_dir, "liana_split_joined.csv"), row.names = FALSE)


# -------------------------
# threshold sensitivity
# -------------------------

thresholds <- c(0.001, 0.005, 0.01, 0.05, 0.1)

threshold_summary <- lapply(thresholds, function(rank_threshold) {
  
  liana_split %>%
    mutate(
      detected_PLA = !is.na(aggregate_rank_PLA) & aggregate_rank_PLA < rank_threshold,
      detected_pf  = !is.na(aggregate_rank_pf) & aggregate_rank_pf < rank_threshold,
      
      direction = case_when(
        detected_PLA & !detected_pf ~ "PLA-only",
        !detected_PLA & detected_pf ~ "platelet-free-only",
        detected_PLA & detected_pf  ~ "shared",
        TRUE ~ "not detected"
      )
    ) %>%
    filter(direction != "not detected") %>%
    count(direction) %>%
    mutate(rank_threshold = rank_threshold)
  
}) %>%
  bind_rows()

write.csv(threshold_summary, file.path(table_dir, "threshold_summary.csv"), row.names = FALSE)

p_threshold_summary <- ggplot(
  threshold_summary,
  aes(
    x = factor(rank_threshold),
    y = n,
    fill = direction
  )
) +
  geom_col(position = "dodge") +
  theme_bw() +
  scale_fill_manual(
    values = c(
      "PLA-only" = "#F8766D",
      "platelet-free-only" = "#00BFC4",
      "shared" = "#7CAE00"
    )
  ) +
  labs(
    title = plot_title("LIANA split sensitivity to aggregate-rank threshold"),
    x = "Aggregate-rank threshold",
    y = "# LR interactions",
    fill = "Category"
  )

save_plot(
  p_threshold_summary,
  "01_threshold_sensitivity.png",
  width = 7,
  height = 5
)


# -------------------------
# fixed threshold
# -------------------------

rank_threshold <- 0.05

liana_split_thresh <- liana_split %>%
  mutate(
    detected_PLA = !is.na(aggregate_rank_PLA) & aggregate_rank_PLA < rank_threshold,
    detected_pf  = !is.na(aggregate_rank_pf) & aggregate_rank_pf < rank_threshold,
    
    direction = case_when(
      detected_PLA & !detected_pf ~ "PLA-only",
      !detected_PLA & detected_pf ~ "platelet-free-only",
      detected_PLA & detected_pf  ~ "shared",
      TRUE ~ "not detected"
    )
  ) %>%
  filter(direction != "not detected")

write.csv(liana_split_thresh, file.path(table_dir, "liana_split_thresholded.csv"), row.names = FALSE)


# -------------------------
# 1. global category counts
# -------------------------

p_liana_split_global <- ggplot(liana_split_thresh, aes(x = direction, fill = direction)) +
  geom_bar() +
  theme_bw() +
  scale_fill_manual(
    values = c(
      "PLA-only" = "#F8766D",
      "platelet-free-only" = "#00BFC4",
      "shared" = "#7CAE00"
    )
  ) +
  labs(
    title = plot_title(paste0("LIANA split interactions at aggregate_rank < ", rank_threshold)),
    x = "Category",
    y = "# LR interactions"
  ) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

save_plot(
  p_liana_split_global,
  "02_global_liana_split_categories.png",
  width = 7,
  height = 5
)


# -------------------------
# 2. condition-specific lineage pairs
# -------------------------

pair_counts_liana_split_specific <- liana_split_thresh %>%
  filter(direction != "shared") %>%
  count(lineage_pair, direction) %>%
  group_by(lineage_pair) %>%
  mutate(total = sum(n)) %>%
  ungroup() %>%
  arrange(desc(total)) %>%
  slice_max(total, n = 25, with_ties = FALSE) %>%
  mutate(lineage_pair = fct_reorder(lineage_pair, total))

write.csv(
  pair_counts_liana_split_specific,
  file.path(table_dir, "pair_counts_liana_split_specific.csv"),
  row.names = FALSE
)

p_pair_counts_liana_split_specific <- ggplot(
  pair_counts_liana_split_specific,
  aes(x = n, y = lineage_pair, fill = direction)
) +
  geom_col() +
  theme_bw() +
  scale_fill_manual(
    values = c(
      "PLA-only" = "#F8766D",
      "platelet-free-only" = "#00BFC4"
    )
  ) +
  labs(
    title = plot_title("Top condition-specific LIANA split lineage pairs"),
    x = "# condition-specific LR interactions",
    y = "Source → target lineage",
    fill = "Category"
  ) +
  theme(
    axis.text.y = element_text(size = 7),
    axis.title = element_text(size = 9),
    plot.title = element_text(size = 11),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )

save_plot(
  p_pair_counts_liana_split_specific,
  "03_top_condition_specific_lineage_pairs.png",
  width = 9,
  height = 7
)


# define this for the later dotplot
top_lineage_pairs <- pair_counts_liana_split_specific %>%
  pull(lineage_pair) %>%
  unique()


# -------------------------
# 3. net direction heatmap
# -------------------------

net_direction_liana_split <- liana_split_thresh %>%
  filter(direction != "shared") %>%
  count(source, target, direction) %>%
  complete(
    source,
    target,
    direction = c("PLA-only", "platelet-free-only"),
    fill = list(n = 0)
  ) %>%
  pivot_wider(
    names_from = direction,
    values_from = n,
    values_fill = 0
  ) %>%
  mutate(
    net = `PLA-only` - `platelet-free-only`
  )

write.csv(
  net_direction_liana_split,
  file.path(table_dir, "net_direction_liana_split.csv"),
  row.names = FALSE
)

p_net_direction_liana_split <- ggplot(
  net_direction_liana_split,
  aes(x = target, y = source, fill = net)
) +
  geom_tile(color = "white") +
  geom_text(aes(label = net), size = 4) +
  theme_bw() +
  labs(
    title = plot_title("Net LIANA split condition-specific interactions"),
    x = "Target lineage",
    y = "Source lineage",
    fill = "PLA-only minus\nplatelet-free-only"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

save_plot(
  p_net_direction_liana_split,
  "04_net_direction_liana_split_heatmap.png",
  width = 8,
  height = 7
)


# -------------------------
# 4. recurrent LR pairs
# -------------------------

top_recurrent_liana_split <- liana_split_thresh %>%
  filter(direction != "shared") %>%
  count(interaction, direction, name = "n_lineage_pairs") %>%
  group_by(interaction) %>%
  mutate(total = sum(n_lineage_pairs)) %>%
  ungroup() %>%
  arrange(desc(total)) %>%
  slice_max(total, n = 25, with_ties = FALSE) %>%
  mutate(interaction = fct_reorder(interaction, total))

write.csv(
  top_recurrent_liana_split,
  file.path(table_dir, "top_recurrent_liana_split.csv"),
  row.names = FALSE
)

p_top_recurrent_liana_split <- ggplot(
  top_recurrent_liana_split,
  aes(x = n_lineage_pairs, y = interaction, fill = direction)
) +
  geom_col() +
  theme_bw() +
  scale_fill_manual(
    values = c(
      "PLA-only" = "#F8766D",
      "platelet-free-only" = "#00BFC4"
    )
  ) +
  labs(
    title = plot_title("Most recurrent LIANA split LR pairs"),
    x = "# source-target lineage pairs",
    y = "Ligand → receptor",
    fill = "Category"
  )

save_plot(
  p_top_recurrent_liana_split,
  "05_top_recurrent_liana_split_lr_pairs.png",
  width = 8,
  height = 8
)


# -------------------------
# 5. top interactions per lineage pair
# -------------------------

top_liana_split_per_pair_specific <- liana_split_thresh %>%
  filter(direction != "shared") %>%
  filter(lineage_pair %in% top_lineage_pairs) %>%
  mutate(
    best_rank = pmin(
      aggregate_rank_PLA,
      aggregate_rank_pf,
      na.rm = TRUE
    ),
    rank_score = -log10(best_rank + 1e-300)
  ) %>%
  group_by(lineage_pair) %>%
  slice_max(rank_score, n = 3, with_ties = FALSE) %>%
  ungroup()

write.csv(
  top_liana_split_per_pair_specific,
  file.path(table_dir, "top_liana_split_per_pair_specific.csv"),
  row.names = FALSE
)

p_top_liana_split_per_pair_specific <- ggplot(
  top_liana_split_per_pair_specific,
  aes(
    x = lineage_pair,
    y = interaction,
    size = rank_score,
    color = direction
  )
) +
  geom_point(alpha = 0.8) +
  theme_bw() +
  scale_color_manual(
    values = c(
      "PLA-only" = "#F8766D",
      "platelet-free-only" = "#00BFC4"
    )
  ) +
  labs(
    title = plot_title("Top condition-specific LIANA split interactions"),
    x = "Source → target lineage",
    y = "Ligand → receptor",
    size = "-log10(best aggregate rank)",
    color = "Category"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y = element_text(size = 7),
    axis.title = element_text(size = 9),
    plot.title = element_text(size = 11),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7)
  )

save_plot(
  p_top_liana_split_per_pair_specific,
  "06_top_condition_specific_liana_split_interactions.png",
  width = 13,
  height = 10
)

message("Finished LIANA split analysis for: ", dataset_name)
message("Saved plots to: ", plot_dir)
message("Saved tables to: ", table_dir)