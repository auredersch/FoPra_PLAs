#!/usr/bin/env Rscript

options(
  error = function() {
    message("R error occurred. Traceback:")
    traceback(2)
    quit(save = "no", status = 1)
  }
)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(nichenetr)
  library(multinichenetr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
  library(circlize)
  library(RColorBrewer)
})

# ============================================================
# 1. arguments and output
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript multinichetr.R <input_rds> <output_dir>")
}

input_file <- args[[1]]
base_output_dir <- args[[2]]

dataset_name <- tools::file_path_sans_ext(
  basename(input_file)
)

dataset_mode <- dplyr::case_when(
  stringr::str_ends(
    dataset_name,
    "_healthyOnly"
  ) ~ "healthyOnly",

  stringr::str_ends(
    dataset_name,
    "_diseasedOnly"
  ) ~ "diseasedOnly",

  TRUE ~ "all"
)

dataset_clean <- dataset_name %>%
  stringr::str_remove(
    "_(healthyOnly|diseasedOnly)$"
  )

plot_title <- function(title) {
  paste0(title, "\n", dataset_clean, " | ", dataset_mode)
}

plot_dir <- file.path(
  base_output_dir,
  dataset_name,
  "plots"
)

dir.create(
  plot_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

save_plot <- function(
    plot,
    filename,
    width = 8,
    height = 6,
    dpi = 250
) {
  output_file <- file.path(
    plot_dir,
    filename
  )

  message("Saving plot: ", output_file)

  ggsave(
    filename = output_file,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )

  if (!file.exists(output_file)) {
    stop("Plot was not created: ", output_file)
  }

  message(
    "Saved plot: ",
    filename,
    " [",
    file.info(output_file)$size
  )

  invisible(output_file)
}

message("Dataset: ", dataset_name)
message("Input: ", input_file)
message("Plots: ", plot_dir)

# ============================================================
# 2. load data and define metadata columns
# ============================================================

seurat_obj <- readRDS(input_file)

condition_col <- "pla_status"
lineage_col <- "lineage"

if (grepl("ImmuneAging", dataset_name)) {
  sample_col <- "donor_id"
  pair_col <- "donor_id"

} else if (grepl("our_dataset", dataset_name)) {
  sample_col <- "sample_ID"
  pair_col <- "patient"

} else {
  sample_col <- "sample"
  pair_col <- "sample"
}

required_cols <- unique(
  c(
    sample_col,
    pair_col,
    condition_col,
    lineage_col
  )
)

missing_cols <- setdiff(
  required_cols,
  colnames(seurat_obj@meta.data)
)

if (length(missing_cols) > 0) {
  stop(
    "Missing metadata columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

keep <- complete.cases(
  seurat_obj@meta.data[
    ,
    required_cols,
    drop = FALSE
  ]
)

seurat_obj <- subset(
  seurat_obj,
  cells = colnames(seurat_obj)[keep]
)

DefaultAssay(seurat_obj) <- "RNA"

required_conditions <- c(
  "PLA",
  "platelet-free"
)

available_conditions <- unique(
  as.character(
    seurat_obj@meta.data[[condition_col]]
  )
)

missing_conditions <- setdiff(
  required_conditions,
  available_conditions
)

if (length(missing_conditions) > 0) {
  stop(
    "Missing conditions: ",
    paste(missing_conditions, collapse = ", ")
  )
}

# One pseudobulk sample per biological sample and condition
seurat_obj$sample_condition <- interaction(
  as.character(
    seurat_obj@meta.data[[sample_col]]
  ),
  as.character(
    seurat_obj@meta.data[[condition_col]]
  ),
  sep = "_",
  drop = TRUE
)

# Pair identifier for the paired model
seurat_obj$pair_id <- as.character(
  seurat_obj@meta.data[[pair_col]]
)

pairing_check <- seurat_obj@meta.data %>%
  dplyr::as_tibble() %>%
  dplyr::transmute(
    pair_id = as.character(.data[["pair_id"]]),
    condition = as.character(.data[[condition_col]])
  ) %>%
  dplyr::distinct(
    pair_id,
    condition
  ) %>%
  dplyr::count(
    pair_id,
    name = "n_conditions"
  )


n_complete_pairs <- sum(
  pairing_check$n_conditions == 2
)

message("Sample column: ", sample_col)
message("Pairing column: ", pair_col)
message("Complete pairs: ", n_complete_pairs)

if (n_complete_pairs < 2) {
  stop(
    "Fewer than two complete PLA/platelet-free pairs are present."
  )
}

# ============================================================
# 3. analysis settings
# ============================================================

organism <- "human"

sample_id <- "sample_condition"
group_id <- condition_col
celltype_id <- lineage_col
covariates <- "pair_id"
batches <- NA

min_cells_config <- c(
  "gated_heart_processed" = 3L,
  "gated_heart_processed_diseasedOnly" = 3L,
  "gated_heart_processed_healthyOnly" = NA_integer_,

  "gated_ImmuneAging" = 8L,
  "gated_ImmuneAging_diseasedOnly" = NA_integer_,
  "gated_ImmuneAging_healthyOnly" = 8L,

  "gated_our_dataset_processed" = 1L,
  "gated_our_dataset_processed_diseasedOnly" = 1L,
  "gated_our_dataset_processed_healthyOnly" = 1L,

  "gated_sepsis_processed" = 20L,
  "gated_sepsis_processed_diseasedOnly" = 20L,
  "gated_sepsis_processed_healthyOnly" = 15L,

  "gated_vaccine_processed" = 5L,
  "gated_vaccine_processed_diseasedOnly" = 5L,
  "gated_vaccine_processed_healthyOnly" = 5L,

  "gated_skin_processed" = 3L,
  "gated_skin_processed_diseasedOnly" = 3L,
  "gated_skin_processed_healthyOnly" = 5L
)

if (!dataset_name %in% names(min_cells_config)) {
  stop(
    "No min_cells configured for dataset: ",
    dataset_name
  )
}

min_cells <- as.numeric(
  unname(
    min_cells_config[[dataset_name]]
  )
)

if (is.na(min_cells)) {
  message(
    "Skipping ",
    dataset_name,
    ": insufficient paired samples."
  )
  quit(save = "no", status = 0)
}

min_sample_prop <- 0.50
fraction_cutoff <- 0.05
logFC_threshold <- 0.50
p_val_threshold <- 0.05
p_val_adj <- TRUE
top_n_target <- 250

options(timeout = 120)

message("min_cells: ", min_cells)

# ============================================================
# 4. load ligand-receptor resources
# ============================================================

lr_network <- readRDS(
  url(
    "https://zenodo.org/record/10229222/files/lr_network_human_allInfo_30112033.rds"
  )
) %>%
  transmute(
    ligand = make.names(
      convert_alias_to_symbols(
        ligand,
        organism = organism
      )
    ),
    receptor = make.names(
      convert_alias_to_symbols(
        receptor,
        organism = organism
      )
    )
  ) %>%
  distinct()

ligand_target_matrix <- readRDS(
  url(
    "https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"
  )
)

colnames(ligand_target_matrix) <- make.names(
  convert_alias_to_symbols(
    colnames(ligand_target_matrix),
    organism = organism
  )
)

rownames(ligand_target_matrix) <- make.names(
  convert_alias_to_symbols(
    rownames(ligand_target_matrix),
    organism = organism
  )
)

lr_network <- lr_network %>%
  filter(
    ligand %in% colnames(
      ligand_target_matrix
    )
  )

ligand_target_matrix <- ligand_target_matrix[
  ,
  unique(lr_network$ligand),
  drop = FALSE
]

# ============================================================
# 5. convert to SingleCellExperiment
# ============================================================

sce <- as.SingleCellExperiment(
  seurat_obj,
  assay = "RNA"
)

sce <- alias_to_symbol_SCE(
  sce,
  organism
) %>%
  makenames_SCE()

for (
  column in c(
    sample_id,
    group_id,
    celltype_id,
    covariates
  )
) {
  colData(sce)[[column]] <- factor(
    make.names(
      as.character(
        colData(sce)[[column]]
      )
    )
  )
}

celltypes_oi <- unique(
  as.character(
    colData(sce)[[celltype_id]]
  )
)

# Colors are named with the same make.names() lineage labels
lineage_colors <- colorRampPalette(
  brewer.pal(11, "Spectral")
)(
  length(sort(celltypes_oi))
)

names(lineage_colors) <- sort(celltypes_oi)

direction_colors <- c(
  "PLA-up" = "#F8766D",
  "platelet-free-up" = "#00BFC4"
)

# ============================================================
# 6. expression and abundance
# ============================================================

abundance_info <- get_abundance_info(
  sce = sce,
  sample_id = sample_id,
  group_id = group_id,
  celltype_id = celltype_id,
  min_cells = min_cells,
  senders_oi = celltypes_oi,
  receivers_oi = celltypes_oi,
  batches = batches
)

frq_list <- get_frac_exprs(
  sce = sce,
  sample_id = sample_id,
  celltype_id = celltype_id,
  group_id = group_id,
  batches = batches,
  min_cells = min_cells,
  fraction_cutoff = fraction_cutoff,
  min_sample_prop = min_sample_prop
)

genes_oi <- frq_list$expressed_df %>%
  filter(expressed) %>%
  pull(gene) %>%
  unique()

if (length(genes_oi) == 0) {
  stop(
    "No expressed genes passed filtering."
  )
}

sce <- sce[
  genes_oi,
  ,
  drop = FALSE
]

abundance_expression_info <-
  process_abundance_expression_info(
    sce = sce,
    sample_id = sample_id,
    group_id = group_id,
    celltype_id = celltype_id,
    min_cells = min_cells,
    senders_oi = celltypes_oi,
    receivers_oi = celltypes_oi,
    lr_network = lr_network,
    batches = batches,
    frq_list = frq_list,
    abundance_info = abundance_info
  )

# ============================================================
# 7. differential expression
# ============================================================

contrasts_oi <- c(
  "'PLA-platelet.free','platelet.free-PLA'"
)

contrast_tbl <- tibble(
  contrast = c(
    "PLA-platelet.free",
    "platelet.free-PLA"
  ),
  group = c(
    "PLA",
    "platelet.free"
  )
)

cat("\n=== DE INPUT CHECK ===\n")

cat("sample_id:", sample_id, "\n")
cat("group_id:", group_id, "\n")
cat("celltype_id:", celltype_id, "\n")
cat("covariates:", covariates, "\n")

cat("\nGroup levels:\n")
print(levels(colData(sce)[[group_id]]))

cat("\nGroup counts:\n")
print(table(colData(sce)[[group_id]]))

cat("\nCell type counts:\n")
print(table(colData(sce)[[celltype_id]]))

cat("\nSample-condition mapping:\n")
print(
  as.data.frame(colData(sce)) %>%
    distinct(
      sample = .data[[sample_id]],
      group = .data[[group_id]],
      pair = .data[[covariates]]
    ) %>%
    head(20)
)

DE_info <- get_DE_info(
  sce = sce,
  sample_id = sample_id,
  group_id = group_id,
  celltype_id = celltype_id,
  batches = batches,
  covariates = covariates,
  contrasts_oi = contrasts_oi,
  min_cells = min_cells,
  expressed_df = frq_list$expressed_df
)

if (
  !is.null(DE_info$hist_pvals) &&
    inherits(
      DE_info$hist_pvals,
      c(
        "gg",
        "ggplot",
        "patchwork"
      )
    )
) {
  save_plot(
    DE_info$hist_pvals,
    "08_deinfo_hist_pvals.png",
    width = 12,
    height = 8
  )
}

celltype_de <-
  DE_info$celltype_de$de_output_tidy

if (
  is.null(celltype_de) ||
    nrow(celltype_de) == 0
) {
  stop(
    "No cell-type DE results were generated."
  )
}

included_celltypes <- unique(
  celltype_de$cluster_id
)

sender_receiver_de <-
  combine_sender_receiver_de(
    sender_de = celltype_de,
    receiver_de = celltype_de,
    senders_oi = included_celltypes,
    receivers_oi = included_celltypes,
    lr_network = lr_network
  )

if (
  is.null(sender_receiver_de) ||
    nrow(sender_receiver_de) == 0
) {
  stop(
    "No ligand-receptor combinations were generated."
  )
}

# ============================================================
# 8. ligand activity
# ============================================================

ligand_activities <- suppressMessages(
  suppressWarnings(
    get_ligand_activities_targets_DEgenes(
      receiver_de = celltype_de,
      receivers_oi = included_celltypes,
      ligand_target_matrix =
        ligand_target_matrix,
      logFC_threshold =
        logFC_threshold,
      p_val_threshold =
        p_val_threshold,
      p_val_adj = p_val_adj,
      top_n_target = top_n_target
    )
  )
)

if (
  is.null(ligand_activities) ||
    NROW(ligand_activities) == 0
) {
  stop(
    "No ligand activities were generated."
  )
}

# ============================================================
# 9. interaction prioritization
# ============================================================

sender_receiver_tbl <-
  sender_receiver_de %>%
  distinct(
    sender,
    receiver
  )

grouping_tbl <-
  as_tibble(
    colData(sce)
  ) %>%
  distinct(
    sample = .data[[sample_id]],
    group = .data[[group_id]]
  )

prioritization_tables <-
  generate_prioritization_tables(
    sender_receiver_info =
      abundance_expression_info$
      sender_receiver_info,
    sender_receiver_de =
      sender_receiver_de,
    ligand_activities_targets_DEgenes =
      ligand_activities,
    contrast_tbl = contrast_tbl,
    sender_receiver_tbl =
      sender_receiver_tbl,
    grouping_tbl = grouping_tbl,
    scenario = "regular",
    fraction_cutoff =
      fraction_cutoff,
    abundance_data_receiver =
      abundance_expression_info$
      abundance_data_receiver,
    abundance_data_sender =
      abundance_expression_info$
      abundance_data_sender,
    ligand_activity_down = FALSE
  )

# ============================================================
# 10. common plotting table
# ============================================================

mn_plot <-
  prioritization_tables$
  group_prioritization_tbl %>%
  as_tibble() %>%
  mutate(
    lineage_pair = paste(
      sender,
      receiver,
      sep = " → "
    ),
    interaction = paste(
      ligand,
      receptor,
      sep = " → "
    ),
    direction = case_when(
      group == "PLA" ~ "PLA-up",
      group == "platelet.free" ~
        "platelet-free-up",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(direction),
    is.finite(
      prioritization_score
    )
  ) %>%
  mutate(
    direction = factor(
      direction,
      levels = c(
        "PLA-up",
        "platelet-free-up"
      )
    )
  )

if (nrow(mn_plot) == 0) {
  stop(
    "No valid prioritized interactions."
  )
}

# ============================================================
# 10b. selected top interactions per condition
# ============================================================

n_selected <- 20

mn_selected <- mn_plot %>%
  dplyr::group_by(direction) %>%
  dplyr::slice_max(
    order_by = prioritization_score,
    n = n_selected,
    with_ties = FALSE
  ) %>%
  ungroup()

message(
  "Selected top ",
  n_selected,
  " MultiNicheNet interactions per direction."
)

print(
  mn_selected %>%
    dplyr::count(direction)
)

# ============================================================
# 11. prioritization-score distribution
# ============================================================

p_score_distribution <- ggplot(
  mn_plot,
  aes(
    x = direction,
    y = prioritization_score,
    fill = direction
  )
) +
  geom_violin(
    trim = TRUE,
    alpha = 0.45
  ) +
  geom_boxplot(
    width = 0.18,
    outlier.size = 0.3
  ) +
  scale_fill_manual(
    values = direction_colors
  ) +
  theme_bw() +
  labs(
    title = plot_title(
      "MultiNicheNet prioritization-score distribution"
    ),
    x = "Prioritized condition",
    y = "Prioritization score"
  ) +
  theme(
    legend.position = "none"
  )

save_plot(
  p_score_distribution,
  "01_multinichenet_score_distribution.png",
  width = 7,
  height = 5
)

# ============================================================
# 12. all prioritized interactions by lineage pair
# ============================================================

pair_counts <- mn_plot %>%
  dplyr::count(
    lineage_pair,
    direction,
    name = "n_interactions"
  ) %>%
  dplyr::group_by(lineage_pair) %>%
  dplyr::mutate(
    total = sum(n_interactions)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::slice_max(
    total,
    n = 25,
    with_ties = FALSE
  ) %>%
  dplyr::mutate(
    lineage_pair = forcats::fct_reorder(
      lineage_pair,
      total
    )
  )

p_pair_counts <- ggplot2::ggplot(
  pair_counts,
  ggplot2::aes(
    x = n_interactions,
    y = lineage_pair,
    fill = direction
  )
) +
  ggplot2::geom_col() +
  ggplot2::scale_fill_manual(
    values = direction_colors
  ) +
  ggplot2::theme_bw() +
  ggplot2::labs(
    title = plot_title(
      "Prioritized MultiNicheNet interactions by lineage pair"
    ),
    x = "# prioritized interactions",
    y = "Sender → receiver lineage",
    fill = "Direction"
  )

save_plot(
  p_pair_counts,
  "02_multinichenet_interactions_by_lineage_pair.png",
  width = 10,
  height = 8
)
# ============================================================
# 13. net mean prioritization score
# ============================================================

net_score <- mn_plot %>%
  group_by(
    sender,
    receiver,
    direction
  ) %>%
  summarise(
    mean_score = mean(
      prioritization_score,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    direction = recode(
      as.character(direction),
      "PLA-up" = "PLA_up",
      "platelet-free-up" =
        "platelet_free_up"
    )
  ) %>%
  pivot_wider(
    names_from = direction,
    values_from = mean_score
  ) %>%
  filter(
    !is.na(PLA_up),
    !is.na(platelet_free_up)
  ) %>%
  mutate(
    net_score =
      PLA_up -
      platelet_free_up
  )

net_limit <- max(
  abs(net_score$net_score),
  na.rm = TRUE
)

if (
  !is.finite(net_limit) ||
    net_limit == 0
) {
  net_limit <- 1
}

p_net_score <- ggplot(
  net_score,
  aes(
    x = receiver,
    y = sender,
    fill = net_score
  )
) +
  geom_tile(
    color = "white"
  ) +
  geom_text(
    aes(
      label = round(
        net_score,
        2
      )
    ),
    size = 3.5
  ) +
  scale_fill_gradient2(
    low = direction_colors[["platelet-free-up"]],
    mid = "white",
    high = direction_colors[["PLA-up"]],
    midpoint = 0,
    limits = c(
      -net_limit,
      net_limit
    )
  ) +
  theme_bw() +
  labs(
    title = plot_title(
      "Net MultiNicheNet prioritization score"
    ),
    x = "Receiver lineage",
    y = "Sender lineage",
    fill =
      "Mean PLA score\nminus mean\nplatelet-free score"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

save_plot(
  p_net_score,
  "03_multinichenet_net_score_heatmap.png",
  width = 9,
  height = 8
)

# ============================================================
# 14. recurrent LR pairs
# ============================================================

recurrent_lr <- mn_selected %>%
  dplyr::distinct(
    interaction,
    lineage_pair,
    direction
  ) %>%
  dplyr::count(
    interaction,
    direction,
    name = "n_lineage_pairs"
  ) %>%
  group_by(direction) %>%
  slice_max(
    n_lineage_pairs,
    n = 15,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  mutate(
    plot_label = paste(
      interaction,
      direction,
      sep = "___"
    ),
    plot_label = fct_reorder(
      plot_label,
      n_lineage_pairs
    )
  )

p_recurrent_lr <- ggplot(
  recurrent_lr,
  aes(
    x = n_lineage_pairs,
    y = plot_label,
    fill = direction
  )
) +
  geom_col() +
  facet_wrap(
    ~ direction,
    scales = "free_y"
  ) +
  scale_y_discrete(
    labels = function(x) {
      sub(
        "___.*$",
        "",
        x
      )
    }
  ) +
  scale_fill_manual(
    values = direction_colors
  ) +
  theme_bw() +
  labs(
    title = plot_title(
      "Most recurrent selected MultiNicheNet LR pairs"
    ),
    x = "# sender-receiver lineage pairs",
    y = "Ligand → receptor"
  ) +
  theme(
    legend.position = "none"
  )

save_plot(
  p_recurrent_lr,
  "04_multinichenet_recurrent_lr_pairs.png",
  width = 12,
  height = 9
)

# ============================================================
# 15. top LR pairs by mean score
# ============================================================

top_lr <- mn_plot %>%
  group_by(
    interaction,
    direction
  ) %>%
  summarise(
    mean_score = mean(
      prioritization_score,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  group_by(direction) %>%
  slice_max(
    mean_score,
    n = 15,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  mutate(
    plot_label = paste(
      interaction,
      direction,
      sep = "___"
    ),
    plot_label = fct_reorder(
      plot_label,
      mean_score
    )
  )

p_top_lr <- ggplot(
  top_lr,
  aes(
    x = mean_score,
    y = plot_label,
    fill = direction
  )
) +
  geom_col() +
  facet_wrap(
    ~ direction,
    scales = "free_y"
  ) +
  scale_y_discrete(
    labels = function(x) {
      sub(
        "___.*$",
        "",
        x
      )
    }
  ) +
  scale_fill_manual(
    values = direction_colors
  ) +
  theme_bw() +
  labs(
    title = plot_title(
      "Top MultiNicheNet LR pairs by mean prioritization score"
    ),
    x = "Mean prioritization score",
    y = "Ligand → receptor"
  ) +
  theme(
    legend.position = "none"
  )

save_plot(
  p_top_lr,
  "05_multinichenet_top_lr_by_mean_score.png",
  width = 12,
  height = 9
)

# ============================================================
# 16. top interactions in top lineage pairs
# ============================================================

n_top_pairs <- 5
n_top_interactions <- 5

# Rank sender-receiver pairs separately for each direction.
#
# n_selected_interactions is counted within mn_selected, i.e.
# among the top n_selected interactions per condition defined
# in section 12.
#
# Mean and maximum prioritization scores are used as tie-breakers.
top_lineage_pairs <- mn_selected %>%
  filter(
    !is.na(lineage_pair),
    !is.na(interaction)
  ) %>%
  group_by(
    direction,
    lineage_pair
  ) %>%
  summarise(
    n_selected_interactions = n(),
    mean_prioritization_score = mean(
      prioritization_score,
      na.rm = TRUE
    ),
    max_prioritization_score = max(
      prioritization_score,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  group_by(direction) %>%
  arrange(
    desc(n_selected_interactions),
    desc(mean_prioritization_score),
    desc(max_prioritization_score),
    .by_group = TRUE
  ) %>%
  slice_head(
    n = n_top_pairs
  ) %>%
  ungroup()

message(
  "Top lineage pairs selected for MultiNicheNet plot 06:"
)
print(top_lineage_pairs)

# Keep only interactions from the selected lineage pairs.
# Within each pair, retain the five highest-prioritization
# interactions.
top_per_pair <- mn_plot %>%
  filter(
    !is.na(lineage_pair),
    !is.na(interaction)
  ) %>%
  semi_join(
    top_lineage_pairs,
    by = c(
      "direction",
      "lineage_pair"
    )
  ) %>%
  group_by(
    direction,
    lineage_pair
  ) %>%
  slice_max(
    order_by = prioritization_score,
    n = n_top_interactions,
    with_ties = FALSE
  ) %>%
  arrange(
    direction,
    lineage_pair,
    desc(prioritization_score)
  ) %>%
  ungroup() %>%
  mutate(
    plot_label = paste(
      interaction,
      lineage_pair,
      direction,
      sep = "___"
    )
  )

# Arrange interaction labels so the largest score appears
# at the top of each panel.
top_per_pair <- top_per_pair %>%
  group_by(
    direction,
    lineage_pair
  ) %>%
  mutate(
    plot_label = factor(
      plot_label,
      levels = rev(
        unique(
          plot_label[
            order(
              prioritization_score,
              decreasing = TRUE
            )
          ]
        )
      )
    )
  ) %>%
  ungroup()

# Preserve the sender-receiver panel order from the ranking table.
pair_levels <- top_lineage_pairs %>%
  arrange(
    direction,
    desc(n_selected_interactions),
    desc(mean_prioritization_score),
    desc(max_prioritization_score)
  ) %>%
  pull(lineage_pair) %>%
  unique()

top_per_pair <- top_per_pair %>%
  mutate(
    lineage_pair = factor(
      lineage_pair,
      levels = pair_levels
    )
  )

p_top_per_pair <- ggplot(
  top_per_pair,
  aes(
    x = prioritization_score,
    y = plot_label,
    color = direction
  )
) +
  geom_point(
    size = 3,
    alpha = 0.85
  ) +
  facet_grid(
    rows = vars(direction),
    cols = vars(lineage_pair),
    scales = "free_y",
    space = "free_x",
    drop = TRUE
  ) +
  scale_y_discrete(
    labels = function(x) {
      sub(
        "___.*$",
        "",
        x
      )
    }
  ) +
  scale_color_manual(
    values = direction_colors
  ) +
  theme_bw(
    base_size = 10
  ) +
  labs(
    title = plot_title(
      "Top interactions in top 5 lineage pairs per direction"
    ),
    subtitle = paste0(
      "Lineage pairs ranked by representation among the top ",
      n_selected,
      " interactions per condition; ",
      "interactions ranked by descending prioritization score"
    ),
    x = "Prioritization score",
    y = "Ligand → receptor",
    color = "Direction"
  ) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(
      size = 8
    ),
    axis.text.y = element_text(
      size = 7
    )
  )

save_plot(
  p_top_per_pair,
  "06_multinichenet_top_interactions_per_lineage_pair.png",
  width = 18,
  height = 10
)

# ============================================================
# 17. top individual interactions
# ============================================================

top_individual <- mn_plot %>%
  group_by(direction) %>%
  slice_max(
    prioritization_score,
    n = 20,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  mutate(
    plot_label = paste0(
      interaction,
      "\n",
      lineage_pair
    ),
    plot_label = paste(
      plot_label,
      direction,
      sep = "___"
    ),
    plot_label = fct_reorder(
      plot_label,
      prioritization_score
    )
  )

p_top_individual <- ggplot(
  top_individual,
  aes(
    x = prioritization_score,
    y = plot_label,
    color = direction
  )
) +
  geom_point(
    size = 3,
    alpha = 0.85
  ) +
  facet_wrap(
    ~ direction,
    scales = "free_y"
  ) +
  scale_y_discrete(
    labels = function(x) {
      sub(
        "___.*$",
        "",
        x
      )
    }
  ) +
  scale_color_manual(
    values = direction_colors
  ) +
  theme_bw() +
  labs(
    title = plot_title(
      "Top individual MultiNicheNet interactions"
    ),
    x = "Prioritization score",
    y =
      "Ligand → receptor / sender → receiver"
  ) +
  theme(
    legend.position = "none"
  )

save_plot(
  p_top_individual,
  "07_multinichenet_top_individual_interactions.png",
  width = 13,
  height = 10
)

# ============================================================
# 18. circos plots
# ============================================================

top_circos <- get_top_n_lr_pairs(
  prioritization_tables,
  top_n = 20,
  rank_per_group = TRUE
)

circos_input <- prioritization_tables$
  group_prioritization_tbl %>%
  filter(
    id %in% top_circos$id
  ) %>%
  distinct(
    id,
    sender,
    receiver,
    ligand,
    receptor,
    group,
    .keep_all = TRUE
  )

circos_input$prioritization_score[
  is.na(circos_input$prioritization_score)
] <- 0

circos_lineages <- union(
  unique(circos_input$sender),
  unique(circos_input$receiver)
)

missing_colors <- setdiff(
  circos_lineages,
  names(lineage_colors)
)

if (length(missing_colors) > 0) {
  stop(
    "Missing lineage colors for: ",
    paste(
      missing_colors,
      collapse = ", "
    )
  )
}

circos_filename <- file.path(
  plot_dir,
  "09_multinichenet_circos_%02d.png"
)

png(
  filename = circos_filename,
  width = 3500,
  height = 3500,
  res = 250
)

tryCatch(
  {
    make_circos_group_comparison(
      circos_input,
      colors_sender =
        lineage_colors,
      colors_receiver =
        lineage_colors
    )
  },
  error = function(e) {
    message(
      "MultiNicheNet circos plot failed: ",
      conditionMessage(e)
    )
  },
  finally = {
    circos.clear()

    if (dev.cur() > 1) {
      dev.off()
    }
  }
)

message(
  "Finished MultiNicheNet analysis for: ",
  dataset_name
)
message(
  "Saved plots to: ",
  plot_dir
)
