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
  library(magrittr)
})

# ============================================================
# 1. arguments
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript multinichetr.R <input_rds> <output_dir>")
}

input_file <- args[[1]]
base_output_dir <- args[[2]]

dataset_name <- tools::file_path_sans_ext(basename(input_file))
dataset_mode <- stringr::str_extract(
  dataset_name,
  "(all|diseasedOnly|healthyOnly)$"
)

dataset_clean <- stringr::str_remove(
  dataset_name,
  "_(all|diseasedOnly|healthyOnly)$"
)

plot_title <- function(title) {
  paste0(title, "\n", dataset_clean, " | ", dataset_mode)
}

out_dir <- file.path(base_output_dir, dataset_name)
table_dir <- file.path(out_dir, "tables")
plot_dir <- file.path(out_dir, "plots")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

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

symmetric_limit <- function(x, fallback = 1) {
  value <- suppressWarnings(
    max(abs(x), na.rm = TRUE)
  )

  if (!is.finite(value) || value <= 0) {
    return(fallback)
  }

  value
}

message("Dataset: ", dataset_name)

# ============================================================
# 2. load and prepare metadata
# ============================================================

seurat_obj <- readRDS(input_file)

condition_col <- "pla_status"
lineage_col <- "lineage"


# ============================================================
# shared plotting configuration
# ============================================================

direction_colors <- c(
  "PLA-up" = "#F8766D",
  "platelet-free-up" = "#00BFC4"
)

direction_levels <- c(
  "PLA-up",
  "platelet-free-up"
)

all_dataset_lineages <- sort(
  unique(
    as.character(
      seurat_obj@meta.data[[lineage_col]]
    )
  )
)

all_dataset_lineages <- all_dataset_lineages[
  !is.na(all_dataset_lineages) &
    all_dataset_lineages != ""
]

if (length(all_dataset_lineages) == 0) {
  stop("No valid lineage values found.")
}

base_palette <- RColorBrewer::brewer.pal(
  n = 11,
  name = "Spectral"
)

global_lineage_colors <- grDevices::colorRampPalette(
  base_palette
)(length(all_dataset_lineages))

names(global_lineage_colors) <- all_dataset_lineages

lineage_color_table <- tibble::tibble(
  lineage = names(global_lineage_colors),
  color = unname(global_lineage_colors)
)

common_theme <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold"
      ),
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1
      ),
      panel.grid.minor = ggplot2::element_blank()
    )
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]

  if (length(x) == 0) {
    return(NA_real_)
  }

  mean(x)
}

# ============================================================
# 2. load and prepare metadata
# ============================================================

seurat_obj <- readRDS(input_file)

condition_col <- "pla_status"
lineage_col <- "lineage"

# ============================================================
# dataset-specific sample and pairing configuration
# ============================================================

if (grepl("ImmuneAging", dataset_name)) {

  # One donor contributes PLA and platelet-free cells
  sample_col <- "donor_id"
  pair_col <- "donor_id"
  paired_analysis <- TRUE

} else if (grepl("our_dataset", dataset_name)) {

  # sample_ID identifies the condition-specific sample,
  # patient links PLA and platelet-free
  sample_col <- "sample_ID"
  pair_col <- "patient"
  paired_analysis <- TRUE

} else {

  # For heart, sepsis, vaccine and skin, the same sample ID
  # occurs in PLA and platelet-free
  sample_col <- "sample"
  pair_col <- "sample"
  paired_analysis <- TRUE
}

message("Using sample column: ", sample_col)
message("Using pairing column: ", pair_col)
message("Paired analysis: ", paired_analysis)

required_cols <- unique(c(
  sample_col,
  pair_col,
  condition_col,
  lineage_col
))

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
  seurat_obj@meta.data[, required_cols, drop = FALSE]
)

seurat_obj <- subset(
  seurat_obj,
  cells = colnames(seurat_obj)[keep]
)

DefaultAssay(seurat_obj) <- "RNA"

required_conditions <- c("PLA", "platelet-free")

available_conditions <- unique(
  as.character(seurat_obj@meta.data[[condition_col]])
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

# ============================================================
# create analysis identifiers
# ============================================================

# Unique pseudobulk unit:
# one sample/patient/donor in one condition
seurat_obj$sample_condition <- interaction(
  as.character(seurat_obj@meta.data[[sample_col]]),
  as.character(seurat_obj@meta.data[[condition_col]]),
  sep = "_",
  drop = TRUE
)

# Biological pairing identifier
seurat_obj$pair_id <- as.character(
  seurat_obj@meta.data[[pair_col]]
)

# ============================================================
# validate identifiers
# ============================================================

# each pseudobulk ID must belong to exactly one condition
sample_group_check <- seurat_obj@meta.data %>%
  dplyr::distinct(
    sample_condition,
    pla_status = .data[[condition_col]]
  ) %>%
  dplyr::count(sample_condition, name = "n_groups")

if (any(sample_group_check$n_groups != 1)) {
  stop(
    "At least one sample_condition is assigned to more than one condition."
  )
}

# check how many pair IDs occur in both conditions
pairing_check <- seurat_obj@meta.data %>%
  dplyr::distinct(
    pair_id,
    condition = .data[[condition_col]]
  ) %>%
  dplyr::count(pair_id, name = "n_conditions")

n_complete_pairs <- sum(pairing_check$n_conditions == 2)
n_incomplete_pairs <- sum(pairing_check$n_conditions < 2)

message("Complete pair IDs: ", n_complete_pairs)
message("Incomplete pair IDs: ", n_incomplete_pairs)

if (paired_analysis && n_complete_pairs < 2) {
  stop(
    "Paired analysis set, but fewer than two complete pairs are present."
  )
}

# ============================================================
# analysis settings
# ============================================================

organism <- "human"

sample_id <- "sample_condition"
group_id <- condition_col
celltype_id <- lineage_col

covariates <- if (paired_analysis) "pair_id" else NA
batches <- NA

min_cells_config <- c(
  "gated_heart_processed_all" = 3L,
  "gated_heart_processed_diseasedOnly" = 3L,
  "gated_heart_processed_healthyOnly" = NA_integer_,

  "gated_ImmuneAging_all" = 8L,
  "gated_ImmuneAging_diseasedOnly" = NA_integer_,
  "gated_ImmuneAging_healthyOnly" = 8L,

  "gated_our_dataset_processed_all" = 8L,
  "gated_our_dataset_processed_diseasedOnly" = 8L,
  "gated_our_dataset_processed_healthyOnly" = 5L,

  "gated_sepsis_processed_all" = 20L,
  "gated_sepsis_processed_diseasedOnly" = 20L,
  "gated_sepsis_processed_healthyOnly" = 15L,

  "gated_vaccine_processed_all" = 5L,
  "gated_vaccine_processed_diseasedOnly" = 5L,
  "gated_vaccine_processed_healthyOnly" = 5L,

  "gated_skin_processed_all" = 3L,
  "gated_skin_processed_diseasedOnly" = 3L,
  "gated_skin_processed_healthyOnly" = 5L
)

if (!dataset_name %in% names(min_cells_config)) {
  stop("No min_cells configured for dataset: ", dataset_name)
}

min_cells <- as.numeric(
  unname(min_cells_config[[dataset_name]])
)
if (is.na(min_cells)) {
  message(
    "Skipping ",
    dataset_name,
    ": insufficient paired samples."
  )
  quit(save = "no", status = 0)
}

message("Using min_cells: ", min_cells)

min_sample_prop <- 0.50 #gene must pass expression criterion (fraction cutoff) in at least 50% of eligible samples for cell type
fraction_cutoff <- 0.05 # gene is considered expressed in particular sample x cell type combination when it has non-0 expression in at least 5% of cells of that cell type in sample

logFC_threshold <- 0.50
p_val_threshold <- 0.05
p_val_adj <- TRUE
top_n_target <- 250

options(timeout = 120)

# ============================================================
# 4. resources
# ============================================================

lr_network <- readRDS(url(
  "https://zenodo.org/record/10229222/files/lr_network_human_allInfo_30112033.rds"
)) %>%
  transmute(
    ligand = make.names(
      convert_alias_to_symbols(ligand, organism = organism)
    ),
    receptor = make.names(
      convert_alias_to_symbols(receptor, organism = organism)
    )
  ) %>%
  distinct()

ligand_target_matrix <- readRDS(url(
  "https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"
))

colnames(ligand_target_matrix) <- colnames(ligand_target_matrix) %>%
  convert_alias_to_symbols(organism = organism) %>%
  make.names()

rownames(ligand_target_matrix) <- rownames(ligand_target_matrix) %>%
  convert_alias_to_symbols(organism = organism) %>%
  make.names()

lr_network <- lr_network %>%
  filter(ligand %in% colnames(ligand_target_matrix))

ligand_target_matrix <- ligand_target_matrix[
  ,
  unique(lr_network$ligand),
  drop = FALSE
]

# ============================================================
# 5. convert to SCE
# ============================================================

sce <- Seurat::as.SingleCellExperiment(
  seurat_obj,
  assay = "RNA"
)

sce <- alias_to_symbol_SCE(sce, organism) %>%
  makenames_SCE()

for (column in c(sample_id, group_id, celltype_id, covariates)) {
  colData(sce)[[column]] <- factor(
    make.names(
      as.character(colData(sce)[[column]])
    )
  )
}

celltypes_oi <- unique(
  as.character(colData(sce)[[celltype_id]])
)

# ============================================================
# 6. abundance and expression
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

# keep genes expressed by at least one cell type
genes_oi <- frq_list$expressed_df %>%
  filter(expressed) %>%
  pull(gene) %>%
  unique()

if (length(genes_oi) == 0) {
  stop("No expressed genes passed filtering.")
}

sce <- sce[genes_oi, ]

abundance_expression_info <- process_abundance_expression_info(
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
    c("gg", "ggplot", "patchwork")
  )
) {
  save_plot(
    DE_info$hist_pvals,
    "14_deinfo_hist_pvals.png",
    width = 12,
    height = 8
  )
}


celltype_de <- DE_info$celltype_de$de_output_tidy

if (is.null(celltype_de) || nrow(celltype_de) == 0) {
  stop("No cell-type DE results were generated.")
}

included_celltypes <- unique(celltype_de$cluster_id)

sender_receiver_de <- combine_sender_receiver_de(
  sender_de = celltype_de,
  receiver_de = celltype_de,
  senders_oi = included_celltypes,
  receivers_oi = included_celltypes,
  lr_network = lr_network
)

if (is.null(sender_receiver_de) || nrow(sender_receiver_de) == 0) {
  stop("No ligand-receptor combinations were generated.")
}

# ============================================================
# 8. ligand activity
# ============================================================

ligand_activities <- suppressMessages(
  suppressWarnings(
    get_ligand_activities_targets_DEgenes(
      receiver_de = celltype_de,
      receivers_oi = included_celltypes,
      ligand_target_matrix = ligand_target_matrix,
      logFC_threshold = logFC_threshold,
      p_val_threshold = p_val_threshold,
      p_val_adj = p_val_adj,
      top_n_target = top_n_target
    )
  )
)

if (is.null(ligand_activities) || NROW(ligand_activities) == 0) {
  stop(
    "No ligand activities were generated. ",
    "Check DE thresholds and included receiver cell types."
  )
}

# ============================================================
# 9. prioritization
# ============================================================

sender_receiver_tbl <- sender_receiver_de %>%
  distinct(sender, receiver)

grouping_tbl <- as_tibble(colData(sce)) %>%
  distinct(
    sample = .data[[sample_id]],
    group = .data[[group_id]]
  )

prioritization_tables <- generate_prioritization_tables(
  sender_receiver_info =
    abundance_expression_info$sender_receiver_info,
  sender_receiver_de = sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities,
  contrast_tbl = contrast_tbl,
  sender_receiver_tbl = sender_receiver_tbl,
  grouping_tbl = grouping_tbl,
  scenario = "regular",
  fraction_cutoff = fraction_cutoff,
  abundance_data_receiver =
    abundance_expression_info$abundance_data_receiver,
  abundance_data_sender =
    abundance_expression_info$abundance_data_sender,
  ligand_activity_down = FALSE
)

# ============================================================
# 10. save
# ============================================================

multinichenet_output <- list(
  dataset_name = dataset_name,
  settings = list(
    sample_col = sample_col,
    pair_col = pair_col,
    min_cells = min_cells,
    min_sample_prop = min_sample_prop,
    fraction_cutoff = fraction_cutoff
  ),
  abundance_info = abundance_info,
  expressed_df = frq_list$expressed_df,
  abundance_expression_info = abundance_expression_info,
  DE_info = DE_info,
  celltype_de = celltype_de,
  sender_receiver_de = sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities,
  prioritization_tables = prioritization_tables,
  grouping_tbl = grouping_tbl
)

settings_table <- tibble(
  dataset = dataset_name,
  sample_col = sample_col,
  pair_col = pair_col,
  min_cells = min_cells,
  min_sample_prop = min_sample_prop,
  fraction_cutoff = fraction_cutoff,
  logFC_threshold = logFC_threshold,
  p_val_threshold = p_val_threshold,
  p_val_adj = p_val_adj
)

# ============================================================
# MultiNicheNet standardized plotting table
# ============================================================

mn_prio_plot <- prioritization_tables$group_prioritization_tbl %>%
  tibble::as_tibble() %>%
  dplyr::mutate(
    lineage_pair = paste(sender, receiver, sep = " → "),
    interaction = paste(ligand, receptor, sep = " → "),

    direction = dplyr::case_when(
      group == "PLA" ~ "PLA-up",
      group == "platelet.free" ~ "platelet-free-up",
      TRUE ~ as.character(group)
    ),

    direction = factor(
      direction,
      levels = direction_levels
    ),

    lfc_score = (
      scaled_lfc_ligand +
        scaled_lfc_receptor
    ) / 2,

    expression_score = (
      scaled_pb_ligand +
        scaled_pb_receptor
    ) / 2,

    activity_score = max_scaled_activity,
    score = prioritization_score
  ) %>%
  dplyr::filter(
    direction %in% direction_levels,
    is.finite(prioritization_score)
  )

# ============================================================
# Select high-priority interactions for count-based plots
# ============================================================

top_n_comparison_per_direction <- 500L

available_per_direction <- mn_prio_plot %>%
  dplyr::count(direction, name = "available_n")

actual_top_n <- min(
  top_n_comparison_per_direction,
  min(available_per_direction$available_n)
)

message(
  "Using top ",
  actual_top_n,
  " MultiNicheNet interactions per direction ",
  "for count-based comparison plots."
)

mn_selected <- mn_prio_plot %>%
  dplyr::group_by(direction) %>%
  dplyr::slice_max(
    order_by = prioritization_score,
    n = actual_top_n,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup()

# ============================================================
# 1. prioritization-score distribution
# ============================================================

p_mn_score_distribution <- ggplot(
  mn_prio_plot,
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
    outlier.size = 0.3,
    alpha = 0.8
  ) +
  scale_fill_manual(
    values = direction_colors,
    drop = FALSE
  ) +
  common_theme() +
  labs(
    title = plot_title(
      "MultiNicheNet prioritization-score distribution"
    ),
    x = "Prioritized condition",
    y = "Prioritization score",
    fill = "Direction"
  ) +
  theme(
    legend.position = "none"
  )

save_plot(
  p_mn_score_distribution,
  "01_multinichenet_score_distribution.png",
  width = 7,
  height = 5
)

# ============================================================
# 2. selected interactions by lineage pair
# ============================================================

pair_counts_mn <- mn_selected %>%
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
    order_by = total,
    n = 25,
    with_ties = FALSE
  ) %>%
  dplyr::mutate(
    lineage_pair = forcats::fct_reorder(
      lineage_pair,
      total
    )
  )

p_pair_counts_mn <- ggplot(
  pair_counts_mn,
  aes(
    x = n_interactions,
    y = lineage_pair,
    fill = direction
  )
) +
  geom_col() +
  scale_fill_manual(
    values = direction_colors,
    drop = FALSE
  ) +
  common_theme() +
  labs(
    title = plot_title(
      paste0(
        "Top ",
        actual_top_n,
        " MultiNicheNet interactions per condition by lineage pair"
      )
    ),
    x = "# selected prioritized interactions",
    y = "Sender → receiver lineage",
    fill = "Direction"
  ) +
  theme(
    axis.text.x = element_text(angle = 0)
  )

save_plot(
  p_pair_counts_mn,
  "02_multinichenet_selected_interactions_by_lineage_pair.png",
  width = 10,
  height = 8
)

# ============================================================
# 3. net mean prioritization-score heatmap
# ============================================================

net_score_mn <- mn_prio_plot %>%
  dplyr::mutate(
    direction_key = dplyr::case_when(
      direction == "PLA-up" ~ "PLA_up",
      direction == "platelet-free-up" ~ "platelet_free_up"
    )
  ) %>%
  dplyr::group_by(
    sender,
    receiver,
    direction_key
  ) %>%
  dplyr::summarise(
    mean_score = safe_mean(prioritization_score),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = direction_key,
    values_from = mean_score
  ) %>%
  dplyr::filter(
    !is.na(PLA_up),
    !is.na(platelet_free_up)
  ) %>%
  dplyr::mutate(
    net_score = PLA_up - platelet_free_up
  )

net_score_limit <- symmetric_limit(
  net_score_mn$net_score
)

p_net_score_mn <- ggplot(
  net_score_mn,
  aes(
    x = receiver,
    y = sender,
    fill = net_score
  )
) +
  geom_tile(color = "white") +
  geom_text(
    aes(label = round(net_score, 2)),
    size = 3.5
  ) +
  scale_fill_gradient2(
    low = direction_colors[["platelet-free-up"]],
    mid = "white",
    high = direction_colors[["PLA-up"]],
    midpoint = 0,
    limits = c(
      -net_score_limit,
      net_score_limit
    )
  ) +
  common_theme() +
  labs(
    title = plot_title(
      "Net MultiNicheNet prioritization score"
    ),
    x = "Receiver lineage",
    y = "Sender lineage",
    fill = paste0(
      "Mean PLA score\nminus mean\n",
      "platelet-free score"
    )
  )

save_plot(
  p_net_score_mn,
  "03_multinichenet_net_mean_score_heatmap.png",
  width = 9,
  height = 8
)

# ============================================================
# 4. recurrent prioritized LR pairs
# ============================================================

recurrent_mn <- mn_selected %>%
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
  dplyr::group_by(direction) %>%
  dplyr::slice_max(
    order_by = n_lineage_pairs,
    n = 15,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    interaction_facet = paste(
      interaction,
      direction,
      sep = "___"
    ),
    interaction_facet = forcats::fct_reorder(
      interaction_facet,
      n_lineage_pairs
    )
  )

p_recurrent_mn <- ggplot(
  recurrent_mn,
  aes(
    x = n_lineage_pairs,
    y = interaction_facet,
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
      sub("___.*$", "", x)
    }
  ) +
  scale_fill_manual(
    values = direction_colors,
    drop = FALSE
  ) +
  common_theme() +
  labs(
    title = plot_title(
      "Most recurrent selected MultiNicheNet LR pairs"
    ),
    x = "# sender–receiver lineage pairs",
    y = "Ligand → receptor",
    fill = "Direction"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 0)
  )

save_plot(
  p_recurrent_mn,
  "04_multinichenet_top_recurrent_lr_pairs.png",
  width = 12,
  height = 9
)

# ============================================================
# 5. top LR pairs by mean prioritization score
# ============================================================

top_lr_by_score_mn <- mn_prio_plot %>%
  dplyr::group_by(
    interaction,
    direction
  ) %>%
  dplyr::summarise(
    n_lineage_pairs = dplyr::n_distinct(lineage_pair),
    mean_score = safe_mean(prioritization_score),
    max_score = max(
      prioritization_score,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  dplyr::group_by(direction) %>%
  dplyr::slice_max(
    order_by = mean_score,
    n = 15,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    interaction_facet = paste(
      interaction,
      direction,
      sep = "___"
    ),
    interaction_facet = forcats::fct_reorder(
      interaction_facet,
      mean_score
    )
  )

p_top_lr_by_score_mn <- ggplot(
  top_lr_by_score_mn,
  aes(
    x = mean_score,
    y = interaction_facet,
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
      sub("___.*$", "", x)
    }
  ) +
  scale_fill_manual(
    values = direction_colors,
    drop = FALSE
  ) +
  common_theme() +
  labs(
    title = plot_title(
      "Top MultiNicheNet LR pairs by mean prioritization score"
    ),
    x = "Mean prioritization score",
    y = "Ligand → receptor",
    fill = "Direction"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 0)
  )

save_plot(
  p_top_lr_by_score_mn,
  "05_multinichenet_top_lr_by_mean_score_faceted.png",
  width = 12,
  height = 9
)

# ============================================================
# 6. top interactions per lineage pair
# ============================================================

top_mn_per_pair <- mn_prio_plot %>%
  dplyr::group_by(
    direction,
    lineage_pair
  ) %>%
  dplyr::slice_max(
    order_by = prioritization_score,
    n = 5,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup()

p_top_mn_per_pair <- ggplot(
  top_mn_per_pair,
  aes(
    x = lineage_pair,
    y = interaction,
    size = prioritization_score,
    color = direction
  )
) +
  geom_point(alpha = 0.8) +
  facet_wrap(
    ~ direction,
    scales = "free"
  ) +
  scale_color_manual(
    values = direction_colors,
    drop = FALSE
  ) +
  common_theme(base_size = 9) +
  labs(
    title = plot_title(
      "Top MultiNicheNet interactions per lineage pair"
    ),
    x = "Sender → receiver lineage",
    y = "Ligand → receptor",
    size = "Prioritization score",
    color = "Direction"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 6
    ),
    axis.text.y = element_text(size = 6)
  )

save_plot(
  p_top_mn_per_pair,
  "06_multinichenet_top_interactions_per_lineage_pair.png",
  width = 16,
  height = 11
)

# ============================================================
# 7. strongest individual interactions
# ============================================================

top_individual_mn <- mn_prio_plot %>%
  dplyr::group_by(direction) %>%
  dplyr::slice_max(
    order_by = prioritization_score,
    n = 20,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    label = paste0(
      interaction,
      "\n",
      lineage_pair
    ),
    label_facet = paste(
      label,
      direction,
      sep = "___"
    ),
    label_facet = forcats::fct_reorder(
      label_facet,
      prioritization_score
    )
  )

p_top_individual_mn <- ggplot(
  top_individual_mn,
  aes(
    x = prioritization_score,
    y = label_facet,
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
      sub("___.*$", "", x)
    }
  ) +
  scale_color_manual(
    values = direction_colors,
    drop = FALSE
  ) +
  common_theme() +
  labs(
    title = plot_title(
      "Top individual MultiNicheNet interactions"
    ),
    x = "Prioritization score",
    y = "Ligand → receptor / sender → receiver",
    color = "Direction"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 0)
  )

save_plot(
  p_top_individual_mn,
  "07_multinichenet_top_individual_interactions.png",
  width = 13,
  height = 10
)

# ============================================================
# 8. MultiNicheNet circos plots
# ============================================================

top_n_circos <- 20L

mn_circos_tbl <- mn_prio_plot %>%
  dplyr::filter(
    !is.na(sender),
    !is.na(receiver),
    !is.na(ligand),
    !is.na(receptor),
    ligand != "",
    receptor != "",
    is.finite(prioritization_score)
  ) %>%
  dplyr::group_by(direction) %>%
  dplyr::slice_max(
    order_by = prioritization_score,
    n = top_n_circos,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(direction) %>%
  dplyr::arrange(
    dplyr::desc(prioritization_score),
    .by_group = TRUE
  ) %>%
  dplyr::mutate(
    rank_within_direction = dplyr::row_number()
  ) %>%
  dplyr::ungroup()

# Convert labels back to the group names expected by
# make_circos_group_comparison().
mn_circos_input <- mn_circos_tbl %>%
  dplyr::mutate(
    group = dplyr::case_when(
      direction == "PLA-up" ~ "PLA",
      direction == "platelet-free-up" ~ "platelet.free",
      TRUE ~ as.character(direction)
    )
  )

circos_filename <- file.path(
  plot_dir,
  "08_multinichenet_circos_%02d.png"
)

grDevices::png(
  filename = circos_filename,
  width = 5000,
  height = 5000,
  res = 300
)

circos_list <- tryCatch(
  {
    multinichenetr::make_circos_group_comparison(
      mn_circos_input,
      colors_sender = global_lineage_colors,
      colors_receiver = global_lineage_colors
    )
  },
  error = function(e) {
    message(
      "MultiNicheNet circos plot failed: ",
      conditionMessage(e)
    )
    NULL
  },
  finally = {
    circlize::circos.clear()

    if (grDevices::dev.cur() > 1) {
      grDevices::dev.off()
    }
  }
)