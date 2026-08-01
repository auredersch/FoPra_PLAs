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

save_table <- function(x, filename) {
  write.csv(
    x,
    file.path(table_dir, filename),
    row.names = FALSE
  )
}

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

message("Dataset: ", dataset_name)

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

saveRDS(
  multinichenet_output,
  file.path(out_dir, "multinichenet_output.rds")
)

save_table(
  prioritization_tables$group_prioritization_tbl,
  "group_prioritization_tbl.csv"
)

save_table(
  celltype_de,
  "celltype_de.csv"
)

save_table(
  sender_receiver_de,
  "sender_receiver_de.csv"
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

save_table(
  settings_table,
  "analysis_settings.csv"
)

# ============================================================
# 12. chord plots
# ============================================================

top_n_circos <- 20

prioritized_tbl_oi_all <- get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables,
  top_n = top_n_circos,
  rank_per_group = TRUE
)

prioritized_tbl_oi <- multinichenet_output$
  prioritization_tables$
  group_prioritization_tbl %>%
  dplyr::filter(
    id %in% prioritized_tbl_oi_all$id
  ) %>%
  dplyr::distinct(
    id,
    sender,
    receiver,
    ligand,
    receptor,
    group
  ) %>%
  dplyr::left_join(
    prioritized_tbl_oi_all,
    by = c(
      "id",
      "sender",
      "receiver",
      "ligand",
      "receptor",
      "group"
    )
  ) %>%
  dplyr::mutate(
    prioritization_score = tidyr::replace_na(
      prioritization_score,
      0
    )
  )

if (nrow(prioritized_tbl_oi) == 0) {

  message(
    "Skipping MultiNicheNet chord plots: ",
    "no prioritized interactions."
  )

} else {

  senders_receivers <- sort(
    union(
      unique(prioritized_tbl_oi$sender),
      unique(prioritized_tbl_oi$receiver)
    )
  )

  base_palette <- RColorBrewer::brewer.pal(
    n = 11,
    name = "Spectral"
  )

  celltype_colors <- grDevices::colorRampPalette(
    base_palette
  )(length(senders_receivers))

  names(celltype_colors) <- senders_receivers

  colors_sender <- celltype_colors
  colors_receiver <- celltype_colors

  # %02d causes separate PNG files to be created for each
  # graphical page produced by make_circos_group_comparison().
  circos_filename <- file.path(
    plot_dir,
    "03_circos_group_comparison_%02d.png"
  )

  message(
    "Saving MultiNicheNet chord plot pages to: ",
    circos_filename
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
        prioritized_tbl_oi,
        colors_sender,
        colors_receiver
      )
    },
    error = function(e) {
      message(
        "MultiNicheNet chord plot failed: ",
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

  if (!is.null(circos_list)) {

    saveRDS(
      circos_list,
      file.path(
        out_dir,
        "circos_group_comparison.rds"
      )
    )

    message(
      "Saved MultiNicheNet chord plot output."
    )

    generated_circos_files <- list.files(
      plot_dir,
      pattern = "^03_circos_group_comparison_[0-9]+\\.png$",
      full.names = TRUE
    )

    message(
      "Generated chord plot files: ",
      paste(
        basename(generated_circos_files),
        collapse = ", "
      )
    )
  }
}
# ============================================================
# 13. ligand-receptor product/activity plot
# ============================================================

group_oi <- "PLA"

prioritized_tbl_oi_PLA_50 <- get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables,
  top_n = 50,
  groups_oi = group_oi
)

save_table(
  prioritized_tbl_oi_PLA_50,
  "top_50_lr_pairs_PLA.csv"
)

if (nrow(prioritized_tbl_oi_PLA_50) > 0) {

  plot_oi <- make_sample_lr_prod_activity_plots(
    multinichenet_output$prioritization_tables,
    prioritized_tbl_oi_PLA_50
  )

  saveRDS(
    plot_oi,
    file.path(out_dir, "sample_lr_prod_activity_plots_PLA.rds")
  )

  if (inherits(plot_oi, c("gg", "ggplot", "patchwork"))) {
    save_plot(
      plot_oi,
      "04_sample_lr_prod_activity_PLA.png",
      width = 12,
      height = 8
    )
  }

} else {
  message("Skipping PLA product/activity plot: no PLA interactions.")
}

# ============================================================
# 14. plotting table
# ============================================================

mn_prio_plot <- prioritization_tables$group_prioritization_tbl %>%
  dplyr::mutate(
    lineage_pair = paste(sender, receiver, sep = " → "),
    interaction = paste(ligand, receptor, sep = " → "),
    direction = dplyr::case_when(
      group == "PLA" ~ "PLA-up",
      group == "platelet.free" ~ "platelet-free-up",
      TRUE ~ as.character(group)
    ),
    lfc_score = (scaled_lfc_ligand + scaled_lfc_receptor) / 2,
    expression_score = (scaled_pb_ligand + scaled_pb_receptor) / 2,
    activity_score = max_scaled_activity,
    score = prioritization_score
  )

# ============================================================
# 15. global prioritized interaction counts
# ============================================================

p_global_mn_prio <- ggplot(mn_prio_plot, aes(x = direction, fill = direction)) +
  geom_bar() +
  theme_bw() +
  scale_fill_manual(
    values = c(
      "PLA-up" = "#F8766D",
      "platelet-free-up" = "#00BFC4"
    )
  ) +
  labs(
    title = plot_title("Global distribution of prioritized MultiNicheNet interactions"),
    x = "Prioritized condition",
    y = "# prioritized LR interactions"
  ) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

save_plot(
  p_global_mn_prio,
  "05_global_distribution_prioritized_multinichenet.png",
  width = 7,
  height = 5
)

# ============================================================
# 16. prioritized interactions by lineage pair
# ============================================================

pair_counts_mn_prio <- mn_prio_plot %>%
  dplyr::count(lineage_pair, direction) %>%
  dplyr::group_by(lineage_pair) %>%
  dplyr::mutate(total = sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::slice_max(total, n = 25, with_ties = FALSE) %>%
  dplyr::mutate(lineage_pair = fct_reorder(lineage_pair, total))

p_pair_counts_mn_prio <- ggplot(
  pair_counts_mn_prio,
  aes(x = n, y = lineage_pair, fill = direction)
) +
  geom_col() +
  theme_bw() +
  scale_fill_manual(
    values = c(
      "PLA-up" = "#F8766D",
      "platelet-free-up" = "#00BFC4"
    )
  ) +
  labs(
    title = plot_title("Prioritized MultiNicheNet interactions by lineage pair"),
    x = "# prioritized LR interactions",
    y = "Sender → receiver lineage",
    fill = "Direction"
  )

save_plot(
  p_pair_counts_mn_prio,
  "06_prioritized_interactions_by_lineage_pair.png",
  width = 9,
  height = 7
)

# ============================================================
# 18. net prioritization score heatmap
# ============================================================

net_score_mn_prio <- mn_prio_plot %>%
  dplyr::group_by(sender, receiver, direction) %>%
  dplyr::summarise(
    mean_score = mean(prioritization_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::complete(
    sender,
    receiver,
    direction = c("PLA-up", "platelet-free-up"),
    fill = list(mean_score = 0)
  ) %>%
  tidyr::pivot_wider(
    names_from = direction,
    values_from = mean_score,
    values_fill = 0
  ) %>%
  dplyr::mutate(
    net_score = `PLA-up` - `platelet-free-up`
  )

p_net_score_mn_prio <- ggplot(
  net_score_mn_prio,
  aes(x = receiver, y = sender, fill = net_score)
) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(net_score, 2)), size = 4) +
  theme_bw() +
  labs(
    title = plot_title("Net prioritization score of MultiNicheNet interactions"),
    x = "Receiver lineage",
    y = "Sender lineage",
    fill = "Mean PLA score minus\nmean platelet-free score"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(
  p_net_score_mn_prio,
  "08_net_score_mn_prio_heatmap.png",
  width = 8,
  height = 7
)

# ============================================================
# 19. top prioritized LR pairs by direction
# ============================================================

top_lr_by_score <- mn_prio_plot %>%
  dplyr::filter(
    direction %in% c(
      "PLA-up",
      "platelet-free-up"
    )
  ) %>%
  dplyr::group_by(
    interaction,
    direction
  ) %>%
  dplyr::summarise(
    n_lineage_pairs = dplyr::n_distinct(lineage_pair),
    mean_score = mean(
      prioritization_score,
      na.rm = TRUE
    ),
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

p_top_lr_by_score <- ggplot(
  top_lr_by_score,
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
    values = c(
      "PLA-up" = "#F8766D",
      "platelet-free-up" = "#00BFC4"
    )
  ) +
  theme_bw() +
  labs(
    title = plot_title(
      "Top prioritized MultiNicheNet LR pairs by mean score"
    ),
    x = "Mean prioritization score",
    y = "Ligand → receptor",
    fill = "Direction"
  ) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 8)
  )

save_plot(
  p_top_lr_by_score,
  "09_top_lr_by_mean_score_faceted.png",
  width = 12,
  height = 9
)

# ============================================================
# 20. top interactions per lineage pair
# ============================================================

top_mn_prio_per_pair <- mn_prio_plot %>%
  dplyr::group_by(lineage_pair) %>%
  dplyr::slice_max(prioritization_score, n = 5, with_ties = FALSE) %>%
  dplyr::ungroup()

p_top_mn_prio_per_pair <- ggplot(
  top_mn_prio_per_pair,
  aes(
    x = lineage_pair,
    y = interaction,
    size = prioritization_score,
    color = direction
  )
) +
  geom_point(alpha = 0.8) +
  theme_bw() +
  scale_color_manual(
    values = c(
      "PLA-up" = "#F8766D",
      "platelet-free-up" = "#00BFC4"
    )
  ) +
  labs(
    title = plot_title("Top prioritized MultiNicheNet interactions per lineage pair"),
    x = "Sender → receiver lineage",
    y = "Ligand → receptor",
    size = "Prioritization score",
    color = "Direction"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    plot.title = element_text(size = 10),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7)
  )

save_plot(
  p_top_mn_prio_per_pair,
  "10_top_prioritized_interactions_per_lineage_pair.png",
  width = 13,
  height = 10
)

# ============================================================
# 21. top prioritized interactions faceted by direction
# ============================================================

top_prioritized_tbl <- mn_prio_plot %>%
  dplyr::filter(direction %in% c("PLA-up", "platelet-free-up")) %>%
  dplyr::group_by(direction) %>%
  dplyr::slice_max(order_by = prioritization_score, n = 20, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    label = paste0(interaction, "\n", lineage_pair),
    label = fct_reorder(label, prioritization_score)
  )

p_top_prioritized <- ggplot(
  top_prioritized_tbl,
  aes(
    x = prioritization_score,
    y = label,
    color = direction
  )
) +
  geom_point(size = 3, alpha = 0.85) +
  facet_wrap(~ direction, scales = "free_y") +
  theme_bw(base_size = 11) +
  scale_color_manual(
    values = c(
      "PLA-up" = "#F8766D",
      "platelet-free-up" = "#00BFC4"
    )
  ) +
  labs(
    x = "MultiNicheNet prioritization score",
    y = "Ligand → receptor / sender → receiver",
    color = "Direction",
    title = plot_title("Top prioritized MultiNicheNet interactions")
  )

save_plot(
  p_top_prioritized,
  "11_top_prioritized_interactions_faceted.png",
  width = 12,
  height = 10
)

# ============================================================
# 22. prioritization components
# ============================================================

mn_prio_top <- mn_prio_plot %>%
  dplyr::slice_max(prioritization_score, n = 2000, with_ties = FALSE)

p_prioritization_components <- ggplot(
  mn_prio_top,
  aes(
    x = lfc_score,
    y = activity_score,
    size = prioritization_score,
    color = direction
  )
) +
  geom_point(alpha = 0.45) +
  theme_bw() +
  scale_color_manual(
    values = c(
      "PLA-up" = "#F8766D",
      "platelet-free-up" = "#00BFC4"
    )
  ) +
  labs(
    title = plot_title("MultiNicheNet prioritization components"),
    x = "Mean scaled ligand/receptor logFC",
    y = "Max scaled ligand activity",
    size = "Prioritization score",
    color = "Direction"
  )

save_plot(
  p_prioritization_components,
  "12_prioritization_components.png",
  width = 8,
  height = 6
)

# ============================================================
# 23. prioritization score distribution
# ============================================================

p_score_distribution <- ggplot(
  prioritization_tables$group_prioritization_tbl,
  aes(x = group, y = prioritization_score, fill = group)
) +
  geom_boxplot(outlier.size = 0.5) +
  theme_bw() +
  labs(
    title = plot_title("Prioritization score distribution by group"),
    x = "Prioritized group",
    y = "Prioritization score"
  )

save_plot(
  p_score_distribution,
  "13_prioritization_score_distribution_by_group.png",
  width = 7,
  height = 5
)

# ============================================================
# 24. finish
# ============================================================

message("Finished MultiNicheNet analysis for: ", dataset_name)
message("Saved plots to: ", plot_dir)
message("Saved tables to: ", table_dir)