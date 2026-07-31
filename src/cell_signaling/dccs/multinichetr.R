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
  library(tidyr)
  library(tibble)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(forcats)
  library(circlize)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(nichenetr)
  library(multinichenetr)
  library(RColorBrewer)
  library(magrittr)
})

# ============================================================
# 1. arguments and output folders
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

dataset_clean <- dataset_name %>%
  stringr::str_remove("_(all|diseasedOnly|healthyOnly)$")

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
message("Table directory: ", table_dir)

save_plot <- function(plot, filename, width = 8, height = 6, dpi = 300) {
  out_file <- file.path(plot_dir, filename)

  ggsave(
    filename = out_file,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )

  message("Saved plot: ", out_file)
}

save_table <- function(x, filename) {
  out_file <- file.path(table_dir, filename)
  write.csv(x, out_file, row.names = FALSE)
  message("Saved table: ", out_file)
}

# ============================================================
# 2. load seurat object and metadata settings
# ============================================================

seurat_obj <- readRDS(input_file)
seurat_obj <- subset(
  seurat_obj,
  subset = !is.na(lineage) & !is.na(pla_status) & !is.na(sample)
)
DefaultAssay(seurat_obj)

sample_col <- "sample"
donor_col <- "donor_id"
celltype_col <- "celltype"
celltype_full <- "celltype_full"
condition_col <- "pla_status"
lineage_col <- "lineage"

required_conditions <- c("PLA", "platelet-free")
available_conditions <- unique(as.character(seurat_obj$pla_status))

missing_conditions <- setdiff(
  required_conditions,
  available_conditions
)

message(
  "Available pla_status values: ",
  paste(sort(available_conditions), collapse = ", ")
)

if (length(missing_conditions) > 0) {
  stop(
    "Cannot perform PLA versus platelet-free comparison. ",
    "Missing condition(s): ",
    paste(missing_conditions, collapse = ", ")
  )
}

seurat_obj$sample_condition <- paste(
  seurat_obj[[sample_col]][, 1],
  seurat_obj[[condition_col]][, 1],
  sep = "_"
)

metadata_summary <- seurat_obj@meta.data %>%
  dplyr::count(
    sample = .data[[sample_col]],
    condition = .data[[condition_col]],
    lineage = .data[[lineage_col]],
    name = "n_cells"
  )

save_table(metadata_summary, "metadata_sample_condition_lineage_counts.csv")

sample_condition_check <- seurat_obj@meta.data %>%
  dplyr::distinct(
    sample = .data[[sample_col]],
    condition = .data[[condition_col]]
  ) %>%
  dplyr::count(sample, name = "n_conditions") %>%
  dplyr::filter(n_conditions > 1)

save_table(sample_condition_check, "sample_condition_check.csv")

# ============================================================
# 3. analysis settings
# ============================================================

organism <- "human"

sample_id <- "sample_condition"
group_id <- "pla_status"
celltype_id <- "lineage"

covariates <- "sample"
batches <- NA

min_cells <- 10
min_sample_prop <- 0.50
fraction_cutoff <- 0.05

logFC_threshold <- 0.50
p_val_threshold <- 0.05
p_val_adj <- TRUE

top_n_target <- 250
ligand_activity_down <- FALSE
analyse_condition_specific_celltypes <- TRUE

options(timeout = 120)

# ============================================================
# 4. load ligand-receptor and ligand-target resources
# ============================================================

lr_network_all <- readRDS(url(
  "https://zenodo.org/record/10229222/files/lr_network_human_allInfo_30112033.rds"
)) %>%
  dplyr::mutate(
    ligand = convert_alias_to_symbols(ligand, organism = organism),
    receptor = convert_alias_to_symbols(receptor, organism = organism),
    ligand = make.names(ligand),
    receptor = make.names(receptor)
  )

lr_network <- lr_network_all %>%
  dplyr::distinct(ligand, receptor)

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
  dplyr::filter(ligand %in% colnames(ligand_target_matrix))

ligand_target_matrix <- ligand_target_matrix[, unique(lr_network$ligand)]

save_table(lr_network, "lr_network_filtered.csv")

# ============================================================
# 5. convert to SingleCellExperiment and clean identifiers
# ============================================================

sce <- Seurat::as.SingleCellExperiment(seurat_obj, assay = "RNA")
sce <- alias_to_symbol_SCE(sce, organism) %>%
  makenames_SCE()

colData(sce)[, sample_id] <- colData(sce)[, sample_id] %>%
  as.character() %>%
  make.names() %>%
  factor()

colData(sce)[, group_id] <- colData(sce)[, group_id] %>%
  as.character() %>%
  make.names() %>%
  factor()

colData(sce)[, celltype_id] <- colData(sce)[, celltype_id] %>%
  as.character() %>%
  make.names() %>%
  factor()

senders_oi <- colData(sce)[, celltype_id] %>%
  unique() %>%
  as.character()

receivers_oi <- colData(sce)[, celltype_id] %>%
  unique() %>%
  as.character()

sce <- sce[, colData(sce)[, celltype_id] %in% union(senders_oi, receivers_oi)]

# ============================================================
# 6. abundance filtering
# ============================================================

abundance_info <- get_abundance_info(
  sce = sce,
  sample_id = sample_id,
  group_id = group_id,
  celltype_id = celltype_id,
  min_cells = min_cells,
  senders_oi = senders_oi,
  receivers_oi = receivers_oi,
  batches = batches
)

save_table(
  abundance_info$abundance_data,
  "abundance_data_raw.csv"
)

if (!is.null(abundance_info$abund_plot_sample)) {
  save_plot(
    abundance_info$abund_plot_sample,
    "01_abundance_per_sample.png",
    width = 10,
    height = 7
  )
}

sample_group_celltype_df <- abundance_info$abundance_data %>%
  dplyr::filter(n > min_cells) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(sample_id, group_id) %>%
  cross_join(
    abundance_info$abundance_data %>%
      dplyr::ungroup() %>%
      dplyr::distinct(celltype_id)
  ) %>%
  dplyr::arrange(sample_id)

abundance_df <- sample_group_celltype_df %>%
  dplyr::left_join(
    abundance_info$abundance_data %>% dplyr::ungroup(),
    by = c("sample_id", "group_id", "celltype_id")
  )

abundance_df$n[is.na(abundance_df$n)] <- 0
abundance_df$keep[is.na(abundance_df$keep)] <- FALSE

abundance_df_summarized <- abundance_df %>%
  dplyr::mutate(keep = as.logical(keep)) %>%
  dplyr::group_by(group_id, celltype_id) %>%
  dplyr::summarise(
    samples_present = sum(keep),
    .groups = "drop"
  )

celltypes_absent_one_condition <- abundance_df_summarized %>%
  dplyr::filter(samples_present == 0) %>%
  dplyr::pull(celltype_id) %>%
  unique()

celltypes_present_one_condition <- abundance_df_summarized %>%
  dplyr::filter(samples_present >= 2) %>%
  dplyr::pull(celltype_id) %>%
  unique()

condition_specific_celltypes <- intersect(
  celltypes_absent_one_condition,
  celltypes_present_one_condition
)

total_nr_conditions <- colData(sce)[, group_id] %>%
  unique() %>%
  length()

absent_celltypes <- abundance_df_summarized %>%
  dplyr::filter(samples_present < 2) %>%
  dplyr::group_by(celltype_id) %>%
  dplyr::count() %>%
  dplyr::filter(n == total_nr_conditions) %>%
  dplyr::pull(celltype_id)

save_table(abundance_df, "abundance_df_completed.csv")
save_table(abundance_df_summarized, "abundance_df_summarized.csv")

condition_specific_celltypes <- as.character(condition_specific_celltypes)
absent_celltypes <- as.character(absent_celltypes)

writeLines(
  condition_specific_celltypes,
  file.path(table_dir, "condition_specific_celltypes.txt")
)

writeLines(
  absent_celltypes,
  file.path(table_dir, "absent_celltypes.txt")
)

message("Condition-specific celltypes:")
message(paste(condition_specific_celltypes, collapse = ", "))

message("Absent celltypes:")
message(paste(absent_celltypes, collapse = ", "))

if (analyse_condition_specific_celltypes) {
  senders_oi <- setdiff(senders_oi, absent_celltypes)
  receivers_oi <- setdiff(receivers_oi, absent_celltypes)
} else {
  excluded_celltypes <- union(absent_celltypes, condition_specific_celltypes)
  senders_oi <- setdiff(senders_oi, excluded_celltypes)
  receivers_oi <- setdiff(receivers_oi, excluded_celltypes)
}

sce <- sce[, colData(sce)[, celltype_id] %in% union(senders_oi, receivers_oi)]

# ============================================================
# 7. gene filtering
# ============================================================

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
  dplyr::filter(expressed == TRUE) %>%
  dplyr::pull(gene) %>%
  unique()

sce <- sce[genes_oi, ]

save_table(frq_list$expressed_df, "expressed_genes_by_celltype.csv")

# ============================================================
# 8. pseudobulk expression processing
# ============================================================

abundance_expression_info <- process_abundance_expression_info(
  sce = sce,
  sample_id = sample_id,
  group_id = group_id,
  celltype_id = celltype_id,
  min_cells = min_cells,
  senders_oi = senders_oi,
  receivers_oi = receivers_oi,
  lr_network = lr_network,
  batches = batches,
  frq_list = frq_list,
  abundance_info = abundance_info
)

saveRDS(
  abundance_expression_info,
  file.path(out_dir, "abundance_expression_info.rds")
)

# ============================================================
# 9. differential expression
# ============================================================

contrasts_oi <- c("'PLA-platelet.free','platelet.free-PLA'")

contrast_tbl <- tibble(
  contrast = c("PLA-platelet.free", "platelet.free-PLA"),
  group = c("PLA", "platelet.free")
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

celltype_de <- DE_info$celltype_de$de_output_tidy
sender_receiver_de <- DE_info$sender_receiver_de

if (is.null(sender_receiver_de)) {
  stop("sender_receiver_de is NULL. Check names(DE_info) and the output of get_DE_info().")
}

saveRDS(DE_info, file.path(out_dir, "DE_info.rds"))
save_table(celltype_de, "celltype_de.csv")
save_table(sender_receiver_de, "sender_receiver_de.csv")

abundance_kept_summary <- abundance_info$abundance_data %>%
  dplyr::ungroup() %>%
  dplyr::filter(keep == TRUE) %>%
  dplyr::count(group_id, celltype_id)

save_table(abundance_kept_summary, "abundance_kept_summary.csv")

if (!is.null(DE_info$hist_pvals)) {
  save_plot(
    DE_info$hist_pvals,
    "02_deseq_pvalue_histograms.png",
    width = 10,
    height = 7
  )
}

# ============================================================
# 10. geneset and ligand activity analysis
# ============================================================

geneset_assessment <- contrast_tbl$contrast %>%
  lapply(
    process_geneset_data,
    celltype_de,
    logFC_threshold,
    p_val_adj,
    p_val_threshold
  ) %>%
  bind_rows()

save_table(geneset_assessment, "geneset_assessment.csv")

ligand_activities_targets_DEgenes <- suppressMessages(suppressWarnings(
  get_ligand_activities_targets_DEgenes(
    receiver_de = celltype_de,
    receivers_oi = intersect(receivers_oi, unique(celltype_de$cluster_id)),
    ligand_target_matrix = ligand_target_matrix,
    logFC_threshold = logFC_threshold,
    p_val_threshold = p_val_threshold,
    p_val_adj = p_val_adj,
    top_n_target = top_n_target
  )
))

save_table(
  ligand_activities_targets_DEgenes,
  "ligand_activities_targets_DEgenes.csv"
)

# ============================================================
# 11. prioritization
# ============================================================

sender_receiver_tbl <- sender_receiver_de %>%
  dplyr::distinct(sender, receiver)

metadata_combined <- colData(sce) %>%
  dplyr::as_tibble()

if (!is.na(batches)) {
  grouping_tbl <- metadata_combined[, c(sample_id, group_id, batches)] %>%
    dplyr::as_tibble() %>%
    dplyr::distinct()

  colnames(grouping_tbl) <- c("sample", "group", batches)
} else {
  grouping_tbl <- metadata_combined[, c(sample_id, group_id)] %>%
    dplyr::as_tibble() %>%
    dplyr::distinct()

  colnames(grouping_tbl) <- c("sample", "group")
}

prioritization_tables <- suppressMessages(generate_prioritization_tables(
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de = sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  contrast_tbl = contrast_tbl,
  sender_receiver_tbl = sender_receiver_tbl,
  grouping_tbl = grouping_tbl,
  scenario = "regular",
  fraction_cutoff = fraction_cutoff,
  abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
  abundance_data_sender = abundance_expression_info$abundance_data_sender,
  ligand_activity_down = ligand_activity_down
))

lr_target_prior_cor = lr_target_prior_cor_inference(
  receivers_oi = prioritization_tables$group_prioritization_tbl$receiver %>% unique(), 
  abundance_expression_info = abundance_expression_info, 
  celltype_de = celltype_de, 
  grouping_tbl = grouping_tbl, 
  prioritization_tables = prioritization_tables, 
  ligand_target_matrix = ligand_target_matrix, 
  logFC_threshold = logFC_threshold, 
  p_val_threshold = p_val_threshold, 
  p_val_adj = p_val_adj
  )

multinichenet_output <- list(
  celltype_info = abundance_expression_info$celltype_info,
  celltype_de = celltype_de,
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de = sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  prioritization_tables = prioritization_tables,
  grouping_tbl = grouping_tbl,
  lr_target_prior_cor = lr_target_prior_cor
)

multinichenet_output <- make_lite_output(multinichenet_output)

saveRDS(
  multinichenet_output,
  file.path(out_dir, "multinichenet_output.rds")
)

save_table(grouping_tbl, "grouping_tbl.csv")
save_table(
  prioritization_tables$group_prioritization_tbl,
  "group_prioritization_tbl.csv"
)

# ============================================================
# 12. chord plot
# ============================================================

prioritized_tbl_oi_all <- get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables,
  top_n = 50,
  rank_per_group = FALSE
)

save_table(
  prioritized_tbl_oi_all,
  "top_50_lr_pairs_all_groups.csv"
)

prioritized_tbl_oi <- multinichenet_output$prioritization_tables$group_prioritization_tbl %>%
  dplyr::filter(id %in% prioritized_tbl_oi_all$id) %>%
  dplyr::distinct(id, sender, receiver, ligand, receptor, group) %>%
  dplyr::left_join(
    prioritized_tbl_oi_all,
    by = c("id", "sender", "receiver", "ligand", "receptor", "group")
  )

prioritized_tbl_oi$prioritization_score[is.na(prioritized_tbl_oi$prioritization_score)] <- 0

save_table(
  prioritized_tbl_oi,
  "prioritized_tbl_oi_top50.csv"
)

senders_receivers <- union(
  unique(prioritized_tbl_oi$sender),
  unique(prioritized_tbl_oi$receiver)
) %>%
  sort()

base_palette <- RColorBrewer::brewer.pal(
  n = 11,
  name = "Spectral"
)

colors_sender <- colorRampPalette(base_palette)(length(senders_receivers)) %>%
  magrittr::set_names(senders_receivers)

colors_receiver <- colors_sender

png(
  filename = file.path(plot_dir, "03_circos_group_comparison.png"),
  width = 3000,
  height = 3000,
  res = 300
)

circos_list <- make_circos_group_comparison(
  prioritized_tbl_oi,
  colors_sender,
  colors_receiver
)

dev.off()

saveRDS(
  circos_list,
  file.path(out_dir, "circos_group_comparison.rds")
)

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

save_table(
  mn_prio_plot,
  "mn_prio_plot_table.csv"
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

save_table(
  pair_counts_mn_prio,
  "pair_counts_mn_prio.csv"
)

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
# 17. net direction heatmap
# ============================================================

net_direction_mn_prio <- mn_prio_plot %>%
  dplyr::count(sender, receiver, direction) %>%
  tidyr::complete(
    sender,
    receiver,
    direction = c("PLA-up", "platelet-free-up"),
    fill = list(n = 0)
  ) %>%
  tidyr::pivot_wider(
    names_from = direction,
    values_from = n,
    values_fill = 0
  ) %>%
  dplyr::mutate(
    net = `PLA-up` - `platelet-free-up`
  )

save_table(
  net_direction_mn_prio,
  "net_direction_mn_prio.csv"
)

p_net_direction_mn_prio <- ggplot(
  net_direction_mn_prio,
  aes(x = receiver, y = sender, fill = net)
) +
  geom_tile(color = "white") +
  geom_text(aes(label = net), size = 4) +
  theme_bw() +
  labs(
    title = plot_title("Net direction of prioritized MultiNicheNet interactions"),
    x = "Receiver lineage",
    y = "Sender lineage",
    fill = "PLA-up minus\nplatelet-free-up"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(
  p_net_direction_mn_prio,
  "07_net_direction_mn_prio_heatmap.png",
  width = 8,
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

save_table(
  net_score_mn_prio,
  "net_score_mn_prio.csv"
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
# 19. top LR pairs by mean score
# ============================================================

top_lr_by_score <- mn_prio_plot %>%
  dplyr::group_by(interaction, direction) %>%
  dplyr::summarise(
    n_lineage_pairs = dplyr::n_distinct(lineage_pair),
    mean_score = mean(prioritization_score, na.rm = TRUE),
    max_score = max(prioritization_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(mean_score)) %>%
  dplyr::slice_head(n = 30) %>%
  dplyr::mutate(
    interaction = fct_reorder(interaction, mean_score)
  )

save_table(
  top_lr_by_score,
  "top_lr_by_score.csv"
)

p_top_lr_by_score <- ggplot(
  top_lr_by_score,
  aes(x = mean_score, y = interaction, fill = direction)
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
    title = plot_title("Top prioritized MultiNicheNet LR pairs by mean score"),
    x = "Mean prioritization score",
    y = "Ligand → receptor",
    fill = "Direction"
  )

save_plot(
  p_top_lr_by_score,
  "09_top_lr_by_mean_score.png",
  width = 8,
  height = 8
)

# ============================================================
# 20. top interactions per lineage pair
# ============================================================

top_mn_prio_per_pair <- mn_prio_plot %>%
  dplyr::group_by(lineage_pair) %>%
  dplyr::slice_max(prioritization_score, n = 5, with_ties = FALSE) %>%
  dplyr::ungroup()

save_table(
  top_mn_prio_per_pair,
  "top_mn_prio_per_pair.csv"
)

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

save_table(
  top_prioritized_tbl,
  "top_prioritized_interactions_by_direction.csv"
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

save_table(
  mn_prio_top,
  "mn_prio_top_2000.csv"
)

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