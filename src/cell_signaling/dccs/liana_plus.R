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
  library(liana)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(forcats)
  library(circlize)
  library(DESeq2)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript liana_plus.R <input_rds> <output_dir>")
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

# pseudobulking prep
seurat_obj$sample_clean <- gsub("_", "-", seurat_obj$sample)

seurat_obj$pla_status[is.na(seurat_obj$pla_status)] <- "platelet-free"

seurat_obj$pb_group <- paste(
  seurat_obj$sample_clean,
  seurat_obj$lineage,
  seurat_obj$pla_status,
  sep = "|"
)

# pseudobulking
pseudobulk <- AggregateExpression(
    object = seurat_obj,
    group.by = "pb_group",
    assays = "RNA",
    slot = "counts",
    verbose = TRUE
)

counts <- pseudobulk$RNA

# deseq prep
pb_meta <- data.frame(
  pb_group = colnames(counts),
  stringsAsFactors = FALSE
)

pb_meta <- tidyr::separate(
  data = pb_meta,
  col = "pb_group",
  into = c("sample_id", "lineage", "condition"),
  sep = "\\|",
  remove = FALSE
)

rownames(pb_meta) <- pb_meta$pb_group

pb_meta$condition <- factor(
  pb_meta$condition,
  levels = c("platelet-free", "PLA")
)

pb_meta$sample_id <- factor(pb_meta$sample_id)
pb_meta$lineage <- factor(pb_meta$lineage)

stopifnot(all(rownames(pb_meta) == colnames(counts)))

table(pb_meta$lineage, pb_meta$condition)

# deseq
dea_results <- list()

for (lin in unique(pb_meta$lineage)) {
  
  message("Running DESeq2 for lineage: ", lin)
  
  meta_lin <- pb_meta %>%
    dplyr::filter(lineage == lin) %>%
    as.data.frame()
  
  if (nrow(meta_lin) == 0) {
    message("Skipping ", lin, ": no pseudobulk samples")
    next
  }
  
  if (length(unique(meta_lin$condition)) < 2) {
    message("Skipping ", lin, ": only one condition present")
    next
  }

  
  keep_genes <- rowSums(counts_lin) >= 10
  counts_lin <- counts_lin[keep_genes, , drop = FALSE]

  if (nrow(counts_lin) == 0) {
    message("Skipping ", lin, ": no genes left after count filtering")
    next
  }
  
  meta_lin$condition <- factor(
    meta_lin$condition,
    levels = c("platelet-free", "PLA")
  )
  
  meta_lin$sample_id <- factor(meta_lin$sample_id)
  
  stopifnot(all(rownames(meta_lin) == colnames(counts_lin)))

  # DESeq2: PLA vs platelet-free inside this lineage
  # creates deseq2 object
  

  design_formula <- if (all(table(meta_lin$sample_id, meta_lin$condition) > 0)) {
    message("Using paired DESeq2 design for ", lin, ": ~ sample_id + condition")
    ~ sample_id + condition
  } else {
    message("Using unpaired DESeq2 design for ", lin, ": ~ condition")
    ~ condition
  }
  
  dds <- DESeqDataSetFromMatrix(
    countData = round(counts_lin),
    colData = meta_lin, # rows of colData correspond to cols of countData
    design = design_formula # model gene counts as a function of pla_status + sample_id (paired analysis)
    # for every gene deseq asks: does expression differ between PLA & platelet-free
    # first account for baseline differences between samples, then estimate the PLA vs platelet-free effect within those samples.

  )
  
  dds <- DESeq(dds)
  
  res <- results(
    dds,
    contrast = c("condition", "PLA", "platelet-free")
    # log2FoldChange = log2(PLA / platelet-free)
    # positive log2FC = higher in PLA
    # negative log2FC = higher in platelet-free
  )
  
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    mutate(lineage = lin)
  
  dea_results[[as.character(lin)]] <- res_df
}

dea_table <- dplyr::bind_rows(dea_results)

if (nrow(dea_table) == 0) {
  stop("No DESeq2 results generated. All lineages were skipped.")
}

dea_table %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::arrange(padj) %>%
  head()

# dea to ligand-receptor interactions
# 1. clean DE table (keep only genes that deseq actually tested)
dea_table_clean <- dea_table %>%
  select(lineage, gene, log2FoldChange, stat, pvalue, padj) %>%
  filter(!is.na(stat))

# 2. get known ligand-receptor pairs (consensus)
lr_resource_raw <- liana::select_resource("Consensus")

lr_resource <- lr_resource_raw$Consensus %>%
  dplyr::select(
    ligand = source_genesymbol,
    receptor = target_genesymbol
  ) %>%
  dplyr::distinct()

# 3. make all source-target lineage pairs
lineages <- unique(dea_table_clean$lineage)

cell_pairs <- expand.grid(
  source = lineages,
  target = lineages,
  stringsAsFactors = FALSE
)

# 4. prepare ligand DE stats
ligand_stats <- dea_table_clean %>%
  select(
    source = lineage,
    ligand = gene,
    ligand_log2FC = log2FoldChange,
    ligand_stat = stat,
    ligand_pvalue = pvalue,
    ligand_padj = padj
  )

# 5. prepare receptor DE stats
receptor_stats <- dea_table_clean %>%
  select(
    target = lineage,
    receptor = gene,
    receptor_log2FC = log2FoldChange,
    receptor_stat = stat,
    receptor_pvalue = pvalue,
    receptor_padj = padj
  )

# 6. combine source-target pairs with ligand-receptor resource
# contains all possible combinations now
lr_res <- cell_pairs %>%
  crossing(lr_resource)

# 7. join ligand stats from the source lineage
# this 
lr_res <- lr_res %>%
  left_join(ligand_stats, by = c("source", "ligand"))

# 8. join receptor stats from the target lineage
lr_res <- lr_res %>%
  left_join(receptor_stats, by = c("target", "receptor"))

# 9. keep only pairs where both ligand and receptor have DE stats
# keep only interactions where ligand and receptor were both tested by DESeq2

# 10. make interaction scores
  # try stouffers method to combine 2 wald statistics

lr_res <- lr_res %>%
  filter(!is.na(ligand_stat), !is.na(receptor_stat)) %>%
  mutate(
    stouffer_interaction_z = (ligand_stat + receptor_stat) / sqrt(2),
    stouffer_interaction_pvalue = 2 * pnorm(-abs(stouffer_interaction_z)),
    stouffer_interaction_padj = p.adjust(stouffer_interaction_pvalue, method = "BH"),
    
    interaction_mean = (ligand_stat + receptor_stat) / 2,
    interaction_abs_mean = (abs(ligand_stat) + abs(receptor_stat)) / 2,
    
    interaction_max_padj = pmax(ligand_padj, receptor_padj, na.rm = TRUE),
    
    ligand_sig = !is.na(ligand_padj) & ligand_padj < 0.05,
    receptor_sig = !is.na(receptor_padj) & receptor_padj < 0.05,
    both_sig = ligand_sig & receptor_sig,
    
    interaction = paste(ligand, receptor, sep = " → "),
    lineage_pair = paste(source, target, sep = " → "),
    
    direction = case_when(
      ligand_stat > 0 & receptor_stat > 0 ~ "PLA-up",
      ligand_stat < 0 & receptor_stat < 0 ~ "platelet-free-up",
      TRUE ~ "discordant"
    ),
    
    strict_direction = case_when(
      both_sig & ligand_log2FC > 0 & receptor_log2FC > 0 ~ "PLA-up",
      both_sig & ligand_log2FC < 0 & receptor_log2FC < 0 ~ "platelet-free-up",
      both_sig ~ "discordant",
      TRUE ~ "not both significant"
    )
  ) %>%
  arrange(stouffer_interaction_padj)

# plots
# save tables
write.csv(dea_table, file.path(table_dir, "dea_table.csv"), row.names = FALSE)
write.csv(lr_res, file.path(table_dir, "lr_res_differential_lr.csv"), row.names = FALSE)

lr_plot <- lr_res %>%
  mutate(
    mean_signed_score = interaction_mean,
    mean_abs_score = interaction_abs_mean,
    stouffer_signed_score = stouffer_interaction_z,
    stouffer_abs_score = abs(stouffer_interaction_z)
  )

lr_plot_stouffer_sig <- lr_plot %>%
  filter(stouffer_interaction_padj < 0.05)

lr_plot_strict <- lr_plot %>%
  filter(strict_direction %in% c("PLA-up", "platelet-free-up", "discordant")) %>%
  mutate(direction = strict_direction)

write.csv(lr_plot, file.path(table_dir, "lr_plot_all.csv"), row.names = FALSE)
write.csv(lr_plot_stouffer_sig, file.path(table_dir, "lr_plot_stouffer_sig.csv"), row.names = FALSE)
write.csv(lr_plot_strict, file.path(table_dir, "lr_plot_strict.csv"), row.names = FALSE)


# 1. Global strict distribution
plot_df <- lr_plot_strict

p_global_strict <- ggplot(plot_df, aes(x = direction, fill = direction)) +
  geom_bar() +
  theme_bw() +
  labs(
    title = plot_title("Global distribution of differential ligand-receptor interactions"),
    x = "Interaction category",
    y = "# LR pairs"
  ) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

save_plot(
  p_global_strict,
  "01_global_distribution_strict.png",
  width = 7,
  height = 5
)


# 2. Global Stouffer-significant distribution
plot_df <- lr_plot_stouffer_sig

p_global_stouffer <- ggplot(plot_df, aes(x = direction, fill = direction)) +
  geom_bar() +
  theme_bw() +
  labs(
    title = plot_title("Global distribution of Stouffer-significant LR interactions"),
    x = "Interaction category",
    y = "# LR pairs"
  ) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

save_plot(
  p_global_stouffer,
  "02_global_distribution_stouffer_sig.png",
  width = 7,
  height = 5
)


# 3. Direction by lineage pair, strict
plot_df <- lr_plot_strict

pair_counts <- plot_df %>%
  dplyr::count(lineage_pair, direction) %>%
  dplyr::group_by(lineage_pair) %>%
  dplyr::mutate(total = sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(total)) %>%
  dplyr::mutate(lineage_pair = fct_reorder(lineage_pair, total))

write.csv(pair_counts, file.path(table_dir, "pair_counts_strict.csv"), row.names = FALSE)

p_pair_counts <- ggplot(pair_counts, aes(x = n, y = lineage_pair, fill = direction)) +
  geom_col() +
  theme_bw() +
  labs(
    title = plot_title("Direction of differential ligand-receptor interactions by lineage pair"),
    x = "# differential LR pairs",
    y = "Source → target lineage",
    fill = "Direction"
  )

save_plot(
  p_pair_counts,
  "03_direction_by_lineage_pair_strict.png",
  width = 9,
  height = 7
)


# 4. Direction by lineage pair, Stouffer
plot_df <- lr_plot_stouffer_sig

pair_counts_stouffer <- plot_df %>%
  dplyr::count(lineage_pair, direction) %>%
  dplyr::group_by(lineage_pair) %>%
  dplyr::mutate(total = sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(total)) %>%
  dplyr::mutate(lineage_pair = fct_reorder(lineage_pair, total))

write.csv(pair_counts_stouffer, file.path(table_dir, "pair_counts_stouffer_sig.csv"), row.names = FALSE)

p_pair_counts_stouffer <- ggplot(pair_counts_stouffer, aes(x = n, y = lineage_pair, fill = direction)) +
  geom_col() +
  theme_bw() +
  labs(
    title = plot_title("Direction of differential ligand-receptor interactions by lineage pair (Stouffer stats)"),
    x = "# differential LR pairs",
    y = "Source → target lineage",
    fill = "Direction"
  )

save_plot(
  p_pair_counts_stouffer,
  "04_direction_by_lineage_pair_stouffer_sig.png",
  width = 9,
  height = 7
)


# 5. Net direction heatmap, strict
plot_df <- lr_plot_strict

net_direction <- plot_df %>%
  dplyr::filter(direction != "discordant") %>%
  dplyr::count(source, target, direction) %>%
  complete(
    source,
    target,
    direction = c("PLA-up", "platelet-free-up"),
    fill = list(n = 0)
  ) %>%
  pivot_wider(
    names_from = direction,
    values_from = n,
    values_fill = 0
  ) %>%
  mutate(
    net = `PLA-up` - `platelet-free-up`
  )

write.csv(net_direction, file.path(table_dir, "net_direction_strict.csv"), row.names = FALSE)

p_net_direction <- ggplot(net_direction, aes(x = target, y = source, fill = net)) +
  geom_tile(color = "white") +
  geom_text(aes(label = net), size = 4) +
  theme_bw() +
  labs(
    title = plot_title("Net direction of differential ligand-receptor interactions"),
    x = "Target lineage",
    y = "Source lineage",
    fill = "PLA-up minus\nplatelet-free-up"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

save_plot(
  p_net_direction,
  "05_net_direction_heatmap_strict.png",
  width = 8,
  height = 7
)


# 6. Net direction heatmap, Stouffer
plot_df <- lr_plot_stouffer_sig

net_direction_stouffer <- plot_df %>%
  dplyr::filter(direction != "discordant") %>%
  dplyr::count(source, target, direction) %>%
  complete(
    source,
    target,
    direction = c("PLA-up", "platelet-free-up"),
    fill = list(n = 0)
  ) %>%
  pivot_wider(
    names_from = direction,
    values_from = n,
    values_fill = 0
  ) %>%
  mutate(
    net = `PLA-up` - `platelet-free-up`
  )

write.csv(net_direction_stouffer, file.path(table_dir, "net_direction_stouffer_sig.csv"), row.names = FALSE)

p_net_direction_stouffer <- ggplot(net_direction_stouffer, aes(x = target, y = source, fill = net)) +
  geom_tile(color = "white") +
  geom_text(aes(label = net), size = 4) +
  theme_bw() +
  labs(
    title = plot_title("Net direction of Stouffer-significant ligand-receptor interactions"),
    x = "Target lineage",
    y = "Source lineage",
    fill = "PLA-up minus\nplatelet-free-up"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

save_plot(
  p_net_direction_stouffer,
  "06_net_direction_heatmap_stouffer_sig.png",
  width = 8,
  height = 7
)


# 7. Most recurrent LR pairs, strict
plot_df <- lr_plot_strict

top_recurrent <- plot_df %>%
  dplyr::filter(direction != "discordant") %>%
  dplyr::count(interaction, direction, name = "n_lineage_pairs") %>%
  dplyr::group_by(interaction) %>%
  dplyr::mutate(total = sum(n_lineage_pairs)) %>%
  dplyr::ungroup() %>%
  dplyr::slice_max(total, n = 25) %>%
  dplyr::mutate(interaction = fct_reorder(interaction, total))

write.csv(top_recurrent, file.path(table_dir, "top_recurrent_strict.csv"), row.names = FALSE)

p_top_recurrent <- ggplot(top_recurrent, aes(x = n_lineage_pairs, y = interaction, fill = direction)) +
  geom_col() +
  theme_bw() +
  labs(
    title = plot_title("Most recurrent differential ligand-receptor pairs"),
    x = "# source-target lineage pairs",
    y = "Ligand → receptor",
    fill = "Direction"
  )

save_plot(
  p_top_recurrent,
  "07_top_recurrent_lr_pairs_strict.png",
  width = 8,
  height = 8
)


# 8. Most recurrent LR pairs, Stouffer
plot_df <- lr_plot_stouffer_sig

top_recurrent_stouffer <- plot_df %>%
  dplyr::filter(direction != "discordant") %>%
  dplyr::count(interaction, direction, name = "n_lineage_pairs") %>%
  dplyr::group_by(interaction) %>%
  dplyr::mutate(total = sum(n_lineage_pairs)) %>%
  dplyr::ungroup() %>%
  dplyr::slice_max(total, n = 25) %>%
  dplyr::mutate(interaction = fct_reorder(interaction, total))

write.csv(top_recurrent_stouffer, file.path(table_dir, "top_recurrent_stouffer_sig.csv"), row.names = FALSE)

p_top_recurrent_stouffer <- ggplot(top_recurrent_stouffer, aes(x = n_lineage_pairs, y = interaction, fill = direction)) +
  geom_col() +
  theme_bw() +
  labs(
    title = plot_title("Most recurrent Stouffer-significant ligand-receptor pairs"),
    x = "# source-target lineage pairs",
    y = "Ligand → receptor",
    fill = "Direction"
  )

save_plot(
  p_top_recurrent_stouffer,
  "08_top_recurrent_lr_pairs_stouffer_sig.png",
  width = 8,
  height = 8
)


# 9. Top LR per pair, mean Wald score
plot_df <- lr_plot_strict

top_lr_per_pair_mean <- plot_df %>%
  filter(direction != "discordant") %>%
  group_by(lineage_pair) %>%
  slice_max(mean_abs_score, n = 5, with_ties = FALSE) %>%
  ungroup()

write.csv(top_lr_per_pair_mean, file.path(table_dir, "top_lr_per_pair_mean_wald.csv"), row.names = FALSE)

p_top_lr_per_pair_mean <- ggplot(
  top_lr_per_pair_mean,
  aes(
    x = lineage_pair,
    y = interaction,
    size = mean_abs_score,
    color = mean_signed_score
  )
) +
  geom_point(alpha = 0.8) +
  theme_bw() +
  labs(
    title = plot_title("Top differential LR interactions per lineage pair using mean Wald score"),
    x = "Source → target lineage",
    y = "Ligand → receptor",
    size = "Mean |Wald stat|",
    color = "Mean Wald stat"
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
  p_top_lr_per_pair_mean,
  "09_top_lr_per_pair_mean_wald.png",
  width = 13,
  height = 10
)


# 10. Top LR per pair, Stouffer
plot_df <- lr_plot_stouffer_sig

top_lr_per_pair_stouffer <- plot_df %>%
  filter(direction != "discordant") %>%
  group_by(lineage_pair) %>%
  slice_max(stouffer_abs_score, n = 5, with_ties = FALSE) %>%
  ungroup()

write.csv(top_lr_per_pair_stouffer, file.path(table_dir, "top_lr_per_pair_stouffer.csv"), row.names = FALSE)

p_top_lr_per_pair_stouffer <- ggplot(
  top_lr_per_pair_stouffer,
  aes(
    x = lineage_pair,
    y = interaction,
    size = stouffer_abs_score,
    color = stouffer_signed_score
  )
) +
  geom_point(alpha = 0.8) +
  theme_bw() +
  labs(
    title = plot_title("Top differential LR interactions per lineage pair"),
    x = "Source → target lineage",
    y = "Ligand → receptor",
    size = "|Stouffer z|",
    color = "Stouffer z"
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
  p_top_lr_per_pair_stouffer,
  "10_top_lr_per_pair_stouffer.png",
  width = 13,
  height = 10
)


# 11. Mean Wald score vs Stouffer z-score
compare_scores <- lr_plot %>%
  filter(direction != "discordant") %>%
  mutate(
    interaction_label = paste(lineage_pair, interaction, sep = " | ")
  )

write.csv(compare_scores, file.path(table_dir, "compare_scores.csv"), row.names = FALSE)

p_compare_signed <- ggplot(
  compare_scores,
  aes(
    x = mean_signed_score,
    y = stouffer_signed_score,
    color = direction
  )
) +
  geom_point(alpha = 0.4, size = 1.5) +
  theme_bw() +
  labs(
    title = plot_title("Comparison of mean Wald score and Stouffer z-score"),
    x = "Mean Wald score",
    y = "Stouffer z-score",
    color = "Direction"
  ) +
  theme(
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 9),
    plot.title = element_text(size = 10),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7)
  )

save_plot(
  p_compare_signed,
  "11_compare_mean_wald_vs_stouffer_signed.png",
  width = 7,
  height = 6
)


# 12. Absolute Wald score vs absolute Stouffer score
p_compare_abs <- ggplot(
  lr_plot,
  aes(
    x = mean_abs_score,
    y = stouffer_abs_score,
    color = direction
  )
) +
  geom_point(alpha = 0.5, size = 1.5) +
  theme_bw() +
  labs(
    title = plot_title("Mean absolute Wald score vs absolute Stouffer score"),
    x = "Mean |Wald stat|",
    y = "|Stouffer z|",
    color = "Direction"
  )

save_plot(
  p_compare_abs,
  "12_compare_mean_abs_wald_vs_abs_stouffer.png",
  width = 7,
  height = 6
)

message("Finished analysis for: ", dataset_name)
message("Saved plots to: ", plot_dir)
message("Saved tables to: ", table_dir)