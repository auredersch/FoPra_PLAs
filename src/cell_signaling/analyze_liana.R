#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tibble)
  library(purrr)
  library(liana)
  library(ggplot2)
  library(circlize)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript analyze_liana.R <input_liana_rds> <output_dir>")
}

input_file <- args[[1]]
output_dir <- args[[2]]

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dataset_name <- gsub("_liana_results\\.rds$", "", basename(input_file))

message("Analyzing: ", dataset_name)
message("Input: ", input_file)

liana_result <- readRDS(input_file)

analyze_liana <- function(
    liana_result,
    top_n = 50,
    aggregate_rank_threshold = 0.01,
    genes_of_interest = NULL,
    source_groups = NULL,
    target_groups = NULL
) {
  consensus <- liana_result$consensus %>%
    arrange(aggregate_rank)

  if (!is.null(source_groups)) {
    consensus <- consensus %>% filter(source %in% source_groups)
  }

  if (!is.null(target_groups)) {
    consensus <- consensus %>% filter(target %in% target_groups)
  }

  significant_interactions <- consensus %>%
    filter(!is.na(aggregate_rank), aggregate_rank <= aggregate_rank_threshold)

  source_target_summary <- significant_interactions %>%
    group_by(source, target) %>%
    summarise(
      n_interactions = n(),
      n_unique_lr_pairs = n_distinct(ligand.complex, receptor.complex),
      best_aggregate_rank = min(aggregate_rank, na.rm = TRUE),
      median_aggregate_rank = median(aggregate_rank, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(n_interactions), best_aggregate_rank)

  ligand_summary <- significant_interactions %>%
    group_by(ligand.complex) %>%
    summarise(
      n_interactions = n(),
      n_source_types = n_distinct(source),
      n_target_types = n_distinct(target),
      best_aggregate_rank = min(aggregate_rank, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(n_interactions), best_aggregate_rank)

  receptor_summary <- significant_interactions %>%
    group_by(receptor.complex) %>%
    summarise(
      n_interactions = n(),
      n_source_types = n_distinct(source),
      n_target_types = n_distinct(target),
      best_aggregate_rank = min(aggregate_rank, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(n_interactions), best_aggregate_rank)

  lr_pair_summary <- significant_interactions %>%
    group_by(ligand.complex, receptor.complex) %>%
    summarise(
      n_interactions = n(),
      n_source_target_pairs = n_distinct(paste(source, target, sep = " -> ")),
      best_aggregate_rank = min(aggregate_rank, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(n_source_target_pairs), best_aggregate_rank)

  gene_interactions <- NULL
  gene_significant <- NULL

  if (!is.null(genes_of_interest)) {
    gene_pattern <- paste0(
      "(^|[^[:alnum:]])(",
      paste(genes_of_interest, collapse = "|"),
      ")([^[:alnum:]]|$)"
    )

    gene_interactions <- consensus %>%
      filter(
        str_detect(coalesce(as.character(ligand.complex), ""), regex(gene_pattern)) |
          str_detect(coalesce(as.character(receptor.complex), ""), regex(gene_pattern))
      )

    gene_significant <- gene_interactions %>%
      filter(!is.na(aggregate_rank), aggregate_rank <= aggregate_rank_threshold)
  }

  list(
    parameters = list(
      top_n = top_n,
      aggregate_rank_threshold = aggregate_rank_threshold,
      genes_of_interest = genes_of_interest,
      source_groups = source_groups,
      target_groups = target_groups
    ),
    filtered_consensus = consensus,
    top_interactions = slice_head(consensus, n = top_n),
    significant_interactions = significant_interactions,
    source_target_summary = source_target_summary,
    ligand_summary = ligand_summary,
    receptor_summary = receptor_summary,
    lr_pair_summary = lr_pair_summary,
    gene_interactions = gene_interactions,
    gene_significant = gene_significant
  )
}

pla_signaling_genes <- list(
  adhesion = c(
    "SELP", "SELPLG", "ITGAM", "ITGB2", "GP1BA", "GP1BB", "GP9", "GP5",
    "VWF", "ITGA2B", "ITGB3", "ICAM2", "ITGAL", "CD40LG", "CD40",
    "PDPN", "CLEC1B", "FCGR2A"
  ),
  cytokines_chemokines = c(
    "CXCL8", "CXCR1", "CXCR2", "CCL2", "CCR2", "PF4", "CCL5",
    "CCR1", "CCR3", "CCR5", "IL1B", "IL1R1", "IL1RAP"
  ),
  lipid_mediator_pathways = c(
    "PTGS2", "TBXAS1", "TBXA2R", "PTAFR"
  )
)

pla_liana_genes <- unique(c(
  pla_signaling_genes$adhesion,
  pla_signaling_genes$cytokines_chemokines
))

overview <- analyze_liana(
  liana_result = liana_result,
  top_n = 50,
  aggregate_rank_threshold = 0.01
)

pla_overview <- analyze_liana(
  liana_result = liana_result,
  top_n = 50,
  aggregate_rank_threshold = 0.01,
  genes_of_interest = pla_liana_genes
)

plot_dir <- file.path(output_dir, "plots", dataset_name)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(plot_expr, filename, width = 9, height = 7) {
  tryCatch({
    p <- plot_expr

    ggsave(
      filename = filename,
      plot = p,
      width = width,
      height = height,
      dpi = 300
    )

    message("Saved plot: ", filename)
  }, error = function(e) {
    message("Could not save plot: ", filename)
    message("Reason: ", conditionMessage(e))
  })
}

p <- liana_result$consensus %>%
  arrange(aggregate_rank) %>%
  liana_dotplot(
    ntop = 30,
    show_complex = TRUE
  ) +
  theme(
    axis.text.x = element_text(size = 5, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 6),
    strip.text = element_text(size = 7),
    axis.title = element_text(size = 9),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    plot.title = element_text(size = 10)
  )

save_plot(
  p,
  file.path(plot_dir, paste0(dataset_name, "_dotplot_top30.pdf")),
  width = 14,
  height = 10
)

save_heat_safely <- function(data, filename, width = 10, height = 8) {
  tryCatch({
    pdf(filename, width = width, height = height)

    p <- heat_freq(data)
    print(p)

    dev.off()
    message("Saved heatmap: ", filename)

  }, error = function(e) {
    if (dev.cur() != 1) dev.off()
    message("Could not save heatmap: ", filename)
    message("Reason: ", conditionMessage(e))
  })
}

save_heat_safely(
  overview$significant_interactions,
  file.path(plot_dir, paste0(dataset_name, "_heat_freq.pdf"))
)

save_chord_safely <- function(data, filename, width = 15, height = 15) {
  tryCatch({
    pdf(filename, width = width, height = height)

    par(
      mar = c(3, 3, 3, 3),
      xpd = NA
    )

    circlize::circos.clear()
    circlize::circos.par(
      canvas.xlim = c(-1.35, 1.35),
      canvas.ylim = c(-1.35, 1.35),
      gap.degree = 4
    )

    p <- chord_freq(data)
    print(p)

    circlize::circos.clear()
    dev.off()

    message("Saved chord plot: ", filename)

  }, error = function(e) {
    if (dev.cur() != 1) dev.off()
    circlize::circos.clear()

    message("Could not save chord plot: ", filename)
    message("Reason: ", conditionMessage(e))
  })
}

save_chord_safely(
  overview$significant_interactions,
  file.path(plot_dir, paste0(dataset_name, "_chord_freq.pdf"))
)

p_pla <- pla_overview$gene_interactions %>%
  arrange(aggregate_rank) %>%
  liana_dotplot(
    ntop = 30,
    show_complex = TRUE
  ) +
  theme(
    axis.text.x = element_text(size = 5, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 6),
    strip.text = element_text(size = 7),
    axis.title = element_text(size = 9),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8)
  )

save_plot(
  p_pla,
  file.path(plot_dir, paste0(dataset_name, "_pla_dotplot_top30.pdf")),
  width = 14,
  height = 10
)

saveRDS(
  overview,
  file.path(output_dir, paste0(dataset_name, "_overview.rds"))
)

saveRDS(
  pla_overview,
  file.path(output_dir, paste0(dataset_name, "_pla_overview.rds"))
)

message("Saved overview files for: ", dataset_name)