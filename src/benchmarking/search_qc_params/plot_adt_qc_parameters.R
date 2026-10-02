# From project root: Rscript src/benchmarking/search_qc_params/plot_adt_qc_parameters.R [grid_dir]
library(dplyr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args)) args[1] else "results/benchmarking/qc_grid_no_na_full"
output_dir <- file.path(input_dir, "plots")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
data <- read.csv(file.path(input_dir, "combined/grid_summary.csv")) %>%
  #filter(baseline | min_cells_per_status == 5) %>%
  arrange(desc(run_id)) %>%
  distinct(dataset, gene_set, mode, scope,
    min_difference, max_overlap, min_cells_per_status,
    .keep_all = TRUE) %>%
  mutate(
    method = paste0(gene_set, "\n", sub("gmm_dist_", "", mode), " | ", scope),
    # Display the unfiltered baseline alongside the min/status = 5 grid.
    min_status_panel = factor(ifelse(baseline, 20, min_cells_per_status)),
    across(c(min_difference, max_overlap, min_cells_per_status), factor)
  )

for (ds in unique(data$dataset)) {
  for (metric in c("F1", "retained_cell_fraction")) {
    points <- filter(data, dataset == ds) %>%
      mutate(label = if (metric == "F1") sprintf("%.3f", F1) else
        sprintf("%.0f%%", 100 * retained_cell_fraction),
        text_color = ifelse(.data[[metric]] < if (metric == "F1")
          mean(range(F1, na.rm = TRUE)) else 0.5, "white", "black"))
    p <- ggplot(points, aes(max_overlap, min_difference, fill = .data[[metric]])) +
      geom_tile(color = "white") +
      geom_text(aes(label = label, color = text_color), size = 2.8) +
      scale_color_identity() +
      facet_grid(method ~ min_status_panel,
        labeller = labeller(min_status_panel = function(x) paste0("min/status = ", x))) +
      scale_x_discrete(labels = function(x) ifelse(x == "Inf", "No ADT QC", x)) +
      scale_y_discrete(labels = function(x) ifelse(x == "-Inf", "No ADT QC", x)) +
      scale_fill_viridis_c(
        limits = if (metric == "F1") NULL else c(0, 1),
        labels = if (metric == "F1") scales::label_number(accuracy = 0.01) else scales::label_percent(),
        na.value = "grey90"
      ) +
      labs(title = paste(ds, "-", if (metric == "F1") "Cell-level F1" else "Cell retention"),
        x = "max_overlap", y = "min_difference", fill = NULL,
        caption = paste("Retention relative to the frequency-only baseline. Grey: undefined metric.",
          "No ADT QC: baseline without ADT filtering, independent of min/status.")) +
      theme_minimal(base_size = 10) +
      theme(panel.grid = element_blank(), strip.text.y = element_text(angle = 0),
            plot.title = element_text(face = "bold"))
    for (ext in c("png", "pdf"))
      ggsave(file.path(output_dir, paste0("adt_parameters_", metric, "_", ds, ".", ext)),
             p, width = 14, height = 14, dpi = 300, bg = "white")
  }
}
