# Run from the project root: Rscript src/benchmarking/search_qc_params/plot_frequency_qc_lines.R
library(dplyr)
library(ggplot2)

input_dir <- "results/benchmarking/frequency_qc_grid"
output_dir <- file.path(input_dir, "plots")
dir.create(output_dir, showWarnings = FALSE)
data <- bind_rows(
  read.csv(file.path(input_dir, "grid_summary.csv")) %>% filter(min_cells_per_pair == 10),
  read.csv("results/benchmarking/frequency_qc_no_filter/grid_summary.csv")
) %>% mutate(
  method = paste(gene_set, sub("gmm_dist_", "", mode), scope, sep = " | "),
  threshold = factor(ifelse(min_cells_per_pair == 0, "No QC", sprintf("%.2f", max_ci_width)),
                     levels = c("No QC", "0.20", "0.15", "0.10", "0.05"))
)

for (dataset_name in unique(data$dataset)) {
  p <- ggplot(filter(data, dataset == dataset_name),
              aes(threshold, F1, color = method, group = method)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.5) +
    scale_color_brewer(palette = "Dark2", limits = sort(unique(data$method))) +
    labs(title = paste(dataset_name, "- frequency QC (min. 10 cells per pair)"),
         subtitle = "",
         x = "Maximum CI width (stricter filtering to the right)",
         y = "Cell-level F1", color = NULL) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold")) +
    guides(color = guide_legend(ncol = 2))
  for (ext in c("png", "pdf"))
    ggsave(file.path(output_dir, paste0("ci_F1_", dataset_name, ".", ext)),
           p, width = 9, height = 6, dpi = 300, bg = "white")
}
