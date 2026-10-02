# Run from the project root: Rscript src/benchmarking/search_qc_params/plot_frequency_qc_heatmaps.R
library(dplyr)
library(tidyr)
library(ggplot2)

input_dir <- "results/benchmarking/frequency_qc_grid"
output_dir <- file.path(input_dir,  "plots", "no_redundancy")
dir.create(output_dir, showWarnings = FALSE)
data <- bind_rows(
  read.csv(file.path(input_dir, "grid_summary.csv")),
  read.csv("results/benchmarking/frequency_qc_no_filter/grid_summary.csv")
) %>%
  arrange(max_ci_width, min_cells_per_pair) %>%
  filter(min_cells_per_pair == 10 | (max_ci_width == 1 & min_cells_per_pair == 0)) %>%
  mutate(method = paste(gene_set, sub("gmm_dist_", "", mode), scope, sep = " | "),
         setting = ifelse(max_ci_width == 1 & min_cells_per_pair == 0, "No QC",
                          sprintf("%.2f / %d", max_ci_width, min_cells_per_pair))) %>%
  mutate(setting = factor(setting, levels = unique(setting)),
         method = factor(method, levels = rev(sort(unique(method)))))
reference <- data %>% filter(setting == "No QC") %>%
  select(dataset, method, reference_F1 = F1)
deltas <- data %>% left_join(reference, by = c("dataset", "method")) %>%
  mutate(delta_F1 = F1 - reference_F1)
delta_limit <- max(abs(deltas$delta_F1), 0.001, na.rm = TRUE)
save_plot <- function(p, name, width, height) {
  for (ext in c("png", "pdf"))
    ggsave(file.path(output_dir, paste0(name, ".", ext)), p,
           width = width, height = height, dpi = 300, bg = "white")
}
plot_theme <- theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), plot.title = element_text(face = "bold"))

for (dataset_name in unique(data$dataset)) {
  p <- ggplot(filter(data, dataset == dataset_name), aes(setting, method, fill = F1)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.3f", F1), color = F1 < 0.5), size = 2.7) +
    scale_color_manual(values = c("FALSE" = "grey20", "TRUE" = "white"), guide = "none") +
    scale_fill_viridis_c(limits = c(0, 1), na.value = "grey85") +
    labs(title = paste(dataset_name, "- cell-level F1"),
         #subtitle = "Absolute scores; No QC = no frequency filtering",
         x = "Maximum CI width / minimum cells per pair", y = NULL, fill = "F1") +
    plot_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p, paste0("F1_", dataset_name), 15, 5)

  p <- ggplot(filter(deltas, dataset == dataset_name), aes(setting, method, fill = delta_F1)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%+.3f", delta_F1),
                  color = abs(delta_F1) > 0.7 * delta_limit), size = 2.7) +
    scale_color_manual(values = c("FALSE" = "grey20", "TRUE" = "white"), guide = "none") +
    scale_fill_gradient2(low = "#B94040", mid = "#FAFAFA", high = "#2465A0",
                         midpoint = 0, limits = c(-delta_limit, delta_limit), na.value = "grey85") +
    labs(title = paste(dataset_name, "- change in cell-level F1"),
         #subtitle = "Relative to No QC within the same method; blue = improvement",
         x = "Maximum CI width / minimum cells per pair", y = NULL, fill = "F1 change") +
    plot_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p, paste0("delta_F1_", dataset_name), 15, 5)

  retention <- data %>% filter(dataset == dataset_name) %>%
    distinct(setting, retained_cell_fraction, retained_PLA_fraction, retained_pair_fraction) %>%
    pivot_longer(-setting, names_to = "metric", values_to = "fraction") %>%
    mutate(metric = factor(metric,
      levels = c("retained_cell_fraction", "retained_PLA_fraction", "retained_pair_fraction"),
      labels = c("All cells", "PLA cells", "Sample-celltype pairs")))
  p <- ggplot(retention, aes(metric, setting, fill = fraction)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.1f%%", 100 * fraction), color = fraction < 0.5), size = 3) +
    scale_color_manual(values = c("FALSE" = "grey20", "TRUE" = "white"), guide = "none") +
    scale_y_discrete(limits = rev(levels(data$setting))) +
    scale_fill_viridis_c(limits = c(0, 1), labels = scales::label_percent()) +
    labs(title = paste(dataset_name, "- retention"), #subtitle = "Relative to all cells / pairs before frequency QC",
         x = NULL, y = "Maximum CI width / minimum cells per pair", fill = "Retained") +
    plot_theme
  save_plot(p, paste0("retention_", dataset_name), 8, 8)
}
