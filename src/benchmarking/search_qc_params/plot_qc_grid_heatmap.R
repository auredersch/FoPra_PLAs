# Run from the project root; optional argument: results/benchmarking/qc_grid_no_na
library(dplyr)
library(tidyr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args)) args[1] else "results/benchmarking/qc_grid_no_na"
output_dir <- file.path(input_dir, "plots")
dir.create(output_dir, showWarnings = FALSE)

# Reading configs also includes older runs and skips the orphan summary duplicate.
configs <- list.files(input_dir, "^run_config\\.rds$", recursive = TRUE, full.names = TRUE)
baseline <- bind_rows(lapply(configs, function(path) {
  config <- readRDS(path)
  signature <- basename(names(config$provenance$signatures)[1])
  read.csv(file.path(dirname(path), "grid_summary.csv")) %>%
    filter(baseline) %>%
    mutate(
      gene_set = ifelse(grepl("^MANNE", signature), "MANNE_DN", "GOBP"),
      mode = sub("gmm_dist_", "", config$mode),
      scope = config$scope,
      F1 = 2 * TP / (2 * TP + FP + FN),
      precision = TP / (TP + FP),
      recall = TP / (TP + FN)
    )
}))

metric_labels <- c(
  F1 = "Cell-level F1", precision = "Cell-level precision",
  recall = "Cell-level recall", macro_F1 = "Macro F1"
)
plot_data <- baseline %>%
  mutate(method = paste(gene_set, mode, scope, sep = " | ")) %>%
  select(dataset, method, all_of(names(metric_labels))) %>%
  pivot_longer(all_of(names(metric_labels)), names_to = "metric", values_to = "value") %>%
  mutate(
    label = sprintf("%.3f", value),
    method = factor(method, levels = rev(sort(unique(method)))),
    metric = factor(metric, levels = names(metric_labels), labels = metric_labels)
  )

p <- ggplot(plot_data, aes(dataset, method, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label,
                color = ifelse(is.na(value) | value > 0.5, "black", "white")), size = 3) +
  scale_color_identity() +
  scale_fill_viridis_c(limits = c(0, 1), na.value = "grey85") +
  facet_wrap(~metric, ncol = 2) +
  labs(title = "QC grid: baseline method comparison",
       x = NULL, y = NULL, fill = "Value") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(output_dir, "baseline_heatmap.png"), p, width = 13, height = 8, dpi = 300)
ggsave(file.path(output_dir, "baseline_heatmap.pdf"), p, width = 13, height = 8)
