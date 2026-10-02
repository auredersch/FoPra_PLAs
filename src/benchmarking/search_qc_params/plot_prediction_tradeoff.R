# From project root: Rscript src/benchmarking/search_qc_params/plot_prediction_tradeoff.R [grid_dir]
library(dplyr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- file.path(if (length(args)) args[1] else "results/benchmarking/qc_grid",
                       "prediction_comparison")
output_dir <- file.path(input_dir, "plots")
dir.create(output_dir, showWarnings = FALSE)
files <- list.files(input_dir, "^prediction_comparison\\.csv$", recursive = TRUE, full.names = TRUE)
results <- bind_rows(lapply(files, read.csv)) %>%
  filter(is.finite(retained_fraction), is.finite(relative_error_reduction)) %>%
  mutate(method = paste(gene_set, sub("gmm_dist_", "", mode), scope, sep = " | "),
         setting = ifelse(frequency_only, "Frequency QC only", "Frequency + ADT QC"))

# Pareto frontier across methods; no QC (retention = 1, reduction = 0) is also an option.
frontier <- results %>% distinct(dataset, retained_fraction, relative_error_reduction) %>%
  group_by(dataset) %>%
  arrange(desc(retained_fraction), desc(relative_error_reduction), .by_group = TRUE) %>%
  filter(relative_error_reduction > lag(cummax(relative_error_reduction), default = 0),
         relative_error_reduction > 0) %>% ungroup()
best <- results %>% semi_join(frontier,
  by = c("dataset", "retained_fraction", "relative_error_reduction"))
write.csv(best, file.path(output_dir, "prediction_tradeoff_best_parameters.csv"), row.names = FALSE)
# Show one representative setting per frontier coordinate; the CSV contains all ties.
labels <- best %>% arrange(method, desc(frequency_only), min_difference, max_overlap, min_cells_per_status) %>%
  distinct(dataset, retained_fraction, relative_error_reduction, .keep_all = TRUE) %>%
  mutate(label = paste0(method, "\n", ifelse(frequency_only, "Frequency QC only",
    paste0("D=", min_difference, "; O=", max_overlap, "; n=", min_cells_per_status))))

for (ds in unique(results$dataset)) {
  p <- ggplot(filter(results, dataset == ds),
              aes(retained_fraction, relative_error_reduction, color = method)) +
    geom_hline(yintercept = 0, color = "grey55", linetype = "dashed") +
    geom_point(aes(shape = setting), size = 2.3, alpha = 0.65) +
    annotate("point", x = 1, y = 0, shape = 4, size = 3, color = "grey25") +
    ggrepel::geom_label_repel(data = filter(labels, dataset == ds), aes(label = label),
      size = 3, fill = "white", seed = 42, box.padding = 0.6, point.padding = 0.4,
      min.segment.length = 0, max.overlaps = Inf, show.legend = FALSE) +
    scale_x_continuous(labels = scales::label_percent(), limits = c(0, 1),
                       expand = expansion(mult = 0.06)) +
    scale_y_continuous(labels = scales::label_percent(), expand = expansion(mult = c(0.08, 0.25))) +
    scale_color_brewer(palette = "Dark2", limits = sort(unique(results$method))) +
    scale_shape_manual(values = c("Frequency QC only" = 17, "Frequency + ADT QC" = 16)) +
    labs(title = tools::toTitleCase(gsub("_", " ", ds)),
      subtitle = "Error reduction on the same retained, jointly classified cells",
      x = "Cells retained relative to no QC", y = "Relative error reduction vs. own no-QC predictions",
      color = NULL, shape = NULL,
      caption = paste("Labels: Pareto-optimal settings; one representative per coordinate. Cross: no QC.",
        "D = min_difference; O = max_overlap; n = min_cells_per_status. All equivalents in CSV.", sep = "\n")) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", legend.box = "vertical", panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"), plot.caption = element_text(hjust = 0, color = "grey40")) +
    guides(color = guide_legend(ncol = 2, order = 1), shape = guide_legend(order = 2))
  for (ext in c("png", "pdf"))
    ggsave(file.path(output_dir, paste0("prediction_tradeoff_", ds, ".", ext)),
           p, width = 12, height = 9, dpi = 300, bg = "white")
}
