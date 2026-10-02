# Run from the project root: Rscript src/benchmarking/search_qc_params/plot_qc_tradeoff.R
library(dplyr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args)) args[1] else "results/benchmarking/qc_grid"
input_file <- file.path(input_dir, "combined/grid_summary.csv")
output_dir <- file.path(input_dir, "plots")
dir.create(output_dir, showWarnings = FALSE)

results <- read.csv(input_file) %>%
  mutate(
    method = paste(gene_set, sub("gmm_dist_", "", mode), scope, sep = " | "),
    setting = ifelse(baseline, "Baseline", "ADT QC")
  ) %>%
  filter(is.finite(F1))
data <- results %>%
  distinct(dataset, method, setting, retained_cell_fraction, n_retained_cells, F1)

# Across all methods: keep points with no alternative having both higher/equal
# retention and higher/equal F1, with at least one strict improvement.
frontier <- data %>%
  distinct(dataset, retained_cell_fraction, n_retained_cells, F1) %>%
  group_by(dataset) %>%
  arrange(desc(retained_cell_fraction), desc(F1), .by_group = TRUE) %>%
  filter(F1 > lag(cummax(F1), default = -Inf)) %>%
  mutate(point = row_number()) %>%
  ungroup()

for (dataset_name in unique(data$dataset)) {
  points <- filter(frontier, dataset == dataset_name)
  parameters <- results %>%
    inner_join(points, by = c("dataset", "retained_cell_fraction", "n_retained_cells", "F1")) %>%
    select(dataset, point, method, gene_set, mode, scope, baseline,
           min_difference, max_overlap, min_cells_per_status,
           n_retained_cells, retained_cell_fraction, F1, run_id, combination) %>%
    arrange(point, method, desc(baseline), min_difference, max_overlap, min_cells_per_status)
  write.csv(parameters, file.path(output_dir, paste0("qc_tradeoff_parameters_", dataset_name, ".csv")),
            row.names = FALSE)
  p <- ggplot(filter(data, dataset == dataset_name),
              aes(retained_cell_fraction, F1, color = method)) +
    geom_point(aes(shape = setting), size = 2.4, alpha = 0.85) +
    geom_point(data = points, aes(retained_cell_fraction, F1),
               inherit.aes = FALSE, shape = 1, size = 4.3, stroke = 0.5, color = "grey25") +
    ggrepel::geom_label_repel(
      data = points, aes(retained_cell_fraction, F1, label = point),
      inherit.aes = FALSE, color = "grey20", fill = "white", label.size = NA,
      size = 3, fontface = "bold", seed = 42, box.padding = 0.3, point.padding = 0.6,
      segment.color = "grey65", segment.size = 0.3, min.segment.length = 0,
      max.overlaps = Inf
    ) +
    scale_x_continuous(labels = scales::label_percent(), limits = c(0, 1)) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
    scale_color_brewer(palette = "Dark2", limits = sort(unique(data$method))) +
    scale_shape_manual(values = c("ADT QC" = 16, "Baseline" = 17)) +
    labs(title = tools::toTitleCase(gsub("_", " ", dataset_name)),
         subtitle = "Performance and cell retention",
         x = "Cells retained relative to frequency-QC baseline", y = "Cell-level F1",
         color = NULL, shape = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 17),
      plot.subtitle = element_text(color = "grey45"),
      legend.position = "bottom", legend.box = "vertical", legend.text = element_text(size = 9),
      panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "grey93"),
      plot.margin = margin(12, 18, 8, 10)
    ) +
    guides(color = guide_legend(ncol = 2, order = 1), shape = guide_legend(order = 2))

  values <- points %>% transmute(
    Point = point, Cells = format(n_retained_cells, big.mark = ",", trim = TRUE),
    F1 = sprintf("%.3f", F1)
  )
  table <- gridExtra::tableGrob(values, rows = NULL, theme = gridExtra::ttheme_minimal(
    base_size = 10, padding = grid::unit(c(3, 4), "mm"),
    core = list(bg_params = list(fill = c("#F4F6F8", "white"), col = NA)),
    colhead = list(fg_params = list(fontface = "bold"),
                   bg_params = list(fill = "#E8EDF2", col = NA))
  ))
  side <- gridExtra::arrangeGrob(
    grid::textGrob("Pareto-optimal points", gp = grid::gpar(fontsize = 11, fontface = "bold")),
    table, grid::nullGrob(), ncol = 1,
    heights = grid::unit.c(grid::unit(0.5, "in"), sum(table$heights), grid::unit(1, "null"))
  )
  figure <- gridExtra::arrangeGrob(p, side, ncol = 2, widths = c(3.6, 1))
  ggsave(file.path(output_dir, paste0("qc_tradeoff_", dataset_name, ".png")),
         figure, width = 11, height = 7.5, dpi = 300, bg = "white")
  ggsave(file.path(output_dir, paste0("qc_tradeoff_", dataset_name, ".pdf")),
         figure, width = 11, height = 7.5, bg = "white")
}
