library(dplyr)
library(plotly)
library(htmlwidgets)

input_dir <- "results/benchmarking/qc_grid_no_na/prediction_comparison"
output_dir <- file.path(input_dir, "interactive")
dir.create(output_dir, showWarnings = FALSE)

files <- list.files(input_dir, "^prediction_comparison\\.csv$",
                    recursive = TRUE, full.names = TRUE)

results <- bind_rows(lapply(files, read.csv)) %>%
  mutate(
    method = paste(gene_set, sub("gmm_dist_", "", mode), scope, sep = " | "),
    baseline_specificity = (TN_to_TN + TN_to_FP) /
      (TN_to_TN + TN_to_FP + FP_to_TN + FP_to_FP),
    qc_specificity = (TN_to_TN + FP_to_TN) /
      (TN_to_TN + TN_to_FP + FP_to_TN + FP_to_FP),
    baseline_sensitivity = (TP_to_TP + TP_to_FN) /
      (TP_to_TP + TP_to_FN + FN_to_TP + FN_to_FN),
    qc_sensitivity = (TP_to_TP + FN_to_TP) /
      (TP_to_TP + TP_to_FN + FN_to_TP + FN_to_FN),
    parameters = paste0(
      ifelse(frequency_only, "Frequency QC only",
        paste0("D=", min_difference, "; O=", max_overlap,
               "; n=", min_cells_per_status)),
      sprintf(" | retained: %.1f%%", 100 * retained_fraction)
    )
  )

results <- results %>%
  mutate(score_text = sprintf(
    paste0(
      "F1: %.4f → %.4f<br>",
      "Specificity: %.4f → %.4f<br>",
      "Sensitivity: %.4f → %.4f<br>",
      "Relative error reduction: %.1f%%"
    ),
    baseline_F1_matched, qc_F1_matched,
    baseline_specificity, qc_specificity,
    baseline_sensitivity, qc_sensitivity,
    100 * relative_error_reduction
  ))

columns <- list(
  F1 = c("baseline_F1_matched", "qc_F1_matched"),
  Specificity = c("baseline_specificity", "qc_specificity"),
  Sensitivity = c("baseline_sensitivity", "qc_sensitivity")
)

colors <- setNames(
  RColorBrewer::brewer.pal(8, "Dark2"),
  sort(unique(results$method))
)

for (metric in names(columns)) {
  data <- results %>%
    mutate(x = .data[[columns[[metric]][1]]],
           y = .data[[columns[[metric]][2]]]) %>%
    filter(is.finite(x), is.finite(y))

  for (ds in unique(data$dataset)) {
    rows <- filter(data, dataset == ds) %>%
      group_by(method, x, y) %>%
      mutate(point_id = cur_group_id()) %>%
      ungroup()

    points <- rows %>%
      group_by(point_id, method, x, y) %>%
      summarise(
        settings = paste(unique(parameters), collapse = "<br>"),
        .groups = "drop"
      ) %>%
      mutate(hover = paste0(
        "<b>Point ", point_id, "</b><br>", method,
        sprintf("<br>Baseline: %.4f<br>QC: %.4f<br>Difference: %+.4f",
                x, y, y - x),
        "<br><br>", settings
      ))

    limits <- c(0, min(1, max(c(points$x, points$y, 0.1)) * 1.05))

    p <- plot_ly(
      points, x = ~x, y = ~y, color = ~method, colors = colors,
      type = "scatter", mode = "markers",
      text = ~hover, hoverinfo = "text",
      marker = list(size = 9, opacity = 0.8)
    ) %>%
      layout(
        title = paste(ds, metric),
        xaxis = list(title = paste("Baseline", metric), range = limits),
        yaxis = list(title = paste("QC", metric), range = limits,
                     scaleanchor = "x", scaleratio = 1),
        shapes = list(list(
          type = "line",
          x0 = 0, y0 = 0, x1 = limits[2], y1 = limits[2],
          line = list(color = "grey", dash = "dash")
        )),
        hoverlabel = list(align = "left")
      )

    filename <- file.path(output_dir, paste0(tolower(metric), "_", ds))
    saveWidget(p, paste0(filename, ".html"), selfcontained = TRUE)
    write.csv(rows, paste0(filename, ".csv"), row.names = FALSE)
  }
}