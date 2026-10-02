# From project root: Rscript src/benchmarking/search_qc_params/plot_adt_qc_delta_f1.R [grid_dir]
library(dplyr)
library(ggplot2)

prepare_delta_data <- function(summary) {
  # Use the baseline from the same run, including when older runs are present.
  reference <- summary %>%
    filter(baseline) %>%
    distinct(dataset, run_id, gene_set, mode, scope, .keep_all = TRUE) %>%
    select(dataset, run_id, gene_set, mode, scope, baseline_F1 = F1) %>%
    mutate(baseline_found = TRUE)

  data <- summary %>%
    filter(!baseline) %>%
    arrange(desc(run_id)) %>%
    distinct(dataset, gene_set, mode, scope,
      min_difference, max_overlap, min_cells_per_status,
      .keep_all = TRUE) %>%
    #filter(min_cells_per_status == 5) %>%
    left_join(reference, by = c("dataset", "run_id", "gene_set", "mode", "scope"))

  if (any(is.na(data$baseline_found)))
    stop("Missing baseline for at least one plotted run/method.")
  if (!nrow(data)) stop("No grid rows with min_cells_per_status = 5.")

  data %>%
    mutate(
      delta_F1 = F1 - baseline_F1,
      method = paste0(gene_set, "\n", sub("gmm_dist_", "", mode), " | ", scope),
      across(c(min_difference, max_overlap, min_cells_per_status), factor)
    )
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  input_dir <- if (length(args)) args[1] else "results/benchmarking/qc_grid_no_na_full"
  output_dir <- file.path(input_dir, "plots", "filtered")
  data <- prepare_delta_data(read.csv(file.path(input_dir, "combined/grid_summary.csv")))
  delta_limit <- max(abs(data$delta_F1), 0.001, na.rm = TRUE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  for (ds in unique(data$dataset)) {
    points <- filter(data, dataset == ds) %>%
      mutate(
        label = ifelse(is.na(delta_F1), "", sprintf("%+.3f", delta_F1)),
        text_color = ifelse(is.na(delta_F1), "black",
          ifelse(abs(delta_F1) > 0.7 * delta_limit, "white", "black"))
      )
    p <- ggplot(points, aes(max_overlap, min_difference, fill = delta_F1)) +
      geom_tile(color = "white") +
      geom_text(aes(label = label, color = text_color), size = 2.8) +
      scale_color_identity() +
      facet_grid(method ~ min_cells_per_status,
        labeller = labeller(min_cells_per_status = function(x) paste0("min/status = ", x))) +
      scale_fill_gradient2(
        low = "#B94040", mid = "#FAFAFA", high = "#2465A0",
        midpoint = 0, limits = c(-delta_limit, delta_limit),
        labels = scales::label_number(accuracy = 0.01), na.value = "grey90"
      ) +
      labs(title = paste(ds, "- change in cell-level F1"),
        x = "max_overlap", y = "min_difference", fill = "F1 change",
        caption = paste("F1 minus the frequency-only baseline of the same run and method.",
          "Blue: improvement; red: decrease. Grey: undefined metric.")) +
      theme_minimal(base_size = 10) +
      theme(panel.grid = element_blank(), strip.text.y = element_text(angle = 0),
            plot.title = element_text(face = "bold"))
    for (ext in c("png", "pdf"))
      ggsave(file.path(output_dir, paste0("adt_parameters_delta_F1_", ds, ".", ext)),
             p, width = 14, height = 14, dpi = 300, bg = "white")
  }
}

if (sys.nframe() == 0L) main()
