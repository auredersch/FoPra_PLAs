# Run from the project root: Rscript src/benchmarking/search_qc_params/Benchmarking_Frequency_QC_Grid.R
# Uses combine_qc_runs.R output and existing AUCell caches; no additional ADT QC.
source("src/benchmarking/Benchmarking_QC_Grid.R")

runs <- read.csv("results/benchmarking/qc_grid_no_na/combined/run_config.csv")
#grid <- expand.grid(max_ci_width = c(0.05, 0.10, 0.15, 0.20), min_cells_per_pair = c(5L, 10L, 15L, 20L, 50L))
grid <- data.frame(max_ci_width = 1, min_cells_per_pair = 0L)
#output_dir <- "results/benchmarking/frequency_qc_grid"
output_dir <- "results/benchmarking/frequency_qc_no_na"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
summaries <- pair_results <- list()

for (i in seq_len(nrow(runs))) {
  run <- runs[i, ]
  scores <- normalize_reference(
    readRDS(file.path(run$run_dir, "scores.rds"))$cells,
    exclude_na = isTRUE(run$exclude_na_metadata)
  )  
  
  pairs <- read.csv(file.path(run$qc_dir, "01_sample_celltype_summary.csv"),
                    colClasses = c(sample_id = "character")) %>%
    transmute(sample_id, lineage = canonical_lineage(celltype_id), ci_width, n_cells)
  before <- scores %>% group_by(sample_id, lineage) %>% summarise(
    n_input_cells = n(), n_input_PLA = sum(reference_status == "PLA", na.rm = TRUE),
    .groups = "drop"
  )
  message(i, "/", nrow(runs), ": ", run$dataset, " / ", run$gene_set,
          " / ", run$mode, " / ", run$scope)

  for (j in seq_len(nrow(grid))) {
    setting <- grid[j, ]
    kept <- pairs %>% filter(ci_width <= setting$max_ci_width,
                             n_cells >= setting$min_cells_per_pair)
    #cells <- semi_join(scores, kept, by = c("sample_id", "lineage"))
    cells <- if (setting$max_ci_width == 1 && setting$min_cells_per_pair == 0) {
      scores
    } else {
      semi_join(scores, kept, by = c("sample_id", "lineage"))
    }
    cells$Prediction <- predict_gmm(cells, run$mode, run$scope, seed = run$seed)
    evaluated <- evaluate_pairs(cells)
    id <- bind_cols(run %>% select(dataset, gene_set, mode, scope, variant,
                                   source_run_id = run_id), setting)
    per_pair <- before %>%
      left_join(evaluated, by = c("sample_id", "lineage")) %>%
      mutate(n_retained = coalesce(n_retained, 0L), retain_pair = n_retained > 0)
    pair_results[[length(pair_results) + 1L]] <- bind_cols(id, per_pair)

    summary <- evaluated %>% summarise(
      across(c(TP, FP, FN, TN, n_evaluated), sum),
      macro_F1 = mean(F1, na.rm = TRUE), n_pairs_with_F1 = sum(!is.na(F1))
    ) %>% mutate(
      F1 = ratio(2 * TP, 2 * TP + FP + FN),
      precision = ratio(TP, TP + FP), recall = ratio(TP, TP + FN),
      n_input_cells = nrow(scores), n_retained_cells = nrow(cells),
      n_input_pairs = nrow(before), n_retained_pairs = nrow(evaluated),
      n_input_PLA = sum(before$n_input_PLA),
      n_retained_PLA = sum(evaluated$n_reference_positive),
      retained_cell_fraction = ratio(n_retained_cells, n_input_cells),
      retained_pair_fraction = ratio(n_retained_pairs, n_input_pairs),
      retained_PLA_fraction = ratio(n_retained_PLA, n_input_PLA),
      n_unclassified = sum(is.na(cells$Prediction)),
      classification_coverage = ratio(n_retained_cells - n_unclassified, n_retained_cells)
    )
    summaries[[length(summaries) + 1L]] <- bind_cols(id, summary)
  }
}

write.csv(runs, file.path(output_dir, "run_config.csv"), row.names = FALSE)
write.csv(grid, file.path(output_dir, "grid.csv"), row.names = FALSE)
write.csv(bind_rows(summaries), file.path(output_dir, "grid_summary.csv"), row.names = FALSE)
write.csv(bind_rows(pair_results), file.path(output_dir, "pair_metrics.csv"), row.names = FALSE)
message("Results: ", output_dir)
