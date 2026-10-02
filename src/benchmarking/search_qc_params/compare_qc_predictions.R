# From project root: Rscript src/benchmarking/search_qc_params/compare_qc_predictions.R [grid_dir] [run_index]
# Baseline: no frequency/ADT QC, with the source run's metadata exclusion rule.
# Uses combined/run_config.csv, score caches and saved pair decisions; refits GMMs.
source("src/benchmarking/Benchmarking_QC_Grid.R")

compare_predictions <- function(cells, groups = character()) {
  cells %>% group_by(across(all_of(groups))) %>% summarise(
    n_input = n(), n_retained = sum(retained), n_removed = sum(!retained),
    n_baseline_missing = sum(retained & is.na(before)),
    n_qc_missing = sum(retained & is.na(after)),
    n_both_missing = sum(retained & is.na(before) & is.na(after)),
    n_reference_unknown = sum(retained & !reference_status %in% c("PLA", "platelet-free")),
    n_both_classified = sum(retained & !is.na(before) & !is.na(after)),
    n_unchanged = sum(retained & before == after, na.rm = TRUE),
    TP_to_TP = sum(old == "TP" & new == "TP", na.rm = TRUE),
    TP_to_FN = sum(old == "TP" & new == "FN", na.rm = TRUE),
    FN_to_TP = sum(old == "FN" & new == "TP", na.rm = TRUE),
    FN_to_FN = sum(old == "FN" & new == "FN", na.rm = TRUE),
    TN_to_TN = sum(old == "TN" & new == "TN", na.rm = TRUE),
    TN_to_FP = sum(old == "TN" & new == "FP", na.rm = TRUE),
    FP_to_TN = sum(old == "FP" & new == "TN", na.rm = TRUE),
    FP_to_FP = sum(old == "FP" & new == "FP", na.rm = TRUE),
    .groups = "drop"
  ) %>% mutate(
    n_compared = TP_to_TP + TP_to_FN + FN_to_TP + FN_to_FN +
      TN_to_TN + TN_to_FP + FP_to_TN + FP_to_FP,
    retained_fraction = ratio(n_retained, n_input),
    unchanged_fraction = ratio(n_unchanged, n_both_classified),
    n_corrected = FP_to_TN + FN_to_TP,
    n_introduced = TN_to_FP + TP_to_FN,
    baseline_errors = FP_to_TN + FP_to_FP + FN_to_TP + FN_to_FN,
    qc_errors = TN_to_FP + FP_to_FP + TP_to_FN + FN_to_FN,
    relative_error_reduction = ratio(n_corrected - n_introduced, baseline_errors),
    net_accuracy_gain = ratio(n_corrected - n_introduced, n_compared),
    baseline_F1_matched = ratio(2 * (TP_to_TP + TP_to_FN),
      2 * (TP_to_TP + TP_to_FN) + baseline_errors),
    qc_F1_matched = ratio(2 * (TP_to_TP + FN_to_TP),
      2 * (TP_to_TP + FN_to_TP) + qc_errors),
    delta_F1_matched = qc_F1_matched - baseline_F1_matched
  )
}

outcome <- function(prediction, reference) {
  if_else(!is.na(prediction) & reference %in% c("PLA", "platelet-free"),
    paste0(ifelse(prediction == reference, "T", "F"), ifelse(prediction == "PLA", "P", "N")),
    NA_character_)
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  input_dir <- if (length(args)) args[1] else "results/benchmarking/qc_grid"
  runs <- read.csv(file.path(input_dir, "combined/run_config.csv"))
  if (length(args) >= 2) runs <- runs[as.integer(args[2]), , drop = FALSE]

  for (i in seq_len(nrow(runs))) {
    run <- runs[i, ]
    config <- readRDS(file.path(run$run_dir, "run_config.rds"))
    cache <- if (nzchar(config$score_cache)) config$score_cache else file.path(run$run_dir, "scores.rds")
    exclude_na <- isTRUE(config$exclude_na_metadata)
    scores <- normalize_reference(readRDS(cache)$cells, exclude_na = exclude_na)
    grid <- read.csv(file.path(run$run_dir, "grid.csv"))
    pairs <- read.csv(file.path(run$run_dir, "pair_metrics.csv"),
                      colClasses = c(sample_id = "character"))
    output_dir <- file.path(input_dir, "prediction_comparison", run$dataset, run$run_id)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    message(run$dataset, " / ", run$gene_set, " / ", run$mode, " / ", run$scope, ": no-QC baseline")
    scores$before <- predict_gmm(scores, config$mode, config$scope, seed = config$provenance$seed)
    summaries <- details <- list()

    for (j in seq_len(nrow(grid))) {
      setting <- grid[j, ]
      kept <- pairs %>% filter(combination == setting$combination, retain_pair)
      cells <- semi_join(scores, kept, by = c("sample_id", "lineage"))
      cells$after <- predict_gmm(cells, config$mode, config$scope, seed = config$provenance$seed)
      compared <- scores %>%
        left_join(select(cells, cell_id, after), by = "cell_id") %>%
        mutate(retained = cell_id %in% cells$cell_id,
               old = outcome(before, reference_status), new = outcome(after, reference_status))
      id <- bind_cols(
        run %>% select(dataset, run_id, gene_set, mode, scope, variant),
        tibble(exclude_na_metadata = exclude_na, reference_baseline = "no_frequency_or_adt_qc"),
        setting %>% rename(frequency_only = baseline)
      )
      summaries[[j]] <- bind_cols(id, compare_predictions(compared))
      details[[j]] <- bind_cols(id, compare_predictions(compared, c("sample_id", "lineage")))
      message(j, "/", nrow(grid), " ", setting$combination, ": ", nrow(cells), "/", nrow(scores), " cells")
    }
    write.csv(bind_rows(summaries), file.path(output_dir, "prediction_comparison.csv"), row.names = FALSE)
    write.csv(bind_rows(details), file.path(output_dir, "prediction_comparison_by_celltype.csv"), row.names = FALSE)
    message("Results: ", output_dir)
  }
}

if (sys.nframe() == 0L) main()
