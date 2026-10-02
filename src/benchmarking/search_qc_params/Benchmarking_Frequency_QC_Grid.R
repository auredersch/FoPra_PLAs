# From project root: Rscript .../Benchmarking_Frequency_QC_Grid.R
#   [adtnorm_input_dir] [output_dir] [grid_csv] [dataset] [exclude_na_metadata]
# Read ADTnorm objects directly; no saved QC-grid runs or QC tables are used.
# Import AUCell/GMM functions only; the sourced script does not execute its main().
source("src/benchmarking/search_qc_params/Benchmarking_QC_Grid.R")
source("src/gating_automation/config/datasets.R")

full_frequency_input <- function(run) {
  normalizePath(as.character(run$input), mustWork = TRUE)
}

# Define the same 5 datasets x 2 signatures x 2 GMM modes x 2 scopes independently.
frequency_run_config <- function(input_dir, dataset = "", exclude_na = FALSE) {
  datasets <- c("heart", "vaccine", "immune_aging", "sepsis", "skin")
  if (nzchar(dataset)) {
    if (!dataset %in% datasets) stop("Unknown dataset: ", dataset)
    datasets <- dataset
  }
  root <- Sys.getenv("PLA_PROJECT_ROOT", unset = getwd())
  signatures <- c(
    GOBP = "GOBP_REGULATION_OF_PLATELET_ACTIVATION.v2025.1.Hs.csv",
    MANNE_DN = "MANNE_COVID19_COMBINED_COHORT_VS_HEALTHY_DONOR_PLATELETS_DN.v2025.1.Hs.csv"
  )
  runs <- expand.grid(dataset = datasets, gene_set = names(signatures),
    mode = c("gmm_dist_platelet", "gmm_dist_dual"), scope = c("global", "per_celltype"),
    stringsAsFactors = FALSE) %>% arrange(dataset, gene_set, mode, scope)
  runs$input <- file.path(input_dir, runs$dataset, paste0(runs$dataset, "_ADTnorm_seurat.rds"))
  runs$variant <- Sys.getenv("PLA_QC_VARIANT", unset = "adt_tuned")
  runs$platelet_signature <- file.path(root, "data/signatures", unname(signatures[runs$gene_set]))
  runs$immune_signature <- file.path(root, "data/signatures",
    "GOBP_LEUKOCYTE_ACTIVATION_INVOLVED_IN_INFLAMMATORY_RESPONSE.v2025.1.Hs.csv")
  runs$seed <- 42L
  runs$exclude_na_metadata <- exclude_na
  runs$run_id <- paste(runs$dataset, runs$gene_set, runs$mode, runs$scope,
                       if (exclude_na) "no_na" else "all_metadata", sep = "_")
  runs
}

validate_full_frequency_scores <- function(scores) {
  if (anyDuplicated(scores$cell_id)) stop("Duplicate cell IDs in score cache.")
  if (anyNA(scores$sample_id)) stop("Missing sample IDs in unfiltered scores.")
  if (!all(c("lineage_original", "reference_status_original") %in% names(scores)))
    stop("Full score cache must preserve original metadata labels.")
  invisible(TRUE)
}

# Recalculate Wilson CI widths from the cells being analysed, before frequency filtering.
frequency_pairs <- function(scores, conf_level = 0.95) {
  pairs <- scores %>% group_by(sample_id, lineage) %>% summarise(
    n_cells = n(), n_PLA = sum(reference_status == "PLA"),
    n_platelet_free = sum(reference_status == "platelet-free"), .groups = "drop")
  z <- qnorm(1 - (1 - conf_level) / 2)
  p <- pairs$n_PLA / pairs$n_cells
  denominator <- 1 + z^2 / pairs$n_cells
  center <- (p + z^2 / (2 * pairs$n_cells)) / denominator
  half_width <- z * sqrt(p * (1 - p) / pairs$n_cells + z^2 / (4 * pairs$n_cells^2)) / denominator
  pairs$ci_lower <- pmax(0, center - half_width)
  pairs$ci_upper <- pmin(1, center + half_width)
  pairs$ci_width <- pairs$ci_upper - pairs$ci_lower
  pairs
}

full_frequency_scores <- function(run, cache_root) {
  input <- full_frequency_input(run)
  signatures <- c(platelet = as.character(run$platelet_signature),
                  immune = as.character(run$immune_signature))
  signatures <- vapply(signatures, normalizePath, character(1), mustWork = TRUE)
  provenance <- list(input_stage = "before_frequency_qc", input = input,
    input_size = file.info(input)$size, input_mtime = as.numeric(file.info(input)$mtime),
    dataset = run$dataset, sample_col = get_pla_dataset_config(run$dataset)$sample_col,
    signatures = tools::md5sum(signatures), seed = as.integer(run$seed),
    scoring = "AUCell counts; fixed genes; aucMaxRank=ceiling(n_genes*0.05); no extension",
    cell_coverage_check = "all input metadata cells must have scores")
  path <- file.path(cache_root, run$variant, run$dataset, run$gene_set, "scores.rds")
  cached <- if (file.exists(path)) readRDS(path) else NULL
  reuse <- !is.null(cached) && identical(cached$provenance, provenance)
  if (reuse) {
    scores <- cached$cells
    message("Reusing full score cache: ", path)
  } else {
    message("Scoring full input before frequency QC: ", input)
    scores <- prepare_scores(input, provenance$sample_col, signatures, seed = provenance$seed)
  }
  validate_full_frequency_scores(scores)
  if (!reuse) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(list(provenance = provenance, cells = scores), path)
  }
  list(cells = scores, path = normalizePath(path), input = input)
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  arg <- function(i, default) if (length(args) >= i && nzchar(args[i])) args[i] else default
  root <- Sys.getenv("PLA_PROJECT_ROOT", unset = getwd())
  input_dir <- arg(1, Sys.getenv("PLA_ADTNORM_TUNED_INPUT_DIR",
                                unset = file.path(root, "src/PLA_QC_ADTnorm_tuned")))
  metadata_policy <- toupper(arg(5, "FALSE"))
  if (!metadata_policy %in% c("TRUE", "FALSE")) stop("exclude_na_metadata must be TRUE or FALSE.")
  no_na <- metadata_policy == "TRUE"
  runs <- frequency_run_config(input_dir, dataset = arg(4, ""), exclude_na = no_na)
  required_files <- unique(c(runs$input, runs$platelet_signature, runs$immune_signature))
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files)) stop("Missing input files: ", paste(missing_files, collapse = ", "))
  output_dir <- arg(2, if (no_na) "results/benchmarking/frequency_qc_full_no_na" else
                      "results/benchmarking/frequency_qc_full")
  cache_root <- Sys.getenv("PLA_FREQUENCY_SCORE_CACHE_DIR",
                          "results/benchmarking/frequency_qc_score_cache")
  grid <- if (nzchar(arg(3, ""))) read.csv(args[3]) else
    expand.grid(max_ci_width = c(0.05, 0.10, 0.15, 0.20, 1.0),
                min_cells_per_pair = c(0L, 5L, 10L, 15L, 20L, 50L))
  stopifnot(all(c("max_ci_width", "min_cells_per_pair") %in% names(grid)),
            all(is.finite(grid$max_ci_width) & grid$max_ci_width >= 0 & grid$max_ci_width <= 1),
            all(is.finite(grid$min_cells_per_pair) & grid$min_cells_per_pair >= 0 &
                grid$min_cells_per_pair == floor(grid$min_cells_per_pair)))
  grid <- bind_rows(data.frame(max_ci_width = 1, min_cells_per_pair = 0L), grid) %>%
    distinct(max_ci_width, min_cells_per_pair)
  grid$baseline <- grid$max_ci_width == 1 & grid$min_cells_per_pair == 0
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  summaries <- pair_results <- frequency_tables <- list()
  score_key <- NULL

  for (i in seq_len(nrow(runs))) {
    run <- runs[i, ]
    key <- paste(run$dataset, run$variant, run$gene_set, run$platelet_signature,
                 run$immune_signature, run$seed, run$input, sep = "|")
                 
    message(i, "/", nrow(runs), ": ", run$dataset, " / ", run$gene_set, " / ",
            run$mode, " / ", run$scope, "; exclude_na_metadata=", run$exclude_na_metadata)

    if (!identical(key, score_key)) {
      full <- full_frequency_scores(run, cache_root)
      score_key <- key
    }
    scores <- normalize_reference(full$cells, exclude_na = isTRUE(run$exclude_na_metadata))
    if (!nrow(scores)) stop("No cells remain after the metadata filter for ", run$dataset)
    pairs <- frequency_pairs(scores)
    frequency_tables[[i]] <- bind_cols(
      run %>% select(dataset, gene_set, mode, scope, variant, run_id, exclude_na_metadata), pairs)
    runs$full_input[i] <- full$input
    runs$full_score_cache[i] <- full$path
    before <- scores %>% group_by(sample_id, lineage) %>% summarise(
      n_input_cells = n(), n_input_PLA = sum(reference_status == "PLA"), .groups = "drop")
    message(i, "/", nrow(runs), ": ", run$dataset, " / ", run$gene_set, " / ",
            run$mode, " / ", run$scope, "; ", nrow(full$cells), " full cells, ",
            nrow(scores), " cells after metadata policy")

    for (j in seq_len(nrow(grid))) {
      setting <- grid[j, ]
      kept <- pairs %>% filter(ci_width <= setting$max_ci_width,
                               n_cells >= setting$min_cells_per_pair)
      cells <- if (setting$baseline) scores else
        semi_join(scores, kept, by = c("sample_id", "lineage"))
      cells$Prediction <- predict_gmm(cells, run$mode, run$scope, seed = run$seed)
      evaluated <- evaluate_pairs(cells)
      id <- bind_cols(run %>% select(dataset, gene_set, mode, scope, variant,
                                      run_id, exclude_na_metadata),
        tibble(full_input = full$input, full_score_cache = full$path), setting)
      per_pair <- before %>% left_join(evaluated, by = c("sample_id", "lineage")) %>%
        mutate(n_retained = coalesce(n_retained, 0L), retain_pair = n_retained > 0)
      pair_results[[length(pair_results) + 1L]] <- bind_cols(id, per_pair)
      summary <- evaluated %>% summarise(
        across(c(TP, FP, FN, TN, n_evaluated), sum),
        macro_F1 = if (all(is.na(F1))) NA_real_ else mean(F1, na.rm = TRUE),
        n_pairs_with_F1 = sum(!is.na(F1))
      ) %>% mutate(
        F1 = ratio(2 * TP, 2 * TP + FP + FN), precision = ratio(TP, TP + FP), recall = ratio(TP, TP + FN),
        n_full_input_cells = nrow(full$cells), n_removed_metadata = nrow(full$cells) - nrow(scores),
        n_input_cells = nrow(scores), n_retained_cells = nrow(cells),
        n_input_pairs = nrow(before), n_retained_pairs = nrow(evaluated),
        n_input_PLA = sum(before$n_input_PLA), n_retained_PLA = sum(evaluated$n_reference_positive),
        retained_cell_fraction = ratio(n_retained_cells, n_input_cells),
        retained_cell_fraction_full_input = ratio(n_retained_cells, n_full_input_cells),
        retained_pair_fraction = ratio(n_retained_pairs, n_input_pairs),
        retained_PLA_fraction = ratio(n_retained_PLA, n_input_PLA),
        n_unclassified = sum(is.na(cells$Prediction)),
        classification_coverage = ratio(sum(!is.na(cells$Prediction)), n_retained_cells)
      )
      summaries[[length(summaries) + 1L]] <- bind_cols(id, summary)
      message(j, "/", nrow(grid), ": retained ", nrow(cells), "/", nrow(scores), " cells")
    }
  }
  write.csv(runs, file.path(output_dir, "run_config.csv"), row.names = FALSE)
  write.csv(grid, file.path(output_dir, "grid.csv"), row.names = FALSE)
  write.csv(bind_rows(summaries), file.path(output_dir, "grid_summary.csv"), row.names = FALSE)
  write.csv(bind_rows(pair_results), file.path(output_dir, "pair_metrics.csv"), row.names = FALSE)
  write.csv(bind_rows(frequency_tables), file.path(output_dir, "frequency_pairs.csv"), row.names = FALSE)
  message("Results: ", output_dir)
}

if (sys.nframe() == 0L) main()
