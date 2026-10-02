suppressPackageStartupMessages({
  library(dplyr)
  library(mclust)
})

canonical_lineage <- function(x) {
  x <- coalesce(as.character(x), "Unassigned")
  recode(x, "CD4 T" = "CD4 T cells", "PLA neutrohils" = "Neutrophils")
}

# Match the status convention used by the frequency/ADT QC tables.
# Applied to both fresh scores and older caches; preserve the original labels.
normalize_reference <- function(scores, exclude_na = FALSE) {
  if (!"reference_status_original" %in% names(scores))
    scores$reference_status_original <- as.character(scores$reference_status)
  if (exclude_na)
    scores <- filter(scores,
      !is.na(lineage_original), tolower(trimws(lineage_original)) != "unassigned",
      !is.na(reference_status_original), tolower(trimws(reference_status_original)) != "unassigned"
    )
  original <- as.character(scores$reference_status_original)
  scores$reference_status <- if_else(
    is.na(original) | tolower(trimws(original)) == "unassigned",
    "platelet-free", original
  )
  scores
}

prepare_scores <- function(input, sample_col, signature_files, seed = 42L) {
  suppressPackageStartupMessages({ library(Seurat); library(AUCell) })
  pbmc <- readRDS(input)
  meta <- pbmc[[]]
  counts <- GetAssayData(pbmc, assay = "RNA", layer = "counts")
  signatures <- lapply(signature_files, function(path) unique(read.csv(path)$geneName))
  signatures <- lapply(signatures, function(genes) {
    matched <- intersect(genes, rownames(counts))
    if (length(matched) < 3 || length(matched) / length(genes) < 0.5)
      stop("Too few signature genes match RNA counts; check gene identifiers.")
    matched
  })
  set.seed(seed)
  rankings <- AUCell_buildRankings(counts, plotStats = FALSE)
  auc <- getAUC(AUCell_calcAUC(signatures, rankings, aucMaxRank = ceiling(nrow(counts) * 0.05)))
  cells <- colnames(auc)
  meta <- meta[cells, , drop = FALSE]
  lineage_col <- if ("lineage" %in% names(meta)) "lineage" else "celltype_clean"
  reference_col <- intersect(c("pla_status", "pla.status"), names(meta))
  if (!length(reference_col)) stop("Input has no pla_status/pla.status reference column.")
  reference <- as.character(meta[[reference_col[1]]])
  tibble(
    cell_id = cells, sample_id = as.character(meta[[sample_col]]),
    lineage = canonical_lineage(meta[[lineage_col]]),
    lineage_original = as.character(meta[[lineage_col]]),
    reference_status_original = reference,
    reference_status = reference,  # Normalized together with cached labels below.
    Raw_Score = as.numeric(auc["platelet", cells]),
    Immune_Score = as.numeric(auc["immune", cells])
  )
}

# Each level passes if any testable CORE marker passes both thresholds.
# Explicitly untestable groups are retained; missing table rows are not.
qc_keep <- function(metrics, id, min_difference, max_overlap, qc, min_cells_per_status = qc$min_cells_per_status) {
  metrics %>%
    filter(marker_label %in% qc$core_marker_labels) %>%
    mutate(
      testable = n_PLA >= min_cells_per_status &
        n_platelet_free >= min_cells_per_status &
        is.finite(median_difference) & is.finite(platelet_free_above_pla_q25),
      pass = testable & median_difference >= min_difference &
        platelet_free_above_pla_q25 <= max_overlap
    ) %>%
    group_by(across(all_of(id))) %>%
    summarise(n_testable = sum(testable, na.rm = TRUE),
              keep = n_testable == 0 | any(pass, na.rm = TRUE), .groups = "drop")
}

select_pairs <- function(pairs, sample_metrics, lineage_metrics, min_difference, max_overlap, qc, min_cells_per_status = qc$min_cells_per_status) {
  samples <- qc_keep(sample_metrics, "sample_id", min_difference, max_overlap, qc, min_cells_per_status) %>%
    rename(sample_keep = keep, sample_testable_markers = n_testable)
  lineages <- qc_keep(lineage_metrics, "lineage", min_difference, max_overlap, qc, min_cells_per_status) %>%
    rename(lineage_keep = keep, lineage_testable_markers = n_testable)
  pairs %>%
    left_join(samples, by = "sample_id") %>%
    left_join(lineages, by = "lineage") %>%
    mutate(retain_pair = usable_celltype_sample &
             coalesce(sample_keep, FALSE) & coalesce(lineage_keep, FALSE))
}

fit_high <- function(x) {
  high <- rep(NA, length(x))
  valid <- which(is.finite(x))
  if (length(valid) < 3 || length(unique(x[valid])) < 2) return(high)
  fit <- tryCatch(Mclust(x[valid], G = 2, verbose = FALSE), error = function(e) NULL)
  if (!is.null(fit)) high[valid] <- fit$classification == which.max(fit$parameters$mean)
  high
}

predict_gmm <- function(cells, mode, scope, seed = 42L) {
  set.seed(seed)
  positive <- rep(NA, nrow(cells))
  groups <- if (scope == "global") list(seq_len(nrow(cells))) else
    split(seq_len(nrow(cells)), cells$lineage)
  for (idx in groups) {
    if (scope == "per_celltype" && length(idx) < 100) next
    platelet_high <- fit_high(cells$Raw_Score[idx])
    prediction <- platelet_high
    if (mode == "gmm_dist_dual") {
      selected <- which(platelet_high %in% TRUE)
      # Failed immune fits leave platelet-high cells unclassified, not negative.
      prediction[selected] <- if (scope == "per_celltype" && length(selected) < 100)
        rep(NA, length(selected)) else fit_high(cells$Immune_Score[idx[selected]])
    }
    positive[idx] <- prediction
  }
  ifelse(is.na(positive), NA_character_, ifelse(positive, "PLA", "platelet-free"))
}

ratio <- function(a, b) ifelse(b > 0, a / b, NA_real_)

# Scores concern retained, classified cells with known reference labels only.
# Coverage is reported separately; missing predictions must not silently improve F1.
evaluate_pairs <- function(cells) {
  cells %>% group_by(sample_id, lineage) %>% summarise(
    n_retained = n(), n_classified = sum(!is.na(Prediction)),
    n_reference = sum(reference_status %in% c("PLA", "platelet-free")),
    n_reference_positive = sum(reference_status %in% "PLA"),
    TP = sum(Prediction == "PLA" & reference_status == "PLA", na.rm = TRUE),
    FP = sum(Prediction == "PLA" & reference_status == "platelet-free", na.rm = TRUE),
    FN = sum(Prediction == "platelet-free" & reference_status == "PLA", na.rm = TRUE),
    TN = sum(Prediction == "platelet-free" & reference_status == "platelet-free", na.rm = TRUE),
    .groups = "drop"
  ) %>% mutate(
    n_evaluated = TP + FP + FN + TN,
    precision = ratio(TP, TP + FP), recall = ratio(TP, TP + FN),
    F1 = ratio(2 * TP, 2 * TP + FP + FN),
    specificity = ratio(TN, TN + FP),
    classification_coverage = ratio(n_classified, n_retained),
    reference_coverage = ratio(n_evaluated, n_reference)
  )
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) < 2) stop("Usage: Benchmarking_QC_Grid.R dataset input.rds [variant] [grid.csv] [mode] [scope] [signature.csv] [run_id] [score_cache.rds]")
  arg <- function(i, default) if (length(args) >= i && nzchar(args[i])) args[i] else default
  script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
  root <- Sys.getenv("PLA_PROJECT_ROOT", unset = dirname(dirname(dirname(normalizePath(script)))))
  source(file.path(root, "src/gating_automation/config/datasets.R"), local = TRUE)
  source(file.path(root, "src/gating_automation/config/qc_parameters.R"), local = TRUE)
  dataset <- args[1]
  input <- normalizePath(args[2], mustWork = TRUE)
  variant <- arg(3, "adt_tuned")
  mode <- arg(5, "gmm_dist_dual")
  scope <- arg(6, "global")
  stopifnot(mode %in% c("gmm_dist_platelet", "gmm_dist_dual"), scope %in% c("global", "per_celltype"))
  signature <- arg(7, file.path(root, "data/signatures/GOBP_REGULATION_OF_PLATELET_ACTIVATION.v2025.1.Hs.csv"))
  signature_files <- c(platelet = signature, immune = file.path(root,
    "data/signatures/GOBP_LEUKOCYTE_ACTIVATION_INVOLVED_IN_INFLAMMATORY_RESPONSE.v2025.1.Hs.csv"))
  run_id <- arg(8, paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid()))
  qc_dir <- Sys.getenv("PLA_QC_DIR", file.path(root, "results/sample_qc", dataset, variant))
  exclude_na <- toupper(Sys.getenv("PLA_EXCLUDE_NA_METADATA", "FALSE")) == "TRUE"
  output_name <- if (exclude_na) "qc_grid_no_na" else "qc_grid"
  out <- file.path(Sys.getenv("PLA_GRID_OUTPUT_ROOT", file.path(root, "results/benchmarking", output_name)), dataset, run_id)
  if (dir.exists(out)) stop("Run directory already exists; choose a new run_id: ", out)
  config <- get_pla_dataset_config(dataset)

  grid <- if (nzchar(arg(4, ""))) read.csv(args[4]) else
    expand.grid(min_difference = c(0, 0.25, 0.5, 1, 1.5),
                max_overlap = c(0.025, 0.05, 0.10, 0.20, 1.0),
                min_cells_per_status = c(0, 3, 5, 10, 20))
  if (!"min_cells_per_status" %in% names(grid))
    grid$min_cells_per_status <- PLA_QC_PARAMETERS$min_cells_per_status
  
  grid <- grid %>% distinct(min_difference, max_overlap, min_cells_per_status) %>% mutate(baseline = FALSE)
  grid <- bind_rows(tibble(min_difference = -Inf, max_overlap = Inf,
                           min_cells_per_status = 0, baseline = TRUE), grid) %>%
    mutate(combination = sprintf("qc_%03d", row_number()), .before = 1)

  pairs <- read.csv(file.path(qc_dir, "01_sample_celltype_summary.csv"),
    colClasses = c(sample_id = "character", celltype_id = "character")) %>%
    transmute(sample_id = as.character(sample_id), lineage = canonical_lineage(celltype_id),
              usable_celltype_sample)
  sample_metrics <- read.csv(file.path(qc_dir, "12_platelet_marker_qc_by_sample.csv"),
    colClasses = c(group_id = "character", marker_label = "character")) %>%
    rename(sample_id = group_id)
  lineage_metrics <- read.csv(file.path(qc_dir, "13_platelet_marker_qc_by_lineage.csv"),
    colClasses = c(lineage = "character", marker_label = "character")) %>%
    mutate(lineage = canonical_lineage(lineage))
  if (anyDuplicated(pairs[c("sample_id", "lineage")])) stop("Duplicate sample-lineage pairs in frequency QC.")

  provenance <- list(input = input, input_size = file.info(input)$size,
    input_mtime = as.numeric(file.info(input)$mtime), dataset = dataset,
    sample_col = config$sample_col, signatures = tools::md5sum(signature_files), seed = 42L,
    scoring = "AUCell counts; fixed genes; aucMaxRank=ceiling(n_genes*0.05); no extension")
  cache <- arg(9, "")
  if (nzchar(cache)) {
    cached <- readRDS(cache)
    if (!identical(cached$provenance, provenance)) stop("Score cache does not match input/signatures.")
    scores <- cached$cells
  } else {
    scores <- prepare_scores(input, config$sample_col, signature_files)
  }
  # Older caches replaced missing lineages with "Unassigned". Recover exact originals.
  if (exclude_na && !all(c("lineage_original", "reference_status_original") %in% names(scores))) {
    meta <- readRDS(input)@meta.data
    lineage_col <- "lineage"
    reference_col <- intersect(c("pla_status", "pla.status"), names(meta))[1]
    scores$lineage_original <- as.character(meta[scores$cell_id, lineage_col])
    scores$reference_status_original <- as.character(meta[scores$cell_id, reference_col])
  }
  all_scores <- normalize_reference(scores)
  scores <- if (exclude_na) normalize_reference(all_scores, exclude_na = TRUE) else all_scores
  message("Metadata filter removed ", nrow(all_scores) - nrow(scores), "/", nrow(all_scores), " cells")
  reference_mapping <- if (exclude_na)
    "Exclude NA/unassigned lineage or pla_status; original labels retained" else
    "NA/unassigned -> platelet-free; original labels retained"
  base_cells <- semi_join(scores, filter(pairs, usable_celltype_sample), by = c("sample_id", "lineage"))
  if (!nrow(base_cells)) stop("No cells match frequency-QC pairs; check sample and lineage identifiers.")
  dir.create(out, recursive = TRUE)
  if (!nzchar(cache)) saveRDS(list(provenance = provenance, cells = all_scores), file.path(out, "scores.rds"))
    saveRDS(list(provenance = provenance, qc_dir = normalizePath(qc_dir),
      qc_md5 = tools::md5sum(file.path(qc_dir, c("01_sample_celltype_summary.csv",
        "12_platelet_marker_qc_by_sample.csv", "13_platelet_marker_qc_by_lineage.csv"))),
      mode = mode, scope = scope, qc_parameters = PLA_QC_PARAMETERS,
      reference_mapping = reference_mapping, exclude_na_metadata = exclude_na,
      n_cells_before_metadata_filter = nrow(all_scores),
      n_cells_removed_metadata_filter = nrow(all_scores) - nrow(scores),
      score_cache = cache, session = sessionInfo()), file.path(out, "run_config.rds"))
  write.csv(grid, file.path(out, "grid.csv"), row.names = FALSE)
  write.csv(data.frame(
    variant = variant, mode = mode, scope = scope, exclude_na_metadata = exclude_na,
    signature_files = paste(names(signature_files), signature_files, sep = "=", collapse = "; ")
  ), file.path(out, "run_config.csv"), row.names = FALSE)
  baseline_counts <- base_cells %>% dplyr::count(sample_id, lineage, name = "n_before_adt_qc")
  save_predictions <- toupper(Sys.getenv("PLA_SAVE_PREDICTIONS", "FALSE")) == "TRUE"
  summaries <- list()
  pair_results <- list()
  mean_defined <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

  for (i in seq_len(nrow(grid))) {
    message("Grid row ", i, "/", nrow(grid), ": min_difference=", grid$min_difference[i],
      ", max_overlap=", grid$max_overlap[i], ", min_cells_per_status=", grid$min_cells_per_status[i])
    setting <- grid[i, ]
    started <- Sys.time()
    decisions <- select_pairs(pairs, sample_metrics, lineage_metrics,
      setting$min_difference, setting$max_overlap, PLA_QC_PARAMETERS,
      if (setting$baseline) PLA_QC_PARAMETERS$min_cells_per_status else setting$min_cells_per_status)
    if (setting$baseline) decisions$retain_pair <- decisions$usable_celltype_sample
    cells <- semi_join(base_cells, filter(decisions, retain_pair), by = c("sample_id", "lineage"))
    cells$Prediction <- predict_gmm(cells, mode, scope)
    evaluated <- evaluate_pairs(cells)
    per_pair <- baseline_counts %>% left_join(decisions, by = c("sample_id", "lineage")) %>%
      left_join(evaluated, by = c("sample_id", "lineage")) %>%
      mutate(n_retained = coalesce(n_retained, 0L), dataset = dataset,
             run_id = run_id, combination = setting$combination,
             min_cells_per_status = setting$min_cells_per_status)
    pair_results[[i]] <- per_pair
    summaries[[i]] <- bind_cols(tibble(dataset = dataset, run_id = run_id), setting,
      tibble(n_cells_before_metadata_filter = nrow(all_scores),
        n_cells_removed_metadata_filter = nrow(all_scores) - nrow(scores),
        n_input_cells = nrow(scores), n_frequency_qc_cells = nrow(base_cells),
        n_retained_cells = nrow(cells), retained_cell_fraction = nrow(cells) / nrow(base_cells),
        n_input_pairs = nrow(baseline_counts), n_retained_pairs = nrow(evaluated),
        retained_pair_fraction = nrow(evaluated) / nrow(baseline_counts),
        n_unclassified = sum(is.na(cells$Prediction)),
        n_evaluated = sum(evaluated$n_evaluated),
        n_pairs_with_F1 = sum(!is.na(evaluated$F1)),
        macro_F1 = mean_defined(evaluated$F1), macro_recall = mean_defined(evaluated$recall),
        macro_precision = mean_defined(evaluated$precision),
        TP = sum(evaluated$TP), FP = sum(evaluated$FP),
        FN = sum(evaluated$FN), TN = sum(evaluated$TN),
        runtime_seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))))
    if (save_predictions) write.csv(cells %>% mutate(dataset = dataset, run_id = run_id,
        combination = setting$combination),
      file.path(out, paste0(setting$combination, "_predictions.csv")), row.names = FALSE)
    write.csv(bind_rows(summaries), file.path(out, "grid_summary.csv"), row.names = FALSE)
    write.csv(bind_rows(pair_results), file.path(out, "pair_metrics.csv"), row.names = FALSE)
    message(setting$combination, ": retained ", nrow(cells), "/", nrow(base_cells), " cells")
  }
  message("Results: ", out)
  invisible(out)
}

if (sys.nframe() == 0L) main()
