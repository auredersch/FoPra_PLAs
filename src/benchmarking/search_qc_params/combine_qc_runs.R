# Run from the project root; optional argument: results/benchmarking/qc_grid_no_na
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args)) args[1] else "results/benchmarking/qc_grid_no_na_full"
output_dir <- file.path(input_dir, "combined")
dir.create(output_dir, showWarnings = FALSE)
configs <- list.files(input_dir, "^run_config\\.rds$", recursive = TRUE, full.names = TRUE)
gene_sets <- c(
  "GOBP_REGULATION_OF_PLATELET_ACTIVATION.v2025.1.Hs.csv" = "GOBP",
  "MANNE_COVID19_COMBINED_COHORT_VS_HEALTHY_DONOR_PLATELETS_DN.v2025.1.Hs.csv" = "MANNE_DN"
)

runs <- lapply(configs, function(path) {
  config <- readRDS(path)
  provenance <- config$provenance
  summary <- read.csv(file.path(dirname(path), "grid_summary.csv"))
  signatures <- names(provenance$signatures)
  gene_set <- unname(gene_sets[basename(signatures[1])])
  if (is.na(gene_set)) gene_set <- tools::file_path_sans_ext(basename(signatures[1]))

  metadata <- tibble(
    dataset = provenance$dataset, run_id = summary$run_id[1],
    run_dir = dirname(path), variant = basename(config$qc_dir),
    gene_set = gene_set, mode = config$mode, scope = config$scope,
    exclude_na_metadata = isTRUE(config$exclude_na_metadata),
    platelet_signature = signatures[1], immune_signature = signatures[2],
    platelet_signature_md5 = unname(provenance$signatures[1]),
    immune_signature_md5 = unname(provenance$signatures[2]),
    input = provenance$input, input_size = provenance$input_size,
    input_mtime = provenance$input_mtime, seed = provenance$seed,
    qc_dir = config$qc_dir,
    frequency_qc_md5 = unname(config$qc_md5[1]),
    sample_qc_md5 = unname(config$qc_md5[2]),
    lineage_qc_md5 = unname(config$qc_md5[3]),
    reference_mapping = config$reference_mapping
  )

  summary <- summary %>%
    left_join(metadata, by = c("dataset", "run_id")) %>%
    mutate(
      # This key describes the values, independently of the run's qc_00x numbering.
      qc_setting = ifelse(baseline, "baseline",
        paste0("difference=", min_difference, ";overlap=", max_overlap,
               ";min_cells=", min_cells_per_status)),
      F1 = 2 * TP / (2 * TP + FP + FN),
      precision = TP / (TP + FP), recall = TP / (TP + FN)
    )

  settings <- summary %>% select(
    dataset, run_id, combination, gene_set, mode, scope, variant,
    qc_setting, baseline, min_difference, max_overlap, min_cells_per_status
  )
  pairs <- read.csv(file.path(dirname(path), "pair_metrics.csv"),
                    colClasses = c(sample_id = "character")) %>%
    select(-min_cells_per_status) %>%
    # combination is only used together with dataset/run_id to resolve local IDs.
    left_join(settings, by = c("dataset", "run_id", "combination"))

  list(config = metadata, summary = summary, pairs = pairs)
})

write.csv(bind_rows(lapply(runs, `[[`, "config")),
          file.path(output_dir, "run_config.csv"), row.names = FALSE)
write.csv(bind_rows(lapply(runs, `[[`, "summary")),
          file.path(output_dir, "grid_summary.csv"), row.names = FALSE)
write.csv(bind_rows(lapply(runs, `[[`, "pairs")),
          file.path(output_dir, "pair_metrics.csv"), row.names = FALSE)
message("Combined ", length(runs), " runs in ", output_dir)
