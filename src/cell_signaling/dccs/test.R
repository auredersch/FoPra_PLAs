#!/usr/bin/env Rscript

library(Seurat)
library(dplyr)

check_object <- function(path, label) {

  cat("\n============================\n")
  cat(label, "\n")
  cat("============================\n")

  obj <- readRDS(path)
  meta <- obj@meta.data

  print(colnames(meta))

  candidate_cols <- intersect(
    c("sample", "sample_ID", "patient", "donor_id", "orig.ident"),
    colnames(meta)
  )

  print(
    head(
      meta[, intersect(
        c(candidate_cols, "pla_status", "lineage"),
        colnames(meta)
      ), drop = FALSE]
    )
  )

  for (x in candidate_cols) {

    cat("\n====================\n")
    cat("COLUMN:", x, "\n")
    cat("====================\n")

    summary_tbl <- meta %>%
      distinct(
        id = .data[[x]],
        pla_status
      ) %>%
      count(id, name = "n_conditions")

    cat("Total IDs:", nrow(summary_tbl), "\n")

    print(count(summary_tbl, n_conditions))
  }
}

check_object(
  "/nfs/home/students/i.kaciran/FoPra_PLAs/data/datasets/gated_our_dataset_processed.rds",
  "OUR DATASET"
)

check_object(
  "/nfs/home/students/i.kaciran/FoPra_PLAs/data/datasets/gated_ImmuneAging.rds",
  "IMMUNE AGING"
)