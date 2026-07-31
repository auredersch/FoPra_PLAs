#!/usr/bin/env Rscript

library(Seurat)
library(magrittr)

obj <- readRDS("/nfs/home/students/i.kaciran/FoPra_PLAs/data/datasets/gated_skin_processed.rds")

cat("\nStatus:\n")
print(table(obj$Status, useNA = "ifany"))

cat("\nSubject:\n")
print(table(obj$Subject, useNA = "ifany"))

print(
  obj@meta.data %>%
    dplyr::filter(!is.na(Subject)) %>%
    dplyr::distinct(Subject, pla_status) %>%
    dplyr::count(Subject, name = "n_conditions") %>%
    dplyr::count(n_conditions)
)