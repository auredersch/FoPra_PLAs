
args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args) && nzchar(args[1])) args[1] else
  "results/benchmarking/frequency_qc_full/by_dataset"
output_dir <- if (length(args) >= 2 && nzchar(args[2])) args[2] else
    "results/benchmarking/frequency_qc_full/combined"

datasets <- c("heart", "vaccine", "immune_aging", "sepsis", "skin")
dirs <- file.path(input_dir, datasets)
tables <- c("grid_summary", "pair_metrics", "frequency_pairs", "run_config")

read_table <- function(path) {
  columns <- names(read.csv(path, nrows = 0))
  classes <- if ("sample_id" %in% columns) c(sample_id = "character") else NA
  read.csv(path, colClasses = classes, stringsAsFactors = FALSE)
}
grids <- lapply(dirs, function(d) {
  g <- read_table(file.path(d, "grid.csv"))
  g <- g[order(g$max_ci_width, g$min_cells_per_pair), , drop = FALSE]
  rownames(g) <- NULL
  g
})
if (!all(vapply(grids, function(g) identical(g, grids[[1]]), logical(1))))
  stop("Die Datensaetze verwenden unterschiedliche Grids.")

merged <- setNames(lapply(tables, function(name)
  do.call(rbind, lapply(dirs, function(d) read_table(file.path(d, paste0(name, ".csv")))))), tables)
policies <- unique(merged$run_config$exclude_na_metadata)
if (length(policies) != 1L || anyNA(policies))
  stop("Varianten mit und ohne NA/Unassigned muessen getrennt zusammengefuehrt werden.")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
# Keep the original grid order, including the baseline as the first row.
write.csv(read_table(file.path(dirs[1], "grid.csv")), file.path(output_dir, "grid.csv"), row.names = FALSE)
for (name in tables) {
  write.csv(merged[[name]], file.path(output_dir, paste0(name, ".csv")), row.names = FALSE)
  message(name, ": ", nrow(merged[[name]]), " Zeilen")
}
message("Zusammengefuehrte Dateien: ", output_dir)

