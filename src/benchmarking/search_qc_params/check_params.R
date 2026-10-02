library(dplyr)

gene_sets <- c("MANNE_DN", "GOBP")
datasets <- c("heart", "immune_aging", "sepsis", "skin", "vaccine")

params <- expand.grid(
    dataset = datasets,
    gene_set = gene_sets
)

for (i in seq_len(nrow(params))){
    setting <- params[i, ]
    name <- setting[["dataset"]]
    gene_list <- setting[["gene_set"]]
    message("Processing dataset: ", name, " with gene set: ", gene_list)
    summary <- read.csv(file.path("results/benchmarking/qc_grid", name, gene_list, "grid_summary.csv"))
    grid <- read.csv(file.path("results/benchmarking/qc_grid", name, gene_list, "grid.csv"), comment.char = "#")

    summary <- summary %>% 
        #select(combination, min_difference, max_overlap, min_cells_per_status, n_input_cells, n_retained_cells, retained_cell_fraction, TP,FP,FN,TN) %>%
        select(combination, n_input_cells, n_frequency_qc_cells, n_retained_cells, retained_cell_fraction,TP,FP,FN,TN) %>%
        mutate(F1 = 2*TP/(2*TP + FP + FN)) %>%
        arrange(desc(F1)) %>% 
        select(combination, n_input_cells, n_frequency_qc_cells, n_retained_cells, retained_cell_fraction, F1)


    grid_stats <- summary %>% left_join(grid, by = "combination")
    write.csv(grid_stats, file.path("results/benchmarking/qc_grid", name, gene_list, "grid_stats.csv"), row.names = FALSE)
}
