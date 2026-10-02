library(ggplot2)
library(dplyr)

input_dir <- "results/benchmarking/qc_grid_no_na"
output_dir <- file.path(input_dir, "combined", "plots")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
pairs <- read.csv(file.path(input_dir, "combined/pair_metrics.csv")) %>%
    filter(tolower(trimws(lineage)) == "unassigned", baseline) %>%
    group_by(dataset, gene_set, mode, scope) %>%
    summarise(
        FP = sum(FP, na.rm = TRUE),
        TN = sum(TN, na.rm = TRUE),
        FPR = FP / (FP + TN),
        .groups = "drop"
    ) %>%
    mutate(
        FPR = FP / (FP + TN),
        method = paste(gene_set, sub("gmm_dist_", "", mode), scope, sep = " | "),
        label = ifelse(
        is.finite(FPR),
        sprintf("%.1f%%\nFP=%s", 100 * FPR, format(FP, big.mark = ",", trim = TRUE)),
        "NA"
        )
    )



p <- ggplot(pairs, aes(dataset, method, fill = FPR)) +
    geom_tile(color="white") +
    geom_text(
        aes(
            label = label,
            color = ifelse(FPR > 0.5, "black", "white")
        ),  
        size = 3.5, 
    ) +
    scale_color_identity() +
    scale_fill_viridis_c(
        limits = c(0, 1), 
        labels = scales::label_percent(),
        na.value = "grey85"
    ) +
    labs(
        title = "False positive rate for unassigned lineages",
        x = NULL, y = NULL, fill = "FPR"
    ) +
    theme_minimal() +
    theme(
        panel.grid = element_blank()
    )

ggsave(file.path(output_dir, "unassigned_FP_heatmap.png"), p, width = 10, height = 6, dpi = 300)
ggsave(file.path(output_dir, "unassigned_FP_heatmap.pdf"), p, width = 10, height = 6)
