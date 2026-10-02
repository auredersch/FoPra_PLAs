# From project root: Rscript src/benchmarking/search_qc_params/plot_prediction_f1.R [grid_dir]
library(dplyr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- file.path(if (length(args)) args[1] else "results/benchmarking/qc_grid_no_na_full",
                       "prediction_comparison")
output_dir <- file.path(input_dir, "plots")
dir.create(output_dir, showWarnings = FALSE)
files <- list.files(input_dir, "^prediction_comparison\\.csv$", recursive = TRUE, full.names = TRUE)
results <- bind_rows(lapply(files, read.csv)) %>%
  filter(is.finite(baseline_F1_matched), is.finite(qc_F1_matched)) %>%
  mutate(method = paste(gene_set, sub("gmm_dist_", "", mode), scope, sep = " | "),
         setting = ifelse(frequency_only, "Frequency QC only", "Frequency + ADT QC")) %>%
  distinct(dataset, method, setting, baseline_F1_matched, qc_F1_matched)

for (ds in unique(results$dataset)) {
  points <- filter(results, dataset == ds)
  limits <- c(0, min(1, max(c(points$baseline_F1_matched, points$qc_F1_matched, 0.1)) * 1.05))
  p <- ggplot(points, aes(baseline_F1_matched, qc_F1_matched, color = method)) +
    geom_abline(slope = 1, intercept = 0, color = "grey55", linetype = "dashed") +
    geom_point(aes(shape = setting), size = 2.5, alpha = 0.75) +
    coord_equal(xlim = limits, ylim = limits) +
    scale_color_brewer(palette = "Dark2", limits = sort(unique(results$method))) +
    scale_shape_manual(values = c("Frequency QC only" = 17, "Frequency + ADT QC" = 16)) +
    labs(title = tools::toTitleCase(gsub("_", " ", ds)),
      subtitle = "Cell-level F1 on the same retained, jointly classified cells",
      x = "F1 of no-QC predictions on retained cells", y = "F1 after QC and GMM refitting",
      color = NULL, shape = NULL,
      caption = "Above diagonal: improved F1. Below diagonal: reduced F1. Matched cells vary by QC setting.") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", legend.box = "vertical", panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"), plot.caption = element_text(hjust = 0, color = "grey40")) +
    guides(color = guide_legend(ncol = 2, order = 1), shape = guide_legend(order = 2))
  for (ext in c("png", "pdf"))
    ggsave(file.path(output_dir, paste0("prediction_f1_", ds, ".", ext)),
           p, width = 10, height = 10, dpi = 300, bg = "white")
}
