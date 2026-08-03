#!/usr/bin/env Rscript

options(
  error = function() {
    message("R error occurred. Traceback:")
    traceback(2)
    quit(save = "no", status = 1)
  }
)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(forcats)
  library(circlize)
  library(scDiffCom)
  library(future)
  library(grid)
  library(RColorBrewer)
})

# ============================================================
# 1. arguments
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript scDiffCom.R <input_rds> <output_dir>")
}

options(future.globals.maxSize = 2 * 1024^3)
future::plan(future::sequential)

input_file <- args[[1]]
base_output_dir <- args[[2]]

dataset_name <- tools::file_path_sans_ext(basename(input_file))

dataset_mode <- stringr::str_extract(
  dataset_name,
  "(all|diseasedOnly|healthyOnly)$"
)

dataset_clean <- stringr::str_remove(
  dataset_name,
  "_(all|diseasedOnly|healthyOnly)$"
)

if (is.na(dataset_mode)) {
  dataset_mode <- "unknownMode"
}

plot_title <- function(title) {
  paste0(title, "\n", dataset_clean, " | ", dataset_mode)
}

out_dir <- file.path(base_output_dir, dataset_name)
plot_dir <- file.path(out_dir, "plots")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

message("Analyzing: ", dataset_name)
message("Input: ", input_file)
message("Output directory: ", out_dir)
message("Plot directory: ", plot_dir)

save_plot <- function(plot, filename, width = 8, height = 6, dpi = 300) {
  out_file <- file.path(plot_dir, filename)

  ggplot2::ggsave(
    filename = out_file,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )

  message("Saved plot: ", out_file)
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]

  if (length(x) == 0) {
    return(NA_real_)
  }

  mean(x)
}

symmetric_limit <- function(x, fallback = 1) {
  value <- suppressWarnings(
    max(abs(x), na.rm = TRUE)
  )

  if (!is.finite(value) || value <= 0) {
    return(fallback)
  }

  value
}

# ============================================================
# 2. load and prepare data
# ============================================================

seurat_obj <- readRDS(input_file)

condition_col <- "pla_status"
lineage_col <- "lineage"

required_cols <- c(condition_col, lineage_col)

missing_cols <- setdiff(
  required_cols,
  colnames(seurat_obj@meta.data)
)

if (length(missing_cols) > 0) {
  stop(
    "Missing metadata columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

keep <- complete.cases(
  seurat_obj@meta.data[, required_cols, drop = FALSE]
)

seurat_obj <- subset(
  seurat_obj,
  cells = colnames(seurat_obj)[keep]
)

DefaultAssay(seurat_obj) <- "RNA"

required_conditions <- c("PLA", "platelet-free")

available_conditions <- unique(
  as.character(
    seurat_obj@meta.data[[condition_col]]
  )
)

missing_conditions <- setdiff(
  required_conditions,
  available_conditions
)

message(
  "Available pla_status values: ",
  paste(sort(available_conditions), collapse = ", ")
)

if (length(missing_conditions) > 0) {
  stop(
    "Cannot perform PLA versus platelet-free comparison. Missing condition(s): ",
    paste(missing_conditions, collapse = ", ")
  )
}

message("Cell counts by condition:")
print(table(seurat_obj@meta.data[[condition_col]], useNA = "ifany"))

message("Cell counts by condition and lineage:")
print(
  table(
    seurat_obj@meta.data[[condition_col]],
    seurat_obj@meta.data[[lineage_col]]
  )
)

# ============================================================
# 3. shared plotting configuration
# ============================================================

direction_colors <- c(
  "PLA-up" = "#F8766D",
  "platelet-free-up" = "#00BFC4"
)

direction_levels <- c(
  "PLA-up",
  "platelet-free-up"
)

all_dataset_lineages <- sort(
  unique(
    as.character(
      seurat_obj@meta.data[[lineage_col]]
    )
  )
)

all_dataset_lineages <- all_dataset_lineages[
  !is.na(all_dataset_lineages) &
    all_dataset_lineages != ""
]

if (length(all_dataset_lineages) == 0) {
  stop("No valid lineage values found.")
}

base_palette <- RColorBrewer::brewer.pal(
  n = 11,
  name = "Spectral"
)

global_lineage_colors <- grDevices::colorRampPalette(
  base_palette
)(length(all_dataset_lineages))

names(global_lineage_colors) <- all_dataset_lineages

lineage_color_table <- tibble::tibble(
  lineage = names(global_lineage_colors),
  color = unname(global_lineage_colors)
)

common_theme <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1
      ),
      panel.grid.minor = ggplot2::element_blank()
    )
}

# ============================================================
# 4. run scDiffCom
# ============================================================

scdiffcom_object <- scDiffCom::run_interaction_analysis(
  seurat_object = seurat_obj,
  LRI_species = "human",
  seurat_celltype_id = lineage_col,
  seurat_condition_id = list(
    column_name = condition_col,
    cond1_name = "platelet-free",
    cond2_name = "PLA"
  )
)

# ============================================================
# 5. extract results
# ============================================================

CCI_detected <- scDiffCom::GetTableCCI(
  scdiffcom_object,
  type = "detected",
  simplified = TRUE
)

ORA_results <- scDiffCom::GetTableORA(
  scdiffcom_object,
  categories = "all",
  simplified = TRUE
)

cci_regulation_counts <- CCI_detected %>%
  dplyr::count(REGULATION)

# ============================================================
# 6. ORA plots
# ============================================================

p_scdiff_ORA_UP <- tryCatch(
  {
    scDiffCom::PlotORA(
      object = scdiffcom_object,
      category = "LRI",
      regulation = "UP"
    ) +
      ggplot2::theme(
        legend.position = c(0.85, 0.4),
        legend.key.size = grid::unit(0.4, "cm")
      ) +
      ggplot2::labs(
        title = plot_title(
          "scDiffCom ORA: up-regulated LRIs"
        )
      )
  },
  error = function(e) {
    message("PlotORA UP failed: ", conditionMessage(e))

    ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::annotate(
        "text",
        x = 0,
        y = 0,
        label = "No UP ORA results available",
        size = 5
      ) +
      ggplot2::labs(
        title = plot_title(
          "scDiffCom ORA: up-regulated LRIs"
        )
      )
  }
)

save_plot(
  p_scdiff_ORA_UP,
  "02_scdiffcom_ORA_LRI_UP.png",
  width = 9,
  height = 7
)

p_scdiff_ORA_DOWN <- tryCatch(
  {
    scDiffCom::PlotORA(
      object = scdiffcom_object,
      category = "LRI",
      regulation = "DOWN"
    ) +
      ggplot2::theme(
        legend.position = c(0.85, 0.4),
        legend.key.size = grid::unit(0.4, "cm")
      ) +
      ggplot2::labs(
        title = plot_title(
          "scDiffCom ORA: down-regulated LRIs"
        )
      )
  },
  error = function(e) {
    message("PlotORA DOWN failed: ", conditionMessage(e))

    ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::annotate(
        "text",
        x = 0,
        y = 0,
        label = "No DOWN ORA results available",
        size = 5
      ) +
      ggplot2::labs(
        title = plot_title(
          "scDiffCom ORA: down-regulated LRIs"
        )
      )
  }
)

save_plot(
  p_scdiff_ORA_DOWN,
  "03_scdiffcom_ORA_LRI_DOWN.png",
  width = 9,
  height = 7
)

# ============================================================
# 7. robust source-target parsing
# ============================================================

parse_lineage_pair <- function(x, valid_lineages) {
  valid_lineages <- valid_lineages[
    order(nchar(valid_lineages), decreasing = TRUE)
  ]

  parsed <- lapply(
    x,
    function(value) {
      if (is.na(value) || value == "") {
        return(
          tibble::tibble(
            source = NA_character_,
            target = NA_character_
          )
        )
      }

      possible_sources <- valid_lineages[
        vapply(
          valid_lineages,
          function(lineage) {
            startsWith(
              value,
              paste0(lineage, "_")
            )
          },
          logical(1)
        )
      ]

      if (length(possible_sources) == 0) {
        return(
          tibble::tibble(
            source = NA_character_,
            target = NA_character_
          )
        )
      }

      source <- possible_sources[[1]]

      target <- substring(
        value,
        nchar(source) + 2
      )

      if (!target %in% valid_lineages) {
        return(
          tibble::tibble(
            source = NA_character_,
            target = NA_character_
          )
        )
      }

      tibble::tibble(
        source = source,
        target = target
      )
    }
  )

  dplyr::bind_rows(parsed)
}

parsed_lineages <- parse_lineage_pair(
  CCI_detected$ER_CELLTYPES,
  all_dataset_lineages
)

scdiff_plot <- CCI_detected %>%
  tibble::as_tibble() %>%
  dplyr::bind_cols(parsed_lineages) %>%
  tidyr::separate(
    LRI,
    into = c("ligand", "receptor"),
    sep = ":",
    remove = FALSE,
    extra = "merge",
    fill = "right"
  ) %>%
  dplyr::mutate(
    lineage_pair = paste(
      source,
      target,
      sep = " → "
    ),
    interaction = paste(
      ligand,
      receptor,
      sep = " → "
    ),
    direction = dplyr::case_when(
      REGULATION == "UP" ~ "PLA-up",
      REGULATION == "DOWN" ~ "platelet-free-up",
      REGULATION == "FLAT" ~ "flat",
      REGULATION == "NSC" ~ "not significant",
      TRUE ~ as.character(REGULATION)
    ),
    direction = factor(
      direction,
      levels = c(
        direction_levels,
        "flat",
        "not significant"
      )
    ),
    score_signed = LOGFC,
    score_abs = abs(LOGFC),
    padj = BH_P_VALUE_DE
  )

n_unparsed <- sum(
  is.na(scdiff_plot$source) |
    is.na(scdiff_plot$target)
)

if (n_unparsed > 0) {
  warning(
    n_unparsed,
    " scDiffCom rows could not be parsed into source and target lineages."
  )
}

scdiff_sig <- scdiff_plot %>%
  dplyr::filter(
    direction %in% direction_levels,
    is.finite(score_signed),
    is.finite(score_abs),
    is.finite(padj)
  )

# ============================================================
# 8. volcano plot
# ============================================================

p_scdiff_volcano <- ggplot2::ggplot(
  scdiff_plot,
  ggplot2::aes(
    x = score_signed,
    y = -log10(pmax(padj, 1e-300)),
    color = direction
  )
) +
  ggplot2::geom_point(
    alpha = 0.6,
    size = 1.3
  ) +
  ggplot2::scale_color_manual(
    values = c(
      direction_colors,
      "flat" = "#7CAE00",
      "not significant" = "grey70"
    ),
    drop = FALSE
  ) +
  common_theme() +
  ggplot2::labs(
    title = plot_title(
      "scDiffCom differential cell-cell interactions"
    ),
    x = "log fold change: PLA / platelet-free",
    y = expression(-log[10]("adjusted p-value")),
    color = "Regulation"
  ) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 0)
  )

save_plot(
  p_scdiff_volcano,
  "01_scdiffcom_volcano.png",
  width = 8,
  height = 6
)

# ============================================================
# 9. global significant interaction counts
# ============================================================

global_counts_scdiff <- scdiff_sig %>%
  dplyr::count(
    direction,
    name = "n_significant_ccis"
  )

p_scdiff_global <- ggplot2::ggplot(
  global_counts_scdiff,
  ggplot2::aes(
    x = direction,
    y = n_significant_ccis,
    fill = direction
  )
) +
  ggplot2::geom_col() +
  ggplot2::geom_text(
    ggplot2::aes(label = n_significant_ccis),
    vjust = -0.4,
    size = 4
  ) +
  ggplot2::scale_fill_manual(
    values = direction_colors,
    drop = FALSE
  ) +
  common_theme() +
  ggplot2::labs(
    title = plot_title(
      "Global number of significant scDiffCom CCIs"
    ),
    x = "Differential direction",
    y = "# significant CCIs",
    fill = "Direction"
  ) +
  ggplot2::theme(
    legend.position = "none"
  )

if (nrow(global_counts_scdiff) > 0) {
  p_scdiff_global <- p_scdiff_global +
    ggplot2::expand_limits(
      y = max(
        global_counts_scdiff$n_significant_ccis
      ) * 1.08
    )
}

save_plot(
  p_scdiff_global,
  "04_scdiffcom_global_significant_counts.png",
  width = 7,
  height = 5
)

# ============================================================
# 10. effect-size distribution
# ============================================================

p_scdiff_score_distribution <- ggplot2::ggplot(
  scdiff_sig,
  ggplot2::aes(
    x = direction,
    y = score_abs,
    fill = direction
  )
) +
  ggplot2::geom_violin(
    trim = TRUE,
    alpha = 0.45
  ) +
  ggplot2::geom_boxplot(
    width = 0.18,
    outlier.size = 0.3,
    alpha = 0.8
  ) +
  ggplot2::scale_fill_manual(
    values = direction_colors,
    drop = FALSE
  ) +
  common_theme() +
  ggplot2::labs(
    title = plot_title(
      "scDiffCom differential effect-size distribution"
    ),
    x = "Differential direction",
    y = "|log fold change|",
    fill = "Direction"
  ) +
  ggplot2::theme(
    legend.position = "none"
  )

save_plot(
  p_scdiff_score_distribution,
  "05_scdiffcom_absolute_logfc_distribution.png",
  width = 7,
  height = 5
)

# ============================================================
# 11. significant CCIs by lineage pair
# ============================================================

pair_counts_scdiff <- scdiff_sig %>%
  dplyr::filter(
    !is.na(lineage_pair)
  ) %>%
  dplyr::count(
    lineage_pair,
    direction,
    name = "n_interactions"
  ) %>%
  dplyr::group_by(lineage_pair) %>%
  dplyr::mutate(
    total = sum(n_interactions)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::slice_max(
    order_by = total,
    n = 25,
    with_ties = FALSE
  ) %>%
  dplyr::mutate(
    lineage_pair = forcats::fct_reorder(
      lineage_pair,
      total
    )
  )

p_pair_counts_scdiff <- ggplot2::ggplot(
  pair_counts_scdiff,
  ggplot2::aes(
    x = n_interactions,
    y = lineage_pair,
    fill = direction
  )
) +
  ggplot2::geom_col() +
  ggplot2::scale_fill_manual(
    values = direction_colors,
    drop = FALSE
  ) +
  common_theme() +
  ggplot2::labs(
    title = plot_title(
      "Significant scDiffCom CCIs by lineage pair"
    ),
    x = "# significant differential CCIs",
    y = "Source → target lineage",
    fill = "Direction"
  ) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 0)
  )

save_plot(
  p_pair_counts_scdiff,
  "06_scdiffcom_significant_interactions_by_lineage_pair.png",
  width = 10,
  height = 8
)

# ============================================================
# 12. net significant CCI-count heatmap
# ============================================================

net_count_scdiff <- scdiff_sig %>%
  dplyr::filter(
    !is.na(source),
    !is.na(target)
  ) %>%
  dplyr::count(
    source,
    target,
    direction,
    name = "n"
  ) %>%
  tidyr::complete(
    source,
    target,
    direction = factor(
      direction_levels,
      levels = direction_levels
    ),
    fill = list(n = 0)
  ) %>%
  tidyr::pivot_wider(
    names_from = direction,
    values_from = n,
    values_fill = 0
  ) %>%
  dplyr::mutate(
    net_count =
      `PLA-up` -
      `platelet-free-up`
  )

net_count_limit <- symmetric_limit(
  net_count_scdiff$net_count
)

p_net_count_scdiff <- ggplot2::ggplot(
  net_count_scdiff,
  ggplot2::aes(
    x = target,
    y = source,
    fill = net_count
  )
) +
  ggplot2::geom_tile(color = "white") +
  ggplot2::geom_text(
    ggplot2::aes(label = net_count),
    size = 3.5
  ) +
  ggplot2::scale_fill_gradient2(
    low = direction_colors[["platelet-free-up"]],
    mid = "white",
    high = direction_colors[["PLA-up"]],
    midpoint = 0,
    limits = c(
      -net_count_limit,
      net_count_limit
    )
  ) +
  common_theme() +
  ggplot2::labs(
    title = plot_title(
      "Net number of significant scDiffCom CCIs"
    ),
    x = "Target lineage",
    y = "Source lineage",
    fill = paste0(
      "PLA-up count\nminus platelet-",
      "free-up count"
    )
  )

save_plot(
  p_net_count_scdiff,
  "07_scdiffcom_net_significant_count_heatmap.png",
  width = 9,
  height = 8
)

# ============================================================
# 13. mean signed logFC heatmap
# ============================================================

net_logfc_scdiff <- scdiff_sig %>%
  dplyr::filter(
    !is.na(source),
    !is.na(target)
  ) %>%
  dplyr::group_by(
    source,
    target
  ) %>%
  dplyr::summarise(
    mean_signed_logfc = safe_mean(score_signed),
    median_signed_logfc = stats::median(
      score_signed,
      na.rm = TRUE
    ),
    n_significant_ccis = dplyr::n(),
    .groups = "drop"
  )

net_logfc_limit <- symmetric_limit(
  net_logfc_scdiff$mean_signed_logfc
)

p_net_logfc_scdiff <- ggplot2::ggplot(
  net_logfc_scdiff,
  ggplot2::aes(
    x = target,
    y = source,
    fill = mean_signed_logfc
  )
) +
  ggplot2::geom_tile(color = "white") +
  ggplot2::geom_text(
    ggplot2::aes(
      label = round(
        mean_signed_logfc,
        2
      )
    ),
    size = 3.5
  ) +
  ggplot2::scale_fill_gradient2(
    low = direction_colors[["platelet-free-up"]],
    mid = "white",
    high = direction_colors[["PLA-up"]],
    midpoint = 0,
    limits = c(
      -net_logfc_limit,
      net_logfc_limit
    )
  ) +
  common_theme() +
  ggplot2::labs(
    title = plot_title(
      "Mean signed logFC among significant scDiffCom CCIs"
    ),
    x = "Target lineage",
    y = "Source lineage",
    fill = "Mean signed\nlogFC"
  )

save_plot(
  p_net_logfc_scdiff,
  "08_scdiffcom_mean_signed_logfc_heatmap.png",
  width = 9,
  height = 8
)

# ============================================================
# 14. recurrent differential LR pairs
# ============================================================

top_recurrent_scdiff <- scdiff_sig %>%
  dplyr::filter(
    !is.na(interaction),
    !is.na(lineage_pair)
  ) %>%
  dplyr::distinct(
    interaction,
    lineage_pair,
    direction
  ) %>%
  dplyr::count(
    interaction,
    direction,
    name = "n_lineage_pairs"
  ) %>%
  dplyr::group_by(direction) %>%
  dplyr::slice_max(
    order_by = n_lineage_pairs,
    n = 15,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    interaction_facet = paste(
      interaction,
      direction,
      sep = "___"
    ),
    interaction_facet = forcats::fct_reorder(
      interaction_facet,
      n_lineage_pairs
    )
  )

p_top_recurrent_scdiff <- ggplot2::ggplot(
  top_recurrent_scdiff,
  ggplot2::aes(
    x = n_lineage_pairs,
    y = interaction_facet,
    fill = direction
  )
) +
  ggplot2::geom_col() +
  ggplot2::facet_wrap(
    ~ direction,
    scales = "free_y"
  ) +
  ggplot2::scale_y_discrete(
    labels = function(x) {
      sub("___.*$", "", x)
    }
  ) +
  ggplot2::scale_fill_manual(
    values = direction_colors,
    drop = FALSE
  ) +
  common_theme() +
  ggplot2::labs(
    title = plot_title(
      "Most recurrent differential scDiffCom LR pairs"
    ),
    x = "# source-target lineage pairs",
    y = "Ligand → receptor",
    fill = "Direction"
  ) +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(angle = 0)
  )

save_plot(
  p_top_recurrent_scdiff,
  "09_scdiffcom_top_recurrent_lr_pairs.png",
  width = 12,
  height = 9
)

# ============================================================
# 15. top LR pairs by mean absolute logFC
# ============================================================

top_lr_by_score_scdiff <- scdiff_sig %>%
  dplyr::filter(
    !is.na(interaction),
    !is.na(lineage_pair)
  ) %>%
  dplyr::group_by(
    interaction,
    direction
  ) %>%
  dplyr::summarise(
    n_lineage_pairs = dplyr::n_distinct(lineage_pair),
    mean_abs_logfc = safe_mean(score_abs),
    max_abs_logfc = max(
      score_abs,
      na.rm = TRUE
    ),
    min_padj = min(
      padj,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  dplyr::group_by(direction) %>%
  dplyr::slice_max(
    order_by = mean_abs_logfc,
    n = 15,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    interaction_facet = paste(
      interaction,
      direction,
      sep = "___"
    ),
    interaction_facet = forcats::fct_reorder(
      interaction_facet,
      mean_abs_logfc
    )
  )

p_top_lr_by_score_scdiff <- ggplot2::ggplot(
  top_lr_by_score_scdiff,
  ggplot2::aes(
    x = mean_abs_logfc,
    y = interaction_facet,
    fill = direction
  )
) +
  ggplot2::geom_col() +
  ggplot2::facet_wrap(
    ~ direction,
    scales = "free_y"
  ) +
  ggplot2::scale_y_discrete(
    labels = function(x) {
      sub("___.*$", "", x)
    }
  ) +
  ggplot2::scale_fill_manual(
    values = direction_colors,
    drop = FALSE
  ) +
  common_theme() +
  ggplot2::labs(
    title = plot_title(
      "Top scDiffCom LR pairs by mean absolute logFC"
    ),
    x = "Mean |logFC|",
    y = "Ligand → receptor",
    fill = "Direction"
  ) +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(angle = 0)
  )

save_plot(
  p_top_lr_by_score_scdiff,
  "10_scdiffcom_top_lr_by_mean_abs_logfc_faceted.png",
  width = 12,
  height = 9
)

# ============================================================
# 16. top interactions per lineage pair
# ============================================================

top_scdiff_per_pair <- scdiff_sig %>%
  dplyr::filter(
    !is.na(lineage_pair),
    !is.na(interaction)
  ) %>%
  dplyr::group_by(
    direction,
    lineage_pair
  ) %>%
  dplyr::slice_max(
    order_by = score_abs,
    n = 5,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    interaction_facet = paste(
      interaction,
      lineage_pair,
      direction,
      sep = "___"
    ),
    interaction_facet = forcats::fct_reorder(
      interaction_facet,
      score_abs
    )
  )

p_top_scdiff_per_pair <- ggplot2::ggplot(
  top_scdiff_per_pair,
  ggplot2::aes(
    x = score_abs,
    y = interaction_facet,
    color = direction,
    size = -log10(pmax(padj, 1e-300))
  )
) +
  ggplot2::geom_point(alpha = 0.85) +
  ggplot2::facet_wrap(
    ggplot2::vars(
      direction,
      lineage_pair
    ),
    scales = "free_y",
    ncol = 4
  ) +
  ggplot2::scale_y_discrete(
    labels = function(x) {
      sub("___.*$", "", x)
    }
  ) +
  ggplot2::scale_color_manual(
    values = direction_colors,
    drop = FALSE
  ) +
  common_theme(base_size = 9) +
  ggplot2::labs(
    title = plot_title(
      "Top scDiffCom interactions per lineage pair"
    ),
    x = "|logFC|",
    y = "Ligand → receptor",
    color = "Direction",
    size = expression(-log[10]("adjusted p-value"))
  ) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 0),
    legend.position = "bottom"
  )

save_plot(
  p_top_scdiff_per_pair,
  "11_scdiffcom_top_interactions_per_lineage_pair.png",
  width = 16,
  height = 14
)

# ============================================================
# 17. strongest individual differential interactions
# ============================================================

top_individual_scdiff <- scdiff_sig %>%
  dplyr::filter(
    !is.na(interaction),
    !is.na(lineage_pair)
  ) %>%
  dplyr::group_by(direction) %>%
  dplyr::slice_max(
    order_by = score_abs,
    n = 20,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    label = paste0(
      interaction,
      "\n",
      lineage_pair
    ),
    label_facet = paste(
      label,
      direction,
      sep = "___"
    ),
    label_facet = forcats::fct_reorder(
      label_facet,
      score_abs
    )
  )

p_top_individual_scdiff <- ggplot2::ggplot(
  top_individual_scdiff,
  ggplot2::aes(
    x = score_abs,
    y = label_facet,
    color = direction
  )
) +
  ggplot2::geom_point(
    size = 3,
    alpha = 0.85
  ) +
  ggplot2::facet_wrap(
    ~ direction,
    scales = "free_y"
  ) +
  ggplot2::scale_y_discrete(
    labels = function(x) {
      sub("___.*$", "", x)
    }
  ) +
  ggplot2::scale_color_manual(
    values = direction_colors,
    drop = FALSE
  ) +
  common_theme() +
  ggplot2::labs(
    title = plot_title(
      "Top individual scDiffCom interactions"
    ),
    x = "|logFC|",
    y = "Ligand → receptor / source → target",
    color = "Direction"
  ) +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(angle = 0)
  )

save_plot(
  p_top_individual_scdiff,
  "12_scdiffcom_top_individual_interactions.png",
  width = 13,
  height = 10
)

# ============================================================
# 18. gene-level circos plots
# ============================================================

top_n_circos <- 20L

scdiff_circos_tbl <- scdiff_sig %>%
  dplyr::filter(
    is.finite(score_abs),
    !is.na(source),
    !is.na(target),
    !is.na(ligand),
    !is.na(receptor),
    ligand != "",
    receptor != ""
  ) %>%
  dplyr::group_by(direction) %>%
  dplyr::slice_max(
    order_by = score_abs,
    n = top_n_circos,
    with_ties = FALSE
  ) %>%
  dplyr::arrange(
    dplyr::desc(score_abs),
    .by_group = TRUE
  ) %>%
  dplyr::mutate(
    rank_within_direction = dplyr::row_number(),
    ligand_node = paste(
      "L",
      source,
      ligand,
      sep = "::"
    ),
    receptor_node = paste(
      "R",
      target,
      receptor,
      sep = "::"
    )
  ) %>%
  dplyr::ungroup()


plot_scdiff_gene_circos <- function(
    data,
    direction_oi,
    output_file,
    title_text,
    lineage_colors
) {
  plot_data <- data %>%
    dplyr::filter(
      direction == direction_oi
    )

  if (nrow(plot_data) == 0) {
    message(
      "Skipping ",
      direction_oi,
      " circos plot: no interactions."
    )
    return(invisible(NULL))
  }

  links <- plot_data %>%
    dplyr::group_by(
      ligand_node,
      receptor_node,
      source,
      target
    ) %>%
    dplyr::summarise(
      value = max(
        score_abs,
        na.rm = TRUE
      ),
      n_rows = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::transmute(
      from = ligand_node,
      to = receptor_node,
      value = value,
      source = source,
      target = target,
      n_rows = n_rows
    )

  ligand_annotation <- plot_data %>%
    dplyr::distinct(
      node = ligand_node,
      lineage = source
    )

  receptor_annotation <- plot_data %>%
    dplyr::distinct(
      node = receptor_node,
      lineage = target
    )

  node_annotation <- dplyr::bind_rows(
    ligand_annotation,
    receptor_annotation
  ) %>%
    dplyr::distinct(
      node,
      lineage
    )

  sector_colors <- lineage_colors[
    node_annotation$lineage
  ]

  names(sector_colors) <- node_annotation$node

  ligand_order <- plot_data %>%
    dplyr::arrange(
      source,
      dplyr::desc(score_abs),
      ligand
    ) %>%
    dplyr::pull(ligand_node) %>%
    unique()

  receptor_order <- plot_data %>%
    dplyr::arrange(
      target,
      dplyr::desc(score_abs),
      receptor
    ) %>%
    dplyr::pull(receptor_node) %>%
    unique()

  sector_order <- c(
    ligand_order,
    receptor_order
  )

  n_ligand <- length(ligand_order)
  n_receptor <- length(receptor_order)

  gap_after <- c(
    if (n_ligand > 1) {
      rep(1.5, n_ligand - 1)
    } else {
      numeric(0)
    },
    8,
    if (n_receptor > 1) {
      rep(1.5, n_receptor - 1)
    } else {
      numeric(0)
    },
    8
  )

  grDevices::png(
    filename = output_file,
    width = 5000,
    height = 5000,
    res = 300
  )

  graphics::par(
    mar = c(1, 1, 4, 1),
    xpd = NA
  )

  circlize::circos.clear()

  circlize::circos.par(
    start.degree = 90,
    canvas.xlim = c(-1.15, 1.15),
    canvas.ylim = c(-1.15, 1.15),
    gap.after = gap_after,
    track.margin = c(0.005, 0.005),
    points.overflow.warning = FALSE
  )

  circlize::chordDiagram(
    x = links %>%
      dplyr::select(
        from,
        to,
        value
      ),
    order = sector_order,
    grid.col = sector_colors,
    col = sector_colors[links$from],
    transparency = 0.15,
    directional = 1,
    direction.type = c(
      "arrows",
      "diffHeight"
    ),
    diffHeight = -0.03,
    link.arr.type = "big.arrow",
    link.sort = TRUE,
    link.decreasing = TRUE,
    link.border = NA,
    annotationTrack = "grid",
    preAllocateTracks = list(
      track.height = 0.16
    )
  )

  circlize::circos.trackPlotRegion(
    track.index = 1,
    bg.border = NA,
    panel.fun = function(x, y) {
      sector <- circlize::get.cell.meta.data(
        "sector.index"
      )

      xlim <- circlize::get.cell.meta.data(
        "xlim"
      )

      ylim <- circlize::get.cell.meta.data(
        "ylim"
      )

      gene_label <- sub(
        "^[LR]::[^:]+::",
        "",
        sector
      )

      circlize::circos.text(
        x = mean(xlim),
        y = mean(ylim),
        labels = gene_label,
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0, 0.5),
        cex = 0.8
      )
    }
  )

  graphics::title(
    main = title_text,
    cex.main = 1.3,
    line = 1.5
  )

  circlize::circos.clear()
  grDevices::dev.off()
}

save_scdiff_circos_legend <- function(
    data,
    direction_oi,
    output_file,
    lineage_colors
) {
  legend_data <- data %>%
    dplyr::filter(
      direction == direction_oi
    )

  if (nrow(legend_data) == 0) {
    return(invisible(NULL))
  }

  lineages <- sort(
    unique(
      c(
        legend_data$source,
        legend_data$target
      )
    )
  )

  grDevices::png(
    filename = output_file,
    width = 1800,
    height = max(
      900,
      180 + 140 * length(lineages)
    ),
    res = 300
  )

  graphics::par(
    mar = c(1, 1, 1, 1)
  )

  graphics::plot.new()

  graphics::legend(
    "center",
    legend = lineages,
    fill = lineage_colors[lineages],
    border = NA,
    title = paste0(
      "Sender and receiver lineages\n",
      direction_oi
    ),
    cex = 1.2,
    bty = "n",
    ncol = 1
  )

  grDevices::dev.off()
}

plot_scdiff_gene_circos(
  data = scdiff_circos_tbl,
  direction_oi = "PLA-up",
  output_file = file.path(
    plot_dir,
    "13a_scdiffcom_circos_PLA_up.png"
  ),
  title_text = plot_title(
    "scDiffCom gene-level interactions: PLA-up"
  ),
  lineage_colors = global_lineage_colors
)

plot_scdiff_gene_circos(
  data = scdiff_circos_tbl,
  direction_oi = "platelet-free-up",
  output_file = file.path(
    plot_dir,
    "13b_scdiffcom_circos_platelet_free_up.png"
  ),
  title_text = plot_title(
    "scDiffCom gene-level interactions: platelet-free-up"
  ),
  lineage_colors = global_lineage_colors
)

save_scdiff_circos_legend(
  data = scdiff_circos_tbl,
  direction_oi = "PLA-up",
  output_file = file.path(
    plot_dir,
    "13c_scdiffcom_circos_PLA_up_legend.png"
  ),
  lineage_colors = global_lineage_colors
)

save_scdiff_circos_legend(
  data = scdiff_circos_tbl,
  direction_oi = "platelet-free-up",
  output_file = file.path(
    plot_dir,
    "13d_scdiffcom_circos_platelet_free_up_legend.png"
  ),
  lineage_colors = global_lineage_colors
)

message("Finished scDiffCom analysis for: ", dataset_name)
message("Saved plots to: ", plot_dir)
