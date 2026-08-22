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
  library(patchwork)
})

# ============================================================
# 1. arguments and output
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript scDiffCom.R <input_rds> <output_dir>")
}

options(future.globals.maxSize = 2 * 1024^3)
future::plan(future::sequential)

input_file <- args[[1]]
base_output_dir <- args[[2]]

dataset_name <- tools::file_path_sans_ext(
  basename(input_file)
)

dataset_mode <- dplyr::case_when(
  stringr::str_ends(
    dataset_name,
    "_healthyOnly"
  ) ~ "healthyOnly",

  stringr::str_ends(
    dataset_name,
    "_diseasedOnly"
  ) ~ "diseasedOnly",

  TRUE ~ "all"
)

dataset_clean <- dataset_name %>%
  stringr::str_remove(
    "_(healthyOnly|diseasedOnly)$"
  )

plot_title <- function(title) {
  paste0(
    title,
    "\n",
    dataset_clean,
    " | ",
    dataset_mode
  )
}

plot_dir <- file.path(
  base_output_dir,
  dataset_name,
  "plots"
)

dir.create(
  plot_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

save_plot <- function(
    plot,
    filename,
    width = 8,
    height = 6,
    dpi = 250
) {
  output_file <- file.path(
    plot_dir,
    filename
  )

  message("Saving plot: ", output_file)

  ggsave(
    filename = output_file,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )

  if (!file.exists(output_file)) {
    stop("Plot was not created: ", output_file)
  }

  message(
    "Saved plot: ",
    filename,
    " [",
    file.info(output_file)$size
  )

  invisible(output_file)
}

message("Dataset: ", dataset_name)
message("Input: ", input_file)
message("Plots: ", plot_dir)

# ============================================================
# 2. load and validate data
# ============================================================

seurat_obj <- readRDS(input_file)

condition_col <- "pla_status"
lineage_col <- "lineage"

required_cols <- c(
  condition_col,
  lineage_col
)

missing_cols <- setdiff(
  required_cols,
  colnames(seurat_obj@meta.data)
)

if (length(missing_cols) > 0) {
  stop(
    "Missing metadata columns: ",
    paste(
      missing_cols,
      collapse = ", "
    )
  )
}

keep <- complete.cases(
  seurat_obj@meta.data[,required_cols,drop = FALSE]
)

seurat_obj <- subset(
  seurat_obj,
  cells = colnames(seurat_obj)[keep]
)

DefaultAssay(seurat_obj) <- "RNA"

required_conditions <- c(
  "PLA",
  "platelet-free"
)

available_conditions <- unique(
  as.character(
    seurat_obj@meta.data[[condition_col]]
  )
)

missing_conditions <- setdiff(
  required_conditions,
  available_conditions
)

if (length(missing_conditions) > 0) {
  stop(
    "Missing conditions: ",
    paste(
      missing_conditions,
      collapse = ", "
    )
  )
}

message("Cell counts by condition:")
print(
  table(
    seurat_obj@meta.data[[condition_col]]
  )
)

message("Cell counts by condition and lineage:")
print(
  table(
    seurat_obj@meta.data[[condition_col]],
    seurat_obj@meta.data[[lineage_col]]
  )
)

# ============================================================
# 3. colors
# ============================================================

direction_colors <- c(
  "PLA-up" = "#F8766D",
  "platelet-free-up" = "#00BFC4"
)

lineages <- sort(
  unique(
    as.character(
      seurat_obj@meta.data[[lineage_col]]
    )
  )
)

lineages <- lineages[!is.na(lineages) & lineages != ""]

lineage_colors <- colorRampPalette(
  brewer.pal(11, "Spectral")
)(
  length(lineages)
)

names(lineage_colors) <- lineages

# ============================================================
# 4. run scDiffCom
# ============================================================

scdiffcom_object <- run_interaction_analysis(
  seurat_object = seurat_obj,
  LRI_species = "human",
  seurat_celltype_id = lineage_col,
  seurat_condition_id = list(
    column_name = condition_col,
    cond1_name = "platelet-free",
    cond2_name = "PLA"
  )
)

CCI_detected <- GetTableCCI(
  scdiffcom_object,
  type = "detected",
  simplified = TRUE
)

# ============================================================
# 5. ORA plots
# ============================================================

p_ora_up <- tryCatch(
  {
    PlotORA(
      object = scdiffcom_object,
      category = "LRI",
      regulation = "UP"
    ) +
      theme(
        legend.position = c(
          0.85,
          0.4
        ),
        legend.key.size = unit(
          0.4,
          "cm"
        )
      ) +
      labs(
        title = plot_title(
          "scDiffCom ORA: up-regulated LRIs"
        )
      )
  },
  error = function(e) {
    message(
      "UP ORA plot failed: ",
      conditionMessage(e)
    )

    ggplot() +
      theme_void() +
      annotate(
        "text",
        x = 0,
        y = 0,
        label =
          "No UP ORA results available",
        size = 5
      ) +
      labs(
        title = plot_title(
          "scDiffCom ORA: up-regulated LRIs"
        )
      )
  }
)

save_plot(
  p_ora_up,
  "02_scdiffcom_ORA_LRI_UP.png",
  width = 9,
  height = 7
)

p_ora_down <- tryCatch(
  {
    PlotORA(
      object = scdiffcom_object,
      category = "LRI",
      regulation = "DOWN"
    ) +
      theme(
        legend.position = c(
          0.85,
          0.4
        ),
        legend.key.size = unit(
          0.4,
          "cm"
        )
      ) +
      labs(
        title = plot_title(
          "scDiffCom ORA: down-regulated LRIs"
        )
      )
  },
  error = function(e) {
    message(
      "DOWN ORA plot failed: ",
      conditionMessage(e)
    )

    ggplot() +
      theme_void() +
      annotate(
        "text",
        x = 0,
        y = 0,
        label =
          "No DOWN ORA results available",
        size = 5
      ) +
      labs(
        title = plot_title(
          "scDiffCom ORA: down-regulated LRIs"
        )
      )
  }
)

save_plot(
  p_ora_down,
  "03_scdiffcom_ORA_LRI_DOWN.png",
  width = 9,
  height = 7
)

# ============================================================
# 6. prepare plotting table
# ============================================================

# ER_CELLTYPES is stored as source_target.
# This parser also works when lineage names contain underscores.
parse_lineage_pair <- function(
    values,
    valid_lineages
) {
  valid_lineages <- valid_lineages[order(nchar(valid_lineages), decreasing = TRUE)]

  bind_rows(
    lapply(
      values,
      function(value) {
        if (
          is.na(value) ||
            value == ""
        ) {
          return(
            tibble(
              source = NA_character_,
              target = NA_character_
            )
          )
        }

        source <- valid_lineages[startsWith(value, paste0(valid_lineages,"_"))]

        if (length(source) == 0) {
          return(
            tibble(
              source = NA_character_,
              target = NA_character_
            )
          )
        }

        source <- source[[1]]

        target <- substring(
          value,
          nchar(source) + 2
        )

        if (!target %in% valid_lineages) {
          target <- NA_character_
        }

        tibble(
          source = source,
          target = target
        )
      }
    )
  )
}

parsed_lineages <- parse_lineage_pair(
  CCI_detected$ER_CELLTYPES,
  lineages
)

scdiff_plot <- CCI_detected %>%
  as_tibble() %>%
  bind_cols(parsed_lineages) %>%
  separate(
    LRI,
    into = c(
      "ligand",
      "receptor"
    ),
    sep = ":",
    remove = FALSE,
    extra = "merge",
    fill = "right"
  ) %>%
  mutate(
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
    direction = case_when(
      REGULATION == "UP" ~
        "PLA-up",
      REGULATION == "DOWN" ~
        "platelet-free-up",
      REGULATION == "FLAT" ~
        "flat",
      REGULATION == "NSC" ~
        "not significant",
      TRUE ~ NA_character_
    ),
    direction = factor(
      direction,
      levels = c(
        "PLA-up",
        "platelet-free-up",
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
    " rows could not be parsed into source and target lineages."
  )
}

scdiff_sig <- scdiff_plot %>%
  filter(
    direction %in% c(
      "PLA-up",
      "platelet-free-up"
    ),
    is.finite(score_signed),
    is.finite(score_abs),
    is.finite(padj)
  )

if (nrow(scdiff_sig) == 0) {
  stop(
    "No significant UP or DOWN interactions were found."
  )
}

# ============================================================
# 7. volcano plot
# ============================================================
p_volcano <- ggplot(
  scdiff_plot,
  aes(
    x = score_signed,
    y = -log10(padj + 1e-2),
    color = direction
  )
) +
  geom_point(
    alpha = 0.6,
    size = 1.3
  ) +
  scale_color_manual(
    values = c(
      direction_colors,
      "flat" = "#7CAE00",
      "not significant" =
        "grey70"
    ),
    drop = FALSE
  ) +
  theme_bw() +
  labs(
    title = plot_title(
      "scDiffCom differential cell-cell interactions"
    ),
    x =
      "log fold change: PLA / platelet-free",
    y =
      expression(
        -log[10](
          "adjusted p-value"
        )
      ),
    color = "Regulation"
  )

save_plot(
  p_volcano,
  "01_scdiffcom_volcano.png",
  width = 8,
  height = 6
)

# ============================================================
# 8. global significant interaction counts
# ============================================================

global_counts <- scdiff_sig %>%
  count(
    direction,
    name = "n_interactions"
  )

p_global_counts <- ggplot(
  global_counts,
  aes(
    x = direction,
    y = n_interactions,
    fill = direction
  )
) +
  geom_col() +
  geom_text(
    aes(
      label = n_interactions
    ),
    vjust = -0.4
  ) +
  scale_fill_manual(
    values = direction_colors
  ) +
  theme_bw() +
  labs(
    title = plot_title(
      "Global number of significant scDiffCom CCIs"
    ),
    x = "Direction",
    y = "# significant CCIs"
  ) +
  theme(
    legend.position = "none"
  ) +
  expand_limits(
    y = max(
      global_counts$n_interactions
    ) * 1.08
  )

save_plot(
  p_global_counts,
  "04_scdiffcom_global_significant_counts.png",
  width = 7,
  height = 5
)

# ============================================================
# 9. effect-size distribution
# ============================================================

p_score_distribution <- ggplot(
  scdiff_sig,
  aes(
    x = direction,
    y = score_abs,
    fill = direction
  )
) +
  geom_violin(
    trim = TRUE,
    alpha = 0.45
  ) +
  geom_boxplot(
    width = 0.18,
    outlier.size = 0.3
  ) +
  scale_fill_manual(
    values = direction_colors
  ) +
  theme_bw() +
  labs(
    title = plot_title(
      "scDiffCom differential effect-size distribution"
    ),
    x = "Direction",
    y = "|log fold change|"
  ) +
  theme(
    legend.position = "none"
  )

save_plot(
  p_score_distribution,
  "05_scdiffcom_absolute_logfc_distribution.png",
  width = 7,
  height = 5
)

# ============================================================
# 10. significant interactions by lineage pair
# ============================================================

pair_counts <- scdiff_sig %>%
  filter(
    !is.na(source),
    !is.na(target)
  ) %>%
  count(
    lineage_pair,
    direction,
    name = "n_interactions"
  ) %>%
  group_by(lineage_pair) %>%
  mutate(
    total = sum(n_interactions)
  ) %>%
  ungroup() %>%
  slice_max(
    total,
    n = 25,
    with_ties = FALSE
  ) %>%
  mutate(
    lineage_pair = fct_reorder(
      lineage_pair,
      total
    )
  )

p_pair_counts <- ggplot(
  pair_counts,
  aes(
    x = n_interactions,
    y = lineage_pair,
    fill = direction
  )
) +
  geom_col() +
  scale_fill_manual(
    values = direction_colors
  ) +
  theme_bw() +
  labs(
    title = plot_title(
      "Significant scDiffCom CCIs by lineage pair"
    ),
    x = "# significant CCIs",
    y = "Source → target lineage",
    fill = "Direction"
  )

save_plot(
  p_pair_counts,
  "06_scdiffcom_interactions_by_lineage_pair.png",
  width = 10,
  height = 8
)

# ============================================================
# 11. net significant interaction count
# ============================================================

net_count <- scdiff_sig %>%
  filter(
    !is.na(source),
    !is.na(target)
  ) %>%
  count(
    source,
    target,
    direction,
    name = "n"
  ) %>%
  complete(
    source,
    target,
    direction = factor(
      c(
        "PLA-up",
        "platelet-free-up"
      ),
      levels = c(
        "PLA-up",
        "platelet-free-up"
      )
    ),
    fill = list(n = 0)
  ) %>%
  pivot_wider(
    names_from = direction,
    values_from = n,
    values_fill = 0
  ) %>%
  mutate(
    net_count =
      `PLA-up` -
      `platelet-free-up`
  )

net_count_limit <- max(
  abs(net_count$net_count),
  na.rm = TRUE
)

if (
  !is.finite(net_count_limit) ||
    net_count_limit == 0
) {
  net_count_limit <- 1
}

p_net_count <- ggplot(
  net_count,
  aes(
    x = target,
    y = source,
    fill = net_count
  )
) +
  geom_tile(
    color = "white"
  ) +
  geom_text(
    aes(
      label = net_count
    ),
    size = 3.5
  ) +
  scale_fill_gradient2(
    low = direction_colors[["platelet-free-up"]],
    mid = "white",
    high = direction_colors[["PLA-up"]],
    midpoint = 0,
    limits = c(
      -net_count_limit,
      net_count_limit
    )
  ) +
  theme_bw() +
  labs(
    title = plot_title(
      "Net number of significant scDiffCom CCIs"
    ),
    x = "Target lineage",
    y = "Source lineage",
    fill =
      "PLA-up count\nminus platelet-free-up count"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

save_plot(
  p_net_count,
  "07_scdiffcom_net_count_heatmap.png",
  width = 9,
  height = 8
)

# ============================================================
# 12. mean signed logFC
# ============================================================

mean_logfc <- scdiff_sig %>%
  filter(
    !is.na(source),
    !is.na(target)
  ) %>%
  group_by(
    source,
    target
  ) %>%
  summarise(
    mean_logfc = mean(
      score_signed,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

logfc_limit <- max(
  abs(mean_logfc$mean_logfc),
  na.rm = TRUE
)

if (
  !is.finite(logfc_limit) ||
    logfc_limit == 0
) {
  logfc_limit <- 1
}

p_mean_logfc <- ggplot(
  mean_logfc,
  aes(
    x = target,
    y = source,
    fill = mean_logfc
  )
) +
  geom_tile(
    color = "white"
  ) +
  geom_text(
    aes(
      label = round(
        mean_logfc,
        2
      )
    ),
    size = 3.5
  ) +
  scale_fill_gradient2(
    low = direction_colors[["platelet-free-up"]],
    mid = "white",
    high = direction_colors[["PLA-up"]],
    midpoint = 0,
    limits = c(
      -logfc_limit,
      logfc_limit
    )
  ) +
  theme_bw() +
  labs(
    title = plot_title(
      "Mean signed logFC among significant scDiffCom CCIs"
    ),
    x = "Target lineage",
    y = "Source lineage",
    fill = "Mean signed\nlogFC"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

save_plot(
  p_mean_logfc,
  "08_scdiffcom_mean_logfc_heatmap.png",
  width = 9,
  height = 8
)

# ============================================================
# 13. recurrent LR pairs
# ============================================================

recurrent_lr <- scdiff_sig %>%
  filter(
    !is.na(interaction),
    !is.na(lineage_pair)
  ) %>%
  distinct(
    interaction,
    lineage_pair,
    direction
  ) %>%
  count(
    interaction,
    direction,
    name = "n_lineage_pairs"
  ) %>%
  group_by(direction) %>%
  slice_max(
    n_lineage_pairs,
    n = 15,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  mutate(
    plot_label = paste(
      interaction,
      direction,
      sep = "___"
    ),
    plot_label = fct_reorder(
      plot_label,
      n_lineage_pairs
    )
  )

p_recurrent_lr <- ggplot(
  recurrent_lr,
  aes(
    x = n_lineage_pairs,
    y = plot_label,
    fill = direction
  )
) +
  geom_col() +
  facet_wrap(
    ~ direction,
    scales = "free_y"
  ) +
  scale_y_discrete(
    labels = function(x) {
      sub(
        "___.*$",
        "",
        x
      )
    }
  ) +
  scale_fill_manual(
    values = direction_colors
  ) +
  theme_bw() +
  labs(
    title = plot_title(
      "Most recurrent differential scDiffCom LR pairs"
    ),
    x = "# source-target lineage pairs",
    y = "Ligand → receptor"
  ) +
  theme(
    legend.position = "none"
  )

save_plot(
  p_recurrent_lr,
  "09_scdiffcom_recurrent_lr_pairs.png",
  width = 12,
  height = 9
)

# ============================================================
# 14. top LR pairs by mean absolute logFC
# ============================================================

top_lr <- scdiff_sig %>%
  filter(
    !is.na(interaction)
  ) %>%
  group_by(
    interaction,
    direction
  ) %>%
  summarise(
    mean_abs_logfc = mean(
      score_abs,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  group_by(direction) %>%
  slice_max(
    mean_abs_logfc,
    n = 15,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  mutate(
    plot_label = paste(
      interaction,
      direction,
      sep = "___"
    ),
    plot_label = fct_reorder(
      plot_label,
      mean_abs_logfc
    )
  )

p_top_lr <- ggplot(
  top_lr,
  aes(
    x = mean_abs_logfc,
    y = plot_label,
    fill = direction
  )
) +
  geom_col() +
  facet_wrap(
    ~ direction,
    scales = "free_y"
  ) +
  scale_y_discrete(
    labels = function(x) {
      sub(
        "___.*$",
        "",
        x
      )
    }
  ) +
  scale_fill_manual(
    values = direction_colors
  ) +
  theme_bw() +
  labs(
    title = plot_title(
      "Top scDiffCom LR pairs by mean absolute logFC"
    ),
    x = "Mean |logFC|",
    y = "Ligand → receptor"
  ) +
  theme(
    legend.position = "none"
  )

save_plot(
  p_top_lr,
  "10_scdiffcom_top_lr_by_mean_logfc.png",
  width = 12,
  height = 9
)

# ============================================================
# 15. top interactions in top lineage pairs
# ============================================================

n_top_pairs <- 5
n_top_interactions <- 5

# Rank lineage pairs separately within each direction.
# Here, pairs are ranked by the number of significant interactions.
# Mean and maximum absolute logFC are used as tie-breakers.
top_lineage_pairs <- scdiff_sig %>%
  filter(
    !is.na(lineage_pair),
    !is.na(interaction)
  ) %>%
  group_by(
    direction,
    lineage_pair
  ) %>%
  summarise(
    n_significant = n(),
    mean_abs_logfc = mean(
      score_abs,
      na.rm = TRUE
    ),
    max_abs_logfc = max(
      score_abs,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  group_by(direction) %>%
  arrange(
    desc(n_significant),
    desc(mean_abs_logfc),
    desc(max_abs_logfc),
    .by_group = TRUE
  ) %>%
  slice_head(
    n = n_top_pairs
  ) %>%
  ungroup()

message("Top lineage pairs selected for plot 11:")
print(top_lineage_pairs)

# Keep only the selected lineage pairs.
# Within each panel, retain the five interactions with the
# largest absolute logFC.
top_per_pair <- scdiff_sig %>%
  filter(
    !is.na(lineage_pair),
    !is.na(interaction)
  ) %>%
  semi_join(
    top_lineage_pairs,
    by = c(
      "direction",
      "lineage_pair"
    )
  ) %>%
  group_by(
    direction,
    lineage_pair
  ) %>%
  slice_max(
    order_by = score_abs,
    n = n_top_interactions,
    with_ties = FALSE
  ) %>%
  arrange(
    direction,
    lineage_pair,
    desc(score_abs)
  ) %>%
  ungroup() %>%
  mutate(
    plot_label = paste(
      interaction,
      lineage_pair,
      direction,
      sep = "___"
    )
  )

# Reorder interaction labels separately inside every facet.
# The reversed levels put the largest |logFC| at the top.
top_per_pair <- top_per_pair %>%
  group_by(
    direction,
    lineage_pair
  ) %>%
  mutate(
    plot_label = factor(
      plot_label,
      levels = rev(
        unique(
          plot_label[
            order(
              score_abs,
              decreasing = TRUE
            )
          ]
        )
      )
    )
  ) %>%
  ungroup()

# Keep facet order consistent with the lineage-pair ranking.
pair_levels <- top_lineage_pairs %>%
  arrange(
    direction,
    desc(n_significant),
    desc(mean_abs_logfc),
    desc(max_abs_logfc)
  ) %>%
  pull(lineage_pair) %>%
  unique()

top_per_pair <- top_per_pair %>%
  mutate(
    lineage_pair = factor(
      lineage_pair,
      levels = pair_levels
    )
  )

p_top_per_pair <- ggplot(
  top_per_pair,
  aes(
    x = score_abs,
    y = plot_label,
    color = direction,
    size = pmin(
      -log10(
        pmax(
          padj,
          1e-300
        )
      ),
      20
    )
  )
) +
  geom_point(
    alpha = 0.85
  ) +
  facet_grid(
    rows = vars(direction),
    cols = vars(lineage_pair),
    scales = "free_y",
    space = "free_x",
    drop = TRUE
  ) +
  scale_y_discrete(
    labels = function(x) {
      sub(
        "___.*$",
        "",
        x
      )
    }
  ) +
  scale_color_manual(
    values = direction_colors
  ) +
  scale_size_continuous(
    name = expression(
      min(
        -log[10](
          "adjusted p-value"
        ),
        20
      )
    )
  ) +
  theme_bw(
    base_size = 10
  ) +
  labs(
    title = plot_title(
      "Top interactions in top 5 lineage pairs per direction"
    ),
    subtitle = paste0(
      "Lineage pairs ranked by number of significant interactions; ",
      "interactions ranked by descending |logFC|"
    ),
    x = "|logFC|",
    y = "Ligand → receptor",
    color = "Direction"
  ) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(
      size = 8
    ),
    axis.text.y = element_text(
      size = 7
    )
  )

save_plot(
  p_top_per_pair,
  "11_scdiffcom_top_interactions_per_lineage_pair.png",
  width = 18,
  height = 10
)

# ============================================================
# 16. top individual interactions
# ============================================================

top_individual <- scdiff_sig %>%
  filter(
    !is.na(interaction),
    !is.na(lineage_pair)
  ) %>%
  group_by(direction) %>%
  slice_max(
    score_abs,
    n = 20,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  mutate(
    plot_label = paste0(
      interaction,
      "\n",
      lineage_pair
    ),
    plot_label = paste(
      plot_label,
      direction,
      sep = "___"
    ),
    plot_label = fct_reorder(
      plot_label,
      score_abs
    )
  )

p_top_individual <- ggplot(
  top_individual,
  aes(
    x = score_abs,
    y = plot_label,
    color = direction
  )
) +
  geom_point(
    size = 3,
    alpha = 0.85
  ) +
  facet_wrap(
    ~ direction,
    scales = "free_y"
  ) +
  scale_y_discrete(
    labels = function(x) {
      sub(
        "___.*$",
        "",
        x
      )
    }
  ) +
  scale_color_manual(
    values = direction_colors
  ) +
  theme_bw() +
  labs(
    title = plot_title(
      "Top individual scDiffCom interactions"
    ),
    x = "|logFC|",
    y =
      "Ligand → receptor / source → target"
  ) +
  theme(
    legend.position = "none"
  )

save_plot(
  p_top_individual,
  "12_scdiffcom_top_individual_interactions.png",
  width = 13,
  height = 10
)

# ============================================================
# 17. gene-level circos plots
# ============================================================

circos_data <- scdiff_sig %>%
  filter(
    !is.na(source),
    !is.na(target),
    !is.na(ligand),
    !is.na(receptor),
    ligand != "",
    receptor != ""
  ) %>%
  group_by(direction) %>%
  slice_max(
    score_abs,
    n = 20,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  mutate(
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
  )

plot_circos <- function(
    direction_oi,
    output_file
) {
  plot_data <- circos_data %>%
    filter(
      direction == direction_oi
    )

  if (nrow(plot_data) == 0) {
    message(
      "Skipping ",
      direction_oi,
      " circos plot."
    )
    return(invisible(NULL))
  }

  links <- plot_data %>%
    group_by(
      ligand_node,
      receptor_node
    ) %>%
    summarise(
      value = max(
        score_abs,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    transmute(
      from = ligand_node,
      to = receptor_node,
      value
    )

  node_colors <- c(
    setNames(
      lineage_colors[plot_data$source],
      plot_data$ligand_node
    ),
    setNames(
      lineage_colors[plot_data$target],
      plot_data$receptor_node
    )
  )

  node_colors <- node_colors[!duplicated(names(node_colors))]

  sector_order <- c(
    unique(
      plot_data$ligand_node
    ),
    unique(
      plot_data$receptor_node
    )
  )

  gap_after <- c(
    rep(
      1.5,
      max(
        length(
          unique(
            plot_data$ligand_node
          )
        ) - 1,
        0
      )
    ),
    8,
    rep(
      1.5,
      max(
        length(
          unique(
            plot_data$receptor_node
          )
        ) - 1,
        0
      )
    ),
    8
  )

  png(
    filename = output_file,
    width = 3500,
    height = 3500,
    res = 250
  )

  par(
    mar = c(
      1,
      1,
      4,
      1
    ),
    xpd = NA
  )

  circos.clear()

  circos.par(
    start.degree = 90,
    canvas.xlim = c(
      -1.15,
      1.15
    ),
    canvas.ylim = c(
      -1.15,
      1.15
    ),
    gap.after = gap_after,
    points.overflow.warning = FALSE
  )

  chordDiagram(
    links,
    order = sector_order,
    grid.border = "grey35",
    grid.col = node_colors,
    col = node_colors[links$from],
    transparency = 0.10,
    directional = 1,
    direction.type = c(
      "arrows",
      "diffHeight"
    ),
    diffHeight = -0.03,
    link.arr.type = "triangle",
    link.sort = TRUE,
    link.decreasing = TRUE,
    link.border = adjustcolor(
      "grey35",
      alpha.f = 0.55
    ),
    annotationTrack = "grid",
    preAllocateTracks = list(
      track.height = 0.16
    )
  )

  circos.trackPlotRegion(
    track.index = 1,
    bg.border = NA,
    panel.fun = function(x, y) {
      sector <- get.cell.meta.data(
        "sector.index"
      )

      xlim <- get.cell.meta.data(
        "xlim"
      )

      ylim <- get.cell.meta.data(
        "ylim"
      )

      gene <- sub(
        "^[LR]::[^:]+::",
        "",
        sector
      )

      circos.text(
        mean(xlim),
        mean(ylim),
        gene,
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(
          0,
          0.5
        ),
        cex = 1.2
      )
    }
  )

  title(
    main = plot_title(
      paste0(
        "scDiffCom gene-level interactions: ",
        direction_oi
      )
    ),
    cex.main = 2.0,
    line = 0.5
  )

  legend(
    "topleft",
    legend = names(
      lineage_colors
    ),
    fill = lineage_colors,
    border = NA,
    cex = 1.2,
    pt.cex = 1.3,
    x.intersp = 0.8,
    y.intersp = 1.15,
    bty = "n"
  )

  circos.clear()
  dev.off()
}

plot_circos(
  "PLA-up",
  file.path(
    plot_dir,
    "13a_scdiffcom_circos_PLA_up.png"
  )
)

plot_circos(
  "platelet-free-up",
  file.path(
    plot_dir,
    "13b_scdiffcom_circos_platelet_free_up.png"
  )
)

message(
  "Finished scDiffCom analysis for: ",
  dataset_name
)
message(
  "Saved plots to: ",
  plot_dir
)
