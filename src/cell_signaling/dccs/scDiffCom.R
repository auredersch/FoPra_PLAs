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
  library(multinichenetr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript scDiffCom.R <input_rds> <output_dir>")
}

options(future.globals.maxSize = 2 * 1024^3)  # allow up to 2 GB
plan(sequential)

input_file <- args[[1]]
base_output_dir <- args[[2]]

dataset_name <- tools::file_path_sans_ext(basename(input_file))

dataset_mode <- stringr::str_extract(
  dataset_name,
  "(all|diseasedOnly|healthyOnly)$"
)

dataset_clean <- dataset_name %>%
  stringr::str_remove("_(all|diseasedOnly|healthyOnly)$")

if (is.na(dataset_mode)) {
  dataset_mode <- "unknownMode"
}

plot_title <- function(title) {
  paste0(title, "\n", dataset_clean, " | ", dataset_mode)
}

out_dir <- file.path(base_output_dir, dataset_name)
plot_dir <- file.path(out_dir, "plots")
table_dir <- file.path(out_dir, "tables")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

message("Analyzing: ", dataset_name)
message("Input: ", input_file)
message("Output directory: ", out_dir)
message("Plot directory: ", plot_dir)

save_plot <- function(plot, filename, width = 8, height = 6, dpi = 300) {
  out_file <- file.path(plot_dir, filename)

  ggsave(
    filename = out_file,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )

  message("Saved plot: ", out_file)
}

seurat_obj <- readRDS(input_file)
seurat_obj <- subset(
  seurat_obj,
  subset = !is.na(lineage) & !is.na(pla_status)
)
DefaultAssay(seurat_obj)

sample_col    <- "sample"       # technical / sample-level replicate
donor_col     <- "donor_id"     # biological donor
celltype_col  <- "celltype"     # broad cell type
celltype_full <- "celltype_full" # more detailed cell type
condition_col <- "pla_status"   # PLA vs platelet-free
lineage_col   <- "lineage"

table(seurat_obj$pla_status, useNA = "ifany")
table(seurat_obj$pla_status, seurat_obj$lineage)

required_conditions <- c("PLA", "platelet-free")
available_conditions <- unique(as.character(seurat_obj$pla_status))

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
    "Cannot perform PLA versus platelet-free comparison. ",
    "Missing condition(s): ",
    paste(missing_conditions, collapse = ", ")
  )
}

# run default analyis
scdiffcom_object <- run_interaction_analysis(
  seurat_object = seurat_obj,
  LRI_species = "human",
  seurat_celltype_id = "lineage",
  seurat_condition_id = list(
    column_name = "pla_status",
    cond1_name = "platelet-free", #log(score(cond2_name) / score(cond1_name))
    cond2_name = "PLA"
  )
)

saveRDS(
  scdiffcom_object,
  file.path(out_dir, "scdiffcom_object.rds")
)

# -------------------------
# explore and save results
# -------------------------

CCI_detected <- GetTableCCI(
  scdiffcom_object,
  type = "detected",
  simplified = TRUE
)

ORA_results <- GetTableORA(
  scdiffcom_object,
  categories = "all",
  simplified = TRUE
)

write.csv(
  CCI_detected,
  file.path(table_dir, "scdiffcom_CCI_detected.csv"),
  row.names = FALSE
)

saveRDS(
  ORA_results,
  file.path(out_dir, "scdiffcom_ORA_results.rds")
)

cci_regulation_counts <- CCI_detected %>%
  dplyr::count(REGULATION)


# -------------------------
# 1. volcano-like CCI plot
# -------------------------

p_scdiff_volcano <- ggplot(
  CCI_detected,
  aes(
    x = LOGFC,
    y = -log10(BH_P_VALUE_DE + 1E-2),
    colour = REGULATION
  )
) +
  geom_point() +
  scale_colour_manual(
    values = c(
      "UP" = "#F8766D",
      "DOWN" = "#00BFC4",
      "FLAT" = "#7CAE00",
      "NSC" = "grey70"
    )
  ) +
  theme_bw() +
  xlab("log(FC)") +
  ylab("-log10(Adj. p-value)") +
  labs(
    title = plot_title("scDiffCom detected CCIs"),
    colour = "Regulation"
  )

save_plot(
  p_scdiff_volcano,
  "01_scdiffcom_detected_CCI_volcano.png",
  width = 7,
  height = 6
)


# -------------------------
# 2. ORA plot: up-regulated LRIs
# -------------------------

p_scdiff_ORA_UP <- tryCatch(
  {
    PlotORA(
      object = scdiffcom_object,
      category = "LRI",
      regulation = "UP"
    ) +
      theme(
        legend.position = c(0.85, 0.4),
        legend.key.size = unit(0.4, "cm")
      ) +
      labs(title = plot_title("scDiffCom ORA: up-regulated LRIs"))
  },
  error = function(e) {
    message("PlotORA UP failed: ", conditionMessage(e))

    ggplot() +
      theme_void() +
      annotate("text", x = 0, y = 0, label = "No UP ORA results available", size = 5) +
      labs(title = plot_title("scDiffCom ORA: up-regulated LRIs"))
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
    PlotORA(
      object = scdiffcom_object,
      category = "LRI",
      regulation = "DOWN"
    ) +
      theme(
        legend.position = c(0.85, 0.4),
        legend.key.size = unit(0.4, "cm")
      ) +
      labs(title = plot_title("scDiffCom ORA: down-regulated LRIs"))
  },
  error = function(e) {
    message("PlotORA DOWN failed: ", conditionMessage(e))

    ggplot() +
      theme_void() +
      annotate("text", x = 0, y = 0, label = "No DOWN ORA results available", size = 5) +
      labs(title = plot_title("scDiffCom ORA: down-regulated LRIs"))
  }
)

save_plot(
  p_scdiff_ORA_DOWN,
  "03_scdiffcom_ORA_LRI_DOWN.png",
  width = 9,
  height = 7
)


# -------------------------
# prepare plotting table
# -------------------------

scdiff_plot <- CCI_detected %>%
  tibble::as_tibble() %>%
  tidyr::separate(
    ER_CELLTYPES,
    into = c("source", "target"),
    sep = "_",
    remove = FALSE
  ) %>%
  tidyr::separate(
    LRI,
    into = c("ligand", "receptor"),
    sep = ":",
    remove = FALSE
  ) %>%
  dplyr::mutate(
    lineage_pair = paste(source, target, sep = " → "),
    interaction = paste(ligand, receptor, sep = " → "),

    direction = dplyr::case_when(
      REGULATION == "UP" ~ "PLA-up",
      REGULATION == "DOWN" ~ "platelet-free-up",
      REGULATION == "FLAT" ~ "flat",
      REGULATION == "NSC" ~ "not significant",
      TRUE ~ as.character(REGULATION)
    ),

    score_signed = LOGFC,
    score_abs = abs(LOGFC),
    padj = BH_P_VALUE_DE
  )

scdiff_sig <- scdiff_plot %>%
  dplyr::filter(direction %in% c("PLA-up", "platelet-free-up"))

scdiff_direction_counts <- scdiff_plot %>%
  dplyr::count(direction) # hier sind auch flats dabei

scdiff_sig_direction_counts <- scdiff_sig %>%
  dplyr::count(direction)

#this is counting rows, not necessarily unique ligand–receptor pairs
#same LR pair can appear several times in different source–target lineage combinations

write.csv(
  scdiff_plot,
  file.path(table_dir, "scdiffcom_plot_table_all.csv"),
  row.names = FALSE
)

write.csv(
  scdiff_sig,
  file.path(table_dir, "scdiffcom_plot_table_significant.csv"),
  row.names = FALSE
)


# -------------------------
# 4. global significant CCI distribution
# -------------------------

p_scdiff_global <- ggplot(scdiff_sig, aes(x = direction, fill = direction)) +
  geom_bar() +
  theme_bw() +
  scale_fill_manual(
    values = c(
      "PLA-up" = "#F8766D",
      "platelet-free-up" = "#00BFC4"
    )
  ) +
  labs(
    title = plot_title("Global distribution of significant scDiffCom CCIs"),
    x = "Interaction category",
    y = "# significant CCIs"
  ) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

save_plot(
  p_scdiff_global,
  "04_global_distribution_significant_scdiffcom_CCIs.png",
  width = 7,
  height = 5
)


# -------------------------
# 5. direction by lineage pair
# -------------------------

pair_counts_scdiff <- scdiff_sig %>%
  dplyr::count(lineage_pair, direction) %>%
  dplyr::group_by(lineage_pair) %>%
  dplyr::mutate(total = sum(n)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(total)) %>%
  dplyr::mutate(lineage_pair = fct_reorder(lineage_pair, total))

p_pair_counts_scdiff <- ggplot(
  pair_counts_scdiff,
  aes(x = n, y = lineage_pair, fill = direction)
) +
  geom_col() +
  theme_bw() +
  scale_fill_manual(
    values = c(
      "PLA-up" = "#F8766D",
      "platelet-free-up" = "#00BFC4"
    )
  ) +
  labs(
    title = plot_title("Direction of differential scDiffCom CCIs by lineage pair"),
    x = "# differential CCIs",
    y = "Source → target lineage",
    fill = "Direction"
  )

save_plot(
  p_pair_counts_scdiff,
  "05_direction_by_lineage_pair_scdiffcom.png",
  width = 9,
  height = 7
)


# -------------------------
# 6. net direction heatmap
# -------------------------

net_direction_scdiff <- scdiff_sig %>%
  dplyr::count(source, target, direction) %>%
  tidyr::complete(
    source,
    target,
    direction = c("PLA-up", "platelet-free-up"),
    fill = list(n = 0)
  ) %>%
  tidyr::pivot_wider(
    names_from = direction,
    values_from = n,
    values_fill = 0
  ) %>%
  dplyr::mutate(
    net = `PLA-up` - `platelet-free-up`
  )

p_net_direction_scdiff <- ggplot(
  net_direction_scdiff,
  aes(x = target, y = source, fill = net)
) +
  geom_tile(color = "white") +
  geom_text(aes(label = net), size = 4) +
  theme_bw() +
  labs(
    title = plot_title("Net direction of differential scDiffCom CCIs"),
    x = "Target lineage",
    y = "Source lineage",
    fill = "PLA-up minus\nplatelet-free-up"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

save_plot(
  p_net_direction_scdiff,
  "06_net_direction_scdiffcom_heatmap.png",
  width = 8,
  height = 7
)


# -------------------------
# 7. most recurrent LR pairs
# -------------------------

top_recurrent_scdiff <- scdiff_sig %>%
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
  dplyr::ungroup()

p_top_recurrent_scdiff <- ggplot(
  top_recurrent_scdiff,
  aes(
    x = n_lineage_pairs,
    y = forcats::fct_reorder(
      interaction,
      n_lineage_pairs
    ),
    fill = direction
  )
) +
  geom_col() +
  facet_wrap(
    ~ direction,
    scales = "free_y"
  ) +
  theme_bw() +
  scale_fill_manual(
    values = c(
      "PLA-up" = "#F8766D",
      "platelet-free-up" = "#00BFC4"
    )
  ) +
  labs(
    title = plot_title(
      "Most recurrent differential scDiffCom ligand-receptor pairs"
    ),
    x = "# source-target lineage pairs",
    y = "Ligand → receptor",
    fill = "Direction"
  )

save_plot(
  p_top_recurrent_scdiff,
  "07_top_recurrent_scdiffcom_lr_pairs.png",
  width = 8,
  height = 8
)


# -------------------------
# 8. top CCIs per lineage pair
# -------------------------

top_scdiff_per_pair <- scdiff_sig %>%
  dplyr::group_by(lineage_pair) %>%
  dplyr::slice_max(score_abs, n = 5, with_ties = FALSE) %>%
  dplyr::ungroup()

p_top_scdiff_per_pair <- ggplot(
  top_scdiff_per_pair,
  aes(
    x = lineage_pair,
    y = interaction,
    size = score_abs,
    color = score_signed
  )
) +
  geom_point(alpha = 0.8) +
  theme_bw() +
  labs(
    title = plot_title("Top differential scDiffCom CCIs per lineage pair"),
    x = "Source → target lineage",
    y = "Ligand → receptor",
    size = "|logFC|",
    color = "logFC"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    plot.title = element_text(size = 10),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7)
  )

save_plot(
  p_top_scdiff_per_pair,
  "08_top_scdiffcom_CCI_per_lineage_pair.png",
  width = 13,
  height = 10
)

# -------------------------
# 9. gene-level scDiffCom circos plots
# -------------------------

top_n_circos <- 20

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
  dplyr::ungroup() %>%
  dplyr::mutate(
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


all_circos_celltypes <- sort(unique(c(
  scdiff_circos_tbl$source,
  scdiff_circos_tbl$target
)))

base_palette <- RColorBrewer::brewer.pal(
  n = 11,
  name = "Spectral"
)

global_celltype_colors <- grDevices::colorRampPalette(
  base_palette
)(length(all_circos_celltypes))

names(global_celltype_colors) <- all_circos_celltypes

plot_scdiff_gene_circos <- function(
    data,
    direction_oi,
    output_file,
    title_text,
    celltype_colors
) {

  plot_data <- data %>%
    dplyr::filter(direction == direction_oi)

  if (nrow(plot_data) == 0) {
    message("Skipping ", direction_oi, " circos plot: no interactions.")
    return(invisible(NULL))
  }

  # One row per source-ligand-target-receptor interaction
  links <- plot_data %>%
    dplyr::transmute(
      from = ligand_node,
      to = receptor_node,
      value = score_abs, # chord wdith
      source = source,
      target = target
    )

  # Cell-type colour palette
  celltypes <- sort(unique(c(
    plot_data$source,
    plot_data$target
  )))

  direction_colors <- celltype_colors[celltypes]

  # Assign each ligand/receptor sector to its corresponding lineage
  ligand_annotation <- plot_data %>%
    dplyr::distinct(
      node = ligand_node,
      celltype = source
    )

  receptor_annotation <- plot_data %>%
    dplyr::distinct(
      node = receptor_node,
      celltype = target
    )

  node_annotation <- dplyr::bind_rows(
    ligand_annotation,
    receptor_annotation
  ) %>%
    dplyr::distinct(node, celltype)

  sector_colors <- direction_colors[node_annotation$celltype]
  names(sector_colors) <- node_annotation$node

  # Keep ligands together and receptors together
  sector_order <- c(
    unique(plot_data$ligand_node),
    unique(plot_data$receptor_node)
  )

  grDevices::png(
    filename = output_file,
    width = 4200,
    height = 4200,
    res = 300
  )

 graphics::par(
  mar = c(1, 1, 3, 1),
  xpd = NA
)

circlize::circos.par(
  start.degree = 90,
  canvas.xlim = c(-1.08, 1.08),
  canvas.ylim = c(-1.08, 1.08),
  gap.after = c(
    rep(1.5, length(unique(plot_data$ligand_node)) - 1),
    8,
    rep(1.5, length(unique(plot_data$receptor_node)) - 1),
    8
  ),
  track.margin = c(0.005, 0.005),
  points.overflow.warning = FALSE
)

  circlize::chordDiagram(
    x = links %>%
      dplyr::select(from, to, value),
    order = sector_order,
    grid.col = sector_colors,
    col = sector_colors[links$from],
    transparency = 0.35,
    directional = 1,
    direction.type = c("arrows", "diffHeight"),
    diffHeight = -0.03,
    link.arr.type = "big.arrow",
    link.sort = TRUE,
    link.decreasing = TRUE,
    link.border = "grey35",
    link.lwd = 0.5,
    annotationTrack = "grid",
    preAllocateTracks = list(
      track.height = 0.14
    )
  )

  circlize::circos.trackPlotRegion(
    track.index = 1,
    bg.border = NA,
    panel.fun = function(x, y) {

      sector <- circlize::get.cell.meta.data("sector.index")
      xlim <- circlize::get.cell.meta.data("xlim")
      ylim <- circlize::get.cell.meta.data("ylim")

      # Remove internal L:/R: prefix from displayed label
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
        cex = 0.9
      )
    }
  )

  graphics::title(
    main = title_text,
    cex.main = 1.3,
    line = 1
  )

  circlize::circos.clear()
  grDevices::dev.off()
}

save_scdiff_circos_legend <- function(
    data,
    direction_oi,
    output_file,
    celltype_colors
) {

  legend_data <- data %>%
    dplyr::filter(direction == direction_oi)

  if (nrow(legend_data) == 0) {
    return(invisible(NULL))
  }

  celltypes <- sort(unique(c(
    legend_data$source,
    legend_data$target
  )))

  direction_colors <- celltype_colors[celltypes]

  grDevices::png(
    filename = output_file,
    width = 1800,
    height = max(900, 180 + 140 * length(celltypes)),
    res = 300
  )

  graphics::par(mar = c(1, 1, 1, 1))
  graphics::plot.new()

  graphics::legend(
    "center",
    legend = celltypes,
    fill = direction_colors[celltypes],
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
    "09a_scdiffcom_gene_circos_PLA_up.png"
  ),
  title_text = plot_title(
    "scDiffCom gene-level interactions: PLA-up"
  ),
  celltype_colors = global_celltype_colors
)

plot_scdiff_gene_circos(
  data = scdiff_circos_tbl,
  direction_oi = "platelet-free-up",
  output_file = file.path(
    plot_dir,
    "09b_scdiffcom_gene_circos_platelet_free_up.png"
  ),
  title_text = plot_title(
    "scDiffCom gene-level interactions: platelet-free-up"
  ),
  celltype_colors = global_celltype_colors
)

save_scdiff_circos_legend(
  data = scdiff_circos_tbl,
  direction_oi = "PLA-up",
  output_file = file.path(
    plot_dir,
    "09c_scdiffcom_gene_circos_PLA_up_legend.png"
  ),
  celltype_colors = global_celltype_colors
)

save_scdiff_circos_legend(
  data = scdiff_circos_tbl,
  direction_oi = "platelet-free-up",
  output_file = file.path(
    plot_dir,
    "09d_scdiffcom_gene_circos_platelet_free_up_legend.png"
  ),
  celltype_colors = global_celltype_colors
)

message("Finished scDiffCom analysis for: ", dataset_name)
message("Saved plots to: ", plot_dir)
message("Saved tables to: ", table_dir)