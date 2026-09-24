# TODO Violin Plot - COunt Distribution per Sample pre/post filtering
# Vln Mitrochondirial content per sample pre/post filtering
# nach scDblFinder UMAP plot colored singlet/doublet
# Distribution Highly Variable Genes
# res 0.5, 1.0, 1.5, 2.0 Clustree plot (UMAP)
# Annotation by Cell Type plotten
# Standard QC Plots aus Tutorial

# Keine Thresholds, approach aus scanpy tutorial

library(optparse)
library(Seurat)
library(scDblFinder)
library(SingleCellExperiment)
library(SoupX)
library(Matrix)
library(matrixStats)
library(scry)
library(ggplot2)
library(ggrepel)
library(clustree)
library(harmony)

options <- list(
  make_option(c("-p", "--plot-dir"), type = "character", default = "new_plots", required = FALSE,
              help = "Directory to save plots"),
  make_option(c("-d", "--data-path"), type = "character", default = "data/GSM5008737_RNA_3P/", required = FALSE,
              help = "Path to data directory"),
  make_option(c("-s", "--sample-col"), type = "character", default = "sample", required = FALSE,
              help = "Column name in metadata for sample identification"),
  make_option(c("-b", "--batch-col"), type = "character", default = "batch", required = FALSE,
              help = "Column name in metadata for batch identification"),
  make_option("--soupx", action = "store_true", default = FALSE,
              help = "Whether to run SoupX correction"),
  make_option("--harmony", action = "store_true", default = FALSE,
              help = "Whether to run Harmony batch correction"),
  make_option("--harmony-col", type = "character", default = "batch", required = FALSE,
              help = "Column name in metadata for Harmony batch correction"),
  make_option(c("-o", "--output"), type = "character", default = "data/seu_final.rds", required = FALSE,
              help = "Output path for the final Seurat object"),
  make_option(c("-m", "--metadata-path"), type = "character", default = NULL, required = FALSE,
              help = "Path to metadata CSV file (optional)")
)

parser <- OptionParser(option_list = options)

if (interactive()) {
  args <- list(
    plot_dir = "results/rna_qc/stemi_no_qc",
    data_path = "/nfs/home/students/f.mathis/FoPra_PLAs/data/stemi_no_qc",
    sample_col = NULL,
    batch_col = NULL,
    soupx = FALSE,
    harmony = FALSE,
    harmony_col = NULL,
    metadata_path = NULL,
    output = "data/stemi_no_qc_final.rds"
  )
} else {
  args <- parse_args(parser)
}

save_plot <- function(plot, filename, width = 8, height = 6) {
  ggsave(
    filename = file.path(plot_dir, filename),
    plot = plot,
    width = width,
    height = height,
    dpi = 300
  )

    ggsave(
    filename = file.path(plot_dir, paste0(filename, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    dpi = 300
  )
}

load_seurat <- function(input_path, metadata_path = NULL, project_name = "plateletEnrichment", min_cells = 3, min_features = 20) {
    
    if(grepl("\\.rds$", input_path, ignore.case = TRUE)) {
      seu <- readRDS(input_path)
      return(seu)
    } 
    counts <- Read10X(data.dir = input_path, gene.column = 1)
    counts <- as(counts, "dgCMatrix")

    seu <- CreateSeuratObject(
      counts = counts,
      project = project_name,
      min.cells = min_cells,
      min.features = min_features
    )

    if(!is.null(metadata_path)) {
      metadata <- read.csv(metadata_path, row.names = 1)
      seu_ids <- colnames(seu)
      meta_use <- metadata[seu_ids, , drop = FALSE]
      seu <- AddMetaData(seu, metadata = meta_use)

    } else {
      warning("No Metadata path provided when loading from non-RDS input.")
    }

  return(seu)
}


set_metadata_group <- function(seu, columns, target, sep = "_") {
  columns <- unlist(columns)

  meta <- seu@meta.data

  if (!all(columns %in% colnames(meta))) {
    missing_cols <- setdiff(columns, colnames(meta))

    stop(
      "Metadata column(s) not found: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  if (target %in% colnames(meta)) {
    seu[[paste0(target, "_old")]] <- meta[[target]]
  }

  if (length(columns) == 1) {
    seu[[target]] <- meta[[columns]]
  } else {
    seu[[target]] <- do.call(
      paste,
      c(meta[columns], sep = sep)
    )
  }

  return(seu)
}


plot_dir <- args$plot_dir
data_path <- args$data_path
sample_col <- args$sample_col
batch_col <- args$batch_col
run_soupx <- args$soupx
run_harmony <- args$harmony
harmony_col <- args$harmony_col
output <- args$output
metadata_path <- args$metadata_path

seu <- load_seurat(data_path, metadata_path = metadata_path)
## Initial Clustering

seu$lane <- sub(".*_", "", colnames(seu))
 
seu$run_lane <- sub("^[^_]+_", "", colnames(seu))

## END

seu <- set_metadata_group(seu, sample_col, "sample")

batch_values <- if ("batch" %in% colnames(seu[[]])) {
  seu$batch
} else {
  NULL
}

# Violin: counts + mt pre-QC
# statistics

plot_qc_prepost <- function(seu, prefix, group.by = sample_col) {
  # Violin: absolute QC metrics per sample
  p <- VlnPlot(
    seu,
    features = c("nCount_RNA", "nFeature_RNA"),
    group.by = group.by,
    pt.size = 0,
    ncol = 2
  )

  save_plot(p, paste0("QC_", prefix, "_violin_counts_features.png"), width = 10, height = 5)

  # Violin: relative QC metrics per sample
  p <- VlnPlot(
    seu,
    features = c("percent.mt"),
    group.by = group.by,
    pt.size = 0,
    ncol = 1
  )

  save_plot(p, paste0("QC_", prefix, "_violin_percent_mt_ribo_hb.png"), width = 10, height = 5)

  # Scatter: nCount vs mt
  p <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt") +
    ggtitle(paste0(prefix, ": nCount_RNA vs percent.mt"))
  save_plot(p, paste0("QC_", prefix, "_scatter_nCount_vs_percentMT.png"))
}

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

plot_qc_prepost(seu, "pre_QC", NULL)

# saveRDS(seu, file = "data/seu_before_filtering.rds")
# print("Saved seu_before_filtering.rds")
# seu_preQC <- readRDS("data/seu_before_filtering.rds")
seu_preQC <- seu

# qc start

is_outlier_mad <- function(x, nmads = 5, type = c("both", "lower", "upper")){
  type <- match.arg(type)

  med <- median(x, na.rm = TRUE)
  madv <- mad(x, constant = 1, na.rm = TRUE)
  if (madv == 0 || is.na(madv)) {
    return (rep(FALSE, length(x)))
  }

  lower <- med - nmads * madv
  upper <- med + nmads * madv

  if(type == "both") {
    return (x < lower | x > upper)
  } else if(type == "lower") {
    return (x < lower)
  } else {
    return (x > upper)
  }
}

is_outlier_mad_by_group <- function (x, group, nmads = 5, type = c("both", "lower", "upper")) {
  type <- match.arg(type)
  groups <- unique(group)

  out <- rep(FALSE, length(x))

  for(g in groups) {
    idx <- which(group == g & !is.na(x))
    if(length(idx) == 0) {
      next
    }

    med <- median(x[idx], na.rm = TRUE)
    madv <- mad(x[idx], constant = 1, na.rm = TRUE)
    if (madv == 0 || is.na(madv)) {
      next
    }
    lower <- med - nmads * madv
    upper <- med + nmads * madv

    if(type == "both") {
      out[idx] <- x[idx] < lower | x[idx] > upper
    } else if(type == "lower") {
      out[idx] <- x[idx] < lower
    } else {
      out[idx] <- x[idx] > upper
    }
  }

  return(out)
}

is_outlier <- function(x, batch = NULL, nmads = 5, type = c("both", "lower", "upper")) {
  if (!is.null(batch)) {
    return(is_outlier_mad_by_group(x, batch, nmads = nmads, type = type))
  } else {
    return(is_outlier_mad(x, nmads = nmads, type = type))
  }
}



out_low_counts <- is_outlier(log1p(seu$nCount_RNA), batch_values, nmads = 5, type = "lower")
out_low_genes <- is_outlier(log1p(seu$nFeature_RNA), batch_values, nmads = 5, type = "lower")

out_high_counts <- is_outlier(log1p(seu$nCount_RNA), batch_values, nmads = 5, type = "upper")
out_high_genes <- is_outlier(log1p(seu$nFeature_RNA), batch_values, nmads = 5, type = "upper")

out_high_mt <- is_outlier(seu$percent.mt, batch_values, nmads = 5, type = "upper")

seu$qc_outlier <- out_low_counts | out_low_genes | out_high_mt | out_high_counts | out_high_genes

seu <- subset(seu, subset = !qc_outlier)

plot_qc_prepost(seu, "post_QC")
dims.use <- 1:20
# cluster for soup
if (run_soupx) {

  tod <- counts
  tod <- as(tod, "dgCMatrix")
  tod <- tod[Matrix::rowSums(tod) > 0, , drop = FALSE] # parameterize if needed
  tod <- tod[, Matrix::colSums(tod) > 2, drop = FALSE] # parameterize if needed

  flt <- GetAssayData(seu, layer = "counts")
  flt <- as(flt, "dgCMatrix")

  common_genes <- intersect(rownames(flt), rownames(tod))
  flt <- flt[common_genes, , drop = FALSE]
  tod <- tod[common_genes, , drop = FALSE]


  seu <- NormalizeData(seu)

  sce <- as.SingleCellExperiment(seu)
  hvgs <- devianceFeatureSelection(sce)

  dev <- rowData(hvgs)$binomial_deviance
  names(dev) <- rownames(hvgs)

  hvgs_genes <- names(sort(dev, decreasing = TRUE))[1:2000]

  VariableFeatures(seu) <- hvgs_genes

  seu <- ScaleData(seu, features = hvgs_genes)
  seu <- RunPCA(seu, features = hvgs_genes)
  seu <- FindNeighbors(seu, dims = dims.use)
  seu <- FindClusters(seu, resolution = 0.5)

  common_cells <- intersect(colnames(flt), colnames(seu))
  flt <- flt[, common_cells, drop = FALSE]

  soupx_groups <- seu$seurat_clusters
  names(soupx_groups) <- colnames(seu)
  soupx_groups <- soupx_groups[common_cells]

  saveRDS(seu, file = "data/seu_before_soup.rds")
  print("Saved before soup")

  sc <- SoupChannel(
    tod,
    flt,
    calcSoupProfile = FALSE
  )

  gene_counts <- Matrix::rowSums(tod)
  gene_counts <- gene_counts[gene_counts > 0]

  soupProf <- data.frame(
    est = as.numeric(gene_counts / sum(gene_counts)),
    counts = as.numeric(gene_counts),
    row.names = names(gene_counts)
  )

  sc <- setSoupProfile(sc, soupProf)
  sc <- setClusters(sc, soupx_groups)
  sc <- autoEstCont(sc, doPlot = TRUE)

  corrected <- adjustCounts(
    sc,
    roundToInt = TRUE
  )

  seu <- CreateSeuratObject(
    counts = corrected,
    meta.data = seu@meta.data[colnames(corrected), , drop = FALSE]
  )

  batch_values <- if (is.null(batch)) {
    NULL
  } else {
    seu$batch
  }

  rm(tod, flt, soupx_groups, corrected, sc)
  gc()

  cnt <- GetAssayData(
    seu,
    assay = "RNA",
    layer = "counts"
  )

  keep <- Matrix::rowSums(cnt > 0) >= 20
  seu <- seu[keep, ]

  seu[["percent.mt"]] <- PercentageFeatureSet(
    seu,
    pattern = "^MT-"
  )

  out_low_counts <- is_outlier_mad_by_group(
    seu$nCount_RNA,
    batch_values,
    nmads = 5,
    type = "lower"
  )

  out_low_genes <- is_outlier_mad_by_group(
    seu$nFeature_RNA,
    batch_values,
    nmads = 5,
    type = "lower"
  )

  out_high_mt <- is_outlier_mad_by_group(
    seu$percent.mt,
    batch_values,
    nmads = 5,
    type = "upper"
  )

  seu$qc_outlier <-
    out_low_counts |
    out_low_genes |
    out_high_mt

  seu <- subset(
    seu,
    subset = !qc_outlier
  )

  saveRDS(
    seu,
    file = "data/seu_after_soup.rds"
  )

  print("Saved seu_after_soup.rds")
}

rm(counts)
gc()
# clustering

seu <- NormalizeData(seu, normalization.method = "LogNormalize")

sce <- as.SingleCellExperiment(seu)

hvgs_sce <- devianceFeatureSelection(sce)

dev <- rowData(hvgs_sce)$binomial_deviance
names(dev) <- rownames(hvgs_sce)

hvgs_genes <- names(sort(dev, decreasing = TRUE))[1:2000]

VariableFeatures(seu) <- hvgs_genes

seu <- ScaleData(seu, features = hvgs_genes) #, vars.to.regress = c("percent.mt"))

seu <- RunPCA(seu, features = hvgs_genes)

# Preliminary clustering for scDblFinder
seu <- FindNeighbors(seu, reduction = "pca", dims = dims.use)

seu <- FindClusters(seu, resolution = 0.5)

# Doublet detection
sce <- as.SingleCellExperiment(seu)

sce <- scDblFinder(sce, clusters = "seurat_clusters", samples = "run_lane")
seu$scDblFinder_class <- sce$scDblFinder.class
seu$scDblFinder_score <- sce$scDblFinder.score

# UMAP only for doublet visualization
seu <- RunUMAP(seu, reduction = "pca", dims = dims.use)

p <- DimPlot(
  seu,
  group.by = "scDblFinder_class",
  reduction = "umap"
) + ggtitle("scDblFinder: singlet/doublet")

save_plot(p, "UMAP_scDblFinder_class.png")

p <- FeaturePlot(
  seu,
  features = "scDblFinder_score",
  reduction = "umap"
) + ggtitle("scDblFinder score")

save_plot(p, "UMAP_scDblFinder_score.png")

table(seu$scDblFinder_class)
#saveRDS(seu, file = "data/seu_after_dd.rds")
print("Saved seu_after_dd.rds")
#print(
  #DimPlot(seu_sx, group.by = "scDblFinder_class", reduction = "umap")
#)

#VlnPlot(seu_sx, "scDblFinder_score", group.by = "scDblFinder_class", pt.size = 0.1)

#FeaturePlot(seu_sx, "scDblFinder_score", reduction = "umap")

seu <- subset(seu, subset = scDblFinder_class == "singlet")

# final clustering

seu <- NormalizeData(seu, normalization.method = "LogNormalize")

sce <- as.SingleCellExperiment(seu)

hvgs_sce <- devianceFeatureSelection(sce)

dev <- rowData(hvgs_sce)$binomial_deviance
names(dev) <- rownames(hvgs_sce)

hvgs_genes <- names(sort(dev, decreasing = TRUE))[1:2000]

VariableFeatures(seu) <- hvgs_genes

seu <- ScaleData(seu, features = hvgs_genes)

seu <- RunPCA(seu, features = hvgs_genes)
if (run_harmony) {
  seu <- RunHarmony(seu, group.by.vars = harmony_col)
  reduction_use <- "harmony"
} else {
  reduction_use <- "pca"
}

seu <- FindNeighbors(seu, reduction = reduction_use, dims = dims.use)

seu <- FindClusters(seu, resolution = c(0.5, 1.0, 1.5, 2.0))

seu <- RunUMAP(seu, reduction = reduction_use, dims = dims.use)

dev_df <- data.frame(deviance = dev)
p <- ggplot(dev_df, aes(x = deviance)) +
  geom_histogram(bins = 50) +
  ggtitle("Distribution of binomial deviance (all genes)")

save_plot(p, "HVG_deviance_distribution.png")

dev_sorted <- sort(dev, decreasing = TRUE)

df_rank <- data.frame(
  rank = seq_along(dev_sorted),
  deviance = as.numeric(dev_sorted),
  HVG = seq_along(dev_sorted) <= 2000
)

p <- ggplot(df_rank, aes(x = rank, y = deviance)) +
  geom_point(aes(alpha = HVG), size = 0.6) +
  ggtitle("Deviance feature selection: rank vs deviance") +
  xlab("Gene rank (1 = highest deviance)") +
  ylab("Binomial deviance")

save_plot(p, "HVG_deviance_rankplot.png")

top_n <- 10
top_genes <- names(dev_sorted)[1:top_n]

df_rank$gene <- names(dev_sorted)
df_rank$is_top10 <- df_rank$gene %in% top_genes

p <- ggplot(df_rank, aes(x = rank, y = deviance, color = HVG)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = c("grey70", "red")) +
  geom_vline(xintercept = 2000, linetype = "dashed") +
  geom_text_repel(
    data = subset(df_rank, is_top10),
    aes(label = gene),
    size = 3,
    max.overlaps = Inf
  ) +
  ggtitle("HVG selection by deviance (top 2000 highlighted)") +
  xlab("Gene rank") +
  ylab("Binomial deviance")

save_plot(p, "HVG_deviance_rankplot_highlight_top10.png")

top_n <- 20
top_genes <- head(names(sort(dev, decreasing = TRUE)), top_n)

df_top <- data.frame(
  gene = names(dev),
  deviance = as.numeric(dev),
  is_top = names(dev) %in% top_genes
)

p <- ggplot(df_top, aes(x = deviance)) +
  geom_histogram(bins = 80) +
  ggtitle("Deviance distribution (top genes flagged)")

save_plot(p, "HVG_deviance_hist_topflag.png")

p <- clustree(seu, prefix = "RNA_snn_res.")
save_plot(p, "clustree_resolutions.png", width = 8, height = 8)

p <- DimPlot(seu, group.by = "RNA_snn_res.0.5",
             reduction = "umap", label = TRUE) +
  ggtitle("Resolution 0.5")

save_plot(p, "UMAP_res_0.5.png")

p <- DimPlot(seu, group.by = "RNA_snn_res.1",
             reduction = "umap", label = TRUE) +
  ggtitle("Resolution 1.0")

save_plot(p, "UMAP_res_1.png")


p <- DimPlot(seu, group.by = "RNA_snn_res.1.5",
             reduction = "umap", label = TRUE) +
  ggtitle("Resolution 1.5")

save_plot(p, "UMAP_res_1.5.png")

p <- DimPlot(seu, group.by = "RNA_snn_res.2",
             reduction = "umap", label = TRUE) +
  ggtitle("Resolution 2.0")

save_plot(p, "UMAP_res_2.png")

p <- DimPlot(
  seu,
  group.by = "celltype.l2",
  reduction = "umap",
  label = TRUE,
  repel = TRUE
) + ggtitle("Cell type annotation (level 2)")

save_plot(p, "UMAP_celltype_L2.png")

p <- DimPlot(
  seu,
  group.by = "celltype.l1",
  reduction = "umap",
  label = TRUE
) + ggtitle("Cell classes (level 1)")

save_plot(p, "UMAP_celltype_L1.png")

p <- DimPlot(
  seu,
  group.by = "celltype.l3",
  reduction = "umap",
  label = FALSE
) + ggtitle("Cell subtypes (level 3)")

save_plot(p, "UMAP_celltype_L3.png")

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
p <- VlnPlot(
  seu,
  features = c("nCount_RNA", "percent.mt"),
  group.by = sample_col,
  pt.size = 0,
  ncol = 2
) + ggtitle("Post-QC: nCount_RNA & percent_final.mt per sample")

save_plot(p, "QC_post_violin_counts_mt.png", width = 10, height = 5)

p <- VlnPlot(
  seu,
  features = "percent.mt",
  group.by = sample_col,
  pt.size = 0,
  ncol = 1
) + ggtitle("Post-QC: percent_final.mt per sample")

save_plot(p, "QC_post_violin_mt.png", width = 10, height = 5)

seu_preQC[["percent.mt"]] <- PercentageFeatureSet(seu_preQC, pattern = "^MT-")
p <- VlnPlot(
  seu_preQC,
  features = "percent.mt",
  group.by = sample_col,
  pt.size = 0,
  ncol = 1
) + ggtitle("Pre-QC: percent_final.mt per sample")

save_plot(p, "QC_pre_violin_mt_10.png", width = 10, height = 5)

seu_preQC <- subset(
  seu_preQC,
  subset = percent.mt < 10 &
           percent.hb < 5
)

p <- VlnPlot(
  seu_preQC,
  features = "nCount_RNA",
  group.by = sample_col,
  pt.size = 0,
  ncol = 1
) + ggtitle("Pre-QC: nCount_RNA")

save_plot(p, "QC_pre_violin_count.png", width = 10, height = 5)

saveRDS(seu, file = output)
print("Saved seu")
# seu_sx <- readRDS("data/seu_sx_final_new.rds")

