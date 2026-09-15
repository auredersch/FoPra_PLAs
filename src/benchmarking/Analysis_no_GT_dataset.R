library(dplyr)
library(tidyr)
library(Seurat)
library(edgeR)
library(ggplot2)
library(gprofiler2)

 plot_threshold_check <- function(pair_table) {
    threshold_check <- lapply(c(1, 3, 5, 10, 20), function(x) {

        pair_table %>%
            mutate(
                valid_pair =
                    PLA >= x &
                    `platelet-free` >= x
            ) %>%
            group_by(cell_type_lowerres, timepoint.final) %>%
            summarise(
                min_cells = x,
                n_pairs = sum(valid_pair),
                .groups = "drop"
            )
    }) %>%
        bind_rows()

    m <- threshold_check %>%
        filter(min_cells %in% c(3, 5, 10)) %>%
        arrange(
            cell_type_lowerres,
            timepoint.final,
            min_cells
        )

    p <- ggplot(
        threshold_check,
        aes(
            x = min_cells,
            y = n_pairs,
            group = cell_type_lowerres
        )
    ) +
        geom_line() +
        geom_point() +
        facet_wrap(~ cell_type_lowerres, scales = "free_y") +
        labs(
            x = "Minimum cells per group",
            y = "Number of paired donors"
        ) +
        theme_bw()

    ggsave(
    p,
    filename = file.path(
        output_path,
        "paired_donors_per_celltype_min_cells.png"
    ),
    width = 10,
    height = 6,
    dpi = 300
)

    return(p)
}


dataset_root <- "/nfs/home/students/f.mathis/FoPra_PLAs/data/datasets/benchmarked_objects/stemi_raw"
relative_dataset_path <- "AUCell_ExtFALSE_gmm_dist_dual_Scopeglobal/GOBP_REG/pbmc_classified_stemi_raw_AUCell_GOBP_REG_gmm_dist_dual_Scopeglobal.rds"
dataset_name <- "stemi_raw"
output_path <- file.path("/nfs/home/students/f.mathis/FoPra_PLAs/results/benchmarking", dataset_name)

dir.create(output_path, recursive = TRUE, showWarnings = FALSE)


pbmc <- readRDS(file.path(dataset_root, relative_dataset_path))

pb <- AggregateExpression(
    pbmc,
    assays = "RNA",
    group.by = c(
        "assignment.final",
        "timepoint.final",
        "cell_type_lowerres",
        "Prediction"
    ),
    return.seurat = TRUE,
    verbose = FALSE
)

PLA_LABEL <- "PLA"
PF_LABEL <- "platelet-free"

min_cells <- 5

pair_table <- pbmc@meta.data %>%
    filter(
        Prediction %in% c("PLA", "platelet-free"),
        !is.na(cell_type_lowerres),
        !is.na(timepoint.final),
        !is.na(assignment.final)
    ) %>%
    count(
        cell_type_lowerres,
        timepoint.final,
        assignment.final,
        Prediction,
        name = "n_cells"
    ) %>%
    pivot_wider(
        names_from = Prediction,
        values_from = n_cells,
        values_fill = 0
    )

pair_table <- pair_table %>%
    mutate(
        valid_pair =
            PLA >= min_cells &
            `platelet-free` >= min_cells
    )

pair_summary <- pair_table %>%
    group_by(
        cell_type_lowerres,
        timepoint.final
    ) %>%
    summarise(
        n_paired_donors = sum(valid_pair),
        .groups = "drop"
    ) %>%
    arrange(
        timepoint.final,
        desc(n_paired_donors)
    )


p <- plot_threshold_check(pair_table)

run_paired_edger <- function(object, ct, tp, min_cells = 5) {

    meta <- object@meta.data
    meta$cell <- rownames(meta)

    meta <- meta %>%
        dplyr::filter(
            cell_type_lowerres == ct,
            timepoint.final == tp,
            Prediction %in% c("PLA", "platelet-free")
        )

    donors <- meta %>%
        dplyr::count(assignment.final, Prediction) %>%
        tidyr::pivot_wider(
            names_from = Prediction,
            values_from = n,
            values_fill = 0
        ) %>%
        dplyr::filter(
            PLA >= min_cells,
            `platelet-free` >= min_cells
        ) %>%
        dplyr::pull(assignment.final)

    obj <- subset(
        object,
        cells = meta$cell[meta$assignment.final %in% donors]
    )

    obj$pb_id <- paste(
        obj$assignment.final,
        obj$Prediction,
        sep = "."
    )

    counts <- AggregateExpression(
        obj,
        assays = "RNA",
        group.by = "pb_id",
        return.seurat = FALSE
    )$RNA

    pb_meta <- obj@meta.data %>%
        dplyr::distinct(
            pb_id,
            assignment.final,
            Prediction
        )

    idx <- match(
        colnames(counts),
        pb_meta$pb_id
    )

    pb_meta <- pb_meta[idx, ]

    pb_meta$Prediction <- factor(
        pb_meta$Prediction,
        levels = c("platelet-free", "PLA")
    )

    design <- model.matrix(
        ~ assignment.final + Prediction,
        data = pb_meta
    )

    y <- DGEList(counts)
    keep <- filterByExpr(y, design)
    y <- y[keep, , keep.lib.sizes = FALSE]
    y <- normLibSizes(y)

    fit <- glmQLFit(y, design, robust = TRUE)
    test <- glmQLFTest(fit, coef = "PredictionPLA")

    res <- edgeR::topTags(test, n = Inf)$table
    res$gene <- rownames(res)

    res <- res %>%
    dplyr::filter(
        is.finite(F),
        is.finite(logFC)
    ) %>%
    dplyr::mutate(
        rank_score = sign(logFC) * sqrt(F)
    ) %>%
    dplyr::arrange(desc(rank_score))


    res
}

de <- run_paired_edger(
    pbmc,
    ct = "monocyte",
    tp = "UT",
    min_cells = 5
)

head(de)

ranked_genes <- de$gene    

gsea_PLA <- gost(
    query = ranked_genes,
    organism = "hsapiens",
    ordered_query = TRUE,
    sources = c("GO:BP", "REAC", "KEGG")
)

ranked_PF <- de %>%
    arrange(rank_score) %>%
    pull(gene)

gsea_PF <- gost(
    query = ranked_PF,
    organism = "hsapiens",
    ordered_query = TRUE,
    sources = c("GO:BP", "REAC", "KEGG")
)


pla_terms <- gsea_PLA$result %>%
    dplyr::mutate(direction = "PLA")

pf_terms <- gsea_PF$result %>%
    dplyr::mutate(direction = "platelet-free")

gsea_terms <- dplyr::bind_rows(
    pla_terms,
    pf_terms
)




# --------- PLOTS ---------

top_terms <- gsea_terms %>%
    dplyr::group_by(direction) %>%
    dplyr::slice_min(
        order_by = p_value,
        n = 10,
        with_ties = FALSE
    ) %>%
    dplyr::ungroup()

p <- ggplot(
    top_terms,
    aes(
        x = -log10(p_value),
        y = reorder(term_name, -log10(p_value)),
        size = intersection_size
    )
) +
    geom_point() +
    facet_wrap(
        ~ direction,
        scales = "free_y"
    ) +
    labs(
        x = "-log10(p-value)",
        y = NULL,
        size = "Genes"
    ) +
    theme_bw()

ggsave(
    p,
    filename = file.path(output_path, "gsea_top_terms.png"),
    width = 10,
    height = 6,
    dpi = 300
)




## Grid

analysis_grid <- pair_summary %>%
    dplyr::filter(n_paired_donors >= 5)

DE_results <- list()
GSEA_results <- list()

for (i in seq_len(nrow(analysis_grid))) {

    ct <- analysis_grid$cell_type_lowerres[i]
    tp <- analysis_grid$timepoint.final[i]
    n_pairs <- analysis_grid$n_paired_donors[i]

    message("Running: ", ct, " | ", tp)

    # Run paired pseudobulk DE
    de <- tryCatch(
        run_paired_edger(
            pbmc,
            ct = ct,
            tp = tp,
            min_cells = 5
        ),
        error = function(e) {
            message("DE failed: ", conditionMessage(e))
            return(NULL)
        }
    )

    if (is.null(de)) {
        next
    }

    # Keep genes with valid ranking statistics
    de <- de %>%
        dplyr::filter(is.finite(rank_score)) %>%
        dplyr::mutate(
            cell_type_lowerres = ct,
            timepoint.final = tp,
            n_pairs = n_pairs
        )

    key <- paste(ct, tp, sep = "__")

    DE_results[[key]] <- de

    # Generate both ranking directions
    rankings <- list(
        PLA = de %>%
            dplyr::arrange(dplyr::desc(rank_score)) %>%
            dplyr::pull(gene),

        `platelet-free` = de %>%
            dplyr::arrange(rank_score) %>%
            dplyr::pull(gene)
    )

    for (direction in names(rankings)) {

        gp <- tryCatch(
            gprofiler2::gost(
                query = rankings[[direction]],
                organism = "hsapiens",
                ordered_query = TRUE,
                significant = TRUE,
                user_threshold = 0.05,
                correction_method = "g_SCS",
                sources = c("GO:BP", "REAC", "KEGG"),
                highlight = TRUE
            ),
            error = function(e) {
                message("g:Profiler failed: ", conditionMessage(e))
                return(NULL)
            }
        )

        if (!is.null(gp) && !is.null(gp$result)) {

            GSEA_results[[paste(key, direction, sep = "__")]] <-
                gp$result %>%
                dplyr::mutate(
                    cell_type_lowerres = ct,
                    timepoint.final = tp,
                    direction = direction,
                    n_pairs = n_pairs
                )
        }
    }
}

DE_all <- dplyr::bind_rows(DE_results)
GSEA_all <- dplyr::bind_rows(GSEA_results)

pathway_counts <- GSEA_all %>%
    dplyr::filter(
        source %in% c("REAC", "KEGG")
    ) %>%
    dplyr::count(
        cell_type_lowerres,
        timepoint.final,
        direction,
        name = "n_pathways"
    )

pathway_balance <- pathway_counts %>%
    tidyr::pivot_wider(
        names_from = direction,
        values_from = n_pathways,
        values_fill = 0
    ) %>%
    dplyr::mutate(
        difference = PLA - `platelet-free`,
        log2_ratio = log2(
            (PLA + 1) /
            (`platelet-free` + 1)
        )
    ) %>%
    dplyr::arrange(dplyr::desc(difference))

pathway_balance