  library(Seurat)
  library(SeuratObject)
  library(ADTnorm)
  library(Matrix)

  script_file <- sub(
    "^--file=",
    "",
    grep("^--file=", commandArgs(), value = TRUE)[[1]]
  )
  GATING_ROOT <- dirname(dirname(normalizePath(script_file, mustWork = TRUE)))
  source(file.path(GATING_ROOT, "config", "datasets.R"))

  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) < 2) {
    stop(
      "Usage: Rscript 03_normalize_adt.R <dataset_name> <dataset_path> ",
      "[output_root]"
    )
  }

  name <- args[[1]]
  dataset_path <- args[[2]]
  output_root <- if (length(args) >= 3) {
    args[[3]]
  } else {
    "PLA_QC_ADTnorm_tuned"
  }

  OUT_DIR <- file.path(output_root, name)
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)


  extract_assay_layer_as_matrix <- function(
    dataset,
    assay_name = "ADT",
    layer_name = "counts"
  ) {
    # Check assay
    if (!assay_name %in% Assays(dataset)) {
      stop(paste0("Assay not found: ", assay_name))
    }

    assay_obj <- dataset[[assay_name]]

    # Show available layers
    available_layers <- SeuratObject::Layers(assay_obj)

    if (!layer_name %in% available_layers) {
      stop(
        paste0(
          "Layer not found: ", layer_name,
          ". Available layers: ",
          paste(available_layers, collapse = ", ")
        )
      )
    }

    # Extract layer from assay object
    mat <- SeuratObject::LayerData(
      object = assay_obj,
      layer = layer_name
    )

    message("Raw extracted object class: ", paste(class(mat), collapse = ", "))

    mat <- as.matrix(mat)

    message("Converted matrix class: ", paste(class(mat), collapse = ", "))
    message("Converted matrix dim: ", paste(dim(mat), collapse = " x "))

    if (!all(colnames(mat) %in% rownames(dataset@meta.data))) {
      stop("Not all matrix column names are present in metadata rownames.")
    }

    return(mat)
  }


  run_adtnorm_on_seurat <- function(
    dataset,
    sample_col,
    batch_col = NULL,
    assay_in = "ADT",
    assay_out = "ADT_ADTnorm",
    normalization_groups,
    marker_order = NULL,
    out_dir,
    run_name = "ADTnorm_run",
    save_fig = TRUE,
    ...
  ) {

    # Check sample column
    if (!sample_col %in% colnames(dataset@meta.data)) {
      stop(paste0("sample_col not found in metadata: ", sample_col))
    }

    # Check batch column if provided
    if (!is.null(batch_col) && !batch_col %in% colnames(dataset@meta.data)) {
      stop(paste0("batch_col not found in metadata: ", batch_col))
    }

    # Extract ADT counts as plain matrix
    adt_counts <- extract_assay_layer_as_matrix(
      dataset = dataset,
      assay_name = assay_in,
      layer_name = "counts"
    )

    if (nrow(adt_counts) == 0 || ncol(adt_counts) == 0) {
      stop("ADT counts matrix is empty.")
    }

    if (
      length(normalization_groups) == 0 ||
      is.null(names(normalization_groups)) ||
      any(names(normalization_groups) == "")
    ) {
      stop("normalization_groups must be a named, non-empty list.")
    }

    valid_alignment_types <- c(
      "negPeak",
      "negPeak_valley",
      "negPeak_valley_posPeak",
      "valley"
    )
    invalid_alignment_types <- setdiff(
      names(normalization_groups),
      valid_alignment_types
    )

    if (length(invalid_alignment_types) > 0) {
      stop(
        "Unsupported landmark alignment type(s): ",
        paste(invalid_alignment_types, collapse = ", ")
      )
    }

    normalization_groups <- lapply(
      normalization_groups,
      function(markers) unique(as.character(markers))
    )

    assigned_markers <- unlist(
      normalization_groups,
      use.names = FALSE
    )

    duplicated_markers <- unique(
      assigned_markers[duplicated(assigned_markers)]
    )
    if (length(duplicated_markers) > 0) {
      stop(
        "Markers assigned to more than one normalization group: ",
        paste(duplicated_markers, collapse = ", ")
      )
    }

    missing_markers <- setdiff(assigned_markers, rownames(adt_counts))
    if (length(missing_markers) > 0) {
      stop(
        "Planned ADT marker(s) not found in the input assay: ",
        paste(missing_markers, collapse = ", ")
      )
    }

    if (is.null(marker_order)) {
      marker_to_process <- assigned_markers
    } else {
      marker_order <- unique(as.character(marker_order))
      unplanned_markers <- setdiff(marker_order, assigned_markers)
      if (length(unplanned_markers) > 0) {
        stop(
          "Selected marker(s) missing from the normalization plan: ",
          paste(unplanned_markers, collapse = ", ")
        )
      }
      marker_to_process <- marker_order[marker_order %in% assigned_markers]
    }

    if (length(marker_to_process) == 0) {
      stop("The normalization plan does not contain any ADT markers.")
    }

    # ADTnorm expects cells x markers.
    cell_x_adt <- t(adt_counts[marker_to_process, , drop = FALSE])

    # Build ADTnorm feature table
    cell_x_feature_raw <- dataset@meta.data[colnames(adt_counts), , drop = FALSE]

    cell_x_feature <- data.frame(
      sample = as.factor(cell_x_feature_raw[[sample_col]]),
      row.names = rownames(cell_x_feature_raw)
    )

    if (!is.null(batch_col)) {
      cell_x_feature$batch <- as.factor(cell_x_feature_raw[[batch_col]])
    } else {
      # If no batch is provided, use sample as batch
      cell_x_feature$batch <- cell_x_feature$sample
    }

    # Make output directory
    adtnorm_out_dir <- file.path(out_dir, paste0(run_name, "_ADTnorm"))

    # Run ADTnorm separately for every marker group so that each group can
    # use the biologically appropriate landmark alignment.
    normalized_parts <- list()
    group_out_dirs <- list()

    for (alignment_type in names(normalization_groups)) {
      markers_this_run <- normalization_groups[[alignment_type]]
      group_out_dir <- file.path(adtnorm_out_dir, alignment_type)
      group_run_name <- paste(run_name, alignment_type, sep = "_")
      dir.create(group_out_dir, recursive = TRUE, showWarnings = FALSE)

      message(
        "Running ADTnorm with landmark_align_type='",
        alignment_type,
        "' for: ",
        paste(markers_this_run, collapse = ", ")
      )

      group_norm <- ADTnorm::ADTnorm(
        cell_x_adt = cell_x_adt[, markers_this_run, drop = FALSE],
        cell_x_feature = cell_x_feature,
        save_outpath = group_out_dir,
        study_name = group_run_name,
        marker_to_process = markers_this_run,
        save_fig = save_fig,
        landmark_align_type = alignment_type,
        ...
      )

      group_norm <- as.matrix(group_norm)

      if (
        !all(rownames(cell_x_adt) %in% rownames(group_norm)) ||
        !all(markers_this_run %in% colnames(group_norm))
      ) {
        stop(
          "ADTnorm output for alignment group '",
          alignment_type,
          "' has unexpected cell or marker names."
        )
      }

      normalized_parts[[alignment_type]] <- group_norm[
        rownames(cell_x_adt),
        markers_this_run,
        drop = FALSE
      ]
      group_out_dirs[[alignment_type]] <- group_out_dir
    }

    adt_norm <- do.call(cbind, normalized_parts)

    # Ensure cell and marker order
    adt_norm <- adt_norm[colnames(adt_counts), marker_to_process, drop = FALSE]

    # Convert back to Seurat format: features x cells
    adt_norm_t <- t(adt_norm)
    adt_norm_t <- as(adt_norm_t, "dgCMatrix")

    # Store ADTnorm output as normalized data layer in a new assay
    dataset[[assay_out]] <- SeuratObject::CreateAssay5Object(data = adt_norm_t)

    # Store run info
    dataset@misc$ADTnorm <- list(
      assay_in = assay_in,
      assay_out = assay_out,
      sample_col = sample_col,
      batch_col = batch_col,
      marker_to_process = marker_to_process,
      normalization_groups = normalization_groups,
      run_name = run_name,
      out_dir = adtnorm_out_dir,
      group_out_dirs = group_out_dirs
    )

    DefaultAssay(dataset) <- assay_out

    return(dataset)
  }

  dataset_config <- get_pla_dataset_config(name)
  plan <- dataset_config$normalization_groups

  if (is.null(plan)) {
    message(
      "Skipping ADTnorm for '",
      name,
      "' because its current density plots do not contain reliable ",
      "landmarks. No normalized RDS was written."
    )
    quit(save = "no", status = 0)
  }

  markers <- dataset_config$markers

  plan_table <- data.frame(
    marker = unlist(plan, use.names = FALSE),
    landmark_align_type = rep(names(plan), lengths(plan)),
    stringsAsFactors = FALSE
  )
  write.csv(
    plan_table,
    file = file.path(OUT_DIR, "normalization_plan.csv"),
    row.names = FALSE
  )

  message("Reading dataset: ", dataset_path)
  dataset <- readRDS(dataset_path)

  dataset <- run_adtnorm_on_seurat(
    dataset = dataset,
    sample_col = dataset_config$sample_col,
    batch_col = NULL,
    assay_in = "ADT",
    assay_out = paste0("ADT_ADTnorm_", name),
    normalization_groups = plan,
    marker_order = markers,
    out_dir = OUT_DIR,
    run_name = paste0(name, "_ADT"),
    save_fig = TRUE,
    save_landmark = TRUE,
    peak_type = "mode",
    exclude_zeroes = FALSE,
    lower_peak_thres = 0.001,
    quantile_clip = 1
  )

  saveRDS(
    dataset,
    file = file.path(OUT_DIR, paste0(name, "_ADTnorm_seurat.rds"))
  )
