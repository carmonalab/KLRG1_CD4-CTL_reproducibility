# scRNA-seq/utils.R
#
# Utility functions for scRNA-seq analysis scripts.
# Source this file in your Rmd preamble:
#   source(here::here("scRNA-seq", "utils.R"))

# cd4_annotation_palette -------------------------------------------------------
# Named color palette for the five canonical CD4 T cell states plus Unassigned.
# Used as the default fill scale in plot_pbmc_tiln_ratio() and dotplot_seurat_style().

cd4_annotation_palette <- c(
  "Naive-like" = "#A6CEE3",
  "Memory"     = "purple2",
  "Treg"       = "#33A02C",
  "TFH"        = "#FB9A99",
  "Cytotoxic"  = "#FF7F00",
  "Unassigned" = "grey70"
)

# dotplot_seurat_style ---------------------------------------------------------
#
# A custom dotplot function that provides Seurat-style marker visualization
# with optional splitting by a grouping variable (e.g., tissue, patient).
#
# Parameters:
#   object: Seurat object
#   features: character vector of gene names to plot
#   ident_col: metadata column defining cell groups (x-axis)
#   split_by: optional metadata column for visual separation (panels on x)
#   subset_col: optional metadata column for subsetting cells
#   subset_values: values to filter by (used with subset_col)
#   assay: which assay to use (default: "RNA")
#   layer: which layer to use (default: "data")
#   scale_expression: if TRUE (default), color scale is centered and scaled
#   low_color, high_color: gradient colors for expression scale
#   dot_scale: max dot size (default: 7)
#
# Returns: ggplot object

dotplot_seurat_style <- function(
    object,
    features,
    ident_col = "annotation",
    split_by = "tissue",
    subset_col = NULL,
    subset_values = NULL,
    assay = "RNA",
    layer = "data",
    scale_expression = TRUE,
    low_color = "gray80",
    high_color = "red3",
    dot_scale = 7
) {
  stopifnot(inherits(object, "Seurat"))

  md <- object@meta.data
  if (!ident_col %in% colnames(md)) {
    stop(paste0("ident_col not found in metadata: ", ident_col))
  }

  if (!is.null(split_by) && !split_by %in% colnames(md)) {
    stop(paste0("split_by not found in metadata: ", split_by))
  }

  if (!is.null(subset_col) && !subset_col %in% colnames(md)) {
    stop(paste0("subset_col not found in metadata: ", subset_col))
  }

  if (is.null(subset_col) != is.null(subset_values)) {
    stop("subset_col and subset_values must be provided together.")
  }

  if (!is.null(subset_col)) {
    keep_cells <- rownames(md)[md[[subset_col]] %in% subset_values]
    if (length(keep_cells) == 0) {
      stop("No cells found after applying subset_col/subset_values.")
    }
    object <- subset(object, cells = keep_cells)
    md <- object@meta.data
  }

  if (!is.null(assay)) {
    DefaultAssay(object) <- assay
  }

  feat_found <- features[features %in% rownames(object)]
  if (length(feat_found) == 0) {
    stop("None of the requested features were found in the object.")
  }

  if (length(feat_found) < length(features)) {
    missing_feats <- setdiff(features, feat_found)
    message("Skipping missing features: ", paste(missing_feats, collapse = ", "))
  }

  vars_to_fetch <- c(ident_col, split_by, feat_found) %>%
    unique() %>%
    .[!is.na(.)]

  expr_df <- Seurat::FetchData(object = object, vars = vars_to_fetch, layer = layer)

  plot_df <- expr_df %>%
    tibble::rownames_to_column("cell") %>%
    tidyr::pivot_longer(
      cols = tidyr::all_of(feat_found),
      names_to = "feature",
      values_to = "expr"
    ) %>%
    dplyr::mutate(
      ident = .data[[ident_col]],
      split = if (!is.null(split_by)) .data[[split_by]] else "all"
    ) %>%
    dplyr::group_by(ident, split, feature) %>%
    dplyr::summarise(
      pct_exp = mean(expr > 0) * 100,
      avg_exp = mean(expr),
      .groups = "drop"
    ) %>%
    dplyr::group_by(feature) %>%
    dplyr::mutate(
      avg_exp_scaled = as.numeric(scale(avg_exp)),
      avg_exp_scaled = ifelse(is.na(avg_exp_scaled), 0, avg_exp_scaled)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      feature = factor(feature, levels = rev(feat_found)),
      ident = droplevels(factor(ident)),
      split = droplevels(factor(split)),
      color_value = if (scale_expression) avg_exp_scaled else avg_exp
    )

  # Rebuild ident levels from observed values only to avoid spacing artifacts.
  ident_levels <- unique(as.character(plot_df$ident))
  plot_df <- plot_df %>%
    dplyr::mutate(
      ident = factor(as.character(ident), levels = ident_levels),
      ident_idx = as.numeric(ident)
    )

  if (!is.null(split_by)) {
    n_split <- nlevels(plot_df$split)
    split_levels <- levels(plot_df$split)
    split_offsets <- seq_len(n_split) - (n_split + 1) / 2
    names(split_offsets) <- split_levels

    plot_df <- plot_df %>%
      dplyr::mutate(
        split_idx = as.numeric(split),
        x_pos = ident_idx + split_offsets[as.character(split)] * 0.28
      )

    split_breaks <- as.vector(
      sapply(seq_len(nlevels(plot_df$ident)), function(i) i + split_offsets * 0.28)
    )
    split_labels <- rep(split_levels, times = nlevels(plot_df$ident))
  } else {
    plot_df <- plot_df %>%
      dplyr::mutate(x_pos = ident_idx)
  }

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = x_pos, y = feature)) +
    ggplot2::geom_point(
      ggplot2::aes(size = pct_exp, color = color_value)
    ) +
    ggplot2::scale_size(range = c(0, dot_scale), name = "% Percent\nExpressed") +
    ggplot2::scale_color_gradient(
      low = low_color,
      high = high_color,
      name = if (scale_expression) "Average\nExpression\n(scaled)" else "Average\nExpression"
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 18),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (!is.null(split_by)) {
    p <- p +
      ggplot2::geom_vline(
        xintercept = seq_len(nlevels(plot_df$ident)) + 0.5,
        color = "gray90",
        linewidth = 0.3
      ) +
      ggplot2::scale_x_continuous(
        breaks = split_breaks,
        labels = split_labels,
        sec.axis = ggplot2::dup_axis(
          breaks = seq_len(nlevels(plot_df$ident)),
          labels = levels(plot_df$ident),
          name = ident_col
        )
      ) +
      ggplot2::labs(x = split_by) +
      ggplot2::theme(
        axis.title.x.top = ggplot2::element_text(margin = ggplot2::margin(b = 6)),
        axis.text.x.top = ggplot2::element_text(angle = 0, hjust = 0.5),
        axis.text.x = ggplot2::element_text(angle = 45,
                                   hjust = 1,
                                   size = 18)
      )
  } else {
    p <- p +
      ggplot2::scale_x_continuous(
        limits = c(0.5, nlevels(plot_df$ident) + 0.5),
        breaks = seq_len(nlevels(plot_df$ident)),
        labels = levels(plot_df$ident),
        expand = ggplot2::expansion(mult = c(0, 0))
      ) +
      ggplot2::labs(x = ident_col)
  }

  return(p)
}

# compute_pbmc_tiln_ratio ------------------------------------------------------
#
# Compute per-patient log2(tissue1 / tissue2) frequency ratios for every
# cell-type annotation state.
#
# Parameters:
#   metadata:       data frame (e.g. seurat_obj@meta.data)
#   annotation_col: column name holding cell-type labels
#   patient_col:    column name holding patient/donor IDs
#   tissue_col:     column name holding tissue labels
#   tissue_levels:  length-2 vector; ratio = (tissue_levels[1]+1)/(tissue_levels[2]+1)
#   min_cells:      minimum total cells (t1 + t2) for a patient x state
#                   combination to be retained
#
# Returns: data frame with columns for patient, annotation, both tissue counts,
#          ratio, and ratio_log2

compute_pbmc_tiln_ratio <- function(
    metadata,
    annotation_col = "annotation",
    patient_col    = "patient",
    tissue_col     = "tissue",
    tissue_levels  = c("PBMC", "TILN"),
    min_cells      = 5
) {
  stopifnot(length(tissue_levels) == 2)

  ann_levels <- levels(metadata[[annotation_col]])
  if (is.null(ann_levels)) {
    ann_levels <- sort(unique(as.character(metadata[[annotation_col]])))
  }

  # Rename to fixed internal names to keep tidyr::complete() simple
  md <- metadata %>%
    dplyr::select(
      patient = dplyr::all_of(patient_col),
      ann     = dplyr::all_of(annotation_col),
      tissue  = dplyr::all_of(tissue_col)
    ) %>%
    dplyr::filter(tissue %in% tissue_levels) %>%
    dplyr::mutate(ann = factor(ann, levels = ann_levels))

  patients <- unique(md$patient)
  t1 <- tissue_levels[1]
  t2 <- tissue_levels[2]

  rat <- md %>%
    dplyr::count(patient, ann, tissue, name = "ncells") %>%
    tidyr::complete(
      patient = patients,
      ann     = factor(ann_levels, levels = ann_levels),
      tissue  = tissue_levels,
      fill    = list(ncells = 0)
    ) %>%
    tidyr::pivot_wider(
      names_from  = tissue,
      values_from = ncells,
      values_fill = 0
    )

  # Ensure both tissue columns exist even if one had zero cells
  for (tis in tissue_levels) {
    if (!tis %in% colnames(rat)) rat[[tis]] <- 0L
  }

  rat <- rat %>%
    dplyr::filter((.data[[t1]] + .data[[t2]]) >= min_cells) %>%
    dplyr::mutate(
      ratio      = (.data[[t1]] + 1) / (.data[[t2]] + 1),
      ratio_log2 = log2(ratio)
    )

  # Restore original column names
  names(rat)[names(rat) == "ann"]     <- annotation_col
  names(rat)[names(rat) == "patient"] <- patient_col

  return(rat)
}

# plot_pbmc_tiln_ratio ---------------------------------------------------------
#
# Boxplot of log2(tissue1 / tissue2) ratios per annotation state.
#
# Parameters:
#   ratio_df:       output of compute_pbmc_tiln_ratio()
#   annotation_col: column name used as x-axis grouping
#   palette:        optional named character vector of fill colors
#   tissue_levels:  used to build the y-axis label (must match what was used
#                   in compute_pbmc_tiln_ratio)
#   title:          optional plot title
#
# Returns: ggplot object

plot_pbmc_tiln_ratio <- function(
    ratio_df,
    annotation_col = "annotation",
    palette        = NULL,
    tissue_levels  = c("PBMC", "TILN"),
    title          = NULL
) {
  y_label <- paste0(
    "log2 ratio (", tissue_levels[1], " / ", tissue_levels[2], ")\nFrequencies per patient"
  )

  p <- ggplot2::ggplot(
    ratio_df,
    ggplot2::aes(
      x    = .data[[annotation_col]],
      y    = ratio_log2,
      fill = .data[[annotation_col]]
    )
  ) +
    ggplot2::geom_boxplot(
      outlier.color = NA,
      show.legend   = FALSE,
      alpha         = 0.4
    ) +
    ggplot2::geom_jitter(
      color       = "grey22",
      show.legend = FALSE,
      shape       = 21,
      width       = 0.1
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      color      = "grey40",
      linetype   = "dashed"
    ) +
    ggplot2::labs(
      x     = NULL,
      y     = y_label,
      title = title
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1)
    )

  if (!is.null(palette)) {
    p <- p + ggplot2::scale_fill_manual(values = palette)
  }

  return(p)
}
