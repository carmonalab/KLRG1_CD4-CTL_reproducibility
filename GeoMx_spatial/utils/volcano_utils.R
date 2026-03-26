volcano_plot <- function(res_df,
                         fc = 2,
                         fdr = 0.05,
                         title = NULL,
                         max_labels = 13,
                         point_size = 2,
                         alpha = 0.95) {
  stopifnot(is.data.frame(res_df))
  required <- c("gene", "log2FoldChange", "padj")
  missing <- setdiff(required, colnames(res_df))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "))
  }

  col_vol <- c(
    downregulated = "firebrick2",
    upregulated = "darkseagreen2",
    NS = "grey"
  )

  vol <- res_df |>
    dplyr::transmute(
      gene = .data$gene,
      log2FC = .data$log2FoldChange,
      padj = .data$padj,
      sig = -log10(.data$padj)
    ) |>
    dplyr::filter(is.finite(.data$log2FC), is.finite(.data$sig)) |>
    dplyr::mutate(
      sig = pmin(.data$sig, 320),
      S = dplyr::case_when(
        .data$padj < fdr & .data$log2FC >= log2(fc) ~ "upregulated",
        .data$padj < fdr & .data$log2FC <= -log2(fc) ~ "downregulated",
        TRUE ~ "NS"
      ),
      S = factor(.data$S, levels = c("downregulated", "upregulated", "NS"))
    )

  if (is.null(title)) {
    title <- "Volcano plot"
  }

  label_df <- vol |>
    dplyr::filter(.data$padj < fdr, abs(.data$log2FC) >= log2(fc)) |>
    dplyr::arrange(dplyr::desc(.data$sig)) |>
    dplyr::slice_head(n = max_labels)

  ggplot2::ggplot(vol, ggplot2::aes(x = .data$log2FC, y = .data$sig, color = .data$S)) +
    ggplot2::geom_point(size = point_size, alpha = alpha, show.legend = TRUE) +
    ggplot2::geom_vline(
      xintercept = c(-log2(fc), log2(fc)),
      linetype = 2,
      linewidth = 0.3,
      colour = "grey50"
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(fdr),
      linetype = 2,
      linewidth = 0.3,
      colour = "grey50"
    ) +
    ggplot2::scale_color_manual(values = col_vol, drop = FALSE) +
    ggrepel::geom_text_repel(
      data = label_df,
      ggplot2::aes(label = .data$gene),
      color = "black",
      size = 5,
      fontface = "italic",
      max.overlaps = max_labels
    ) +
    ggplot2::xlab("log2 Fold Change") +
    ggplot2::ylab("-log10 adj. p-value") +
    ggplot2::scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
    ggplot2::scale_x_continuous(trans = scales::pseudo_log_trans(base = 10)) +
    ggplot2::ggtitle(title) +
    ggpubr::theme_classic2() +
    ggplot2::theme(
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 15, face = "bold"),
      axis.title.x = ggplot2::element_text(size = 15),
      axis.text.x = ggplot2::element_text(size = 12),
      axis.text.y = ggplot2::element_text(size = 12),
      axis.title.y = ggplot2::element_text(size = 15),
      legend.key.height = grid::unit(2, "line"),
      plot.title = ggplot2::element_text(size = 18, face = "bold")
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 4)))
}
