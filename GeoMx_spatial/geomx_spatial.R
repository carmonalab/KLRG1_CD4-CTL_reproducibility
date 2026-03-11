suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readxl)
  library(ggplot2)
  library(ggrepel)
  library(DESeq2)
  library(apeglm)
  library(UCell)
})

set.seed(22)

# -------------------------------
# Paths and parameters
# -------------------------------
data_dir <- Sys.getenv("GEOMX_DATA_DIR", unset = "data/data_CD4CD20 segments")
signatures_rds <- Sys.getenv(
  "UCELL_SIGNATURES_RDS",
  unset = "../scRNA-seq/cache/ucell_signatures_geomx.rds"
)

results_dir <- "results"
plot_dir <- file.path(results_dir, "plots")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

fdr_threshold <- 0.05
fc_threshold <- 2


# -------------------------------
# Load data
# -------------------------------
required_files <- c(
  "Annotation_B_distance.xlsx",
  "Raw data.xlsx",
  "Annotation_Probe.xlsx"
)

missing_files <- required_files[!file.exists(file.path(data_dir, required_files))]
if (length(missing_files) > 0) {
  stop(
    paste0(
      "Missing input files in ", data_dir, ": ", paste(missing_files, collapse = ", "),
      "\nDownload Zenodo DOI 10.5281/zenodo.18761515 and extract data_CD4CD20_segments.zip"
    )
  )
}

annot <- read_xlsx(file.path(data_dir, "Annotation_B_distance.xlsx"))
probe <- read_xlsx(file.path(data_dir, "Annotation_Probe.xlsx")) %>%
  filter(CodeClass == "Endogenous") %>%
  distinct(RTS_ID, TargetName)
raw <- read_xlsx(file.path(data_dir, "Raw data.xlsx"))


# -------------------------------
# Metadata processing
# -------------------------------
annot <- annot %>%
  mutate(
    Patient_letter = Patient,
    Patient = factor(
      Patient,
      levels = c("A", "B", "C", "D", "E"),
      labels = c("LAU706", "LAU380", "LAU1417", "LAU372", "LAU1397")
    )
  ) %>%
  rename(Spot = dcc)

sample_info <- annot %>%
  column_to_rownames("Spot") %>%
  transmute(
    Type = factor(Type),
    Localisation = factor(Localisation),
    Patient = factor(Patient),
    B_dist = factor(B_dist)
  )

raw_counts <- raw %>%
  rename(RTS_ID = gene) %>%
  left_join(probe, by = "RTS_ID") %>%
  mutate(gene = if_else(is.na(TargetName), RTS_ID, TargetName)) %>%
  select(-RTS_ID, -TargetName)

spot_cols <- intersect(colnames(raw_counts), rownames(sample_info))
if (length(spot_cols) == 0) {
  stop("No overlapping spot IDs between count matrix and metadata")
}

sample_info <- sample_info[spot_cols, , drop = FALSE]

raw_counts <- raw_counts %>%
  select(gene, all_of(spot_cols)) %>%
  mutate(across(-gene, as.numeric)) %>%
  group_by(gene) %>%
  summarise(across(everything(), ~ sum(.x, na.rm = TRUE)), .groups = "drop")

count_mat <- raw_counts %>%
  column_to_rownames("gene") %>%
  as.matrix()

count_mat <- count_mat[, rownames(sample_info), drop = FALSE]

write.csv(sample_info %>% rownames_to_column("Spot"),
          file.path(results_dir, "metadata_processed.csv"),
          row.names = FALSE)
write.csv(count_mat %>% as.data.frame() %>% rownames_to_column("gene"),
          file.path(results_dir, "raw_counts_gene_symbols.csv"),
          row.names = FALSE)


# -------------------------------
# Differential expression (DESeq2) and volcano plots
# -------------------------------
for (cell_type in intersect(c("CD4", "CD20"), unique(as.character(sample_info$Type)))) {
  meta_sub <- sample_info[sample_info$Type == cell_type, , drop = FALSE]

  if (nrow(meta_sub) < 4 || length(unique(meta_sub$Localisation)) < 2) {
    message("Skipping localisation DE for ", cell_type, " (insufficient samples)")
    next
  }

  counts_sub <- count_mat[, rownames(meta_sub), drop = FALSE]
  design_formula <- if (length(unique(meta_sub$Patient)) >= 2) {
    ~ Localisation + Patient
  } else {
    ~ Localisation
  }

  dds <- DESeqDataSetFromMatrix(
    countData = round(counts_sub),
    colData = droplevels(meta_sub),
    design = design_formula
  )
  dds <- dds[rowSums(counts(dds)) > 10, ]

  if (nrow(dds) < 20) {
    message("Skipping localisation DE for ", cell_type, " (too few genes after filtering)")
    next
  }

  dds <- tryCatch(DESeq(dds), error = function(e) NULL)
  if (is.null(dds)) {
    message("Skipping localisation DE for ", cell_type, " (DESeq2 fit failed)")
    next
  }

  coef_names <- resultsNames(dds)
  loc_coef <- coef_names[grepl("^Localisation_", coef_names)]

  for (cf in loc_coef) {
    res <- lfcShrink(dds, coef = cf, type = "apeglm") %>%
      as.data.frame() %>%
      rownames_to_column("gene") %>%
      mutate(
        padj = if_else(is.na(padj), 1, padj),
        neglog10_padj = -log10(padj),
        status = case_when(
          padj < fdr_threshold & log2FoldChange >= log2(fc_threshold) ~ "up",
          padj < fdr_threshold & log2FoldChange <= -log2(fc_threshold) ~ "down",
          TRUE ~ "ns"
        )
      )

    write.csv(
      res,
      file.path(results_dir, paste0("DESeq2_Localisation_", cell_type, "_", cf, ".csv")),
      row.names = FALSE
    )

    labels_df <- res %>%
      filter(status != "ns") %>%
      arrange(padj) %>%
      slice_head(n = 20)

    p <- ggplot(res, aes(x = log2FoldChange, y = neglog10_padj, color = status)) +
      geom_point(alpha = 0.8, size = 1.8) +
      geom_vline(xintercept = c(-log2(fc_threshold), log2(fc_threshold)), linetype = 2, color = "grey50") +
      geom_hline(yintercept = -log10(fdr_threshold), linetype = 2, color = "grey50") +
      ggrepel::geom_text_repel(
        data = labels_df,
        aes(label = gene),
        size = 3,
        max.overlaps = 20,
        box.padding = 0.3,
        segment.alpha = 0.5,
        show.legend = FALSE
      ) +
      scale_color_manual(values = c(down = "firebrick", ns = "grey70", up = "steelblue")) +
      labs(
        title = paste0("Localisation ", cell_type, " | ", gsub("_", " ", cf)),
        x = "log2 fold-change",
        y = "-log10 adjusted p-value",
        color = NULL
      ) +
      theme_bw(base_size = 11)

    ggsave(
      file.path(plot_dir, paste0("Volcano_Localisation_", cell_type, "_", cf, ".png")),
      p,
      width = 7,
      height = 5,
      dpi = 300
    )
  }
}

# B-cell distance DE in CD4 spots
meta_cd4_bdist <- sample_info %>%
  as.data.frame() %>%
  rownames_to_column("Spot") %>%
  filter(Type == "CD4", !is.na(B_dist)) %>%
  column_to_rownames("Spot")

if (nrow(meta_cd4_bdist) >= 4 && length(unique(meta_cd4_bdist$B_dist)) >= 2) {
  counts_cd4_bdist <- count_mat[, rownames(meta_cd4_bdist), drop = FALSE]
  design_formula <- if (length(unique(meta_cd4_bdist$Patient)) >= 2) {
    ~ B_dist + Patient
  } else {
    ~ B_dist
  }

  dds <- DESeqDataSetFromMatrix(
    countData = round(counts_cd4_bdist),
    colData = droplevels(meta_cd4_bdist),
    design = design_formula
  )
  dds <- dds[rowSums(counts(dds)) > 10, ]

  if (nrow(dds) >= 20) {
    dds <- tryCatch(DESeq(dds), error = function(e) NULL)

    if (!is.null(dds)) {
      coef_names <- resultsNames(dds)
      bdist_coef <- coef_names[grepl("^B_dist_", coef_names)]

      for (cf in bdist_coef) {
        res <- lfcShrink(dds, coef = cf, type = "apeglm") %>%
          as.data.frame() %>%
          rownames_to_column("gene") %>%
          mutate(
            padj = if_else(is.na(padj), 1, padj),
            neglog10_padj = -log10(padj),
            status = case_when(
              padj < fdr_threshold & log2FoldChange >= log2(fc_threshold) ~ "up",
              padj < fdr_threshold & log2FoldChange <= -log2(fc_threshold) ~ "down",
              TRUE ~ "ns"
            )
          )

        write.csv(
          res,
          file.path(results_dir, paste0("DESeq2_B_dist_CD4_", cf, ".csv")),
          row.names = FALSE
        )

        labels_df <- res %>%
          filter(status != "ns") %>%
          arrange(padj) %>%
          slice_head(n = 20)

        p <- ggplot(res, aes(x = log2FoldChange, y = neglog10_padj, color = status)) +
          geom_point(alpha = 0.8, size = 1.8) +
          geom_vline(xintercept = c(-log2(fc_threshold), log2(fc_threshold)), linetype = 2, color = "grey50") +
          geom_hline(yintercept = -log10(fdr_threshold), linetype = 2, color = "grey50") +
          ggrepel::geom_text_repel(
            data = labels_df,
            aes(label = gene),
            size = 3,
            max.overlaps = 20,
            box.padding = 0.3,
            segment.alpha = 0.5,
            show.legend = FALSE
          ) +
          scale_color_manual(values = c(down = "firebrick", ns = "grey70", up = "steelblue")) +
          labs(
            title = paste0("B_dist CD4 | ", gsub("_", " ", cf)),
            x = "log2 fold-change",
            y = "-log10 adjusted p-value",
            color = NULL
          ) +
          theme_bw(base_size = 11)

        ggsave(
          file.path(plot_dir, paste0("Volcano_B_dist_CD4_", cf, ".png")),
          p,
          width = 7,
          height = 5,
          dpi = 300
        )
      }
    }
  }
}


# -------------------------------
# UCell deconvolution
# -------------------------------
if (!file.exists(signatures_rds)) {
  stop(
    paste0(
      "UCell signatures file not found: ", signatures_rds,
      "\nPlease place it in scRNA-seq/cache (or set UCELL_SIGNATURES_RDS)."
    )
  )
}

sig_obj <- readRDS(signatures_rds)

# Minimal supported signature formats from RDS:
# 1) named list of character vectors
# 2) data.frame with signature/celltype/cluster + gene/marker/feature columns
# 3) list with $signatures following 1 or 2
if (is.list(sig_obj) && !is.null(names(sig_obj)) && all(vapply(sig_obj, is.character, logical(1)))) {
  signatures <- sig_obj
} else if (is.list(sig_obj) && "signatures" %in% names(sig_obj)) {
  tmp <- sig_obj$signatures
  if (is.list(tmp) && !is.null(names(tmp)) && all(vapply(tmp, is.character, logical(1)))) {
    signatures <- tmp
  } else if (inherits(tmp, "data.frame")) {
    nms <- tolower(colnames(tmp))
    sig_col <- colnames(tmp)[match(c("signature", "celltype", "cluster", "name"), nms, nomatch = 0)]
    gene_col <- colnames(tmp)[match(c("gene", "marker", "feature"), nms, nomatch = 0)]
    sig_col <- sig_col[sig_col != ""]
    gene_col <- gene_col[gene_col != ""]
    if (length(sig_col) == 0 || length(gene_col) == 0) {
      stop("Could not identify signature/gene columns in signatures data.frame")
    }
    tmp2 <- tmp %>%
      transmute(signature = as.character(.data[[sig_col[1]]]), gene = as.character(.data[[gene_col[1]]])) %>%
      filter(!is.na(signature), !is.na(gene), signature != "", gene != "")
    signatures <- split(tmp2$gene, tmp2$signature)
  } else {
    stop("Unsupported signatures object in $signatures")
  }
} else if (inherits(sig_obj, "data.frame")) {
  nms <- tolower(colnames(sig_obj))
  sig_col <- colnames(sig_obj)[match(c("signature", "celltype", "cluster", "name"), nms, nomatch = 0)]
  gene_col <- colnames(sig_obj)[match(c("gene", "marker", "feature"), nms, nomatch = 0)]
  sig_col <- sig_col[sig_col != ""]
  gene_col <- gene_col[gene_col != ""]
  if (length(sig_col) == 0 || length(gene_col) == 0) {
    stop("Could not identify signature/gene columns in signatures data.frame")
  }
  tmp2 <- sig_obj %>%
    transmute(signature = as.character(.data[[sig_col[1]]]), gene = as.character(.data[[gene_col[1]]])) %>%
    filter(!is.na(signature), !is.na(gene), signature != "", gene != "")
  signatures <- split(tmp2$gene, tmp2$signature)
} else {
  stop("Unsupported signatures RDS format")
}

signatures <- lapply(signatures, unique)
signatures <- lapply(signatures, function(g) intersect(g, rownames(count_mat)))
signatures <- signatures[vapply(signatures, length, integer(1)) >= 3]

if (length(signatures) == 0) {
  stop("No signature has at least 3 overlapping genes in GeoMx count matrix")
}

u_scores <- UCell::ScoreSignatures_UCell(
  matrix = count_mat,
  features = signatures,
  ncores = max(1, parallel::detectCores() - 1)
)

u_long <- as.data.frame(u_scores) %>%
  rownames_to_column("Spot") %>%
  pivot_longer(-Spot, names_to = "Celltype", values_to = "Score")

u_long <- sample_info %>%
  rownames_to_column("Spot") %>%
  right_join(u_long, by = "Spot")

write.csv(u_long, file.path(results_dir, "UCell_scores_long.csv"), row.names = FALSE)

p_type <- ggplot(u_long, aes(x = Type, y = Score, fill = Type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.35, size = 0.9) +
  facet_wrap(~ Celltype, scales = "free_y") +
  theme_bw(base_size = 10) +
  labs(title = "UCell scores by Type")
ggsave(file.path(plot_dir, "UCell_by_Type.png"), p_type, width = 12, height = 8, dpi = 300)

p_loc <- ggplot(u_long, aes(x = Localisation, y = Score, fill = Localisation)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.35, size = 0.9) +
  facet_wrap(~ Celltype, scales = "free_y") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "UCell scores by Localisation")
ggsave(file.path(plot_dir, "UCell_by_Localisation.png"), p_loc, width = 12, height = 8, dpi = 300)

u_bdist <- u_long %>% filter(!is.na(B_dist))
if (nrow(u_bdist) > 0) {
  p_bdist <- ggplot(u_bdist, aes(x = B_dist, y = Score, fill = B_dist)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.15, alpha = 0.35, size = 0.9) +
    facet_wrap(~ Celltype, scales = "free_y") +
    theme_bw(base_size = 10) +
    labs(title = "UCell scores by B-cell distance")
  ggsave(file.path(plot_dir, "UCell_by_B_dist.png"), p_bdist, width = 12, height = 8, dpi = 300)
}

message("Done. Outputs written to: ", normalizePath(results_dir, winslash = "/", mustWork = FALSE))