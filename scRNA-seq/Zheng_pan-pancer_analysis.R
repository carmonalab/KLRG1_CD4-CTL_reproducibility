suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(Seurat)
  library(tibble)
  library(ProjecTILs)
  library(scGate)
  library(ggpubr)
  library(data.table)
})

set.seed(22)

# ============================================================
# Paths and config
# ============================================================
# Zheng source object can be downloaded from:
# http://cancer-pku.cn:3838/PanC_T/
cache_dir <- "cache"
plot_dir <- "plots"
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

zheng_rds <- Sys.getenv("ZHENG_RDS", unset = file.path(cache_dir, "Zheng_2021_ 34914499_whole.rds"))
cd4_ref <- Sys.getenv("CD4_REF_RDS", unset = file.path(cache_dir, "CD4T_human_ref_v2.rds"))
cd8_ref <- Sys.getenv("CD8_REF_RDS", unset = file.path(cache_dir, "CD8T_human_ref_v1.rds"))

min_cells_per_subset <- 10

if (!file.exists(zheng_rds)) {
  stop(
    "Zheng dataset not found: ", zheng_rds,
    "\nDownload object from http://cancer-pku.cn:3838/PanC_T/ and place it in cache/ or set ZHENG_RDS"
  )
}
if (!file.exists(cd4_ref)) {
  stop("CD4 reference not found: ", cd4_ref)
}
if (!file.exists(cd8_ref)) {
  stop("CD8 reference not found: ", cd8_ref)
}

# ============================================================
# Processing and ProjecTILs annotation
# ============================================================
zheng <- readRDS(zheng_rds)
zheng$authors_Tcell <- ifelse(grepl("^CD4", zheng$meta.cluster), "CD4T", "CD8T")

refcd4 <- load.reference.map(cd4_ref)
refcd8 <- load.reference.map(cd8_ref)

levs <- unique(zheng@meta.data[["Sample"]])
tcell <- list()

for (subset in levs) {
  s <- zheng[, zheng$Sample == subset]

  if ("CD4T" %in% s$authors_Tcell) {
    c4 <- s[, s$authors_Tcell == "CD4T"]
    if (ncol(c4) >= min_cells_per_subset) {
      tcell[[paste0(subset, "_CD4T")]] <- ProjecTILs.classifier(
        query = c4,
        ref = refcd4,
        filter.cells = FALSE,
        ncores = 1
      )
    }
  }

  if ("CD8T" %in% s$authors_Tcell) {
    c8 <- s[, s$authors_Tcell == "CD8T"]
    if (ncol(c8) >= min_cells_per_subset) {
      tcell[[paste0(subset, "_CD8T")]] <- ProjecTILs.classifier(
        query = c8,
        ref = refcd8,
        filter.cells = FALSE,
        ncores = 1
      )
    }
  }
}

if (length(tcell) == 0) {
  stop("No subsets passed the minimum cell threshold for annotation")
}

annotation_md <- lapply(tcell, function(x) {
  x@meta.data %>% rownames_to_column("cellid")
}) %>% data.table::rbindlist()

saveRDS(annotation_md, file.path(cache_dir, "zheng_tcell_annotation_metadata.rds"))

bas <- CreateSeuratObject(
  counts = zheng@assays$RNA$counts[, annotation_md$cellid],
  meta.data = annotation_md %>% column_to_rownames("cellid")
)
bas <- NormalizeData(bas)
rm(zheng)
gc()

# ============================================================
# Proportions by tissue + Wilcoxon tests
# Keep ProjecTILs labels as-is (no Th17/CTL regrouping)
# ============================================================
cd4_md <- bas@meta.data %>%
  as.data.frame() %>%
  filter(authors_Tcell == "CD4T")

s <- "patient"
t <- "Tissue"

md_join <- cd4_md %>%
  distinct(patient, Tissue, cancerType)

compute_props <- function(df, cluster_col) {
  out <- df %>%
    filter(.data[[t]] %in% c("P", "T", "N")) %>%
    group_by(.data[[s]], .data[[t]], .data[[cluster_col]]) %>%
    summarize(n = n(), .groups = "drop") %>%
    complete(.data[[s]], .data[[t]], .data[[cluster_col]], fill = list(n = 0)) %>%
    group_by(.data[[s]], .data[[t]]) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    left_join(md_join, by = c(s, t)) %>%
    mutate(
      Tissue = factor(
        .data[[t]],
        levels = c("T", "N", "P", "L"),
        labels = c("Tumor", "Normal", "Blood", "LN")
      )
    ) %>%
    filter(!is.na(.data$prop))
  out
}

props_functional <- compute_props(cd4_md, "functional.cluster")
props_authors <- compute_props(cd4_md, "meta.cluster")

write.csv(props_functional, file.path(cache_dir, "zheng_cd4_props_functional_cluster.csv"), row.names = FALSE)
write.csv(props_authors, file.path(cache_dir, "zheng_cd4_props_meta_cluster.csv"), row.names = FALSE)

comps <- list(c("Tumor", "Normal"), c("Tumor", "Blood"), c("Normal", "Blood"))

p_prop_functional <- props_functional %>%
  ggplot(aes(Tissue, prop)) +
  geom_boxplot(outlier.colour = NA) +
  geom_point(aes(color = cancerType), alpha = 0.5) +
  geom_line(aes(group = patient, color = cancerType), alpha = 0.6) +
  stat_compare_means(
    comparisons = comps,
    method = "wilcox.test",
    paired = FALSE,
    label = "p.format",
    p.adjust.method = "BH",
    size = 3.2,
    tip.length = 0,
    bracket.size = 0.3
  ) +
  facet_wrap(~functional.cluster, scales = "free_y", ncol = 3) +
  labs(y = "Proportion of cells", title = "CD4 ProjecTILs cell-type proportions", color = "Cancer type") +
  ggpubr::theme_classic2()

ggsave(file.path(plot_dir, "Zheng_CD4_Proportions_ProjecTILs.pdf"), p_prop_functional, width = 12, height = 9)

auth_col <- "meta.cluster"
p_prop_authors <- props_authors %>%
  ggplot(aes(Tissue, prop)) +
  geom_boxplot(outlier.colour = NA) +
  geom_point(aes(color = cancerType), alpha = 0.5) +
  geom_line(aes(group = patient, color = cancerType), alpha = 0.6) +
  stat_compare_means(
    comparisons = comps,
    method = "wilcox.test",
    paired = FALSE,
    label = "p.format",
    p.adjust.method = "BH",
    size = 3.2,
    tip.length = 0,
    bracket.size = 0.3
  ) +
  facet_wrap(as.formula(paste0("~", auth_col)), scales = "free_y", ncol = 3) +
  labs(y = "Proportion of cells", title = "CD4 authors cell-type proportions", color = "Cancer type") +
  ggpubr::theme_classic2()

ggsave(file.path(plot_dir, "Zheng_CD4_Proportions_Authors.pdf"), p_prop_authors, width = 12, height = 9)

# ============================================================
# KLRG1+ vs KLRG1- ratio in CD4 cells using scGate
# ============================================================
my_scGate_model <- gating_model(name = "KLRG1", signature = "KLRG1")
cd4 <- bas[, bas$authors_Tcell == "CD4T"]

levs <- unique(cd4@meta.data[["Sample"]])
cd4_klrg1 <- lapply(levs, function(subset) {
  s_obj <- cd4[, cd4$Sample == subset]
  if (ncol(s_obj) < min_cells_per_subset) {
    return(NULL)
  }

  d <- min(30, ncol(s_obj) - 1)
  if (d < 2) {
    return(NULL)
  }

  s_obj <- scGate(
    s_obj,
    model = my_scGate_model,
    ncores = 1,
    min.cells = d,
    pca.dim = d,
    k.param = d,
    verbose = FALSE
  )
  s_obj$KLRG1_expression <- ifelse(s_obj$is.pure == "Pure", "KLRG1+", "KLRG1-")
  s_obj
})

cd4_klrg1 <- Filter(Negate(is.null), cd4_klrg1)
if (length(cd4_klrg1) == 0) {
  stop("No CD4 subsets passed scGate KLRG1 step")
}

klrg1_md <- lapply(cd4_klrg1, function(x) {
  x@meta.data %>% rownames_to_column("cellid")
}) %>% data.table::rbindlist()

saveRDS(klrg1_md, file.path(cache_dir, "zheng_cd4_klrg1_metadata.rds"))

klrg1_counts <- klrg1_md %>%
  filter(Tissue %in% c("P", "T", "N")) %>%
  group_by(patient, Tissue, cancerType, KLRG1_expression) %>%
  summarize(n = n(), .groups = "drop") %>%
  complete(patient, Tissue, cancerType, KLRG1_expression = c("KLRG1+", "KLRG1-"), fill = list(n = 0))

klrg1_ratio <- klrg1_counts %>%
  group_by(patient, Tissue, cancerType) %>%
  summarize(
    ratio = (n[KLRG1_expression == "KLRG1+"] + 1) / (n[KLRG1_expression == "KLRG1-"] + 1),
    log2_ratio = log2(ratio),
    .groups = "drop"
  ) %>%
  mutate(
    Tissue = factor(
      Tissue,
      levels = c("T", "N", "P", "L"),
      labels = c("Tumor", "Normal", "Blood", "LN")
    )
  )

write.csv(klrg1_ratio, file.path(cache_dir, "zheng_cd4_klrg1_ratio.csv"), row.names = FALSE)

p_klrg1 <- klrg1_ratio %>%
  ggplot(aes(Tissue, log2_ratio)) +
  geom_boxplot(outlier.colour = NA) +
  geom_point(aes(color = cancerType), alpha = 0.5) +
  geom_line(aes(group = patient, color = cancerType), alpha = 0.6) +
  stat_compare_means(
    comparisons = comps,
    method = "wilcox.test",
    paired = FALSE,
    label = "p.format",
    p.adjust.method = "BH",
    size = 3.2,
    tip.length = 0,
    bracket.size = 0.3
  ) +
  labs(
    y = "Log2 ratio of KLRG1+ vs KLRG1- CD4 T cells",
    title = "KLRG1 composition in CD4 compartment"
  ) +
  ggpubr::theme_classic2() +
  theme(aspect.ratio = 1)

ggsave(file.path(plot_dir, "Zheng_CD4_KLRG1_ratio.pdf"), p_klrg1, width = 8, height = 6)

message("Done.")
message("Outputs written in cache and plots directories.")
