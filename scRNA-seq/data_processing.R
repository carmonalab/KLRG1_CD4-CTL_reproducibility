suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(tibble)
  library(ProjecTILs)
  library(STACAS)
})

set.seed(22)

# ============================================================
# Paths and settings
# ============================================================
# Zenodo DOI for this input data: 10.5281/zenodo.18683366
# Expected files in data_dir:
# - C*_sample_feature_bc_matrix.zip
# - metadata.xlsx

data_dir <- Sys.getenv("SCRNA_DATA_DIR", unset = "data/zenodo_18683366")
metadata_file <- Sys.getenv("SCRNA_METADATA_XLSX", unset = file.path(data_dir, "metadata.xlsx"))

cache_dir <- "cache"
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

ref_cd4_path <- Sys.getenv(
  "CD4_REF_RDS",
  unset = file.path(cache_dir, "CD4T_human_ref_v2.rds")
)
ref_cd8_path <- Sys.getenv(
  "CD8_REF_RDS",
  unset = file.path(cache_dir, "CD8T_human_ref_v1.rds")
)

ncores <- as.integer(Sys.getenv("NCORES", unset = "8"))
if (is.na(ncores) || ncores < 1) ncores <- 1

# ============================================================
# Input checks and unzip matrices
# ============================================================
if (!dir.exists(data_dir)) {
  stop("Data directory not found: ", data_dir)
}
if (!file.exists(metadata_file)) {
  stop("metadata.xlsx not found: ", metadata_file)
}

zip_files <- list.files(data_dir, pattern = "^C[0-9]+_sample_feature_bc_matrix\\.zip$", full.names = TRUE)
if (length(zip_files) == 0) {
  stop("No C*_sample_feature_bc_matrix.zip files found in: ", data_dir)
}

extract_root <- file.path(data_dir, "extracted")
dir.create(extract_root, recursive = TRUE, showWarnings = FALSE)

for (zf in zip_files) {
  cap <- sub("_sample_feature_bc_matrix\\.zip$", "", basename(zf))
  out_dir <- file.path(extract_root, cap)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  unzip(zf, exdir = out_dir)
}

captures <- sort(unique(sub("_sample_feature_bc_matrix\\.zip$", "", basename(zip_files))))

# ============================================================
# Read metadata (multiplex map)
# ============================================================
meta <- readxl::read_xlsx(metadata_file)
meta <- as.data.frame(meta, stringsAsFactors = FALSE)
colnames(meta) <- gsub("[[:space:]]+", "_", colnames(meta))

# Flexible column matching to support variants in metadata.xlsx
nms <- tolower(colnames(meta))

capture_col <- colnames(meta)[match(c("capture", "sample", "run", "pool"), nms, nomatch = 0)]
capture_col <- capture_col[capture_col != ""]

tag_col <- colnames(meta)[match(c("multi_id", "hashtag", "htag", "hto", "tag", "hash_id"), nms, nomatch = 0)]
tag_col <- tag_col[tag_col != ""]

patient_col <- colnames(meta)[match(c("patient", "patient_id"), nms, nomatch = 0)]
patient_col <- patient_col[patient_col != ""]

hla_col <- colnames(meta)[match(c("hla", "hla_type"), nms, nomatch = 0)]
hla_col <- hla_col[hla_col != ""]

antigen_col <- colnames(meta)[match(c("antigen", "target_antigen"), nms, nomatch = 0)]
antigen_col <- antigen_col[antigen_col != ""]

spec_col <- colnames(meta)[match(c("specificity", "tumor_specificity", "tumor_spe"), nms, nomatch = 0)]
spec_col <- spec_col[spec_col != ""]

subset_col <- colnames(meta)[match(c("t_subset", "tcell_subset", "subset", "cd4_cd8"), nms, nomatch = 0)]
subset_col <- subset_col[subset_col != ""]

tissue_col <- colnames(meta)[match(c("tissue", "source_tissue", "compartment"), nms, nomatch = 0)]
tissue_col <- tissue_col[tissue_col != ""]

if (length(tag_col) == 0) {
  stop("metadata.xlsx must include a hashtag/multi_id column")
}

tag_col <- tag_col[1]

if (length(capture_col) > 0) {
  capture_col <- capture_col[1]
}

if (length(capture_col) > 0) {
  meta <- meta %>%
    mutate(
      capture = as.character(.data[[capture_col]]),
      multi_id = as.character(.data[[tag_col]])
    )
} else {
  meta <- meta %>% mutate(multi_id = as.character(.data[[tag_col]]))
  # If capture is absent, derive it from multi_id formatted like C1-HT1-...
  meta$capture <- sub("^([Cc][0-9]+)-.*$", "\\1", meta$multi_id)
  meta$capture[!grepl("^[Cc][0-9]+-", meta$multi_id)] <- NA_character_
}

if (length(patient_col) > 0) meta$patient <- as.character(meta[[patient_col[1]]]) else meta$patient <- NA_character_
if (length(hla_col) > 0) meta$HLA <- as.character(meta[[hla_col[1]]]) else meta$HLA <- NA_character_
if (length(antigen_col) > 0) meta$antigen <- as.character(meta[[antigen_col[1]]]) else meta$antigen <- NA_character_
if (length(spec_col) > 0) meta$specificity <- as.character(meta[[spec_col[1]]]) else meta$specificity <- NA_character_
if (length(subset_col) > 0) meta$T_subset <- as.character(meta[[subset_col[1]]]) else meta$T_subset <- NA_character_
if (length(tissue_col) > 0) meta$tissue <- as.character(meta[[tissue_col[1]]]) else meta$tissue <- NA_character_

meta <- meta %>%
  dplyr::select(capture, multi_id, patient, HLA, antigen, specificity, T_subset, tissue) %>%
  mutate(
    capture = trimws(capture),
    multi_id = trimws(multi_id),
    multi_id_short = sub("^[Cc][0-9]+-", "", multi_id)
  ) %>%
  filter(!is.na(multi_id), multi_id != "")

# ============================================================
# Load each capture, QC, demultiplex
# ============================================================
samples <- list()
for (cap in captures) {
  message("Processing capture ", cap)

  cap_root <- file.path(extract_root, cap)

  # Find a folder containing matrix.mtx(.gz), features.tsv(.gz), barcodes.tsv(.gz)
  all_dirs <- unique(c(cap_root, list.dirs(cap_root, recursive = TRUE, full.names = TRUE)))
  matrix_dir <- NULL
  for (d in all_dirs) {
    has_m <- any(file.exists(file.path(d, c("matrix.mtx", "matrix.mtx.gz"))))
    has_f <- any(file.exists(file.path(d, c("features.tsv", "features.tsv.gz", "genes.tsv", "genes.tsv.gz"))))
    has_b <- any(file.exists(file.path(d, c("barcodes.tsv", "barcodes.tsv.gz"))))
    if (has_m && has_f && has_b) {
      matrix_dir <- d
      break
    }
  }

  if (is.null(matrix_dir)) {
    stop("Could not locate matrix directory for capture ", cap)
  }

  mat <- Seurat::Read10X(matrix_dir)

  if (!"Gene Expression" %in% names(mat)) {
    stop("Gene Expression assay missing in ", cap)
  }
  if (!"Antibody Capture" %in% names(mat)) {
    stop("Antibody Capture assay missing in ", cap, " (required for demultiplexing)")
  }

  seu <- CreateSeuratObject(counts = mat$`Gene Expression`, project = cap)

  # Prefix HTO names with capture to create globally unique MULTI_ID values.
  hto <- mat$`Antibody Capture`
  if (!all(grepl(paste0("^", cap, "-"), rownames(hto)))) {
    rownames(hto) <- paste0(cap, "-", rownames(hto))
  }
  seu[["HTG"]] <- CreateAssayObject(counts = hto)

  # QC
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  seu <- subset(seu, subset = nFeature_RNA > 300 & percent.mt < 8)

  # Demultiplex
  seu <- NormalizeData(seu, assay = "HTG", normalization.method = "CLR")
  seu <- MULTIseqDemux(seu, assay = "HTG", autoThresh = TRUE, maxiter = 20)
  seu$MULTI_ID <- factor(seu$MULTI_ID,
                         levels = c(rownames(seu[["HTG"]]), "Negative", "Doublet"))

  # Keep singlets only
  seu <- subset(seu, subset = MULTI_ID %in% c("Negative", "Doublet"), invert = TRUE)

  # Add capture metadata and merge multiplex metadata
  seu$capture <- cap
  m <- seu@meta.data %>%
    rownames_to_column("cell") %>%
    mutate(
      MULTI_ID = as.character(MULTI_ID),
      MULTI_ID_SHORT = sub("^[Cc][0-9]+-", "", MULTI_ID),
      capture = cap
    )

  cap_meta <- meta %>%
    filter(is.na(capture) | capture == cap) %>%
    mutate(capture = cap)

  # First try exact MULTI_ID match, then fallback to short hashtag id.
  m <- m %>%
    left_join(
      cap_meta %>% select(capture, multi_id, patient, HLA, antigen, specificity, T_subset, tissue),
      by = c("capture", "MULTI_ID" = "multi_id")
    )

  miss <- which(is.na(m$patient) & is.na(m$HLA) & is.na(m$antigen) & is.na(m$T_subset))
  if (length(miss) > 0) {
    m2 <- m[miss, ] %>%
      select(cell, capture, MULTI_ID, MULTI_ID_SHORT)
    m2 <- m2 %>%
      left_join(
        cap_meta %>% select(capture, multi_id_short, patient, HLA, antigen, specificity, T_subset, tissue),
        by = c("capture", "MULTI_ID_SHORT" = "multi_id_short")
      )

    m$patient[miss] <- m2$patient
    m$HLA[miss] <- m2$HLA
    m$antigen[miss] <- m2$antigen
    m$specificity[miss] <- m2$specificity
    m$T_subset[miss] <- m2$T_subset
    m$tissue[miss] <- m2$tissue
  }

  # Fallback derivations if metadata columns are missing
  if (all(is.na(m$T_subset))) {
    m$T_subset <- ifelse(grepl("CD8", m$MULTI_ID, ignore.case = TRUE), "CD8",
                         ifelse(grepl("CD4", m$MULTI_ID, ignore.case = TRUE), "CD4", NA))
  }
  if (all(is.na(m$tissue))) {
    m$tissue <- ifelse(grepl("PBMC", m$MULTI_ID, ignore.case = TRUE), "PBMC",
                       ifelse(grepl("TIL", m$MULTI_ID, ignore.case = TRUE), "TIL", NA))
  }
  if (all(is.na(m$specificity))) {
    m$specificity <- ifelse(is.na(m$antigen), NA,
                            ifelse(toupper(m$antigen) %in% c("TT", "HA"), "virus_specific", "tumor_specific"))
  }

  # Write back selected metadata columns
  seu$patient <- m$patient[match(colnames(seu), m$cell)]
  seu$HLA <- m$HLA[match(colnames(seu), m$cell)]
  seu$antigen <- m$antigen[match(colnames(seu), m$cell)]
  seu$specificity <- m$specificity[match(colnames(seu), m$cell)]
  seu$T_subset <- m$T_subset[match(colnames(seu), m$cell)]
  seu$tissue <- m$tissue[match(colnames(seu), m$cell)]

  samples[[cap]] <- seu
}

if (length(samples) == 0) {
  stop("No captures were processed")
}

# ============================================================
# Merge all captures
# ============================================================
merged <- samples[[1]]
if (length(samples) > 1) {
  merged <- merge(samples[[1]], y = samples[2:length(samples)])
}
merged <- NormalizeData(merged)

# Keep hashtags with at least 2 cells
tbl <- table(merged$MULTI_ID)
keep_multi <- names(tbl)[tbl > 1]
merged <- merged[, merged$MULTI_ID %in% keep_multi]

# ============================================================
# Load ProjecTILs references
# ============================================================
if (!file.exists(ref_cd4_path)) {
  stop("CD4 reference not found: ", ref_cd4_path,
       "\nDownload CD4T_human_ref_v2 from https://doi.org/10.6084/m9.figshare.24886611.v1")
}
if (!file.exists(ref_cd8_path)) {
  stop("CD8 reference not found: ", ref_cd8_path,
       "\nDownload CD8T_human_ref_v1 from https://doi.org/10.6084/m9.figshare.23608308.v1")
}

ref_maps <- list(
  CD4 = ProjecTILs::load.reference.map(ref_cd4_path),
  CD8 = ProjecTILs::load.reference.map(ref_cd8_path)
)

# ============================================================
# Annotation + integration for CD4 and CD8 separately
# ============================================================

for (ct in c("CD4", "CD8")) {
  message("Running annotation/integration for ", ct)

  sub <- merged[, merged$T_subset == ct]
  if (ncol(sub) < 50) {
    warning("Skipping ", ct, " (too few cells after filtering): ", ncol(sub))
    next
  }

  pred <- ProjecTILs::ProjecTILs.classifier(
    query = sub,
    ref = ref_maps[[ct]],
    split.by = "MULTI_ID",
    ncores = ncores,
    filter.cells = FALSE,
    skip.normalize = TRUE
  )

  obj_list <- SplitObject(pred, split.by = "orig.ident")

  # Remove very small batches that cannot be integrated reliably
  obj_list <- obj_list[vapply(obj_list, ncol, integer(1)) >= 20]
  if (length(obj_list) < 2) {
    warning("Skipping STACAS for ", ct, " (need at least 2 batches with >=20 cells)")
    saveRDS(pred, file.path(cache_dir, paste0("integrated_annotated_", ct, ".rds")))
    next
  }

  int <- Run.STACAS(
    object.list = obj_list,
    dims = 1:30,
    cell.labels = "functional.cluster"
  )
  int <- RunUMAP(int, dims = 1:30, verbose = FALSE)
  saveRDS(int, file.path(cache_dir, paste0("integrated_annotated_", ct, ".rds")))
}

message("Done.")
message("Outputs in ", normalizePath(cache_dir, winslash = "/", mustWork = FALSE))
