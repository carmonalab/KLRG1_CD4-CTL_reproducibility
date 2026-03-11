# scRNA-seq processing

Reproducible preprocessing pipeline for Zenodo data (DOI: 10.5281/zenodo.18683366).

This pipeline includes:
- demultiplexing (MULTIseqDemux)
- quality control
- ProjecTILs annotation
- STACAS integration

This pipeline does not include TCR processing.

## Input files

Place Zenodo files in a data folder (default: `scRNA-seq/data/zenodo_18683366`):
- C0_sample_feature_bc_matrix.zip
- C1_sample_feature_bc_matrix.zip
- ...
- C8_sample_feature_bc_matrix.zip
- metadata.xlsx

## Reference maps

Place reference RDS files in `scRNA-seq/cache` (or pass env vars):
- CD4T_human_ref_v2.rds
- CD8T_human_ref_v1.rds

References:
- CD4: https://doi.org/10.6084/m9.figshare.24886611.v1
- CD8: https://doi.org/10.6084/m9.figshare.23608308.v1

## Run

From `scRNA-seq`:

```bash
Rscript data_processing.R
```

Optional environment variables:

```bash
SCRNA_DATA_DIR="/absolute/path/to/zenodo_data" \
SCRNA_METADATA_XLSX="/absolute/path/to/metadata.xlsx" \
CD4_REF_RDS="/absolute/path/to/CD4T_human_ref_v2.rds" \
CD8_REF_RDS="/absolute/path/to/CD8T_human_ref_v1.rds" \
NCORES=8 \
Rscript data_processing.R
```

## Outputs

Saved to `scRNA-seq/cache`:
- integrated_annotated_CD4.rds
- integrated_annotated_CD8.rds

## Zheng pan-cancer meta-analysis

Script:
- `Zheng_pan-pancer_analysis.R`

This script summarizes processing from Zheng_correlation_BestHit and includes:
- processing and ProjecTILs annotation of Zheng CD4/CD8 T cells
- proportions of CD4 cell types across tissues with Wilcoxon tests
- KLRG1+ vs KLRG1- classification in CD4 using scGate
- log2 ratio analysis of KLRG1+ vs KLRG1- across tissues with Wilcoxon tests

It does not include TCR processing and does not apply Th17/CD4.CTL relabeling filters.

### Inputs

- Zheng object RDS (default expected path: `scRNA-seq/cache/Zheng_2021_ 34914499_whole.rds`)
- Download source for Zheng object: http://cancer-pku.cn:3838/PanC_T/
- CD4 reference: `scRNA-seq/cache/CD4T_human_ref_v2.rds`
- CD8 reference: `scRNA-seq/cache/CD8T_human_ref_v1.rds`

### Run

From `scRNA-seq`:

```bash
Rscript Zheng_pan-pancer_analysis.R
```

Optional environment variables:

```bash
ZHENG_RDS="/absolute/path/to/Zheng_object.rds" \
CD4_REF_RDS="/absolute/path/to/CD4T_human_ref_v2.rds" \
CD8_REF_RDS="/absolute/path/to/CD8T_human_ref_v1.rds" \
Rscript Zheng_pan-pancer_analysis.R
```

### Outputs

Saved to `scRNA-seq/cache`:
- zheng_tcell_annotation_metadata.rds
- zheng_cd4_props_functional_cluster.csv
- zheng_cd4_props_meta_cluster.csv
- zheng_cd4_klrg1_metadata.rds
- zheng_cd4_klrg1_ratio.csv

Saved to `scRNA-seq/plots`:
- Zheng_CD4_Proportions_ProjecTILs.pdf
- Zheng_CD4_Proportions_Authors.pdf
- Zheng_CD4_KLRG1_ratio.pdf
