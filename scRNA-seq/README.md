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
