# GeoMx spatial reproducibility (CD4/CD20 segments)

Concise workflow to reproduce the GeoMx analysis components needed for the paper:
- Data loading
- Metadata processing
- Differential expression with DESeq2
- Volcano plot rendering
- Deconvolution using UCell only

## Data source

Raw files are in Zenodo DOI: 10.5281/zenodo.18761515

Use the archive:
- data_CD4CD20_segments.zip

Extract it so this folder exists:
- GeoMx_spatial/data/data_CD4CD20 segments

Expected files inside:
- Annotation_B_distance.xlsx
- Raw data.xlsx
- Annotation_Probe.xlsx

## UCell signatures input

The deconvolution step expects a signatures RDS file generated from scRNA-seq analysis.

Default path expected by the script:
- scRNA-seq/cache/ucell_signatures_geomx.rds

You can also override this path with an environment variable (see below).

Accepted RDS structures:
- Named list where each element is a character vector of genes
- Data frame with columns similar to signature/celltype/cluster and gene/marker/feature
- List containing a signatures field with one of the two formats above

## Run

From this directory:

```bash
Rscript geomx_spatial.R
```

Optional custom paths:

```bash
GEOMX_DATA_DIR="/absolute/path/to/data_CD4CD20 segments" \
UCELL_SIGNATURES_RDS="/absolute/path/to/ucell_signatures_geomx.rds" \
Rscript geomx_spatial.R
```

## Outputs

All outputs are written to:
- GeoMx_spatial/results

Main files:
- metadata_processed.csv
- raw_counts_gene_symbols.csv
- DESeq2_*.csv
- UCell_scores_long.csv

Plots are written to:
- GeoMx_spatial/results/plots

Main plots:
- Volcano_*.png
- UCell_by_Type.png
- UCell_by_Localisation.png
- UCell_by_B_dist.png (only when B_dist is available)
