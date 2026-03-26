# KLRG1 identifies circulating cytotoxic CD4 T cells with selective anti-tumor function in human cancer

Reproducibility repository for all computational analyses supporting the manuscript.

---

## Overview

This repository contains the code required to reproduce the main bioinformatic results demonstrating that **KLRG1 marks a circulating cytotoxic CD4 T-cell program** in human cancer.

The analyses integrate:

* **Single-cell RNA-seq (scRNA-seq)**
  Identification of CD4 T-cell states (including cytotoxic CD4 and TFH-like populations), compositional analysis across compartments (blood vs tumor), and differential expression and pathway enrichment.

* **GeoMx spatial transcriptomics**
  Spatial characterization of CD4 T-cell programs in relation to tumor architecture (e.g., tumor proximity, B-cell niches), including marker visualization and deconvolution.

Each notebook explicitly documents:

* **Inputs**
* **Outputs**
* **Associated manuscript figure panels**

---

## Reproducibility

This project uses [`renv`](https://rstudio.github.io/renv/) to ensure a fully reproducible R environment.

To restore the environment:

```r
install.packages("renv")
renv::restore(prompt = FALSE)
```

---

## Repository structure and execution order

Analyses are organized into two main modules:

---

### 1. scRNA-seq analysis

#### Preprocessing and reference construction

1. **Data processing**
   `scRNA-seq/1.1.Data_processing.Rmd`

   * Download, quality control, and merging of datasets
   * Output: merged Seurat object used in downstream analyses

2. **Integration and annotation**
   `scRNA-seq/1.2.Integration_annotation.Rmd`

   * Dataset integration and CD4 T-cell annotation
   * Marker inspection to support manual state annotation
   * Compositional analyses across compartments and identification of cytotoxic CD4 markers (KLRG1-associated program)
   * Output: integrated CD4 reference object

#### External dataset analysis

3. **Zheng cohort processing**
   `scRNA-seq/2.1.Zheng_processing.Rmd`

   * Processing of an independent multi-cancer dataset

4. **Multi-cancer annotation**
   `scRNA-seq/2.2.Annotation_multicancer_dataset.Rmd`

   * Annotation using ProjecTILs in an independent multi-cancer cohort
   * Marker/state mapping and compositional readouts to validate the cytotoxic CD4 program in external data

#### Downstream analyses

5. **Differential expression and enrichment**
   `scRNA-seq/3.DEG_enrichment.Rmd`

   * Pseudobulk differential expression
   * Pathway enrichment analysis
   * STRING network visualization

6. **Custom gene expression analyses**
   `scRNA-seq/4.Custom_Genes_Expression.Rmd`

   * Targeted analysis of selected gene programs (e.g., IFN/MX1)
   * Pseudobulk-based quantification

---

### 2. GeoMx spatial transcriptomics

1. **Data processing**
   `GeoMx_spatial/1.Data_processing.Rmd`

   * Download from Zenodo
   * Construction of cached GeoMx object

2. **Differential expression analysis**
   `GeoMx_spatial/2.DEG_analysis.Rmd`

   * DESeq2-based comparisons
   * Volcano plots corresponding to manuscript figures

3. **Marker visualization**
   `GeoMx_spatial/3.Heatmap_markers.Rmd`

   * TFH marker heatmap reproduction

4. **CD4 spot deconvolution / signature scoring**
   `GeoMx_spatial/4.CD4_spots_Deconvolution.Rmd`

   * Projection of scRNA-seq-derived CD4 signatures into spatial data

---

## Data availability

All raw and processed data used in this study are publicly available:

* **scRNA-seq datasets**
  Zenodo: `10.5281/zenodo.18683366`
  (Downloaded within relevant notebooks)

* **GeoMx spatial transcriptomics data**
  Zenodo: `10.5281/zenodo.18761515`
  (Automatically retrieved in GeoMx workflows)

---

## Notes

* All intermediate objects are cached where possible to reduce recomputation time.
* Notebooks are designed to be run sequentially within each module.

---

## Contact

For questions regarding the code or reproducibility, please open an issue or contact the authors of the manuscript.
