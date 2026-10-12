# Comparative Analysis of Immune Cell Composition in Nasopharyngeal Carcinoma Using scRNA-seq

## Project Overview

This MSc Biotechnology dissertation investigates differences in cellular composition between normal and tumour tissue using single-cell RNA sequencing (scRNA-seq) data from nasopharyngeal carcinoma (NPC).

The project focuses on exploring annotated cell populations, comparing their representation between normal and tumour groups, and documenting a reproducible computational analysis workflow.

## Dataset

- **Source:** NCBI Gene Expression Omnibus (GEO)
- **Accession:** GSE150430
- **Data:** Processed single-cell RNA-seq expression matrix with cell-type annotations
- **Primary tools:** R, Seurat, dplyr, ggplot2

## Dataset Summary

The saved Seurat object currently contains:

| Metric | Value |
|---|---:|
| Cells retained in the object | 47,086 |
| Normal cells | 2,518 |
| Tumour cells | 44,568 |
| Annotated cell populations | 26 |
| Clusters | 17 |
| Stored dimensionality reductions | PCA and UMAP |

These figures describe the saved object and should not be interpreted as independent validation of every downstream result.

## Analysis Scope

The intended analysis workflow covers:

1. Inspection and preparation of the processed expression matrix
2. Cell metadata and quality-control assessment
3. Dimensionality reduction and clustering
4. Visualisation of clusters and annotated cell populations
5. Comparison of cell-type composition between normal and tumour groups
6. Review of marker-gene results and summary tables

## Reproducibility and Current Validation

The source matrix contains substantially more gene rows than the saved Seurat object. The saved object contains 199 features across 47,086 cells, and its layer named `counts` contains fractional values.

These findings require further investigation before the saved PCA, UMAP, and marker-gene results can be treated as a verified genome-wide analysis. The original data and existing outputs are being preserved while the import and analysis workflow is reviewed.

## Repository Structure

- `R/` — R analysis scripts
- `scripts/` — supporting scripts
- `figures/` — selected figures, added after review
- `results/tables/` — reviewed aggregate summary tables
- `docs/` — methods and supplementary documentation
- `data/` — input data location; large source files are not committed

## Tools and Skills

- R programming
- Seurat
- Single-cell RNA-seq data analysis
- Cell-type composition analysis
- Data visualisation with ggplot2
- Reproducible research and Git version control

## Data Source

GEO accession: [GSE150430](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150430)

## Project Status

The repository is being organised for reproducibility and portfolio presentation. Figures, tables, and scripts will be added after their contents and analytical provenance have been reviewed.

