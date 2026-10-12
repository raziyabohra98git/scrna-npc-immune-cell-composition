# Submitted work and post-submission validation status

## Purpose

This note preserves the scope and results reported in the submitted M.Sc. Biotechnology dissertation presentation while separating them from subsequent source-metadata checks. It does not revise the submitted thesis or claim that missing analytical steps have been reproduced.

## Summary reported in the submitted presentation

The submitted presentation reports:

- Dataset: GEO GSE150430, a nasopharyngeal carcinoma single-cell RNA-seq dataset.
- Saved cells: 47,086.
- Normal cells: 2,518.
- Tumor cells: 44,568.
- Features in the saved object: 199.
- Clusters: 17.
- Annotated cell populations: 26.
- Methods described: data import, quality assessment, normalization, dimensionality reduction, clustering, UMAP, cell-type annotation, marker-gene identification, and Normal/Tumor comparison.

These are values and methods reported by the submitted work. The exact original commands for creating the Seurat object and obtaining its PCA, clustering, UMAP, and marker results have not been recovered.

## Separate source-metadata audit

The source metadata audit was run after submission using the barcode and cell-type metadata at the top of the original compressed file. It reported:

| Check | Result |
|---|---:|
| Barcodes in source metadata | 48,584 |
| Normal barcodes in source metadata | 2,583 |
| Tumor barcodes in source metadata | 46,001 |
| Cell types represented in source metadata | 26 |
| Barcodes in saved Seurat object | 47,086 |
| Source barcodes not present in saved object | 1,498 |

All saved-object barcodes were found in the source file. The 1,498 additional source barcodes were not found in the saved object. The reason for the difference has not been established.

The source metadata count table and percentage table were checked: group totals matched the source metadata and percentages summed to 100% within each group.

## Feature-row check

The original source file contains approximately 24,720 gene rows. The saved Seurat object has only 199 features. A separate audit matched each of the 199 saved feature rows to a source row and found that all 199 expression-value rows matched the source values across the 47,086 saved cells within tolerance `1e-8`.

The retained feature set corresponds to source rows 1–200 with source row 132 (`RPL22`) absent and `CASZ1` from source row 200 present. This establishes alignment for the retained rows, not that the analysis used the full transcriptome. The reason for the 199-feature selection has not been recovered.

Additional observations:

- The saved layer named `counts` contains fractional values, so it should not be described as raw integer counts without further evidence.
- The saved object's `percent.mt` values were all zero. With only 199 retained features, that metric should not be interpreted as a validated mitochondrial QC assessment.
- The original `analysis.R` was zero bytes in the local project listing.
- The recovered `analysis_reconstructed.R` begins by loading `seurat_thesis_final.rds`; it does not reconstruct the original object from the full source matrix and has syntax/undefined-object issues. It is not a verified end-to-end pipeline and should not be run as-is.

## Interpretation for this repository

The submitted figures and counts are retained as a historical record of what was presented. They are not claimed to have been independently reproduced from the full gene matrix. The existing `R/audit_source_metadata.R` script documents a separate, reproducible metadata-level check and writes aggregate source counts and percentages.

The source-to-object discrepancy remains unresolved. The repo should not claim a validated genome-wide single-cell analysis until the input representation, feature selection, omitted cells, and original Seurat workflow are resolved.
