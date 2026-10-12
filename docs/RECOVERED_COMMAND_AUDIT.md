# Recovered Dissertation Command Audit

**Project:** Comparative Analysis of Immune Cell Composition Between Normal and Tumor Tissue Using scRNA-seq
**Dataset:** GSE150430 (NPC single-cell RNA-seq)
**Purpose:** Separate commands recovered from the available saved history/logs from steps that are only described in the dissertation or inferred as a possible workflow.

## Scope and limitation

I searched the ChatGPT Library materials available from this conversation. I can search saved files and pasted terminal outputs, but I do not have a tool that searches every conversation in the user's ChatGPT sidebar or opens arbitrary past chat threads. The recovered materials include pasted outputs from the old `.Rhistory`, a file listing, and the reconstructed analysis script. These are useful evidence, but they are not the complete conversations.

**Important distinction:** `.Rhistory` records submitted commands, not whether every command finished successfully. A command's presence alone is not proof that it ran without error.

## 1. Commands present in the recovered R history

### Environment and file setup

The history includes commands such as:

```r
install.packages("Seurat")
library(Seurat)
sessionInfo()
install.packages("tidyverse")
install.packages("patchwork")
library(patchwork)

dir.create("data")
dir.create("scripts")
dir.create("results")
dir.create("figures")
dir.create("docs")
```

### Source-file inspection and import attempts

The history records attempts to read the compressed source file and inspect its structure:

```r
data_file <- "data/GSE150430_npc_scRNA_hg19_processed_data.txt.gz"

raw_data <- read.table(
  gzfile(data_file),
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

test_data <- read.table(
  gzfile(data_file),
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE,
  nrows = 5
)

con <- gzfile(data_file)
lines <- readLines(con, n = 5)
close(con)

count.fields(gzfile(data_file), sep = "\t", n = 2)

data_preview <- read.table(
  gzfile(data_file),
  header = FALSE,
  sep = "\t",
  nrows = 20,
  fill = TRUE
)

con <- gzfile(data_file)
gene_lines <- readLines(con, n = 100)
close(con)
```

Later history also records another full import attempt and a `data.table::fread()` attempt:

```r
file_path <- "data/GSE150430_npc_scRNA_hg19_processed_data.txt.gz"
expr_data <- read.table(
  gzfile(file_path),
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

library(data.table)
expr_data <- fread(cmd = paste("gunzip -c", data_file))
```

The recovered history does **not** establish that either full import completed successfully or produced a correctly structured full expression matrix.

### Metadata and normal/tumor composition commands

The history shows that cell names, cell-type labels, and group labels were derived from the first two lines of the processed file:

```r
cell_names <- strsplit(lines[1], "\t")[[1]]
cell_types <- strsplit(lines[2], "\t")[[1]]

sample_group <- ifelse(
  substr(cell_names, 1, 1) == "N",
  "Normal",
  "Tumor"
)
```

The history then records several attempts to remove leading elements and resize vectors to equal lengths. The exact final sequence is potentially error-prone and should not be copied into a new pipeline without validation. A metadata table and grouped percentages were later created:

```r
metadata <- data.frame(
  Cell = cell_names,
  Group = sample_group,
  CellType = cell_types
)

metadata$CellType <- gsub('"', "", metadata$CellType)

library(dplyr)
prop_table <- metadata %>%
  group_by(Group, CellType) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Group) %>%
  mutate(percent = n / sum(n) * 100)

heatmap_data <- prop_table %>%
  select(Group, CellType, percent) %>%
  tidyr::pivot_wider(
    names_from = Group,
    values_from = percent,
    values_fill = 0
  )
```

A cell-class recoding step was also recorded for broad categories:

```r
metadata$CellClass <- metadata$CellType
metadata$CellClass <- gsub(".*B cell.*", "B cells", metadata$CellClass)
metadata$CellClass <- gsub(".*CD8.*", "T cells (CD8+)", metadata$CellClass)
metadata$CellClass <- gsub(".*CD4.*", "T cells (CD4+)", metadata$CellClass)
metadata$CellClass <- gsub(".*Tregs.*", "T cells (Tregs)", metadata$CellClass)
metadata$CellClass <- gsub(".*Monocyte.*|.*Macrophage.*", "Myeloid cells", metadata$CellClass)
metadata$CellClass <- gsub("CAFs", "Stromal (CAFs)", metadata$CellClass)
metadata$CellClass <- gsub("Malignant", "Tumor cells", metadata$CellClass)
metadata$CellClass <- gsub("Non-malignant epithelial", "Epithelial", metadata$CellClass)

prop_class <- metadata %>%
  group_by(Group, CellClass) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Group) %>%
  mutate(percent = n / sum(n))

table_counts <- table(prop_class$Group, prop_class$CellClass)
chisq_res <- chisq.test(table_counts)
chisq_res$stdres
```

### Figures and saved outputs present in the history

The recovered history contains commands to make the two early composition figures and save summary tables/workspace:

```r
png(
  "Figure1_Tumor_Cell_Composition.png",
  width = 1800,
  height = 1200,
  res = 300
)
barplot(
  percent_table["Tumor", ],
  las = 2,
  cex.names = 0.7,
  main = "Cell Type Composition in Tumor Samples (%)"
)
dev.off()

png(
  "Figure2_Normal_vs_Tumor_Comparison.png",
  width = 2200,
  height = 1400,
  res = 300
)
barplot(
  t(percent_table),
  beside = TRUE,
  las = 2,
  cex.names = 0.7,
  legend = TRUE,
  main = "Comparison of Cell Types in Normal and Tumor Samples (%)"
)
dev.off()

write.csv(percent_table, "Cell_Type_Percentages.csv")
save.image("Thesis_Workspace.RData")
```

The history also contains `ggplot2` composition plots (including `geom_tile()` and `geom_bar()`), and a normalized stacked-bar plot using `position = "fill"` and `scales::percent_format()`.

## 2. What was not recovered from these snippets

I did not find the original commands that built the Seurat object from the complete expression matrix. In the available history/logs, I could not verify commands/results for:

- `CreateSeuratObject()`
- `NormalizeData()`
- `FindVariableFeatures()`
- `ScaleData()`
- `RunPCA()`
- `FindNeighbors()`
- `FindClusters()`
- `RunUMAP()`

The retrieved `analysis.R` was reported as a zero-byte file. The later `analysis_reconstructed.R` starts by loading the already-created object with `readRDS("seurat_thesis_final.rds")`; it is not the original processing pipeline. In that reconstructed script, `cell_compare` is referenced without being created in the displayed code, and one line joins `save.image(...)` and `markers <- FindAllMarkers(...)` without a separator. Therefore, it should not be run as-is or treated as verified original code.

The dissertation document describes a workflow involving normalization, variable feature selection, scaling, PCA, clustering, UMAP, marker identification, and cell-type annotation. That is a **methods description**, not recovered proof of the exact commands, parameters, or successful execution.

## 3. Current data-validation evidence from this audit

The direct source audit performed separately in RStudio reported:

| Item | Source file | Saved Seurat object |
|---|---:|---:|
| Total cells | 48,584 | 47,086 |
| Normal cells | 2,583 | 2,518 |
| Tumor cells | 46,001 | 44,568 |
| Features | approximately 24,720 gene rows in source | 199 |

All 47,086 saved cell barcodes were found in the source. A comparison of every one of the 199 features in the saved object against the corresponding source row reported 199/199 value matches within tolerance 1e-8. The saved features correspond to source rows 1–200 except for source row 132 (`RPL22`); `CASZ1` from source row 200 is present. This supports alignment of the retained rows and cells, but it does not validate an analysis as representing the complete transcriptome.

The source metadata audit produced 2,583 Normal and 46,001 Tumor cells, with cell-type percentages summing to 100% within each group. The smaller `Thesis_Workspace.RData` had 48,582 cells (2,582 Normal and 46,000 Tumor), so it is not the authoritative source for final counts.

## 4. Recommended next action

Do not spend more time loading 820 MB backups blindly. The fastest defensible route is to write a new, clearly labelled **reconstructed pipeline** that imports and validates the original file format, checks all cell and gene identifiers, and only then performs the Seurat analysis using explicitly documented parameters. Preserve old plots and results as historical outputs, not validated final figures, until that pipeline reproduces them from the full source data.
