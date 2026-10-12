# MSc Biotechnology Dissertation
# Source metadata audit for GSE150430
#
# This script reads the first two metadata rows of the processed GEO file.
# It does not load the full expression matrix or treat processed expression
# values as raw UMI counts.
#
# Before running:
# Place GSE150430_npc_scRNA_hg19_processed_data.txt.gz in the local data/
# directory. The dataset is excluded from Git by .gitignore.

input_file <- Sys.getenv(
  "GSE150430_FILE",
  unset = file.path("data", "GSE150430_npc_scRNA_hg19_processed_data.txt.gz")
)

if (!file.exists(input_file)) {
  stop(
    "Input file not found: ", input_file,
    "\nPlace the GSE150430 compressed file in the local data/ directory."
  )
}

con <- gzfile(input_file, open = "rt")

barcode_line <- readLines(con, n = 1, warn = FALSE)
celltype_line <- readLines(con, n = 1, warn = FALSE)

close(con)

if (length(barcode_line) != 1 || length(celltype_line) != 1) {
  stop("Could not read both expected metadata rows.")
}

barcodes <- strsplit(barcode_line, "\t", fixed = TRUE)[[1]]
celltype_fields <- strsplit(celltype_line, "\t", fixed = TRUE)[[1]]

if (celltype_fields[1] != "Cell_type") {
  stop("Unexpected second-row label: ", celltype_fields[1])
}

cell_types <- gsub('"', "", celltype_fields[-1], fixed = TRUE)

if (length(barcodes) != length(cell_types)) {
  stop(
    "Barcode and cell-type lengths do not match: ",
    length(barcodes), " versus ", length(cell_types)
  )
}

# The source dataset uses N-prefixed barcodes for normal samples.
# Confirm prefixes before assigning the remaining cells to tumour.
prefix_counts <- table(substr(barcodes, 1, 1))
print(prefix_counts)

if (!setequal(names(prefix_counts), c("N", "P"))) {
  stop("Unexpected barcode prefixes; review group assignment.")
}

groups <- ifelse(substr(barcodes, 1, 1) == "N", "Normal", "Tumor")

metadata <- data.frame(
  Cell = barcodes,
  Group = groups,
  CellType = cell_types,
  stringsAsFactors = FALSE
)

# Aggregate summaries only; no cell-level metadata is exported.
group_counts <- as.data.frame(
  table(Group = metadata$Group, CellType = metadata$CellType)
)

group_counts <- group_counts[group_counts$Freq > 0, ]

group_percentages <- aggregate(
  cbind(Percentage = rep(1, nrow(metadata))),
  by = list(Group = metadata$Group, CellType = metadata$CellType),
  FUN = sum
)


group_percentages$Percentage <- ave(
  group_percentages$Percentage,
  group_percentages$Group,
  FUN = function(x) 100 * x / sum(x)
)

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

write.csv(
  group_counts,
  "results/tables/source_group_celltype_counts.csv",
  row.names = FALSE
)

write.csv(
  group_percentages,
  "results/tables/source_group_celltype_percentages.csv",
  row.names = FALSE
)

cat("Source cells:", nrow(metadata), "\n")
cat("Distinct cell types:", length(unique(metadata$CellType)), "\n")
cat("Group counts:\n")
print(table(metadata$Group))
cat("Wrote aggregate summary tables to results/tables/.\n")
