suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("usage: Rscript export_pseudobulk_source.R FINAL_GLOBAL.rds FINAL_UMAP.csv.gz OUTPUT_DIR")
}

input_rds <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
final_umap_path <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(args[[3]], winslash = "/", mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(input_rds)
counts <- LayerData(obj, assay = "RNA", layer = "counts")
metadata <- obj[[]]

if (nrow(counts) != 23740L || ncol(counts) != 102036L) {
  stop(sprintf("Unexpected final object dimensions: %d x %d", nrow(counts), ncol(counts)))
}
if (!identical(colnames(counts), rownames(metadata))) {
  stop("Count columns and metadata rows are not in the same order")
}
if (!"sample" %in% colnames(metadata)) {
  stop("Final object metadata does not contain sample")
}

public_cell_id <- function(values) {
  values <- sub("^CON1_", "COM1_", values)
  sub("^CON2_", "COM2_", values)
}
public_sample <- function(values) {
  values <- sub("^CON1$", "COM1", values)
  sub("^CON2$", "COM2", values)
}

final_umap <- read.csv(gzfile(final_umap_path), check.names = FALSE, stringsAsFactors = FALSE)
object_ids <- public_cell_id(colnames(counts))
if (length(object_ids) != nrow(final_umap) || !setequal(object_ids, final_umap$cell_id)) {
  stop("Final RDS cell IDs do not close to the deposited final UMAP table")
}

sample_order <- c("CTL1", "CTL2", "Llac1", "Llac2", "Slim1", "Slim2", "COM1", "COM2")
sample_plot <- c("CTL_1", "CTL_2", "Llac_1", "Llac_2", "Slim_1", "Slim_2", "COM_1", "COM_2")
samples <- public_sample(as.character(metadata$sample))
if (!setequal(unique(samples), sample_order)) {
  stop(sprintf("Unexpected samples: %s", paste(sort(unique(samples)), collapse = ", ")))
}

sample_factor <- factor(samples, levels = sample_order)
design <- sparse.model.matrix(~ 0 + sample_factor)
colnames(design) <- sample_plot
pseudobulk_counts <- counts %*% design
library_sizes <- Matrix::colSums(pseudobulk_counts)
keep <- Matrix::rowSums(pseudobulk_counts) > 0
pseudobulk_counts <- as.matrix(pseudobulk_counts[keep, , drop = FALSE])
log_cpm <- log2(sweep(pseudobulk_counts, 2, library_sizes, "/") * 1e6 + 1)
correlation <- cor(log_cpm, method = "pearson")

log_cpm_out <- data.frame(feature_id = rownames(log_cpm), log_cpm, check.names = FALSE)
correlation_out <- data.frame(sample_plot = rownames(correlation), correlation, check.names = FALSE)
parameters_out <- data.frame(
  parameter = c(
    "input_features",
    "input_nuclei",
    "retained_expressed_features",
    "aggregation",
    "normalization",
    "transformation",
    "correlation"
  ),
  value = c(
    nrow(counts),
    ncol(counts),
    nrow(log_cpm),
    "sum raw RNA counts by capture library",
    "counts per million using total pseudobulk library size",
    "log2(CPM + 1)",
    "Pearson correlation across retained expressed features"
  ),
  stringsAsFactors = FALSE
)

write.csv(log_cpm_out, file.path(output_dir, "SuppFigureS3_pseudobulk_logCPM.csv"), row.names = FALSE)
write.csv(correlation_out, file.path(output_dir, "SuppFigureS3_pseudobulk_correlation.csv"), row.names = FALSE)
write.csv(parameters_out, file.path(output_dir, "SuppFigureS3_pseudobulk_parameters.csv"), row.names = FALSE)

cat(sprintf("Final object: %d features x %d nuclei\n", nrow(counts), ncol(counts)))
cat(sprintf("Pseudobulk matrix: %d expressed features x %d libraries\n", nrow(log_cpm), ncol(log_cpm)))
cat(sprintf("Correlation range off diagonal: %.6f to %.6f\n",
            min(correlation[row(correlation) != col(correlation)]),
            max(correlation[row(correlation) != col(correlation)])))
