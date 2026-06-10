# Export final public-facing metadata for the humpback grouper intestine snRNA-seq atlas.
#
# Inputs:
#   1. Main Seurat object with the final 27-cluster global annotation.
#   2. Immune-cell reclustering Seurat object with immune subtype labels.
#
# Outputs:
#   - cell_metadata.csv
#   - umap_coordinates.csv
#   - sample_metadata_core.csv
#   - cluster_annotation_counts.csv
#   - celltype_counts_by_library.csv
#   - metadata_export_summary.txt
#
# The script is written for repository reuse. By default it expects the companion
# data repository folder to sit next to this code folder. Set GROUPER_DATA_DIR,
# GROUPER_MAIN_RDS, GROUPER_IMMUNE_RDS or GROUPER_OUTPUT_DIR to override paths.

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(Seurat)
})

# Small helper used throughout the workflow. Failing early with a clear file
# path is more useful than a later readRDS/read.csv error.
require_file <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " is missing: ", path)
  }
  invisible(path)
}

# Resolve portable paths instead of using local absolute paths.
repo_dir <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
data_dir <- normalizePath(
  Sys.getenv("GROUPER_DATA_DIR", file.path(repo_dir, "..", "humpback_grouper_intestine_snRNAseq_processed_data")),
  winslash = "/",
  mustWork = FALSE
)
main_rds <- Sys.getenv(
  "GROUPER_MAIN_RDS",
  file.path(data_dir, "integrated_object", "grouper_intestine_snRNAseq_seurat.rds")
)
immune_rds <- Sys.getenv(
  "GROUPER_IMMUNE_RDS",
  file.path(data_dir, "integrated_object", "grouper_intestine_immune_reclustering_seurat.rds")
)
out_dir <- normalizePath(
  Sys.getenv("GROUPER_OUTPUT_DIR", file.path(data_dir, "analysis_outputs_27cluster")),
  winslash = "/",
  mustWork = FALSE
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
require_file(main_rds, "Main Seurat RDS")
require_file(immune_rds, "Immune reclustering Seurat RDS")

# Map historical sample codes from the original Seurat object to public-facing
# library IDs and GEO sample accessions. The original object used CON1/CON2 for
# the combined treatment; the manuscript and GEO records use COM_1/COM_2.
sample_map <- data.frame(
  original_sample_code = c("CON1", "CON2", "CTL1", "CTL2", "Llac1", "Llac2", "Slim1", "Slim2"),
  library_id = c("COM_1", "COM_2", "CTL_1", "CTL_2", "Llac_1", "Llac_2", "Slim_1", "Slim_2"),
  geo_sample_accession = c(
    "GSM9627310", "GSM9627311", "GSM9627312", "GSM9627313",
    "GSM9627314", "GSM9627315", "GSM9627316", "GSM9627317"
  ),
  diet_group = c("COM", "COM", "CTL", "CTL", "Llac", "Llac", "Slim", "Slim"),
  biological_replicate = c(1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L),
  pooled_fish_per_library = 3L,
  stringsAsFactors = FALSE
)

message("Reading main object: ", main_rds)
obj <- readRDS(main_rds)
md <- obj@meta.data

# Fail early if a deposited object does not contain the metadata columns used
# by the manuscript tables and figures.
required <- c(
  "sample", "group", "nCount_RNA", "nFeature_RNA", "percent.mt.true",
  "cluster_res04", "RNA_snn_res.0.3", "RNA_snn_res.0.4", "RNA_snn_res.0.6",
  "celltype_final", "celltype_clean", "celltype_detail"
)
missing <- setdiff(required, colnames(md))
if (length(missing) > 0) {
  stop("Main object is missing required metadata columns: ", paste(missing, collapse = ", "))
}
if (!"umap" %in% names(obj@reductions)) {
  stop("Main object lacks a UMAP reduction named 'umap'.")
}

message("Reading immune object: ", immune_rds)
immune <- readRDS(immune_rds)
imd <- immune@meta.data

# Immune annotations were generated in a separate immune-focused reclustering
# object. Accept either available immune annotation column for compatibility.
immune_col <- if ("immune_anno_final" %in% colnames(imd)) {
  "immune_anno_final"
} else if ("immune_anno" %in% colnames(imd)) {
  "immune_anno"
} else {
  stop("Immune object lacks immune_anno_final/immune_anno metadata.")
}
immune_subtype_by_cell <- as.character(imd[[immune_col]])
names(immune_subtype_by_cell) <- rownames(imd)

map_idx <- match(as.character(md$sample), sample_map$original_sample_code)
if (anyNA(map_idx)) {
  stop("Unmapped sample codes: ", paste(sort(unique(as.character(md$sample)[is.na(map_idx)])), collapse = ", "))
}

cell_barcode <- rownames(md)
barcode_10x <- sub("^[^_]+_", "", cell_barcode)

# Pull the UMAP embedding from the final global object and align it explicitly
# to the metadata row order.
umap <- as.data.frame(Embeddings(obj, "umap"))
umap <- umap[cell_barcode, , drop = FALSE]
colnames(umap)[seq_len(2)] <- c("umap_1", "umap_2")

# Match immune subtype labels back to the global object. Non-immune nuclei are
# left blank so downstream users can distinguish broad labels from immune labels.
immune_subtype <- immune_subtype_by_cell[cell_barcode]
immune_subtype[is.na(immune_subtype)] <- ""

# The pancreatic acinar-like cluster is retained but flagged for transparency
# because possible peri-intestinal tissue carryover should be considered.
tissue_carryover_flag <- as.character(md$celltype_clean) == "Acinar_like" |
  as.character(md$celltype_final) == "Acinar_like.CEL" |
  as.character(md$cluster_res04) == "14"
tissue_carryover_note <- ifelse(
  tissue_carryover_flag,
  "possible_peri_intestinal_carryover_acinar_like",
  ""
)

# Build a level-2 label that keeps broad cell identities for non-immune nuclei
# and immune subtype identities for leukocyte nuclei.
celltype_lvl2 <- as.character(md$celltype_clean)
is_immune <- immune_subtype != ""
celltype_lvl2[is_immune] <- paste0("Leukocyte_", immune_subtype[is_immune])
celltype_lvl2 <- gsub("[^A-Za-z0-9_+.-]+", "_", celltype_lvl2)

cell_metadata <- data.frame(
  cell_barcode = cell_barcode,
  barcode_10x = barcode_10x,
  library_id = sample_map$library_id[map_idx],
  original_sample_code = as.character(md$sample),
  geo_sample_accession = sample_map$geo_sample_accession[map_idx],
  diet_group = sample_map$diet_group[map_idx],
  original_group_code = as.character(md$group),
  biological_replicate = sample_map$biological_replicate[map_idx],
  pooled_fish_per_library = sample_map$pooled_fish_per_library[map_idx],
  nFeature_RNA = md$nFeature_RNA,
  nCount_RNA = md$nCount_RNA,
  percent_mt = md[["percent.mt.true"]],
  pass_qc = "True",
  global_cluster = as.character(md$cluster_res04),
  cluster_res03 = as.character(md[["RNA_snn_res.0.3"]]),
  cluster_res04 = as.character(md[["RNA_snn_res.0.4"]]),
  cluster_res06 = as.character(md[["RNA_snn_res.0.6"]]),
  celltype_final = as.character(md$celltype_final),
  celltype_clean = as.character(md$celltype_clean),
  tissue_carryover_flag = ifelse(tissue_carryover_flag, "True", "False"),
  celltype_detail = as.character(md$celltype_detail),
  immune_subtype = immune_subtype,
  celltype_lvl2 = celltype_lvl2,
  umap_1 = umap$umap_1,
  umap_2 = umap$umap_2,
  tissue_carryover_note = tissue_carryover_note,
  stringsAsFactors = FALSE
)

if (nrow(cell_metadata) != ncol(obj)) {
  stop("Metadata row count does not match the number of nuclei in the main Seurat object.")
}
if (anyNA(cell_metadata$umap_1) || anyNA(cell_metadata$umap_2)) {
  stop("UMAP coordinates contain missing values after alignment to metadata rows.")
}
if (!identical(rownames(md), cell_metadata$cell_barcode)) {
  stop("Cell barcode order changed unexpectedly during metadata export.")
}

# Full per-nucleus metadata table used as the main reusable metadata resource.
write.csv(cell_metadata, file.path(out_dir, "cell_metadata.csv"), row.names = FALSE, na = "")

# A compact UMAP table is useful for plotting without loading the full Seurat object.
write.csv(
  cell_metadata[, c("cell_barcode", "library_id", "diet_group", "global_cluster", "celltype_final", "celltype_clean", "tissue_carryover_flag", "immune_subtype", "umap_1", "umap_2", "tissue_carryover_note")],
  file.path(out_dir, "umap_coordinates.csv"),
  row.names = FALSE,
  na = ""
)

write.csv(sample_map, file.path(out_dir, "sample_metadata_core.csv"), row.names = FALSE)

# Summarize the final 27 global clusters for manuscript Table 4 and reuse.
cluster_counts <- aggregate(cell_barcode ~ global_cluster + celltype_final + celltype_clean, cell_metadata, length)
colnames(cluster_counts)[ncol(cluster_counts)] <- "n_nuclei"
cluster_counts$global_cluster_num <- as.integer(cluster_counts$global_cluster)
cluster_counts <- cluster_counts[order(cluster_counts$global_cluster_num, -cluster_counts$n_nuclei), ]
write.csv(cluster_counts[, c("global_cluster", "celltype_final", "celltype_clean", "n_nuclei")],
          file.path(out_dir, "cluster_annotation_counts.csv"), row.names = FALSE)

# Library-level abundance summaries are descriptive and should not be treated as
# formal diet-effect tests without replicate-aware modeling.
library_counts <- aggregate(cell_barcode ~ library_id + diet_group + biological_replicate + celltype_clean, cell_metadata, length)
colnames(library_counts)[ncol(library_counts)] <- "n_nuclei"
write.csv(library_counts, file.path(out_dir, "celltype_counts_by_library.csv"), row.names = FALSE)

# Small text summary for repository audits.
n_tissue_carryover_flag <- sum(cell_metadata$tissue_carryover_flag == "True")
cluster14_n <- sum(cell_metadata$global_cluster == "14")
summary_lines <- c(
  paste0("main_rds=", main_rds),
  paste0("immune_rds=", immune_rds),
  paste0("output_dir=", out_dir),
  paste0("n_nuclei=", nrow(cell_metadata)),
  paste0("n_global_clusters=", length(unique(cell_metadata$global_cluster))),
  paste0("global_clusters=", paste(sort(as.integer(unique(cell_metadata$global_cluster))), collapse = ",")),
  paste0("n_immune_subtyped=", sum(cell_metadata$immune_subtype != "")),
  paste0("n_celltype_clean=", length(unique(cell_metadata$celltype_clean))),
  paste0("n_tissue_carryover_flag=", n_tissue_carryover_flag),
  paste0("global_cluster_14_n=", cluster14_n)
)
writeLines(summary_lines, file.path(out_dir, "metadata_export_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
