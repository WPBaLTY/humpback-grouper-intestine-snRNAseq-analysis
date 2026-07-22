options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
})

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", script_args, value = TRUE)[1])
script_dir <- if (length(script_file) == 0 || is.na(script_file)) {
  getwd()
} else {
  dirname(normalizePath(script_file, winslash = "/", mustWork = FALSE))
}
repo_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = FALSE)

path_env <- function(name, default) {
  gsub("\\", "/", Sys.getenv(name, unset = default), fixed = TRUE)
}

tenx_root <- path_env("GROUPER_10X_ROOT", file.path(repo_root, "external_inputs", "10x_matrices"))
sample_sheet_path <- path_env(
  "GROUPER_10X_SAMPLE_SHEET",
  file.path(repo_root, "metadata", "tenx_sample_sheet.csv")
)
final_rds <- path_env(
  "GROUPER_FINAL_GLOBAL_RDS",
  file.path(repo_root, "external_inputs", "seurat_FINAL_celltype_fixedMeta.rds")
)
audit_out <- path_env(
  "GROUPER_ANNOTATION_AUDIT_OUT",
  file.path(repo_root, "audit", "annotation")
)

dir.create(audit_out, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(sample_sheet_path), file.exists(final_rds))

safe_join <- function(obj, assay = "RNA") {
  if (!(assay %in% Assays(obj))) return(obj)
  DefaultAssay(obj) <- assay
  layers <- tryCatch(Layers(obj[[assay]]), error = function(e) character(0))
  if (length(layers) > 1) obj <- JoinLayers(obj, assay = assay)
  obj
}

resolve_10x_dir <- function(root, primary, alias = NA_character_, delivery = NA_character_) {
  candidates <- c(primary, alias, delivery)
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  suffixes <- c(
    "",
    "filtered_feature_bc_matrix",
    file.path("outs", "filtered_feature_bc_matrix"),
    "raw_feature_bc_matrix",
    file.path("outs", "raw_feature_bc_matrix")
  )
  for (candidate in candidates) {
    for (suffix in suffixes) {
      path <- file.path(root, candidate, suffix)
      matrix_files <- file.path(path, c("matrix.mtx", "matrix.mtx.gz"))
      if (dir.exists(path) && any(file.exists(matrix_files))) {
        return(gsub("\\", "/", path, fixed = TRUE))
      }
    }
  }
  stop("Could not find a 10x matrix directory for sample ", primary)
}

sample_info <- read.csv(sample_sheet_path)
if ("analysis_order" %in% colnames(sample_info)) {
  sample_info <- sample_info[order(sample_info$analysis_order), , drop = FALSE]
  rownames(sample_info) <- NULL
}
sample_info$folder_alias <- if ("folder_alias" %in% colnames(sample_info)) {
  sample_info$folder_alias
} else {
  NA_character_
}
sample_info$tenx_dir <- mapply(
  resolve_10x_dir,
  root = tenx_root,
  primary = sample_info$folder_primary,
  alias = sample_info$folder_alias,
  delivery = sample_info$final_label
)

obj_list <- lapply(seq_len(nrow(sample_info)), function(i) {
  counts <- Read10X(data.dir = sample_info$tenx_dir[i])
  obj <- CreateSeuratObject(
    counts = counts,
    project = sample_info$treatment[i],
    min.cells = 3,
    min.features = 200
  )
  obj$sample <- sample_info$sample[i]
  obj$sample_plot <- sample_info$final_label[i]
  obj$group <- sample_info$treatment[i]
  obj$replicate <- sample_info$replicate[i]
  RenameCells(obj, add.cell.id = sample_info$sample[i])
})

reconstructed <- merge(obj_list[[1]], y = obj_list[-1])
rm(obj_list)
gc()
reconstructed <- safe_join(reconstructed, "RNA")

genes <- rownames(reconstructed)
mito_true <- unique(c(
  grep("^(MT-|mt-)", genes, value = TRUE),
  grep(
    "^(COX1|COX2|COX3|CYTB|ATP6|ATP8|ND1|ND2|ND3|ND4L|ND4|ND5|ND6|12S|16S)$",
    genes,
    value = TRUE,
    ignore.case = TRUE
  )
))
reconstructed[["percent.mt.true"]] <- PercentageFeatureSet(reconstructed, features = mito_true)

min_feat <- 300
max_feat <- as.numeric(quantile(reconstructed$nFeature_RNA, 0.995))
max_cnt <- as.numeric(quantile(reconstructed$nCount_RNA, 0.995))
max_mt <- 2
reconstructed <- subset(
  reconstructed,
  subset = nFeature_RNA >= min_feat &
    nFeature_RNA <= max_feat &
    nCount_RNA <= max_cnt &
    percent.mt.true <= max_mt
)
reconstructed <- safe_join(reconstructed, "RNA")

final_obj <- safe_join(readRDS(final_rds), "RNA")
final_counts <- GetAssayData(final_obj, assay = "RNA", layer = "counts")
reconstructed_counts <- GetAssayData(reconstructed, assay = "RNA", layer = "counts")

final_public_cells <- sub("^CON([0-9]+)_", "COM\\1_", colnames(final_counts))
cell_set_match <- setequal(final_public_cells, colnames(reconstructed_counts))
feature_order_match <- identical(rownames(final_counts), rownames(reconstructed_counts))

if (!cell_set_match) stop("Filtered cell sets do not match between raw reconstruction and final object")
if (!feature_order_match) stop("Feature order does not match between raw reconstruction and final object")

column_order <- match(final_public_cells, colnames(reconstructed_counts))
reconstructed_counts <- reconstructed_counts[, column_order, drop = FALSE]

count_slots_match <-
  identical(final_counts@Dim, reconstructed_counts@Dim) &&
  identical(final_counts@i, reconstructed_counts@i) &&
  identical(final_counts@p, reconstructed_counts@p) &&
  isTRUE(all.equal(final_counts@x, reconstructed_counts@x, tolerance = 0))

sample_public_final <- sub("^CON", "COM", as.character(final_obj$sample))
final_sample_counts <- as.data.frame(table(sample = sample_public_final), stringsAsFactors = FALSE)
reconstructed_sample_counts <- as.data.frame(
  table(sample = reconstructed$sample),
  stringsAsFactors = FALSE
)
colnames(final_sample_counts)[2] <- "final_object_nuclei"
colnames(reconstructed_sample_counts)[2] <- "raw_reconstruction_nuclei"
sample_counts <- full_join(final_sample_counts, reconstructed_sample_counts, by = "sample") %>%
  mutate(sample_nuclei_match = final_object_nuclei == raw_reconstruction_nuclei)

summary_rows <- data.frame(
  check = c(
    "feature_order",
    "filtered_cell_set",
    "sparse_count_matrix_slots",
    "all_sample_nuclei_counts",
    "overall_nuclei_count"
  ),
  status = c(
    ifelse(feature_order_match, "pass", "fail"),
    ifelse(cell_set_match, "pass", "fail"),
    ifelse(count_slots_match, "pass", "fail"),
    ifelse(all(sample_counts$sample_nuclei_match), "pass", "fail"),
    ifelse(ncol(final_obj) == ncol(reconstructed), "pass", "fail")
  ),
  details = c(
    paste0(nrow(final_counts), " features in identical order"),
    paste0(ncol(final_counts), " filtered nuclei with identical public cell IDs"),
    paste0("Dim/i/p/x slots exact: ", count_slots_match),
    paste(sample_counts$sample, sample_counts$final_object_nuclei, sep = "=", collapse = "; "),
    paste0("final=", ncol(final_obj), "; reconstructed=", ncol(reconstructed))
  ),
  stringsAsFactors = FALSE
)

write.csv(
  summary_rows,
  file.path(audit_out, "raw_10x_to_final_object_counts_validation.csv"),
  row.names = FALSE,
  quote = TRUE
)
write.csv(
  sample_counts,
  file.path(audit_out, "raw_10x_to_final_object_sample_counts.csv"),
  row.names = FALSE
)

if (!count_slots_match) stop("Sparse RNA count matrices are not exact")
message("PASS: raw 10x reconstruction and final annotated object have identical filtered RNA count matrices")
