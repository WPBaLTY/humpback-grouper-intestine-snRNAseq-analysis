options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

set.seed(1234)

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

input_rds <- path_env(
  "GROUPER_FINAL_GLOBAL_RDS",
  file.path(repo_root, "external_inputs", "seurat_FINAL_celltype_fixedMeta.rds")
)
annotation_path <- path_env(
  "GROUPER_GLOBAL_CLUSTER_ANNOTATION_MAP",
  file.path(repo_root, "metadata", "global_cluster_annotation_map.csv")
)
source_out <- path_env(
  "GROUPER_SOURCE_DATA_OUT",
  file.path(repo_root, "source_data", "figures")
)
audit_out <- path_env(
  "GROUPER_FINAL_GLOBAL_AUDIT_OUT",
  file.path(repo_root, "outputs", "global_atlas_final_export")
)

dir.create(source_out, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_out, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(input_rds), file.exists(annotation_path))
obj <- readRDS(input_rds)

required_meta <- c("sample", "group", "replicate", "cluster_res04", "celltype_clean")
missing_meta <- setdiff(required_meta, colnames(obj@meta.data))
if (length(missing_meta) > 0) {
  stop("Final global object is missing metadata columns: ", paste(missing_meta, collapse = ", "))
}
if (!("umap" %in% names(obj@reductions))) {
  stop("Final global object does not contain a UMAP reduction")
}

annotation <- read.csv(annotation_path, stringsAsFactors = FALSE)
required_annotation <- c("cell_type_raw", "global_cell_type_public")
missing_annotation <- setdiff(required_annotation, colnames(annotation))
if (length(missing_annotation) > 0) {
  stop("Annotation map is missing columns: ", paste(missing_annotation, collapse = ", "))
}

raw_map <- unique(annotation[, required_annotation])
duplicate_raw <- raw_map$cell_type_raw[duplicated(raw_map$cell_type_raw)]
if (length(duplicate_raw) > 0) {
  stop("Conflicting duplicated raw annotations: ", paste(unique(duplicate_raw), collapse = ", "))
}
raw_to_public <- setNames(raw_map$global_cell_type_public, raw_map$cell_type_raw)

meta <- obj@meta.data
meta$sample_public <- sub("^CON", "COM", as.character(meta$sample))
meta$group_public <- sub("^CON$", "COM", as.character(meta$group))
meta$cell_type_original <- as.character(meta$celltype_clean)
missing_cell_types <- setdiff(unique(meta$cell_type_original), names(raw_to_public))
if (length(missing_cell_types) > 0) {
  stop("Unmapped final cell types: ", paste(missing_cell_types, collapse = ", "))
}
meta$global_cell_type <- unname(raw_to_public[meta$cell_type_original])

umap <- as.data.frame(Embeddings(obj, "umap"))
colnames(umap)[1:2] <- c("UMAP_1", "UMAP_2")
cell_metadata <- data.frame(
  cell_id = sub("^CON([0-9]+)_", "COM\\1_", rownames(umap)),
  sample = meta$sample_public,
  sample_plot = sub("([A-Za-z]+)([0-9]+)$", "\\1_\\2", meta$sample_public),
  group = meta$group_public,
  replicate = as.character(meta$replicate),
  cluster_res04 = as.character(meta$cluster_res04),
  cell_type_original = meta$cell_type_original,
  global_cell_type = meta$global_cell_type,
  UMAP_1 = umap$UMAP_1,
  UMAP_2 = umap$UMAP_2,
  stringsAsFactors = FALSE
)

metadata_path <- file.path(source_out, "Fig3_global_umap_coordinates_final_annotation.csv.gz")
metadata_connection <- gzfile(metadata_path, open = "wt")
write.csv(cell_metadata, metadata_connection, row.names = FALSE, quote = TRUE)
close(metadata_connection)

composition <- cell_metadata %>%
  count(sample, group, replicate, global_cell_type, name = "n") %>%
  group_by(sample) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(
    sample_plot = sub("([A-Za-z]+)([0-9]+)$", "\\1_\\2", sample),
    group_plot = group
  ) %>%
  rename(cell_type = global_cell_type) %>%
  select(sample, group, replicate, cell_type, n, sample_plot, group_plot, prop)

write.csv(
  composition,
  file.path(source_out, "Fig5a_global_celltype_composition_leukocytes_merged.csv"),
  row.names = FALSE,
  quote = FALSE
)

cluster_sample_counts <- cell_metadata %>%
  count(sample, group, replicate, cluster_res04, cell_type_original, global_cell_type, name = "n")
write.csv(
  cluster_sample_counts,
  file.path(audit_out, "final_global_annotation_counts_by_cluster_sample.csv"),
  row.names = FALSE
)

obj$sample_public <- meta$sample_public
obj$global_cell_type_public <- factor(
  meta$global_cell_type,
  levels = unique(annotation$global_cell_type_public)
)
p_cluster <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "cluster_res04",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.15
) + ggtitle("UMAP by transcriptional cluster")
p_global <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "global_cell_type_public",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.15
) + ggtitle("UMAP by broad annotated category")

ggsave(
  file.path(source_out, "Fig3A_transcriptional_cluster_umap_source.png"),
  p_cluster,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(source_out, "Fig3B_broad_category_umap_source.png"),
  p_global,
  width = 9,
  height = 6,
  dpi = 300,
  bg = "white"
)

writeLines(capture.output(sessionInfo()), file.path(audit_out, "sessionInfo.txt"))
write.csv(
  data.frame(
    metric = c("n_features", "n_nuclei", "n_samples", "n_clusters_res0.4", "n_global_cell_types"),
    value = c(
      nrow(obj),
      ncol(obj),
      length(unique(cell_metadata$sample)),
      length(unique(cell_metadata$cluster_res04)),
      length(unique(cell_metadata$global_cell_type))
    )
  ),
  file.path(audit_out, "final_global_export_summary.csv"),
  row.names = FALSE
)

message("Final global per-cell source data: ", metadata_path)
message("Final Figure 5a composition and Figure 3 UMAP sources written to: ", source_out)
