options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(harmony)
  library(dplyr)
  library(ggplot2)
})

expected_core_versions <- c(
  Seurat = "5.3.0",
  SeuratObject = "5.2.0",
  harmony = "1.2.3"
)
actual_core_versions <- vapply(
  names(expected_core_versions),
  function(package) as.character(packageVersion(package)),
  character(1)
)
version_mismatches <- names(expected_core_versions)[
  actual_core_versions != expected_core_versions
]
allow_version_drift <- tolower(Sys.getenv("GROUPER_ALLOW_VERSION_DRIFT", unset = "0")) %in%
  c("1", "true", "yes")
if (length(version_mismatches) > 0 && !allow_version_drift) {
  mismatch_text <- paste0(
    version_mismatches,
    " expected ",
    expected_core_versions[version_mismatches],
    " but found ",
    actual_core_versions[version_mismatches],
    collapse = "; "
  )
  stop(
    "Core R package versions do not match the validated raw-10x reconstruction: ",
    mismatch_text,
    ". Run environment/install_exact_r_packages.R and use that library, or set ",
    "GROUPER_ALLOW_VERSION_DRIFT=1 only for a non-exact exploratory rerun."
  )
}

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

tenx_root <- path_env("GROUPER_10X_ROOT", file.path(repo_root, "external_inputs", "10x_matrices"))
sample_sheet_path <- path_env("GROUPER_10X_SAMPLE_SHEET", file.path(repo_root, "metadata", "tenx_sample_sheet.csv"))
cluster_annotation_path <- path_env("GROUPER_GLOBAL_CLUSTER_ANNOTATION_MAP", file.path(repo_root, "metadata", "global_cluster_annotation_map.csv"))
harmony_batch_variable <- Sys.getenv("GROUPER_HARMONY_BATCH_VAR", unset = "orig.ident")
if (harmony_batch_variable != "orig.ident" && !allow_version_drift) {
  stop(
    "Exact reconstruction requires GROUPER_HARMONY_BATCH_VAR=orig.ident. ",
    "Set GROUPER_ALLOW_VERSION_DRIFT=1 only for an explicitly non-exact rerun."
  )
}
out_dir <- path_env("GROUPER_GLOBAL_RECON_OUT", file.path(repo_root, "outputs", "global_atlas_reconstruction"))
source_out <- path_env("GROUPER_SOURCE_DATA_OUT", file.path(repo_root, "source_data", "figures"))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_out, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "00_reconstruct_global_atlas_from_10x.log")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

message("Repository root: ", repo_root)
message("10x matrix root: ", tenx_root)
message("Sample sheet: ", sample_sheet_path)
message("Global cluster annotation map: ", cluster_annotation_path)
message("Harmony batch variable: ", harmony_batch_variable)
message("Output directory: ", out_dir)
message("Source-data directory: ", source_out)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

safe_join <- function(obj, assay = "RNA") {
  if (!(assay %in% Assays(obj))) return(obj)
  DefaultAssay(obj) <- assay
  layers <- tryCatch(Layers(obj[[assay]]), error = function(e) character(0))
  if (length(layers) > 1) {
    obj <- JoinLayers(obj, assay = assay)
  }
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
      files <- file.path(path, c("matrix.mtx", "matrix.mtx.gz"))
      if (dir.exists(path) && any(file.exists(files))) return(gsub("\\", "/", path, fixed = TRUE))
    }
  }
  stop("Could not find a 10x matrix directory for sample ", primary,
       ". Checked folder names: ", paste(candidates, collapse = ", "),
       " under ", root)
}

resolve_features <- function(obj, genes) {
  features_all <- rownames(obj[["RNA"]])
  out <- character(0)
  for (gene in genes) {
    if (gene %in% features_all) {
      out <- c(out, gene)
      next
    }
    idx <- which(tolower(features_all) == tolower(gene))
    if (length(idx) > 0) {
      out <- c(out, features_all[idx][which.min(nchar(features_all[idx]))])
      next
    }
    gene_base <- sub("\\.\\d+$", "", gene)
    hit <- grep(paste0("^", gene_base, "(\\.|$)"), features_all, value = TRUE, ignore.case = TRUE, perl = TRUE)
    if (length(hit) > 0) out <- c(out, hit[which.min(nchar(hit))])
  }
  unique(out)
}

export_dotplot_source <- function(obj, labels, feature_order, filename) {
  features <- resolve_features(obj, feature_order)
  expr <- GetAssayData(obj, assay = "RNA", layer = "data")
  groups <- factor(labels, levels = unique(labels))
  rows <- list()
  for (feature in features) {
    values <- as.numeric(expr[feature, ])
    avg_by_group <- tapply(values, groups, mean)
    pct_by_group <- tapply(values > 0, groups, mean) * 100
    avg_by_group <- avg_by_group[levels(groups)]
    pct_by_group <- pct_by_group[levels(groups)]
    scaled <- as.numeric(scale(avg_by_group))
    scaled[is.na(scaled)] <- 0
    rows[[feature]] <- data.frame(
      cell_type = levels(groups),
      marker_gene = feature,
      percent_expressed = as.numeric(pct_by_group),
      average_expression = as.numeric(avg_by_group),
      average_expression_scaled = scaled,
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  write.csv(out, file.path(source_out, filename), row.names = FALSE, quote = TRUE)
  write.csv(data.frame(feature = features), file.path(source_out, sub("_source\\.csv$", "_features_used.csv", filename)), row.names = FALSE)
  out
}

stopifnot(file.exists(sample_sheet_path))
sample_info <- read.csv(sample_sheet_path)
required_cols <- c("sample", "folder_primary", "treatment", "replicate", "final_label")
missing_cols <- setdiff(required_cols, colnames(sample_info))
if (length(missing_cols) > 0) stop("Missing columns in sample sheet: ", paste(missing_cols, collapse = ", "))

sample_info$folder_alias <- if ("folder_alias" %in% colnames(sample_info)) sample_info$folder_alias else NA_character_
if ("analysis_order" %in% colnames(sample_info)) {
  sample_info <- sample_info[order(sample_info$analysis_order), , drop = FALSE]
  rownames(sample_info) <- NULL
}
sample_info$tenx_dir <- mapply(
  resolve_10x_dir,
  root = tenx_root,
  primary = sample_info$folder_primary,
  alias = sample_info$folder_alias,
  delivery = sample_info$final_label
)
write.csv(sample_info, file.path(out_dir, "resolved_10x_sample_sheet.csv"), row.names = FALSE, quote = TRUE)
print(sample_info[, c("sample", "treatment", "replicate", "final_label", "tenx_dir")])

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
names(obj_list) <- sample_info$sample

obj <- merge(obj_list[[1]], y = obj_list[-1])
rm(obj_list)
gc()
obj <- safe_join(obj, "RNA")

genes <- rownames(obj)
mito_true <- unique(c(
  grep("^(MT-|mt-)", genes, value = TRUE),
  grep("^(COX1|COX2|COX3|CYTB|ATP6|ATP8|ND1|ND2|ND3|ND4L|ND4|ND5|ND6|12S|16S)$",
       genes, value = TRUE, ignore.case = TRUE)
))
obj[["percent.mt.true"]] <- PercentageFeatureSet(obj, features = mito_true)
write.csv(data.frame(feature = mito_true), file.path(out_dir, "mitochondrial_features_used.csv"), row.names = FALSE)

qc_before <- obj@meta.data %>%
  group_by(sample, sample_plot, group, replicate) %>%
  summarise(
    nuclei_before_filtering = n(),
    median_nFeature_RNA_before = median(nFeature_RNA),
    median_nCount_RNA_before = median(nCount_RNA),
    median_percent_mt_true_before = median(percent.mt.true),
    .groups = "drop"
  )

min_feat <- 300
max_feat <- as.numeric(quantile(obj$nFeature_RNA, 0.995))
max_cnt <- as.numeric(quantile(obj$nCount_RNA, 0.995))
max_mt <- 2

thresholds <- data.frame(
  min_nFeature_RNA = min_feat,
  max_nFeature_RNA_0_995 = max_feat,
  max_nCount_RNA_0_995 = max_cnt,
  max_percent_mt_true = max_mt
)
write.csv(thresholds, file.path(out_dir, "qc_thresholds.csv"), row.names = FALSE)

obj <- subset(
  obj,
  subset = nFeature_RNA >= min_feat &
    nFeature_RNA <= max_feat &
    nCount_RNA <= max_cnt &
    percent.mt.true <= max_mt
)

qc_after <- obj@meta.data %>%
  group_by(sample, sample_plot, group, replicate) %>%
  summarise(
    nuclei_after_filtering = n(),
    median_nFeature_RNA_after = median(nFeature_RNA),
    median_nCount_RNA_after = median(nCount_RNA),
    median_percent_mt_true_after = median(percent.mt.true),
    .groups = "drop"
  )
qc_summary <- left_join(qc_before, qc_after, by = c("sample", "sample_plot", "group", "replicate"))
write.csv(qc_summary, file.path(out_dir, "qc_summary_by_library.csv"), row.names = FALSE)

DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
obj <- RunPCA(obj, npcs = 50, verbose = FALSE)
obj <- safe_join(obj, "RNA")
if (!(harmony_batch_variable %in% colnames(obj@meta.data))) {
  stop("Harmony batch variable is absent from object metadata: ", harmony_batch_variable)
}
obj <- RunHarmony(obj, group.by.vars = harmony_batch_variable, plot_convergence = TRUE)

dims_use <- 1:20
cluster_resolutions <- c(0.1, 0.2, 0.3, 0.4, 0.6)
analysis_parameters <- data.frame(
  parameter = c(
    "random_seed", "harmony_batch_variable", "pca_npcs", "harmony_umap_dims",
    "cluster_resolutions", "primary_cluster_resolution"
  ),
  value = c(
    "1234", harmony_batch_variable, "50", paste(dims_use, collapse = ","),
    paste(cluster_resolutions, collapse = ","), "0.4"
  ),
  stringsAsFactors = FALSE
)
write.csv(analysis_parameters, file.path(out_dir, "analysis_parameters.csv"), row.names = FALSE)
obj <- FindNeighbors(obj, reduction = "harmony", dims = dims_use, verbose = FALSE)
obj <- FindClusters(obj, resolution = cluster_resolutions, verbose = FALSE)
obj$cluster_res04 <- obj[["RNA_snn_res.0.4"]][, 1]
Idents(obj) <- "cluster_res04"
obj <- RunUMAP(obj, reduction = "harmony", dims = dims_use, umap.method = "uwot", metric = "cosine", n.threads = 1)

cluster_counts <- as.data.frame(table(cluster_res04 = obj$cluster_res04), stringsAsFactors = FALSE)
colnames(cluster_counts) <- c("cluster_res04", "n_cells")
write.csv(cluster_counts, file.path(out_dir, "cluster_res0.4_counts.csv"), row.names = FALSE)

markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
if (is.null(markers) || nrow(markers) == 0) {
  markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.10)
}
write.csv(markers, file.path(out_dir, "markers_res0.4_all.csv"), row.names = FALSE)
write.csv(markers, file.path(source_out, "Fig3_global_markers_res0.4_all.csv"), row.names = FALSE)

stopifnot(file.exists(cluster_annotation_path))
cluster_annotation <- read.csv(cluster_annotation_path, stringsAsFactors = FALSE)
required_anno_cols <- c("cluster", "global_cell_type_public")
missing_anno_cols <- setdiff(required_anno_cols, colnames(cluster_annotation))
if (length(missing_anno_cols) > 0) {
  stop("Missing columns in cluster annotation map: ", paste(missing_anno_cols, collapse = ", "))
}
cluster_to_cell_type <- setNames(
  cluster_annotation$global_cell_type_public,
  as.character(cluster_annotation$cluster)
)
missing_clusters <- setdiff(sort(unique(as.character(obj$cluster_res04))), names(cluster_to_cell_type))
if (length(missing_clusters) > 0) {
  stop("Unmapped clusters in cluster_to_cell_type: ", paste(missing_clusters, collapse = ", "))
}
obj$global_cell_type <- unname(cluster_to_cell_type[as.character(obj$cluster_res04)])
public_cell_type_levels <- unique(cluster_annotation$global_cell_type_public)
obj$global_cell_type <- factor(obj$global_cell_type, levels = public_cell_type_levels)

umap <- as.data.frame(Embeddings(obj, "umap"))
colnames(umap)[1:2] <- c("UMAP_1", "UMAP_2")
umap$cell_id <- rownames(umap)
umap$sample <- obj$sample
umap$sample_plot <- obj$sample_plot
umap$group <- obj$group
umap$cluster_res04 <- as.character(obj$cluster_res04)
umap$global_cell_type <- as.character(obj$global_cell_type)
write.csv(umap, file.path(source_out, "Fig3_global_umap_coordinates.csv"), row.names = FALSE, quote = TRUE)

composition <- obj@meta.data %>%
  count(sample, group, replicate, global_cell_type, name = "n") %>%
  group_by(sample) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(
    sample_plot = gsub("([A-Za-z]+)([0-9]+)$", "\\1_\\2", sample),
    sample_plot = sub("^COM_", "COM_", sample_plot),
    group_plot = group
  ) %>%
  rename(cell_type = global_cell_type)
write.csv(composition, file.path(out_dir, "composition_celltype_bySample.csv"), row.names = FALSE)
write.csv(composition, file.path(source_out, "Fig5a_global_celltype_composition_leukocytes_merged.csv"), row.names = FALSE)

global_marker_genes <- c(
  "fabp2", "cd36", "SI",
  "slc10a2", "CUBN", "LRP2",
  "muc2", "spdef", "FER1L6",
  "pou2f3", "avil", "Pik3ap1",
  "best4", "cftr", "slc20a1a",
  "neurod1", "scgn", "ISL1",
  "syt1", "elavl3", "phox2a",
  "ptprc", "lcp1", "BCL11B", "SATB1", "FCER1G", "XCR1", "Axl",
  "col1a1", "col1a2", "dcn",
  "pecam1", "cdh5", "kdrl",
  "tagln", "CNTNAP5", "Pld5.1",
  "cel", "Cela1.1", "cpa2"
)
dotplot_source <- export_dotplot_source(
  obj,
  as.character(obj$global_cell_type),
  global_marker_genes,
  "Fig3_global_marker_dotplot_source.csv"
)

p_cluster <- DimPlot(obj, reduction = "umap", group.by = "cluster_res04", label = TRUE, repel = TRUE, pt.size = 0.15) +
  ggtitle("UMAP by transcriptional cluster")
p_global <- DimPlot(obj, reduction = "umap", group.by = "global_cell_type", label = TRUE, repel = TRUE, pt.size = 0.15) +
  ggtitle("UMAP by broad annotated category")
ggsave(file.path(out_dir, "Fig3A_transcriptional_cluster_umap_source.png"), p_cluster, width = 8, height = 6, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "Fig3B_broad_category_umap_source.png"), p_global, width = 9, height = 6, dpi = 300, bg = "white")
ggsave(file.path(source_out, "Fig3A_transcriptional_cluster_umap_source.png"), p_cluster, width = 8, height = 6, dpi = 300, bg = "white")
ggsave(file.path(source_out, "Fig3B_broad_category_umap_source.png"), p_global, width = 9, height = 6, dpi = 300, bg = "white")

saveRDS(obj, file.path(out_dir, "GLOBAL_ATLAS_harmony_res0.4.rds"))

leukocyte_obj <- subset(obj, subset = global_cell_type == "Leukocytes")
saveRDS(leukocyte_obj, file.path(out_dir, "GLOBAL_LEUKOCYTES_subset_for_reclustering.rds"))

message("Done. Global atlas outputs written to: ", out_dir)
message("Small figure source-data files written to: ", source_out)
