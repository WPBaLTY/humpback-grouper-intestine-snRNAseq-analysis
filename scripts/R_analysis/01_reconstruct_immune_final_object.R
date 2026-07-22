options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(harmony)
  library(dplyr)
  library(ggplot2)
})

set.seed(1234)

required_core_versions <- c(
  Seurat = "5.3.0",
  SeuratObject = "5.2.0",
  harmony = "1.2.3"
)
installed_core_versions <- vapply(
  names(required_core_versions),
  function(package) as.character(packageVersion(package)),
  character(1)
)
version_mismatch <- installed_core_versions != required_core_versions
allow_version_drift <- identical(Sys.getenv("GROUPER_ALLOW_VERSION_DRIFT", unset = "0"), "1")
if (any(version_mismatch) && !allow_version_drift) {
  mismatch_text <- paste(
    names(required_core_versions)[version_mismatch],
    installed_core_versions[version_mismatch],
    "!=",
    required_core_versions[version_mismatch],
    collapse = "; "
  )
  stop(
    "Unvalidated core package versions: ", mismatch_text,
    ". Run environment/install_exact_r_packages.R and set R_LIBS_USER to the resulting .r-library. ",
    "Set GROUPER_ALLOW_VERSION_DRIFT=1 only for an explicitly non-exact exploratory run."
  )
}

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", script_args, value = TRUE)[1])
script_dir <- if (length(script_file) == 0 || is.na(script_file)) getwd() else dirname(normalizePath(script_file, winslash = "/", mustWork = FALSE))
repo_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = FALSE)
path_env <- function(name, default) gsub("\\\\", "/", Sys.getenv(name, unset = default))

portable_path <- function(path, env_name, kind = "input") {
  normalized <- normalizePath(path, winslash = "/", mustWork = FALSE)
  root_prefix <- paste0(tolower(sub("/+$", "", repo_root)), "/")
  if (startsWith(tolower(normalized), root_prefix)) {
    return(substring(normalized, nchar(root_prefix) + 1L))
  }
  paste0("external ", kind, " via ", env_name, ": ", basename(normalized))
}

in_rds <- path_env(
  "GROUPER_IMMUNE_INPUT_RDS",
  file.path(repo_root, "outputs", "global_atlas_reconstruction", "GLOBAL_LEUKOCYTES_subset_for_reclustering.rds")
)
out_dir <- path_env("GROUPER_IMMUNE_RECON_OUT", file.path(repo_root, "outputs", "immune_final_reconstruction"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

input_stage <- tolower(Sys.getenv("GROUPER_IMMUNE_INPUT_STAGE", unset = "auto"))
if (!(input_stage %in% c("auto", "global_leukocytes", "preclustered"))) {
  stop("GROUPER_IMMUNE_INPUT_STAGE must be auto, global_leukocytes, or preclustered")
}

safe_join <- function(obj, assay = "RNA") {
  if (!(assay %in% Assays(obj))) return(obj)
  DefaultAssay(obj) <- assay
  ly <- tryCatch(Layers(obj[[assay]]), error = function(e) character(0))
  if (length(ly) > 1) {
    obj <- JoinLayers(obj, assay = assay)
  }
  obj
}

rerun_basic <- function(obj, dims_use = 1:10, resolution_use = 0.2) {
  obj <- safe_join(obj, "RNA")
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = 3000, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 50, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = dims_use, verbose = FALSE)
  obj <- FindClusters(obj, resolution = resolution_use, verbose = FALSE)
  obj <- RunUMAP(obj, dims = dims_use, verbose = FALSE)
  obj
}

prepare_initial_immune <- function(obj) {
  required_meta <- c("sample", "group")
  missing_meta <- setdiff(required_meta, colnames(obj@meta.data))
  if (length(missing_meta) > 0) {
    stop("Global leukocyte input is missing metadata: ", paste(missing_meta, collapse = ", "))
  }

  obj <- safe_join(obj, "RNA")
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = 3000, verbose = FALSE)
  obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
  obj <- RunPCA(obj, npcs = 50, verbose = FALSE)
  obj <- RunHarmony(obj, group.by.vars = "sample", plot_convergence = FALSE, verbose = FALSE)
  obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:10, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.2, verbose = FALSE)
  obj$imm_cluster <- as.character(obj[["RNA_snn_res.0.2"]][, 1])
  Idents(obj) <- "imm_cluster"
  obj <- RunUMAP(
    obj,
    reduction = "harmony",
    dims = 1:10,
    umap.method = "uwot",
    metric = "cosine",
    n.threads = 1,
    verbose = FALSE
  )
  obj
}

write_cluster_counts <- function(obj, tag) {
  df <- as.data.frame(table(as.character(Idents(obj))), stringsAsFactors = FALSE)
  colnames(df) <- c("cluster", "n_cells")
  df$tag <- tag
  write.csv(df, file.path(out_dir, paste0(tag, "_cluster_counts.csv")), row.names = FALSE)
  df
}

pick_gene <- function(obj, pattern) {
  genes <- rownames(obj[["RNA"]])
  hit <- grep(pattern, genes, value = TRUE, ignore.case = TRUE, perl = TRUE)
  if (length(hit) == 0) return(NA_character_)
  hit[which.min(nchar(hit))]
}

base_cairo_png <- function(filename, width, height, units, res, ...) {
  grDevices::png(
    filename = filename,
    width = width,
    height = height,
    units = units,
    res = res,
    type = "cairo",
    ...
  )
}

save_umap <- function(obj, group_by, stem, width = 8.8, height = 6.8) {
  p <- DimPlot(obj, group.by = group_by, label = TRUE, repel = TRUE) + NoLegend()
  ggsave(
    file.path(out_dir, paste0(stem, ".png")), p,
    width = width, height = height, dpi = 300, device = base_cairo_png
  )
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), p, width = width, height = height)
}

stopifnot(file.exists(in_rds))
obj0 <- readRDS(in_rds)
obj0 <- safe_join(obj0, "RNA")

is_preclustered_input <- "imm_cluster" %in% colnames(obj0@meta.data)
run_initial_clustering <- input_stage == "global_leukocytes" ||
  (input_stage == "auto" && !is_preclustered_input)
if (run_initial_clustering) {
  obj0 <- prepare_initial_immune(obj0)
} else if (!is_preclustered_input) {
  stop("Preclustered immune input must contain the imm_cluster metadata column")
}
saveRDS(obj0, file.path(out_dir, "IMMUNE_INITIAL_harmony_res0.2.rds"))
Idents(obj0) <- "seurat_clusters"

resolved_input_stage <- if (run_initial_clustering) "global_leukocytes" else "preclustered"
write.csv(
  data.frame(
    parameter = c(
      "input_rds",
      "input_stage",
      "seed",
      "initial_harmony_batch_variable",
      "initial_dims",
      "initial_resolution",
      "cleanup_dims",
      "cleanup_resolution",
      names(required_core_versions)
    ),
    value = c(
      portable_path(in_rds, "GROUPER_IMMUNE_INPUT_RDS"),
      resolved_input_stage,
      "1234",
      "sample",
      "1:10",
      "0.2",
      "1:10",
      "0.2",
      unname(installed_core_versions)
    ),
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "immune_analysis_parameters.csv"),
  row.names = FALSE
)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "immune_rerun_sessionInfo.txt"))

write_cluster_counts(obj0, "step0_original")

obj1 <- subset(obj0, idents = c("0", "13"), invert = TRUE)
obj1 <- rerun_basic(obj1)
Idents(obj1) <- "seurat_clusters"
write_cluster_counts(obj1, "step1_drop0_13_recluster")

obj2 <- subset(obj1, idents = "1", invert = TRUE)
obj2 <- rerun_basic(obj2)
Idents(obj2) <- "seurat_clusters"
write_cluster_counts(obj2, "step2_drop1_recluster")

obj3 <- subset(obj2, idents = c("8", "11"), invert = TRUE)
obj3 <- rerun_basic(obj3)
Idents(obj3) <- "seurat_clusters"
write_cluster_counts(obj3, "step3_drop8_11_recluster")

obj4 <- subset(obj3, idents = "3", invert = TRUE)
Idents(obj4) <- "seurat_clusters"
write_cluster_counts(obj4, "step4_drop3")

base_map <- c(
  "0" = "T_cell_CCR7_like",
  "1" = "NK_like_cytotoxic",
  "2" = "Activated_lymphoid_CCL20_high",
  "4" = "Macrophage_Monocyte_Axl_Csf1r1",
  "5" = "Activated_T_RORA_stress",
  "6" = "Granulocyte_like_EPX_Ncf4_CYBB",
  "7" = "cDC1_XCR1_ZNF366",
  "8" = "MoDC_like_Cd209d_FN1_AOC3",
  "9" = "Cycling_immune_G2M"
)

missing_base <- setdiff(sort(unique(as.character(Idents(obj4)))), names(base_map))
if (length(missing_base) > 0) {
  stop("Unmapped immune clusters after reconstruction: ", paste(missing_base, collapse = ", "))
}

obj4$immune_anno_final <- unname(base_map[as.character(obj4$seurat_clusters)])

sub6 <- subset(obj4, idents = "6")
sub6 <- rerun_basic(sub6)
Idents(sub6) <- "seurat_clusters"
write_cluster_counts(sub6, "step5_cluster6_recluster")

g_cd79 <- pick_gene(sub6, "^cd79a(\\.|$)")
g_pax5 <- pick_gene(sub6, "^pax5(\\.|$)")
g_ebf1 <- pick_gene(sub6, "^ebf1(\\.|$)")
g_blnk <- pick_gene(sub6, "^blnk(\\.|$)")

gene_map <- data.frame(
  label = c("CD79A", "PAX5", "EBF1", "BLNK"),
  feature = c(g_cd79, g_pax5, g_ebf1, g_blnk),
  stringsAsFactors = FALSE
)
write.csv(gene_map, file.path(out_dir, "cluster6_Bcell_marker_feature_map.csv"), row.names = FALSE)

if (is.na(g_cd79) || is.na(g_pax5)) {
  stop("Could not find both CD79A and PAX5 features in cluster6 subclustering object.")
}

mat <- GetAssayData(sub6, assay = "RNA", layer = "data")
clu <- as.character(Idents(sub6))

prop_positive <- function(gene_name) {
  out <- tapply(mat[gene_name, ] > 0, clu, mean)
  out[order(as.numeric(names(out)))]
}

p_cd79 <- prop_positive(g_cd79)
p_pax5 <- prop_positive(g_pax5)

sub6_diag <- data.frame(
  subcluster = names(p_cd79),
  cd79_positive_prop = as.numeric(p_cd79),
  pax5_positive_prop = as.numeric(p_pax5[names(p_cd79)]),
  stringsAsFactors = FALSE
)
write.csv(sub6_diag, file.path(out_dir, "cluster6_Bcell_subcluster_diagnostics.csv"), row.names = FALSE)

preferred_b <- c("2", "3")
if (all(preferred_b %in% levels(Idents(sub6)))) {
  b_subclusters <- preferred_b
} else {
  b_subclusters <- intersect(
    names(p_cd79)[p_cd79 >= 0.15],
    names(p_pax5)[p_pax5 >= 0.15]
  )
  if (length(b_subclusters) == 0) {
    cells_cd79 <- colnames(mat)[mat[g_cd79, ] > 0]
    cells_pax5 <- colnames(mat)[mat[g_pax5, ] > 0]
    cells_double <- intersect(cells_cd79, cells_pax5)
    tab <- sort(table(clu[cells_double]), decreasing = TRUE)
    b_subclusters <- names(tab)[tab == max(tab)]
  }
}

writeLines(
  paste("Selected B-cell subclusters:", paste(b_subclusters, collapse = ", ")),
  con = file.path(out_dir, "cluster6_selected_B_subclusters.txt")
)

cells_B <- WhichCells(sub6, idents = b_subclusters)
cells_cluster6 <- WhichCells(obj4, idents = "6")

obj4$immune_anno_final[cells_B] <- "B_cell"
obj4$immune_anno_final[setdiff(cells_cluster6, cells_B)] <- "Granulocyte_like_EPX_Ncf4_CYBB"

final_levels <- c(
  "T_cell_CCR7_like",
  "Activated_T_RORA_stress",
  "Activated_lymphoid_CCL20_high",
  "NK_like_cytotoxic",
  "B_cell",
  "cDC1_XCR1_ZNF366",
  "Macrophage_Monocyte_Axl_Csf1r1",
  "MoDC_like_Cd209d_FN1_AOC3",
  "Granulocyte_like_EPX_Ncf4_CYBB",
  "Cycling_immune_G2M"
)

obj4$immune_anno_final <- factor(obj4$immune_anno_final, levels = final_levels)
Idents(obj4) <- "immune_anno_final"

final_counts <- obj4@meta.data %>%
  count(immune_anno_final, name = "n_cells") %>%
  arrange(desc(n_cells))
write.csv(final_counts, file.path(out_dir, "immune_final_counts.csv"), row.names = FALSE)

if (all(c("group", "immune_anno_final") %in% colnames(obj4@meta.data))) {
  final_by_group <- obj4@meta.data %>%
    mutate(group_plot = ifelse(group == "CON", "COM", as.character(group))) %>%
    count(group_plot, immune_anno_final, name = "n_cells") %>%
    group_by(group_plot) %>%
    mutate(prop = n_cells / sum(n_cells)) %>%
    ungroup()
  write.csv(final_by_group, file.path(out_dir, "immune_final_counts_by_group.csv"), row.names = FALSE)
}

saveRDS(obj4, file.path(out_dir, "IMMUNE_FINAL_annotated_with_Bcell.rds"))
saveRDS(sub6, file.path(out_dir, "IMMUNE_cluster6_subclustered.rds"))

save_umap(obj4, "immune_anno_final", "Immune_final_UMAP_by_annotation")
Idents(obj4) <- "seurat_clusters"
save_umap(obj4, "seurat_clusters", "Immune_final_UMAP_by_cluster")
Idents(obj4) <- "immune_anno_final"
save_umap(sub6, "seurat_clusters", "Immune_cluster6_subclusters", width = 6.5, height = 5.5)

marker_features <- c(g_cd79, g_pax5, g_ebf1, g_blnk)
marker_features <- unique(marker_features[!is.na(marker_features)])
if (length(marker_features) > 0) {
  p_feat <- FeaturePlot(sub6, features = marker_features, ncol = min(2, length(marker_features)))
  ggsave(
    file.path(out_dir, "Immune_cluster6_B_marker_featureplot.png"), p_feat,
    width = 8, height = 6, dpi = 300, device = base_cairo_png
  )
  ggsave(file.path(out_dir, "Immune_cluster6_B_marker_featureplot.pdf"), p_feat, width = 8, height = 6)
}

cat("Saved reconstructed immune final object to:\n", out_dir, "\n")
