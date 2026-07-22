options(stringsAsFactors = FALSE)

path_env <- function(name, default) {
  gsub("\\\\", "/", Sys.getenv(name, unset = default))
}

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", script_args, value = TRUE)[1])
script_dir <- if (length(script_file) == 0 || is.na(script_file)) {
  getwd()
} else {
  dirname(gsub("\\\\", "/", script_file))
}
repo_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = FALSE)

sibling_reconstruction <- file.path(script_dir, "01_reconstruct_immune_final_object.R")
default_reconstruction <- if (file.exists(sibling_reconstruction)) {
  sibling_reconstruction
} else {
  file.path(repo_root, "scripts", "R_analysis", "01_reconstruct_immune_final_object.R")
}
reconstruction_script <- path_env("GROUPER_IMMUNE_RECON_SCRIPT", default_reconstruction)
out_dir <- path_env(
  "GROUPER_IMMUNE_RECON_OUT",
  file.path(repo_root, "audit", "immune_filtering")
)
doublet_calls <- path_env(
  "GROUPER_DOUBLET_CALLS",
  file.path(
    repo_root, "audit", "doublet_assessment",
    "doublet_assessment_per_cell.csv.gz"
  )
)
final_immune_csv <- path_env(
  "GROUPER_FINAL_IMMUNE_CSV",
  file.path(
    repo_root, "source_data", "figures", "Fig4_immune_umap_coordinates.csv"
  )
)

if (!file.exists(reconstruction_script)) {
  stop("Missing immune reconstruction script: ", reconstruction_script)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

Sys.setenv(GROUPER_IMMUNE_RECON_OUT = out_dir)

message("Running the exact immune reconstruction before auditing exclusions")
sys.source(reconstruction_script, envir = globalenv(), keep.source = FALSE)

required_objects <- c("obj0", "obj1", "obj2", "obj3", "obj4")
missing_objects <- required_objects[!vapply(
  required_objects,
  function(object_name) exists(object_name, envir = globalenv(), inherits = FALSE),
  logical(1)
)]
if (length(missing_objects) > 0) {
  stop("Reconstruction did not expose expected stage objects: ", paste(missing_objects, collapse = ", "))
}
if (ncol(obj0) != 10440L || ncol(obj4) != 5061L) {
  stop("Expected immune stage sizes 10,440 and 5,061; observed ", ncol(obj0), " and ", ncol(obj4))
}

stage_specs <- list(
  list(
    stage = "step1_drop_initial_clusters_0_13",
    before = obj0,
    drop_clusters = c("0", "13"),
    after = obj1
  ),
  list(
    stage = "step2_drop_reclustered_cluster_1",
    before = obj1,
    drop_clusters = "1",
    after = obj2
  ),
  list(
    stage = "step3_drop_reclustered_clusters_8_11",
    before = obj2,
    drop_clusters = c("8", "11"),
    after = obj3
  ),
  list(
    stage = "step4_drop_reclustered_cluster_3",
    before = obj3,
    drop_clusters = "3",
    after = obj4
  )
)

cluster_vector <- function(obj) {
  setNames(as.character(Idents(obj)), colnames(obj))
}

transition_rows <- list()
removed_rows <- list()
all_cluster_qc <- list()
sample_rows <- list()
umap_rows <- list()

numeric_summary <- function(x, prefix) {
  x <- as.numeric(x)
  stats <- c(
    mean = mean(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    min = min(x, na.rm = TRUE),
    q25 = unname(quantile(x, 0.25, na.rm = TRUE)),
    median = median(x, na.rm = TRUE),
    q75 = unname(quantile(x, 0.75, na.rm = TRUE)),
    max = max(x, na.rm = TRUE)
  )
  setNames(as.list(stats), paste0(prefix, "_", names(stats)))
}

qc_row <- function(meta, stage, cluster, disposition) {
  data.frame(
    stage = stage,
    cluster = cluster,
    disposition = disposition,
    n_cells = nrow(meta),
    numeric_summary(meta[["nCount_RNA"]], "nCount_RNA"),
    numeric_summary(meta[["nFeature_RNA"]], "nFeature_RNA"),
    numeric_summary(meta[["percent.mt.true"]], "percent_mt"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

for (i in seq_along(stage_specs)) {
  spec <- stage_specs[[i]]
  before <- spec$before
  after <- spec$after
  before_clusters <- cluster_vector(before)
  drop_ids <- names(before_clusters)[before_clusters %in% spec$drop_clusters]
  retained_ids <- setdiff(colnames(before), drop_ids)

  if (!setequal(retained_ids, colnames(after))) {
    stop("Cell set mismatch after ", spec$stage)
  }

  transition_rows[[i]] <- data.frame(
    step = i,
    stage = spec$stage,
    input_nuclei = ncol(before),
    excluded_clusters = paste(spec$drop_clusters, collapse = ";"),
    excluded_nuclei = length(drop_ids),
    retained_nuclei = length(retained_ids),
    post_reclustering_nuclei = ncol(after),
    cumulative_excluded_nuclei = 10440L - ncol(after),
    stringsAsFactors = FALSE
  )

  meta <- before@meta.data
  meta$cell_id <- rownames(meta)
  removed_meta <- meta[drop_ids, , drop = FALSE]
  removed_rows[[i]] <- data.frame(
    cell_id = drop_ids,
    public_cell_id = sub("^CON", "COM", drop_ids),
    sample = as.character(removed_meta[["sample"]]),
    public_sample = sub("^CON", "COM", as.character(removed_meta[["sample"]])),
    group = as.character(removed_meta[["group"]]),
    public_group = sub("^CON$", "COM", as.character(removed_meta[["group"]])),
    removal_step = i,
    removal_stage = spec$stage,
    removal_cluster = unname(before_clusters[drop_ids]),
    nCount_RNA = as.numeric(removed_meta[["nCount_RNA"]]),
    nFeature_RNA = as.numeric(removed_meta[["nFeature_RNA"]]),
    percent_mt = as.numeric(removed_meta[["percent.mt.true"]]),
    stringsAsFactors = FALSE
  )

  for (cluster_id in sort(unique(before_clusters))) {
    ids <- names(before_clusters)[before_clusters == cluster_id]
    disposition <- if (cluster_id %in% spec$drop_clusters) "excluded" else "retained_at_step"
    all_cluster_qc[[length(all_cluster_qc) + 1L]] <- qc_row(
      meta[ids, , drop = FALSE], spec$stage, cluster_id, disposition
    )
  }
  all_cluster_qc[[length(all_cluster_qc) + 1L]] <- qc_row(
    meta[retained_ids, , drop = FALSE], spec$stage, "ALL_RETAINED", "retained_at_step"
  )

  sample_tab <- as.data.frame(
    table(
      sample = as.character(meta[drop_ids, "sample"]),
      cluster = unname(before_clusters[drop_ids])
    ),
    stringsAsFactors = FALSE
  )
  sample_tab <- sample_tab[sample_tab$Freq > 0, , drop = FALSE]
  sample_tab$public_sample <- sub("^CON", "COM", sample_tab$sample)
  sample_tab$stage <- spec$stage
  colnames(sample_tab)[colnames(sample_tab) == "Freq"] <- "excluded_nuclei"
  sample_rows[[i]] <- sample_tab[, c(
    "stage", "cluster", "sample", "public_sample", "excluded_nuclei"
  )]

  emb <- Embeddings(before, reduction = "umap")
  umap_rows[[i]] <- data.frame(
    stage = spec$stage,
    cell_id = rownames(emb),
    UMAP_1 = emb[, 1],
    UMAP_2 = emb[, 2],
    cluster = unname(before_clusters[rownames(emb)]),
    disposition = ifelse(rownames(emb) %in% drop_ids, "excluded", "retained_at_step"),
    stringsAsFactors = FALSE
  )
}

transition_ledger <- dplyr::bind_rows(transition_rows)
removed_cells <- dplyr::bind_rows(removed_rows)
cluster_qc <- dplyr::bind_rows(all_cluster_qc)
excluded_by_sample <- dplyr::bind_rows(sample_rows)
umap_source <- dplyr::bind_rows(umap_rows)

if (nrow(removed_cells) != 5379L || anyDuplicated(removed_cells$cell_id)) {
  stop("Removal ledger must contain 5,379 unique cell IDs")
}
if (!setequal(c(removed_cells$cell_id, colnames(obj4)), colnames(obj0))) {
  stop("Removed and final cell IDs do not partition the 10,440-cell input")
}

stage_cluster_maps <- list(
  step0_cluster = cluster_vector(obj0),
  step1_cluster = cluster_vector(obj1),
  step2_cluster = cluster_vector(obj2),
  step3_cluster = cluster_vector(obj3),
  final_cluster = cluster_vector(obj4)
)
initial_meta <- obj0@meta.data
fate <- data.frame(
  cell_id = colnames(obj0),
  public_cell_id = sub("^CON", "COM", colnames(obj0)),
  sample = as.character(initial_meta[colnames(obj0), "sample"]),
  public_sample = sub("^CON", "COM", as.character(initial_meta[colnames(obj0), "sample"])),
  group = as.character(initial_meta[colnames(obj0), "group"]),
  public_group = sub("^CON$", "COM", as.character(initial_meta[colnames(obj0), "group"])),
  nCount_RNA = as.numeric(initial_meta[colnames(obj0), "nCount_RNA"]),
  nFeature_RNA = as.numeric(initial_meta[colnames(obj0), "nFeature_RNA"]),
  percent_mt = as.numeric(initial_meta[colnames(obj0), "percent.mt.true"]),
  stringsAsFactors = FALSE
)
for (column_name in names(stage_cluster_maps)) {
  fate[[column_name]] <- unname(stage_cluster_maps[[column_name]][fate$cell_id])
}
removal_match <- match(fate$cell_id, removed_cells$cell_id)
fate$removal_step <- removed_cells$removal_step[removal_match]
fate$removal_stage <- removed_cells$removal_stage[removal_match]
fate$removal_cluster <- removed_cells$removal_cluster[removal_match]
fate$disposition <- ifelse(is.na(removal_match), "retained_final_5061", "excluded")

final_annotations <- setNames(as.character(obj4$immune_anno_final), colnames(obj4))
fate$final_immune_annotation <- unname(final_annotations[fate$cell_id])

if (file.exists(final_immune_csv)) {
  final_public <- read.csv(final_immune_csv, check.names = FALSE)
  final_public_ids <- sub("^COM", "CON", final_public$cell_id)
  if (!setequal(final_public_ids, colnames(obj4))) {
    stop("Reconstructed final immune cells differ from the Figure 4 per-cell source table")
  }
  final_match <- match(colnames(obj4), final_public_ids)
  reconstructed_umap <- Embeddings(obj4, reduction = "umap")[, 1:2, drop = FALSE]
  public_umap <- as.matrix(final_public[final_match, c("UMAP_1", "UMAP_2")])
  max_umap_delta <- max(abs(reconstructed_umap - public_umap))
  annotation_mismatches <- sum(
    as.character(obj4$immune_anno_final) != final_public$immune_subtype_internal[final_match]
  )
  if (max_umap_delta > 1e-12 || annotation_mismatches != 0L) {
    stop(
      "Final immune reconstruction differs from Figure 4 source data: max UMAP delta=",
      format(max_umap_delta, scientific = TRUE),
      "; annotation mismatches=", annotation_mismatches
    )
  }
  write.csv(
    data.frame(
      check = c(
        "final_nuclei", "cell_set_match", "cell_order_after_matching",
        "annotation_mismatches", "maximum_absolute_umap_delta"
      ),
      value = c(ncol(obj4), TRUE, identical(colnames(obj4), final_public_ids[final_match]),
                annotation_mismatches, format(max_umap_delta, scientific = TRUE)),
      stringsAsFactors = FALSE
    ),
    file.path(out_dir, "immune_filtering_final_object_validation.csv"),
    row.names = FALSE
  )
}

if (file.exists(doublet_calls)) {
  message("Merging retrospective scDblFinder calls into the filtering audit")
  dbl <- read.csv(doublet_calls, check.names = FALSE)
  required_doublet_cols <- c(
    "cell_id", "cluster_aware_score", "cluster_aware_class", "random_score", "random_class"
  )
  missing_doublet_cols <- setdiff(required_doublet_cols, colnames(dbl))
  if (length(missing_doublet_cols) > 0) {
    stop("Doublet table is missing: ", paste(missing_doublet_cols, collapse = ", "))
  }
  dbl_match <- match(fate$cell_id, dbl$cell_id)
  if (anyNA(dbl_match)) stop("Some immune cells are absent from the doublet table")
  for (column_name in setdiff(required_doublet_cols, "cell_id")) {
    fate[[column_name]] <- dbl[[column_name]][dbl_match]
  }
}

write.csv(
  transition_ledger,
  file.path(out_dir, "immune_filtering_transition_ledger.csv"),
  row.names = FALSE
)
write.csv(
  excluded_by_sample,
  file.path(out_dir, "immune_filtering_excluded_counts_by_stage_cluster_sample.csv"),
  row.names = FALSE
)
write.csv(
  cluster_qc,
  file.path(out_dir, "immune_filtering_qc_by_stage_cluster.csv"),
  row.names = FALSE
)
write.csv(
  removed_cells,
  file.path(out_dir, "immune_filtering_excluded_5379_cells.csv"),
  row.names = FALSE
)

fate_con <- gzfile(file.path(out_dir, "immune_filtering_fate_10440_cells.csv.gz"), open = "wt")
write.csv(fate, fate_con, row.names = FALSE, quote = TRUE)
close(fate_con)

umap_con <- gzfile(file.path(out_dir, "immune_filtering_stage_umap_source.csv.gz"), open = "wt")
write.csv(umap_source, umap_con, row.names = FALSE, quote = TRUE)
close(umap_con)

if (all(c("cluster_aware_class", "random_class") %in% colnames(fate))) {
  doublet_summary <- dplyr::bind_rows(
    fate %>%
      dplyr::mutate(audit_group = ifelse(disposition == "excluded", removal_stage, disposition)) %>%
      dplyr::group_by(audit_group) %>%
      dplyr::summarise(
        n_cells = dplyr::n(),
        predicted_doublets = sum(cluster_aware_class == "doublet"),
        predicted_doublet_rate = mean(cluster_aware_class == "doublet"),
        median_score = median(cluster_aware_score),
        .groups = "drop"
      ) %>%
      dplyr::mutate(mode = "cluster_aware", .before = 1),
    fate %>%
      dplyr::mutate(audit_group = ifelse(disposition == "excluded", removal_stage, disposition)) %>%
      dplyr::group_by(audit_group) %>%
      dplyr::summarise(
        n_cells = dplyr::n(),
        predicted_doublets = sum(random_class == "doublet"),
        predicted_doublet_rate = mean(random_class == "doublet"),
        median_score = median(random_score),
        .groups = "drop"
      ) %>%
      dplyr::mutate(mode = "random", .before = 1)
  )
  write.csv(
    doublet_summary,
    file.path(out_dir, "immune_filtering_doublet_summary_by_removal_stage.csv"),
    row.names = FALSE
  )
}

marker_definitions <- data.frame(
  marker = c(
    "PTPRC", "LCP1", "BCL11B", "LCK", "ZAP70", "ITK",
    "PRF1", "GZMB", "TYROBP", "FCER1G",
    "CD79A", "PAX5", "EBF1", "BLNK",
    "CSF1R1", "CSF1R2", "SPI1", "MRC1", "C1QA", "AXL",
    "XCR1", "ZNF366", "WDFY4", "IRF8", "CIITA", "CD74",
    "CLEC10A", "CD209", "IRF4",
    "LRP4", "CPNE9", "LAMA3"
  ),
  lineage_panel = c(
    "pan_leukocyte", "pan_leukocyte",
    rep("T_cell", 4),
    rep("cytotoxic", 4),
    rep("B_cell", 4),
    rep("macrophage_monocyte", 6),
    rep("cDC1", 6),
    rep("DC2_APC", 3),
    rep("historical_contaminant_check", 3)
  ),
  pattern = paste0(
    "^",
    c(
      "ptprc", "lcp1", "bcl11b", "lck", "zap70", "itk",
      "prf1", "gzmb", "tyrobp", "fcer1g",
      "cd79a", "pax5", "ebf1", "blnk",
      "csf1r1", "csf1r2", "spi1", "mrc1", "c1qa", "axl",
      "xcr1", "znf366", "wdfy4", "irf8", "ciita", "cd74",
      "clec10a", "cd209", "irf4",
      "lrp4", "cpne9", "lama3"
    ),
    "(\\.|$)"
  ),
  stringsAsFactors = FALSE
)

all_features <- rownames(obj0[["RNA"]])
marker_definitions$feature <- vapply(
  marker_definitions$pattern,
  function(pattern) {
    hits <- grep(pattern, all_features, value = TRUE, ignore.case = TRUE, perl = TRUE)
    if (length(hits) == 0) return(NA_character_)
    hits[which.min(nchar(hits))]
  },
  character(1)
)
write.csv(
  marker_definitions,
  file.path(out_dir, "immune_filtering_marker_panel_feature_map.csv"),
  row.names = FALSE
)

marker_panel_rows <- list()
features_present <- marker_definitions$feature[!is.na(marker_definitions$feature)]
for (spec in stage_specs) {
  panel_obj <- spec$before
  panel_clusters <- cluster_vector(panel_obj)
  panel_matrix <- GetAssayData(panel_obj, assay = "RNA", layer = "data")
  for (cluster_id in sort(unique(panel_clusters))) {
    ids <- names(panel_clusters)[panel_clusters == cluster_id]
    values <- panel_matrix[features_present, ids, drop = FALSE]
    marker_panel_rows[[length(marker_panel_rows) + 1L]] <- data.frame(
      stage = spec$stage,
      cluster = cluster_id,
      disposition = ifelse(cluster_id %in% spec$drop_clusters, "excluded", "retained_at_step"),
      n_cells = length(ids),
      marker = marker_definitions$marker[match(features_present, marker_definitions$feature)],
      lineage_panel = marker_definitions$lineage_panel[match(
        features_present, marker_definitions$feature
      )],
      feature = features_present,
      fraction_detected = as.numeric(Matrix::rowMeans(values > 0)),
      average_log_normalized_expression = as.numeric(Matrix::rowMeans(values)),
      stringsAsFactors = FALSE
    )
  }
}
write.csv(
  dplyr::bind_rows(marker_panel_rows),
  file.path(out_dir, "immune_filtering_marker_panel_by_stage_cluster.csv"),
  row.names = FALSE
)

message("Calculating removed-cluster marker evidence after reconstruction is complete")
marker_tables <- list()
for (spec in stage_specs) {
  marker_obj <- spec$before
  DefaultAssay(marker_obj) <- "RNA"
  Idents(marker_obj) <- "seurat_clusters"
  for (cluster_id in spec$drop_clusters) {
    marker_result <- tryCatch(
      FindMarkers(
        marker_obj,
        ident.1 = cluster_id,
        ident.2 = setdiff(levels(Idents(marker_obj)), cluster_id),
        only.pos = TRUE,
        min.pct = 0.10,
        logfc.threshold = 0.10,
        verbose = FALSE
      ),
      error = function(e) {
        warning("FindMarkers failed for ", spec$stage, " cluster ", cluster_id, ": ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(marker_result) && nrow(marker_result) > 0) {
      marker_result$feature <- rownames(marker_result)
      effect_col <- intersect(c("avg_log2FC", "avg_logFC"), colnames(marker_result))[1]
      marker_result <- marker_result[order(-marker_result[[effect_col]], marker_result$p_val_adj), , drop = FALSE]
      marker_result$rank_by_effect <- seq_len(nrow(marker_result))
      marker_result$stage <- spec$stage
      marker_result$cluster <- cluster_id
      marker_tables[[length(marker_tables) + 1L]] <- marker_result
    }
  }
}

removed_markers <- dplyr::bind_rows(marker_tables)
if (nrow(removed_markers) > 0) {
  marker_cols <- c(
    "stage", "cluster", "rank_by_effect", "feature",
    setdiff(colnames(removed_markers), c("stage", "cluster", "rank_by_effect", "feature"))
  )
  removed_markers <- removed_markers[, marker_cols, drop = FALSE]
  write.csv(
    removed_markers,
    file.path(out_dir, "immune_filtering_excluded_cluster_markers_all.csv"),
    row.names = FALSE
  )
  top_markers <- removed_markers %>%
    dplyr::group_by(stage, cluster) %>%
    dplyr::slice_min(rank_by_effect, n = 25, with_ties = FALSE) %>%
    dplyr::ungroup()
  write.csv(
    top_markers,
    file.path(out_dir, "immune_filtering_excluded_cluster_top25_markers.csv"),
    row.names = FALSE
  )
}

flow_nodes <- data.frame(
  y = 5:1,
  label = c(
    "Initial global leukocyte subset\nn = 10,440",
    "After excluding clusters 0 and 13\nn = 6,782 (excluded 3,658)",
    "After excluding reclustered cluster 1\nn = 6,006 (excluded 776)",
    "After excluding reclustered clusters 8 and 11\nn = 5,701 (excluded 305)",
    "Final immune atlas after excluding cluster 3\nn = 5,061 (excluded 640)"
  ),
  stringsAsFactors = FALSE
)
flow_plot <- ggplot(flow_nodes, aes(x = 0, y = y, label = label)) +
  geom_label(fill = "white", color = "black", linewidth = 0.35, size = 3.2) +
  geom_segment(
    data = data.frame(y = c(4.72, 3.72, 2.72, 1.72), yend = c(4.28, 3.28, 2.28, 1.28)),
    aes(x = 0, xend = 0, y = y, yend = yend),
    inherit.aes = FALSE,
    arrow = grid::arrow(length = grid::unit(0.12, "inches")),
    linewidth = 0.45
  ) +
  xlim(-1, 1) +
  ylim(0.7, 5.3) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(10, 10, 10, 10)
  )

base_cairo_png_audit <- function(filename, width, height, units, res, ...) {
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

ggsave(
  file.path(out_dir, "immune_filtering_flowchart.png"),
  flow_plot, width = 7.2, height = 8, dpi = 300,
  device = base_cairo_png_audit, bg = "white"
)
ggsave(
  file.path(out_dir, "immune_filtering_flowchart.pdf"),
  flow_plot, width = 7.2, height = 8
)

umap_plot_source <- umap_source
umap_plot_source$disposition <- factor(
  umap_plot_source$disposition,
  levels = c("retained_at_step", "excluded")
)
umap_plot <- ggplot(umap_plot_source, aes(UMAP_1, UMAP_2, color = disposition)) +
  geom_point(size = 0.08, alpha = 0.75) +
  facet_wrap(~stage, ncol = 2, scales = "fixed") +
  scale_color_manual(values = c(retained_at_step = "#B8BDC5", excluded = "#C43C35")) +
  coord_equal() +
  theme_void() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 8),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(8, 8, 8, 8)
  )
ggsave(
  file.path(out_dir, "immune_filtering_excluded_clusters_umap.png"),
  umap_plot, width = 10, height = 8, dpi = 300,
  device = base_cairo_png_audit, bg = "white"
)
ggsave(
  file.path(out_dir, "immune_filtering_excluded_clusters_umap.pdf"),
  umap_plot, width = 10, height = 8
)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "immune_filtering_audit_sessionInfo.txt"))
message("Immune filtering audit completed: ", out_dir)
