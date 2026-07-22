options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(SingleCellExperiment)
  library(S4Vectors)
  library(scDblFinder)
  library(BiocParallel)
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
  gsub("\\\\", "/", Sys.getenv(name, unset = default))
}

portable_path <- function(path, env_name, kind = "input") {
  normalized <- normalizePath(path, winslash = "/", mustWork = FALSE)
  root_prefix <- paste0(tolower(sub("/+$", "", repo_root)), "/")
  if (startsWith(tolower(normalized), root_prefix)) {
    return(substring(normalized, nchar(root_prefix) + 1L))
  }
  paste0("external ", kind, " via ", env_name, ": ", basename(normalized))
}

global_rds <- path_env(
  "GROUPER_GLOBAL_FINAL_RDS",
  file.path(repo_root, "outputs", "global_atlas_reconstruction", "GLOBAL_ATLAS_harmony_res0.4.rds")
)
initial_immune_rds <- path_env(
  "GROUPER_INITIAL_IMMUNE_RDS",
  file.path(repo_root, "outputs", "immune_final_reconstruction", "IMMUNE_INITIAL_harmony_res0.2.rds")
)
final_immune_csv <- path_env(
  "GROUPER_FINAL_IMMUNE_CSV",
  file.path(
    repo_root, "source_data", "figures", "Fig4_immune_umap_coordinates.csv"
  )
)
out_dir <- path_env(
  "GROUPER_DOUBLET_OUT",
  file.path(repo_root, "audit", "doublet_assessment")
)

seed_primary <- as.integer(Sys.getenv("GROUPER_DOUBLET_SEED", unset = "20260711"))
seed_random <- seed_primary + 1L
dbr_per_1k <- as.numeric(Sys.getenv("GROUPER_DOUBLET_DBR_PER_1K", unset = "0.008"))
reuse_primary <- identical(Sys.getenv("GROUPER_REUSE_PRIMARY_DOUBLET_CALLS", unset = "0"), "1")

required_r_version <- "4.6.1"
required_package_versions <- c(
  Seurat = "5.5.1",
  SeuratObject = "5.4.0",
  SingleCellExperiment = "1.34.0",
  scDblFinder = "1.26.7",
  BiocParallel = "1.46.0"
)
installed_package_versions <- vapply(
  names(required_package_versions),
  function(package) as.character(packageVersion(package)),
  character(1)
)
allow_version_drift <- identical(
  Sys.getenv("GROUPER_ALLOW_DOUBLET_VERSION_DRIFT", unset = "0"),
  "1"
)
version_mismatch <- installed_package_versions != required_package_versions
if (
  (!identical(as.character(getRversion()), required_r_version) || any(version_mismatch)) &&
    !allow_version_drift
) {
  mismatch_text <- paste(
    names(required_package_versions)[version_mismatch],
    installed_package_versions[version_mismatch],
    "!=",
    required_package_versions[version_mismatch],
    collapse = "; "
  )
  stop(
    "Unvalidated doublet-assessment environment. Expected R ", required_r_version,
    " and package versions recorded in environment/install_doublet_r_packages.R. ",
    if (nzchar(mismatch_text)) paste0("Mismatches: ", mismatch_text, ". ") else "",
    "Set GROUPER_ALLOW_DOUBLET_VERSION_DRIFT=1 only for an explicitly non-exact exploratory run."
  )
}

required_files <- c(global_rds, initial_immune_rds, final_immune_csv)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required input(s): ", paste(missing_files, collapse = "; "))
}
if (is.na(seed_primary) || is.na(dbr_per_1k) || dbr_per_1k <= 0) {
  stop("GROUPER_DOUBLET_SEED and GROUPER_DOUBLET_DBR_PER_1K must be valid numeric values")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
run_started <- Sys.time()

message("Reading the 10,440-cell initial immune object")
initial_obj <- readRDS(initial_immune_rds)
initial_ids <- colnames(initial_obj)
if (length(initial_ids) != 10440L || anyDuplicated(initial_ids)) {
  stop("Initial immune object must contain 10,440 unique cell IDs")
}
rm(initial_obj)
invisible(gc())

message("Reading the public 5,061-cell Figure 4 coordinate table")
final_immune <- read.csv(final_immune_csv, check.names = FALSE)
required_final_cols <- c("cell_id", "immune_subtype_internal", "immune_subtype")
missing_final_cols <- setdiff(required_final_cols, colnames(final_immune))
if (length(missing_final_cols) > 0) {
  stop("Final immune table is missing: ", paste(missing_final_cols, collapse = ", "))
}

# Public display labels use COM1/COM2, whereas the underlying barcodes retain CON1/CON2.
final_immune$canonical_cell_id <- sub("^COM", "CON", final_immune$cell_id)
final_ids <- final_immune$canonical_cell_id
if (length(final_ids) != 5061L || anyDuplicated(final_ids)) {
  stop("Final immune table must contain 5,061 unique cell IDs")
}
if (!all(final_ids %in% initial_ids)) {
  stop("Final immune IDs not found in the initial immune object: ", sum(!(final_ids %in% initial_ids)))
}
if (sum(!(initial_ids %in% final_ids)) != 5379L) {
  stop("Expected exactly 5,379 initial immune cells to be absent from the final immune set")
}

final_subtype_internal <- setNames(final_immune$immune_subtype_internal, final_ids)
final_subtype_display <- setNames(final_immune$immune_subtype, final_ids)
rm(final_immune)
invisible(gc())

message("Reading the 102,036-cell final global object")
global_obj <- readRDS(global_rds)
if (!all(dim(global_obj) == c(23740, 102036))) {
  stop("Unexpected global object dimensions: ", paste(dim(global_obj), collapse = " x "))
}

required_meta <- c(
  "sample", "group", "nCount_RNA", "nFeature_RNA", "percent.mt.true",
  "RNA_snn_res.0.4", "cluster_res04", "celltype_final", "major_final"
)
missing_meta <- setdiff(required_meta, colnames(global_obj@meta.data))
if (length(missing_meta) > 0) {
  stop("Global object is missing metadata: ", paste(missing_meta, collapse = ", "))
}

cell_ids <- colnames(global_obj)
if (anyDuplicated(cell_ids)) stop("Global object contains duplicated cell IDs")
if (!all(initial_ids %in% cell_ids)) stop("Not all initial immune IDs occur in the global object")

global_meta <- global_obj@meta.data
sample_id <- as.character(global_meta[["sample"]])
expected_samples <- c("CON1", "CON2", "CTL1", "CTL2", "Llac1", "Llac2", "Slim1", "Slim2")
if (!setequal(unique(sample_id), expected_samples)) {
  stop("Unexpected sample labels: ", paste(sort(unique(sample_id)), collapse = ", "))
}

cluster_res04 <- as.character(global_meta[["cluster_res04"]])
cluster_direct <- as.character(global_meta[["RNA_snn_res.0.4"]])
if (!identical(cluster_res04, cluster_direct)) {
  stop("cluster_res04 and RNA_snn_res.0.4 are not cell-wise identical")
}
if (length(unique(cluster_res04)) != 27L) {
  stop("Expected 27 resolution-0.4 global clusters")
}

analysis_set <- ifelse(
  cell_ids %in% final_ids,
  "final_immune_5061",
  ifelse(cell_ids %in% initial_ids, "removed_from_initial_immune_5379", "outside_initial_immune")
)

base_meta <- DataFrame(
  cell_id = cell_ids,
  sample = sample_id,
  group = as.character(global_meta[["group"]]),
  global_cluster_res0.4 = cluster_res04,
  global_celltype = as.character(global_meta[["celltype_final"]]),
  global_major_type = as.character(global_meta[["major_final"]]),
  analysis_set = analysis_set,
  immune_subtype_internal = unname(final_subtype_internal[cell_ids]),
  immune_subtype = unname(final_subtype_display[cell_ids]),
  nCount_RNA = as.numeric(global_meta[["nCount_RNA"]]),
  nFeature_RNA = as.numeric(global_meta[["nFeature_RNA"]]),
  percent_mt = as.numeric(global_meta[["percent.mt.true"]]),
  row.names = cell_ids
)

counts_matrix <- GetAssayData(global_obj, assay = "RNA", layer = "counts")
if (!identical(colnames(counts_matrix), cell_ids)) {
  stop("Count matrix columns do not match global object cell IDs")
}

sce_base <- SingleCellExperiment(
  assays = list(counts = counts_matrix),
  colData = base_meta
)
rm(global_obj, global_meta, counts_matrix, base_meta, initial_ids, final_ids)
invisible(gc())

extract_scdblfinder <- function(sce, prefix) {
  cd <- as.data.frame(colData(sce))
  keep <- grep("^scDblFinder\\.", colnames(cd), value = TRUE)
  if (!all(c("scDblFinder.score", "scDblFinder.class") %in% keep)) {
    stop("scDblFinder output is missing score or class columns")
  }
  out <- cd[, keep, drop = FALSE]
  colnames(out) <- paste0(prefix, sub("^scDblFinder\\.", "_", colnames(out)))
  out
}

write_csv_gz <- function(x, path) {
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  write.csv(x, con, row.names = FALSE, quote = TRUE)
}

primary_cache <- file.path(out_dir, "doublet_calls_cluster_aware_per_cell.csv.gz")
if (reuse_primary && file.exists(primary_cache)) {
  message("Reusing the completed cluster-aware per-cell checkpoint")
  primary_checkpoint <- read.csv(primary_cache, check.names = FALSE)
  expected_ids <- as.character(colData(sce_base)$cell_id)
  if (
    nrow(primary_checkpoint) != ncol(sce_base) ||
      !identical(primary_checkpoint$cell_id, expected_ids)
  ) {
    stop("Cluster-aware checkpoint does not match the current global object")
  }
  primary_result <- primary_checkpoint[, grep(
    "^cluster_aware_", colnames(primary_checkpoint), value = TRUE
  ), drop = FALSE]
  rm(primary_checkpoint, expected_ids)
} else {
  message("Running cluster-aware scDblFinder by capture/library")
  set.seed(seed_primary)
  sce_primary <- scDblFinder(
    sce_base,
    clusters = "global_cluster_res0.4",
    samples = "sample",
    dbr.per1k = dbr_per_1k,
    multiSampleMode = "split",
    BPPARAM = SerialParam(RNGseed = seed_primary, progressbar = TRUE)
  )
  primary_result <- extract_scdblfinder(sce_primary, "cluster_aware")
  write_csv_gz(
    cbind(as.data.frame(colData(sce_base)), primary_result),
    primary_cache
  )
  rm(sce_primary)
  invisible(gc())
}

message("Running random-mode scDblFinder sensitivity analysis by capture/library")
set.seed(seed_random)
sce_random <- scDblFinder(
  sce_base,
  clusters = FALSE,
  samples = "sample",
  dbr.per1k = dbr_per_1k,
  multiSampleMode = "split",
  BPPARAM = SerialParam(RNGseed = seed_random, progressbar = TRUE)
)
random_result <- extract_scdblfinder(sce_random, "random")
write_csv_gz(
  cbind(as.data.frame(colData(sce_base)), random_result),
  file.path(out_dir, "doublet_calls_random_per_cell.csv.gz")
)
rm(sce_random)
invisible(gc())

cell_results <- cbind(
  as.data.frame(colData(sce_base)),
  primary_result,
  random_result
)
rm(sce_base, primary_result, random_result)
invisible(gc())

required_result_cols <- c(
  "cluster_aware_score", "cluster_aware_class", "random_score", "random_class"
)
if (!all(required_result_cols %in% colnames(cell_results))) {
  stop("Unexpected renamed result columns: ", paste(colnames(cell_results), collapse = ", "))
}

cell_results$cluster_aware_class <- as.character(cell_results$cluster_aware_class)
cell_results$random_class <- as.character(cell_results$random_class)
write_csv_gz(cell_results, file.path(out_dir, "doublet_assessment_per_cell.csv.gz"))

summarize_calls <- function(data, group_cols, class_col, score_col, mode_label) {
  data %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      n_cells = n(),
      predicted_doublets = sum(.data[[class_col]] == "doublet", na.rm = TRUE),
      predicted_doublet_rate = mean(.data[[class_col]] == "doublet", na.rm = TRUE),
      mean_score = mean(.data[[score_col]], na.rm = TRUE),
      median_score = median(.data[[score_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(mode = mode_label, .before = 1)
}

summaries_for <- function(group_cols) {
  bind_rows(
    summarize_calls(
      cell_results, group_cols,
      "cluster_aware_class", "cluster_aware_score", "cluster_aware"
    ),
    summarize_calls(cell_results, group_cols, "random_class", "random_score", "random")
  )
}

library_summary <- summaries_for("sample") %>%
  mutate(
    expected_rate_from_default = dbr_per_1k * n_cells / 1000,
    expected_doublets_from_default = expected_rate_from_default * n_cells
  )
write.csv(library_summary, file.path(out_dir, "doublet_summary_by_library.csv"), row.names = FALSE)
write.csv(
  summaries_for("global_cluster_res0.4"),
  file.path(out_dir, "doublet_summary_by_global_cluster.csv"),
  row.names = FALSE
)
write.csv(
  summaries_for("global_celltype"),
  file.path(out_dir, "doublet_summary_by_global_celltype.csv"),
  row.names = FALSE
)
write.csv(
  summaries_for("analysis_set"),
  file.path(out_dir, "doublet_summary_by_analysis_set.csv"),
  row.names = FALSE
)
write.csv(
  summaries_for(c("analysis_set", "sample")),
  file.path(out_dir, "doublet_summary_by_analysis_set_and_library.csv"),
  row.names = FALSE
)
write.csv(
  summaries_for("immune_subtype") %>% filter(!is.na(immune_subtype)),
  file.path(out_dir, "doublet_summary_by_final_immune_subtype.csv"),
  row.names = FALSE
)

concordance <- as.data.frame(
  table(
    cluster_aware = cell_results$cluster_aware_class,
    random = cell_results$random_class,
    useNA = "ifany"
  ),
  stringsAsFactors = FALSE
)
write.csv(concordance, file.path(out_dir, "doublet_mode_concordance.csv"), row.names = FALSE)

qc_by_call <- bind_rows(
  cell_results %>%
    group_by(predicted_class = cluster_aware_class) %>%
    summarise(
      n_cells = n(),
      median_nCount_RNA = median(nCount_RNA),
      median_nFeature_RNA = median(nFeature_RNA),
      median_percent_mt = median(percent_mt),
      .groups = "drop"
    ) %>% mutate(mode = "cluster_aware", .before = 1),
  cell_results %>%
    group_by(predicted_class = random_class) %>%
    summarise(
      n_cells = n(),
      median_nCount_RNA = median(nCount_RNA),
      median_nFeature_RNA = median(nFeature_RNA),
      median_percent_mt = median(percent_mt),
      .groups = "drop"
    ) %>% mutate(mode = "random", .before = 1)
)
write.csv(qc_by_call, file.path(out_dir, "doublet_qc_summary_by_call.csv"), row.names = FALSE)

compare_removed_to_final <- function(class_col, mode_label) {
  keep <- cell_results$analysis_set %in% c(
    "removed_from_initial_immune_5379", "final_immune_5061"
  )
  d <- cell_results[keep, , drop = FALSE]
  removed_doublet <- sum(
    d[[class_col]] == "doublet" & d$analysis_set == "removed_from_initial_immune_5379"
  )
  removed_singlet <- sum(
    d[[class_col]] == "singlet" & d$analysis_set == "removed_from_initial_immune_5379"
  )
  final_doublet <- sum(d[[class_col]] == "doublet" & d$analysis_set == "final_immune_5061")
  final_singlet <- sum(d[[class_col]] == "singlet" & d$analysis_set == "final_immune_5061")
  tab <- matrix(
    c(removed_doublet, removed_singlet, final_doublet, final_singlet),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("removed_5379", "final_5061"), c("doublet", "singlet"))
  )
  ft <- fisher.test(tab)
  data.frame(
    mode = mode_label,
    removed_doublets = removed_doublet,
    removed_total = removed_doublet + removed_singlet,
    removed_doublet_rate = removed_doublet / (removed_doublet + removed_singlet),
    final_doublets = final_doublet,
    final_total = final_doublet + final_singlet,
    final_doublet_rate = final_doublet / (final_doublet + final_singlet),
    odds_ratio_removed_vs_final = unname(ft$estimate),
    fisher_p_value = ft$p.value,
    stringsAsFactors = FALSE
  )
}

enrichment <- bind_rows(
  compare_removed_to_final("cluster_aware_class", "cluster_aware"),
  compare_removed_to_final("random_class", "random")
)
write.csv(
  enrichment,
  file.path(out_dir, "doublet_enrichment_removed_5379_vs_final_5061.csv"),
  row.names = FALSE
)

run_finished <- Sys.time()
versions <- c(
  R = paste(R.version$major, R.version$minor, sep = "."),
  Seurat = as.character(packageVersion("Seurat")),
  SeuratObject = as.character(packageVersion("SeuratObject")),
  SingleCellExperiment = as.character(packageVersion("SingleCellExperiment")),
  scDblFinder = as.character(packageVersion("scDblFinder")),
  BiocParallel = as.character(packageVersion("BiocParallel"))
)
parameters <- data.frame(
  parameter = c(
    "global_rds", "initial_immune_rds", "final_immune_csv", "output_directory",
    "primary_mode", "sensitivity_mode", "sample_column", "cluster_column",
    "multiSampleMode", "dbr.per1k", "primary_seed", "random_seed",
    "parallel_backend", "cluster_aware_checkpoint_reused", "run_started", "run_finished",
    names(versions)
  ),
  value = c(
    portable_path(global_rds, "GROUPER_GLOBAL_FINAL_RDS"),
    portable_path(initial_immune_rds, "GROUPER_INITIAL_IMMUNE_RDS"),
    portable_path(final_immune_csv, "GROUPER_FINAL_IMMUNE_CSV"),
    portable_path(out_dir, "GROUPER_DOUBLET_OUT", "output"),
    "cluster-aware", "random", "sample", "global_cluster_res0.4", "split",
    format(dbr_per_1k, scientific = FALSE), seed_primary, seed_random,
    "BiocParallel::SerialParam", as.character(reuse_primary),
    format(run_started, tz = "Asia/Shanghai"),
    format(run_finished, tz = "Asia/Shanghai"), unname(versions)
  ),
  stringsAsFactors = FALSE
)
write.csv(parameters, file.path(out_dir, "doublet_analysis_parameters.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "doublet_sessionInfo.txt"))

message("Doublet assessment completed: ", out_dir)
