# One-command analysis workflow driver for the humpback grouper intestine
# snRNA-seq Data Descriptor.
#
# Run from a terminal:
#   Rscript 00_run_analysis_workflow.R
#
# Optional environment variables:
#   GROUPER_DATA_DIR       Path to the companion data repository.
#   GROUPER_ANALYSIS_OUT_DIR  Directory for all generated outputs from this run.
#
# This driver does not rerun Cell Ranger. It verifies the deposited Seurat
# objects, exports final metadata, regenerates Figures 2-4, and writes a
# compact run summary.

options(stringsAsFactors = FALSE)

command_args <- commandArgs(trailingOnly = FALSE)
file_arg <- command_args[grepl("^--file=", command_args)]
repo_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
setwd(repo_dir)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
run_root <- normalizePath(
  Sys.getenv(
    "GROUPER_ANALYSIS_OUT_DIR",
    file.path(repo_dir, paste0("analysis_run_", timestamp))
  ),
  winslash = "/",
  mustWork = FALSE
)
dir.create(run_root, recursive = TRUE, showWarnings = FALSE)

output_log <- file(file.path(run_root, "run_output.log"), open = "wt")
message_log <- file(file.path(run_root, "run_messages.log"), open = "wt")
sink(output_log, split = TRUE)
sink(message_log, type = "message")
on.exit({
  while (sink.number(type = "message") > 0) sink(type = "message")
  while (sink.number(type = "output") > 0) sink(type = "output")
  close(output_log)
  close(message_log)
}, add = TRUE)

cat("Humpback grouper intestine snRNA-seq analysis workflow\n")
cat("Started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("Code repository: ", repo_dir, "\n", sep = "")
cat("Run output directory: ", run_root, "\n\n", sep = "")

required_packages <- c("Seurat", "ggplot2", "patchwork", "scales", "harmony")
package_available <- vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
if (!all(package_available)) {
  missing <- names(package_available)[!package_available]
  stop(
    "Missing required R packages: ", paste(missing, collapse = ", "),
    "\nInstall them with: Rscript 00_setup/00_install_required_packages.R"
  )
}

cat("R version: ", R.version.string, "\n", sep = "")
cat("Package versions:\n")
print(vapply(required_packages, function(pkg) as.character(utils::packageVersion(pkg)), character(1)))
cat("\n")

data_dir <- normalizePath(
  Sys.getenv("GROUPER_DATA_DIR", file.path(repo_dir, "..", "humpback_grouper_intestine_snRNAseq_processed_data")),
  winslash = "/",
  mustWork = FALSE
)
main_rds <- file.path(data_dir, "integrated_object", "grouper_intestine_snRNAseq_seurat.rds")
immune_rds <- file.path(data_dir, "integrated_object", "grouper_intestine_immune_reclustering_seurat.rds")

required_files <- c(main_rds, immune_rds)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(
    "Required input files are missing:\n",
    paste(paste0("  - ", missing_files), collapse = "\n"),
    "\nSet GROUPER_DATA_DIR to the unpacked companion data repository if needed."
  )
}

metadata_out <- file.path(run_root, "metadata_27cluster")
figure_out <- file.path(run_root, "figures")
inventory_out <- file.path(run_root, "rds_inventory")

Sys.setenv(
  GROUPER_DATA_DIR = data_dir,
  GROUPER_MAIN_RDS = main_rds,
  GROUPER_IMMUNE_RDS = immune_rds,
  GROUPER_OUTPUT_DIR = metadata_out,
  GROUPER_CELL_METADATA_CSV = file.path(metadata_out, "cell_metadata.csv"),
  GROUPER_FIGURE_OUT_DIR = figure_out,
  GROUPER_INVENTORY_OUT_DIR = inventory_out,
  GROUPER_SESSION_INFO_FILE = file.path(run_root, "R_sessionInfo.txt")
)

run_step <- function(label, script) {
  cat("\n--- ", label, " ---\n", sep = "")
  start_time <- Sys.time()
  source(script, local = new.env(parent = globalenv()))
  elapsed <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1)
  cat("Completed ", label, " in ", elapsed, " seconds.\n", sep = "")
}

run_step("Write R session information", "04_environment/01_write_session_info.R")
run_step("Inspect deposited Seurat objects", "01_object_inventory/01_inspect_deposited_rds_objects.R")
run_step("Export final metadata tables", "02_export_metadata/01_export_submission_metadata.R")
run_step("Regenerate Figures 2-4", "03_make_figures/01_make_submission_figures.R")

count_csv_rows <- function(path) {
  if (!file.exists(path)) stop("Expected output is missing: ", path)
  length(readLines(path, warn = FALSE)) - 1L
}

summary_file <- file.path(metadata_out, "metadata_export_summary.txt")
summary_lines <- readLines(summary_file, warn = FALSE)
cell_metadata_file <- file.path(metadata_out, "cell_metadata.csv")
umap_file <- file.path(metadata_out, "umap_coordinates.csv")
cluster_file <- file.path(metadata_out, "cluster_annotation_counts.csv")

cell_rows <- count_csv_rows(cell_metadata_file)
umap_rows <- count_csv_rows(umap_file)
cluster_counts <- utils::read.csv(cluster_file, check.names = FALSE)
c14 <- cluster_counts[as.character(cluster_counts$global_cluster) == "14", , drop = FALSE]

expected <- list(
  n_nuclei = 102036L,
  n_umap_rows = 102036L,
  n_global_clusters = 27L,
  n_acinar_like_cluster14 = 1347L
)

if (cell_rows != expected$n_nuclei) stop("Unexpected cell_metadata row count: ", cell_rows)
if (umap_rows != expected$n_umap_rows) stop("Unexpected UMAP row count: ", umap_rows)
if (nrow(cluster_counts) != expected$n_global_clusters) {
  stop("Unexpected number of global clusters: ", nrow(cluster_counts))
}
if (nrow(c14) != 1L || as.integer(c14$n_nuclei) != expected$n_acinar_like_cluster14) {
  stop("Unexpected global cluster 14 count.")
}

expected_outputs <- c(
  file.path(inventory_out, "rds_object_summary.csv"),
  cell_metadata_file,
  umap_file,
  cluster_file,
  file.path(metadata_out, "celltype_counts_by_library.csv"),
  file.path(figure_out, "Figure2_global_annotation_final.png"),
  file.path(figure_out, "Figure3_immune_reclustering_final.png"),
  file.path(figure_out, "Figure4_composition_final.png"),
  file.path(figure_out, "figure_source_data", "figure4_major_celltype_composition.csv"),
  file.path(figure_out, "figure_source_data", "figure4_immune_subtype_composition.csv")
)
missing_outputs <- expected_outputs[!file.exists(expected_outputs) | file.info(expected_outputs)$size == 0]
if (length(missing_outputs) > 0) {
  stop("Expected outputs are missing or empty:\n", paste(paste0("  - ", missing_outputs), collapse = "\n"))
}

run_summary <- c(
  "Humpback grouper intestine snRNA-seq analysis workflow",
  paste0("completed=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("repo_dir=", repo_dir),
  paste0("data_dir=", data_dir),
  paste0("run_root=", run_root),
  "",
  "Metadata audit:",
  paste0("cell_metadata_rows=", cell_rows),
  paste0("umap_coordinate_rows=", umap_rows),
  paste0("global_cluster_rows=", nrow(cluster_counts)),
  paste0("global_cluster_14_celltype=", c14$celltype_final),
  paste0("global_cluster_14_n_nuclei=", c14$n_nuclei),
  "",
  "metadata_export_summary.txt:",
  summary_lines,
  "",
  "Output files:",
  paste0("  - ", expected_outputs)
)
writeLines(run_summary, file.path(run_root, "RUN_SUMMARY.txt"))

cat("\nAll workflow validation checks passed.\n")
cat("Run summary: ", file.path(run_root, "RUN_SUMMARY.txt"), "\n", sep = "")
