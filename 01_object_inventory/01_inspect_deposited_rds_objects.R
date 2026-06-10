# Inspect the deposited Seurat RDS objects and write lightweight metadata
# summaries for the analysis record.

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(Seurat)
})

repo_dir <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
data_dir <- normalizePath(
  Sys.getenv("GROUPER_DATA_DIR", file.path(repo_dir, "..", "humpback_grouper_intestine_snRNAseq_processed_data")),
  winslash = "/",
  mustWork = FALSE
)
input_dir <- file.path(data_dir, "integrated_object")
out_dir <- normalizePath(
  Sys.getenv("GROUPER_INVENTORY_OUT_DIR", file.path(repo_dir, "rds_inventory")),
  winslash = "/",
  mustWork = FALSE
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rds_files <- c(
  "grouper_intestine_snRNAseq_seurat.rds",
  "grouper_intestine_immune_reclustering_seurat.rds"
)

summaries <- list()

for (file in rds_files) {
  path <- file.path(input_dir, file)
  if (!file.exists(path)) {
    stop("Expected RDS file is missing: ", path)
  }
  message("Reading ", path)
  obj <- readRDS(path)
  md <- obj@meta.data
  reductions <- names(obj@reductions)
  assays <- names(obj@assays)
  dims <- paste(nrow(obj), ncol(obj), sep = " x ")

  metadata_summary_cols <- grep(
    "cluster|res|celltype|annot|immune|sample|group",
    colnames(md),
    ignore.case = TRUE,
    value = TRUE
  )

  col_summaries <- lapply(metadata_summary_cols, function(col) {
    vals <- md[[col]]
    uniq <- unique(as.character(vals))
    data.frame(
      file = file,
      column = col,
      n_unique = length(uniq),
      values_preview = paste(head(sort(uniq), 80), collapse = " | "),
      stringsAsFactors = FALSE
    )
  })
  col_summaries <- if (length(col_summaries)) do.call(rbind, col_summaries) else data.frame()
  write.csv(col_summaries, file.path(out_dir, paste0(file, ".metadata_columns.csv")), row.names = FALSE)

  summaries[[file]] <- data.frame(
    file = file,
    class = paste(class(obj), collapse = ";"),
    dimensions = dims,
    assays = paste(assays, collapse = ";"),
    reductions = paste(reductions, collapse = ";"),
    metadata_columns = ncol(md),
    cells = nrow(md),
    stringsAsFactors = FALSE
  )

  write.csv(head(md, 20), file.path(out_dir, paste0(file, ".metadata_head.csv")), row.names = TRUE)
}

summary_df <- do.call(rbind, summaries)
write.csv(summary_df, file.path(out_dir, "rds_object_summary.csv"), row.names = FALSE)
print(summary_df)
