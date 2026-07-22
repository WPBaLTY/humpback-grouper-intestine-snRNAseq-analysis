options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})

required_core_versions <- c(Seurat = "5.3.0", SeuratObject = "5.2.0")
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
    ". Run environment/install_exact_r_packages.R and set R_LIBS_USER to the resulting .r-library."
  )
}

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", script_args, value = TRUE)[1])
script_dir <- if (length(script_file) == 0 || is.na(script_file)) {
  getwd()
} else {
  dirname(normalizePath(script_file, winslash = "/", mustWork = FALSE))
}
repo_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = FALSE)

path_env <- function(name, default) {
  normalizePath(Sys.getenv(name, unset = default), winslash = "/", mustWork = FALSE)
}

immune_rds <- path_env(
  "GROUPER_IMMUNE_FINAL_RDS",
  file.path(repo_root, "outputs", "immune_final_reconstruction", "IMMUNE_FINAL_annotated_with_Bcell.rds")
)
source_out <- path_env(
  "GROUPER_SOURCE_DATA_OUT",
  file.path(repo_root, "source_data", "figures")
)

dir.create(source_out, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(immune_rds))

safe_join <- function(obj, assay = "RNA") {
  if (!(assay %in% Assays(obj))) return(obj)
  DefaultAssay(obj) <- assay
  layers <- tryCatch(Layers(obj[[assay]]), error = function(e) character(0))
  if (length(layers) > 1) {
    obj <- JoinLayers(obj, assay = assay)
  }
  obj
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
    hit <- grep(
      paste0("^", gene_base, "(\\.|$)"),
      features_all,
      value = TRUE,
      ignore.case = TRUE,
      perl = TRUE
    )
    if (length(hit) > 0) {
      out <- c(out, hit[which.min(nchar(hit))])
    }
  }
  unique(out)
}

display_labels <- c(
  "T_cell_CCR7_like" = "T cell (CCR7+)",
  "Activated_T_RORA_stress" = "Activated T (RORA+)",
  "Activated_lymphoid_CCL20_high" = "Activated lymphoid (CCL20hi)",
  "NK_like_cytotoxic" = "NK-like cytotoxic",
  "B_cell" = "B cells",
  "cDC1_XCR1_ZNF366" = "cDC1-like (XCR1+)",
  "Macrophage_Monocyte_Axl_Csf1r1" = "Monocytes/macrophages",
  "MoDC_like_Cd209d_FN1_AOC3" = "MoDC-like (CD209d+)",
  "Granulocyte_like_EPX_Ncf4_CYBB" = "Granulocyte-like",
  "Cycling_immune_G2M" = "Cycling (G2/M)"
)

subtype_order <- unname(display_labels)

marker_genes <- c(
  "CCR7", "BCL11B", "LCK",
  "roraa", "FKBP5", "ddit4.1",
  "Ccl20", "CCR6",
  "Prf1", "Gzmb.1", "TYROBP", "FCER1G",
  "XCR1", "ZNF366",
  "csf1r1", "Axl", "CMKLR1",
  "Cd209d.1", "FN1", "AOC3",
  "EPX", "Ncf4", "CYBB.1",
  "CD79A", "EBF1", "BLNK"
)

obj <- readRDS(immune_rds)
obj <- safe_join(obj, "RNA")
DefaultAssay(obj) <- "RNA"

stopifnot("immune_anno_final" %in% colnames(obj@meta.data))
stopifnot("umap" %in% Reductions(obj))

anno_internal <- as.character(obj[["immune_anno_final"]][, 1])
anno_display <- unname(display_labels[anno_internal])
if (any(is.na(anno_display))) {
  stop(
    "Unmapped immune annotations: ",
    paste(sort(unique(anno_internal[is.na(anno_display)])), collapse = ", ")
  )
}

emb <- as.data.frame(Embeddings(obj, "umap"))
emb$cell_id <- rownames(emb)
emb$cell_id <- sub("^CON", "COM", emb$cell_id)
colnames(emb)[1:2] <- c("UMAP_1", "UMAP_2")

meta <- obj@meta.data
source_umap <- data.frame(
  cell_id = emb$cell_id,
  UMAP_1 = emb$UMAP_1,
  UMAP_2 = emb$UMAP_2,
  group = ifelse(as.character(meta$group) == "CON", "COM", as.character(meta$group)),
  sample = sub("^CON", "COM", as.character(meta$sample)),
  immune_subtype_internal = anno_internal,
  immune_subtype = anno_display,
  stringsAsFactors = FALSE
)

source_umap$immune_subtype <- factor(source_umap$immune_subtype, levels = subtype_order)
source_umap <- source_umap[order(source_umap$immune_subtype, source_umap$cell_id), ]
write.csv(
  source_umap,
  file.path(source_out, "Fig4_immune_umap_coordinates.csv"),
  row.names = FALSE,
  quote = TRUE
)

features <- resolve_features(obj, marker_genes)
missing_features <- setdiff(tolower(marker_genes), tolower(features))
if (length(missing_features) > 0) {
  warning("Marker genes not found in RNA assay: ", paste(missing_features, collapse = ", "))
}

expr <- GetAssayData(obj, assay = "RNA", layer = "data")
groups <- factor(anno_display, levels = subtype_order)
dot_rows <- list()

for (feature in features) {
  values <- as.numeric(expr[feature, ])
  avg_by_group <- tapply(values, groups, mean)
  pct_by_group <- tapply(values > 0, groups, mean) * 100
  avg_by_group <- avg_by_group[subtype_order]
  pct_by_group <- pct_by_group[subtype_order]

  scaled <- as.numeric(scale(avg_by_group))
  scaled[is.na(scaled)] <- 0
  scaled <- pmax(pmin(scaled, 2.5), -1.5)

  dot_rows[[feature]] <- data.frame(
    immune_subtype = subtype_order,
    marker_gene = feature,
    percent_expressed = as.numeric(pct_by_group),
    average_expression = as.numeric(avg_by_group),
    average_expression_scaled = scaled,
    stringsAsFactors = FALSE
  )
}

source_dot <- do.call(rbind, dot_rows)
source_dot$immune_subtype <- factor(source_dot$immune_subtype, levels = subtype_order)
source_dot$marker_gene <- factor(source_dot$marker_gene, levels = features)
source_dot <- source_dot[order(source_dot$immune_subtype, source_dot$marker_gene), ]
write.csv(
  source_dot,
  file.path(source_out, "Fig4_immune_marker_dotplot_source.csv"),
  row.names = FALSE,
  quote = TRUE
)

write.csv(
  data.frame(feature = features, stringsAsFactors = FALSE),
  file.path(source_out, "Fig4_immune_marker_features_used.csv"),
  row.names = FALSE,
  quote = TRUE
)

sample_plot <- gsub("^CON", "COM", as.character(meta$sample))
sample_plot <- gsub("([A-Za-z]+)([0-9]+)$", "\\1_\\2", sample_plot)
immune_comp <- aggregate(
  x = data.frame(n = rep(1L, length(anno_display))),
  by = list(
    group = ifelse(as.character(meta$group) == "CON", "COM", as.character(meta$group)),
    sample = sub("^CON", "COM", as.character(meta$sample)),
    subtype = anno_display,
    sample_plot = sample_plot
  ),
  FUN = sum
)
immune_comp$prop <- ave(immune_comp$n, immune_comp$sample, FUN = function(x) x / sum(x))
immune_comp <- immune_comp[, c("group", "sample", "subtype", "n", "prop", "sample_plot")]
write.csv(
  immune_comp,
  file.path(source_out, "Fig5b_immune_subtype_composition.csv"),
  row.names = FALSE,
  quote = TRUE
)

message("Wrote Scientific Data Figure 4 source data to: ", source_out)
