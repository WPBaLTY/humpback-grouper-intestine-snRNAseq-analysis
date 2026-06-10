# Regenerate submission-oriented figures for the humpback grouper intestine
# snRNA-seq atlas.
#
# Inputs:
#   - Main final global Seurat object.
#   - Immune reclustering Seurat object.
#   - cell_metadata.csv from the companion data repository.
#
# Outputs:
#   - Figure2_global_annotation_final.png
#   - Figure3_immune_reclustering_final.png
#   - Figure4_composition_final.png
#   - figure_source_data/figure4_major_celltype_composition.csv
#   - figure_source_data/figure4_immune_subtype_composition.csv
#
# Figure 1 is generated from sequencing/alignment metric tables and Figure 5 is
# retained as a final annotation-support heatmap in the companion data package.

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

# Fail early with clear messages when users run the script outside the
# expected folder layout.
require_file <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " is missing: ", path)
  }
  invisible(path)
}

# Resolve paths portably. Set GROUPER_DATA_DIR if the processed-data repository
# is not located next to this code repository.
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
  Sys.getenv("GROUPER_FIGURE_OUT_DIR", file.path(repo_dir, "analysis_figures")),
  winslash = "/",
  mustWork = FALSE
)
source_dir <- data_dir
metadata_path <- Sys.getenv("GROUPER_CELL_METADATA_CSV", file.path(source_dir, "cell_metadata.csv"))
figure_source_out <- file.path(out_dir, "figure_source_data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_source_out, recursive = TRUE, showWarnings = FALSE)
require_file(main_rds, "Main Seurat RDS")
require_file(immune_rds, "Immune reclustering Seurat RDS")
require_file(metadata_path, "Cell metadata CSV")

# Convert the historical combined-treatment prefix CON to the public-facing COM.
sample_label <- function(x) {
  gsub("^CON", "COM", x)
}

# Save figures at manuscript-friendly resolution.
save_png <- function(plot, filename, width, height) {
  ggsave(
    filename = file.path(out_dir, filename),
    plot = plot,
    width = width,
    height = height,
    dpi = 300,
    bg = "white",
    limitsize = FALSE
  )
}

# Keep marker lists robust to annotation/genome-model differences. The deposited
# gene model contains mixed naming conventions, so exact matching is tried first
# and then case-insensitive matching is used before a marker is omitted.
safe_features <- function(obj, features, context) {
  object_features <- rownames(obj)
  object_features_upper <- toupper(object_features)
  selected <- character(0)
  recovered_pairs <- character(0)
  missing <- character(0)

  for (feature in features) {
    if (feature %in% object_features) {
      selected <- c(selected, feature)
      next
    }
    idx <- match(toupper(feature), object_features_upper)
    if (!is.na(idx)) {
      selected <- c(selected, object_features[[idx]])
      recovered_pairs <- c(recovered_pairs, paste0(feature, "->", object_features[[idx]]))
    } else {
      missing <- c(missing, feature)
    }
  }

  selected <- unique(selected)
  if (length(recovered_pairs) > 0) {
    message(context, ": markers recovered by case-insensitive matching: ", paste(recovered_pairs, collapse = ", "))
  }
  if (length(missing) > 0) {
    message(context, ": marker genes omitted because they are absent from the deposited object: ", paste(missing, collapse = ", "))
  }
  if (length(selected) == 0) {
    stop(context, ": no requested marker genes were found in the deposited object.")
  }
  selected
}

# Consistent visual style for submission figures.
theme_submission <- theme_classic(base_family = "Arial", base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )

message("Reading main RDS: ", main_rds)
obj <- readRDS(main_rds)
DefaultAssay(obj) <- "RNA"

# The final manuscript uses 27 global transcriptional clusters labelled 0-26.
obj$cluster_res04 <- factor(as.character(obj$cluster_res04), levels = as.character(0:26))
obj$celltype_clean <- factor(as.character(obj$celltype_clean))

# Representative markers used to support broad intestinal cell-type annotation.
major_markers <- c(
  "fabp2", "cd36", "apoa1", "slc10a2", "fabp6", "muc2", "spdef",
  "best4", "otop2", "pou2f3", "avil", "neurod1", "scgn",
  "pecam1", "cdh5", "col1a1", "col1a2", "tagln", "myh11",
  "ptprc", "lcp1", "xcr1", "fcer1g", "cel", "cela1.1", "cpa2",
  "syt1", "elavl3"
)
major_markers <- safe_features(obj, major_markers, "Major cell-type dot plot")

# Figure 2a: global cluster UMAP.
p2a <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "cluster_res04",
  label = TRUE,
  repel = TRUE,
  raster = TRUE,
  pt.size = 0.05
) +
  labs(title = "a  UMAP by transcriptional cluster", color = "Cluster") +
  theme_submission

# Figure 2b: broad cell-type annotation UMAP.
p2b <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "celltype_clean",
  label = TRUE,
  repel = TRUE,
  raster = TRUE,
  pt.size = 0.05
) +
  labs(title = "b  UMAP by annotated cell type", color = "Cell type") +
  theme_submission

# Figure 2c: marker-gene dot plot for major annotated populations.
p2c <- DotPlot(obj, features = major_markers, group.by = "celltype_clean", dot.scale = 5) +
  labs(title = "c  Representative marker-gene dot plot", x = "Marker genes", y = "Cell type") +
  theme_submission +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
    axis.text.y = element_text(size = 8),
    legend.position = "right"
  )

fig2 <- (p2a + p2b) / p2c + plot_layout(heights = c(1.0, 0.85))
save_png(fig2, "Figure2_global_annotation_final.png", width = 11.8, height = 10.2)

message("Reading immune RDS: ", immune_rds)
immune <- readRDS(immune_rds)
DefaultAssay(immune) <- "RNA"

# Accept either available immune annotation column.
immune_group <- if ("immune_anno_final" %in% colnames(immune@meta.data)) "immune_anno_final" else "immune_anno"

# Representative markers for immune subtype support.
immune_markers <- c(
  "CCR7", "BCL11B", "LCK", "roraa", "FKBP5", "ddit4.1", "CCL20", "CCR6",
  "Prf1", "Gzmb.1", "TYROBP", "FCER1G", "XCR1", "ZNF366", "csf1r1", "Axl",
  "CMKLR1", "Cd209d.1", "FN1", "AOC3", "EPX", "Ncf4", "CYBB.1", "CD79A",
  "EBF1", "BLNK"
)
immune_markers <- safe_features(immune, immune_markers, "Immune subtype dot plot")

# Figure 3a: immune-focused UMAP.
p3a <- DimPlot(
  immune,
  reduction = "umap",
  group.by = immune_group,
  label = TRUE,
  repel = TRUE,
  raster = TRUE,
  pt.size = 0.25
) +
  labs(title = "a  Immune-cell subclustering", color = "Immune subtype") +
  theme_submission

# Figure 3b: immune-subtype marker-gene dot plot.
p3b <- DotPlot(immune, features = immune_markers, group.by = immune_group, dot.scale = 5) +
  labs(title = "b  Immune-subtype marker-gene dot plot", x = "Marker genes", y = "Immune subtype") +
  theme_submission +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
    axis.text.y = element_text(size = 8),
    legend.position = "right"
  )

fig3 <- p3a / p3b + plot_layout(heights = c(1.0, 0.95))
save_png(fig3, "Figure3_immune_reclustering_final.png", width = 11.5, height = 9.3)

cell_md <- read.csv(metadata_path, check.names = FALSE)
required_metadata_cols <- c("library_id", "diet_group", "celltype_clean", "immune_subtype")
missing_metadata_cols <- setdiff(required_metadata_cols, colnames(cell_md))
if (length(missing_metadata_cols) > 0) {
  stop("Cell metadata CSV is missing required columns: ", paste(missing_metadata_cols, collapse = ", "))
}
cell_md$library_id <- factor(cell_md$library_id, levels = c("CTL_1", "CTL_2", "Llac_1", "Llac_2", "Slim_1", "Slim_2", "COM_1", "COM_2"))

# Descriptive library-level composition for major intestinal lineages.
major_counts <- as.data.frame(table(cell_md$library_id, cell_md$diet_group, cell_md$celltype_clean))
colnames(major_counts) <- c("library_id", "diet_group", "cell_type", "n")
major_counts <- major_counts[major_counts$n > 0, ]
major_counts$prop <- ave(major_counts$n, major_counts$library_id, FUN = function(x) x / sum(x))
write.csv(major_counts, file.path(figure_source_out, "figure4_major_celltype_composition.csv"), row.names = FALSE)

# Descriptive library-level composition for immune subtypes only.
immune_md <- cell_md[cell_md$immune_subtype != "" & !is.na(cell_md$immune_subtype), ]
immune_counts <- as.data.frame(table(immune_md$library_id, immune_md$diet_group, immune_md$immune_subtype))
colnames(immune_counts) <- c("library_id", "diet_group", "immune_subtype", "n")
immune_counts <- immune_counts[immune_counts$n > 0, ]
immune_counts$prop <- ave(immune_counts$n, immune_counts$library_id, FUN = function(x) x / sum(x))
write.csv(immune_counts, file.path(figure_source_out, "figure4_immune_subtype_composition.csv"), row.names = FALSE)

p4a <- ggplot(major_counts, aes(x = library_id, y = prop, fill = cell_type)) +
  geom_col(width = 0.82) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(title = "a  Major intestinal cell types", x = "Library", y = "Relative abundance", fill = "Cell type") +
  theme_submission +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")

p4b <- ggplot(immune_counts, aes(x = library_id, y = prop, fill = immune_subtype)) +
  geom_col(width = 0.82) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(title = "b  Immune subtypes", x = "Library", y = "Relative abundance", fill = "Immune subtype") +
  theme_submission +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")

fig4 <- p4a + p4b + plot_layout(widths = c(1, 1))
save_png(fig4, "Figure4_composition_final.png", width = 12.0, height = 6.8)

message("Wrote final figures to ", out_dir)
