# Regenerate submission-oriented figures for the humpback grouper intestine
# snRNA-seq atlas.
#
# Inputs:
#   - Main final global Seurat object.
#   - Immune reclustering Seurat object.
#   - cell_metadata.csv from the companion data repository.
#
# Outputs:
#   - Figure3_global_annotation_final.png/.tiff/.pdf
#   - Figure4_immune_reclustering_final.png/.tiff/.pdf
#   - Figure5_composition_final.png/.tiff/.pdf
#   - figure_source_data/figure5_major_celltype_composition.csv
#   - figure_source_data/figure5_immune_subtype_composition.csv
#
# Figure 1 is the experimental workflow figure and Figure 2 is generated from
# sequencing/alignment metric tables. This script regenerates Figures 3-5.

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

# Save standalone, production-quality figure files. PNG files are retained for
# manuscript preview, TIFF files are suitable for journal upload, and PDF files
# provide a vector-friendly version where the plotting device supports it.
save_publication_figure <- function(plot, filename_base, width, height) {
  ggsave(
    filename = file.path(out_dir, paste0(filename_base, ".png")),
    plot = plot,
    width = width,
    height = height,
    dpi = 600,
    bg = "white",
    limitsize = FALSE
  )
  ggsave(
    filename = file.path(out_dir, paste0(filename_base, ".tiff")),
    plot = plot,
    width = width,
    height = height,
    dpi = 600,
    bg = "white",
    compression = "lzw",
    limitsize = FALSE
  )
  ggsave(
    filename = file.path(out_dir, paste0(filename_base, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    bg = "white",
    device = cairo_pdf,
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
  labels <- character(0)
  recovered_pairs <- character(0)
  missing <- character(0)

  for (feature in features) {
    if (feature %in% object_features) {
      selected <- c(selected, feature)
      labels <- c(labels, feature)
      next
    }
    idx <- match(toupper(feature), object_features_upper)
    if (!is.na(idx)) {
      selected <- c(selected, object_features[[idx]])
      labels <- c(labels, feature)
      recovered_pairs <- c(recovered_pairs, paste0(feature, "->", object_features[[idx]]))
    } else {
      missing <- c(missing, feature)
    }
  }

  keep <- !duplicated(selected)
  selected <- selected[keep]
  labels <- labels[keep]
  if (length(recovered_pairs) > 0) {
    message(context, ": markers recovered by case-insensitive matching: ", paste(recovered_pairs, collapse = ", "))
  }
  if (length(missing) > 0) {
    message(context, ": marker genes omitted because they are absent from the deposited object: ", paste(missing, collapse = ", "))
  }
  if (length(selected) == 0) {
    stop(context, ": no requested marker genes were found in the deposited object.")
  }
  names(labels) <- selected
  labels
}

# Consistent visual style for submission figures.
theme_submission <- theme_classic(base_family = "Arial", base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11, hjust = 0, margin = margin(l = 22, b = 6)),
    plot.title.position = "plot",
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    plot.margin = margin(6, 6, 6, 6)
  )

theme_panel_tags <- theme(
  plot.tag = element_text(face = "bold", size = 15),
  plot.tag.position = c(0, 1)
)

label_values <- function(values, mapping) {
  unname(ifelse(values %in% names(mapping), mapping[values], values))
}

celltype_label_map <- c(
  "Acinar_like" = "Acinar-like",
  "Best4_cell" = "Best4+ epithelial cells",
  "Dc_cell" = "Dendritic cell-like",
  "EEC" = "Enteroendocrine cells",
  "Endothelial_cell" = "Endothelial cells",
  "Enterocyte" = "Enterocytes",
  "Fibroblasts" = "Fibroblasts",
  "Goblet_cell" = "Goblet cells",
  "LRE" = "LREs",
  "Macrophages_monocytes" = "Macrophages/monocytes",
  "Nerve_cell" = "Neuronal cells",
  "NK_cell" = "NK-like cells",
  "Smooth_muscle" = "Smooth muscle",
  "T_cell" = "T cells",
  "Tuft_cell" = "Putative tuft-like epithelial cells"
)

immune_label_map <- c(
  "T_cell_CCR7_like" = "T cell (CCR7+)",
  "Activated_T_RORA_stress" = "Activated T (RORA+)",
  "Activated_lymphoid_CCL20_high" = "Activated lymphoid (CCL20hi)",
  "NK_like_cytotoxic" = "NK-like cytotoxic",
  "B_cell" = "B cells",
  "cDC1_XCR1_ZNF366" = "cDC1-like cells",
  "Macrophage_Monocyte_Axl_Csf1r1" = "Monocytes/macrophages",
  "MoDC_like_Cd209d_FN1_AOC3" = "MoDC-like (CD209d+)",
  "Granulocyte_like_EPX_Ncf4_CYBB" = "Granulocyte-like",
  "Cycling_immune_G2M" = "Cycling (G2/M)"
)

celltype_order <- c(
  "Enterocytes", "LREs", "Goblet cells", "Putative tuft-like epithelial cells", "Best4+ epithelial cells",
  "Enteroendocrine cells", "Dendritic cell-like", "Macrophages/monocytes",
  "T cells", "NK-like cells", "Fibroblasts", "Endothelial cells",
  "Smooth muscle", "Neuronal cells", "Acinar-like"
)
celltype_palette <- c(
  "Enterocytes" = "#1F77B4",
  "LREs" = "#4E9CD3",
  "Goblet cells" = "#FF7F0E",
  "Putative tuft-like epithelial cells" = "#E28E2C",
  "Best4+ epithelial cells" = "#33A02C",
  "Enteroendocrine cells" = "#6FAE45",
  "Dendritic cell-like" = "#7A6A00",
  "Macrophages/monocytes" = "#B15928",
  "T cells" = "#D62728",
  "NK-like cells" = "#6A3D9A",
  "Fibroblasts" = "#D65F5F",
  "Endothelial cells" = "#9467BD",
  "Smooth muscle" = "#9B78BD",
  "Neuronal cells" = "#8C564B",
  "Acinar-like" = "#D45A9D"
)

immune_order <- c(
  "T cell (CCR7+)", "Activated T (RORA+)", "Activated lymphoid (CCL20hi)",
  "NK-like cytotoxic", "B cells", "cDC1-like cells",
  "Monocytes/macrophages", "MoDC-like (CD209d+)", "Granulocyte-like",
  "Cycling (G2/M)"
)
immune_palette <- c(
  "T cell (CCR7+)" = "#4DB6AC",
  "Activated T (RORA+)" = "#D6B400",
  "Activated lymphoid (CCL20hi)" = "#8E7CC3",
  "NK-like cytotoxic" = "#E64B35",
  "B cells" = "#3C8DBC",
  "cDC1-like cells" = "#F28E2B",
  "Monocytes/macrophages" = "#79A832",
  "MoDC-like (CD209d+)" = "#D66BA0",
  "Granulocyte-like" = "#8C8C8C",
  "Cycling (G2/M)" = "#9C4F9F"
)

message("Reading main RDS: ", main_rds)
obj <- readRDS(main_rds)
DefaultAssay(obj) <- "RNA"

# The final manuscript uses 27 global transcriptional clusters labelled 0-26.
obj$cluster_res04 <- factor(as.character(obj$cluster_res04), levels = as.character(0:26))
obj$celltype_clean <- factor(as.character(obj$celltype_clean))
obj$celltype_plot <- factor(
  label_values(as.character(obj$celltype_clean), celltype_label_map),
  levels = celltype_order[celltype_order %in% label_values(levels(obj$celltype_clean), celltype_label_map)]
)

# Representative markers used to support broad intestinal cell-type annotation.
major_markers <- c(
  "fabp2", "cd36", "apoa1", "slc10a2", "fabp6", "muc2", "spdef",
  "best4", "otop2", "pou2f3", "avil", "neurod1", "scgn",
  "pecam1", "cdh5", "col1a1", "col1a2", "tagln", "myh11",
  "ptprc", "lcp1", "xcr1", "fcer1g", "cel", "cela1.1", "cpa2",
  "syt1", "elavl3"
)
major_marker_labels <- safe_features(obj, major_markers, "Major cell-type dot plot")
major_markers <- names(major_marker_labels)

# Figure 3a: global cluster UMAP.
p2a <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "cluster_res04",
  label = TRUE,
  repel = TRUE,
  raster = TRUE,
  raster.dpi = c(600, 600),
  label.size = 4,
  pt.size = 0.65,
  cols = scales::hue_pal(c = 130, l = 45)(27)
) +
  labs(title = "UMAP by transcriptional cluster", color = "Cluster", x = "UMAP 1", y = "UMAP 2") +
  theme_submission +
  theme(legend.position = "right")

# Figure 3b: broad cell-type annotation UMAP.
p2b <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "celltype_plot",
  label = TRUE,
  repel = TRUE,
  raster = TRUE,
  raster.dpi = c(600, 600),
  label.size = 3.3,
  pt.size = 0.65,
  cols = unname(celltype_palette[levels(obj$celltype_plot)])
) +
  labs(title = "UMAP by annotated cell type", color = "Cell type", x = "UMAP 1", y = "UMAP 2") +
  theme_submission +
  theme(legend.position = "right")

# Figure 3c: marker-gene dot plot for major annotated populations.
p2c <- DotPlot(obj, features = major_markers, group.by = "celltype_plot", dot.scale = 4.2) +
  labs(title = "Representative marker-gene dot plot", x = "Marker genes", y = "Cell type") +
  scale_x_discrete(labels = function(x) label_values(x, major_marker_labels)) +
  theme_submission +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7, face = "italic"),
    axis.text.y = element_text(size = 8),
    legend.position = "right"
  )

fig2 <- ((p2a + p2b) / p2c) +
  plot_layout(heights = c(1.10, 0.90)) +
  plot_annotation(tag_levels = "a", tag_suffix = "  ") &
  theme_panel_tags
save_publication_figure(fig2, "Figure3_global_annotation_final", width = 11.8, height = 10.2)

message("Reading immune RDS: ", immune_rds)
immune <- readRDS(immune_rds)
DefaultAssay(immune) <- "RNA"

# Accept either available immune annotation column.
immune_group <- if ("immune_anno_final" %in% colnames(immune@meta.data)) "immune_anno_final" else "immune_anno"
immune$immune_plot <- factor(
  label_values(as.character(immune@meta.data[[immune_group]]), immune_label_map),
  levels = immune_order[immune_order %in% label_values(unique(as.character(immune@meta.data[[immune_group]])), immune_label_map)]
)

# Representative markers for immune subtype support.
immune_markers <- c(
  "CCR7", "BCL11B", "LCK", "roraa", "FKBP5", "ddit4.1", "CCL20", "CCR6",
  "Prf1", "Gzmb.1", "TYROBP", "FCER1G", "XCR1", "ZNF366", "csf1r1", "Axl",
  "CMKLR1", "Cd209d.1", "FN1", "AOC3", "EPX", "Ncf4", "CYBB.1", "CD79A",
  "EBF1", "BLNK"
)
immune_marker_labels <- safe_features(immune, immune_markers, "Immune subtype dot plot")
immune_markers <- names(immune_marker_labels)

# Figure 4a: immune-focused UMAP.
p3a <- DimPlot(
  immune,
  reduction = "umap",
  group.by = "immune_plot",
  label = TRUE,
  repel = TRUE,
  raster = TRUE,
  raster.dpi = c(600, 600),
  label.size = 3.2,
  pt.size = 1.55,
  cols = unname(immune_palette[levels(immune$immune_plot)])
) +
  labs(title = "Immune-lineage UMAP", color = "Immune subtype", x = "UMAP 1", y = "UMAP 2") +
  theme_submission +
  theme(legend.position = "right")

# Figure 4b: immune-subtype marker-gene dot plot.
p3b <- DotPlot(immune, features = immune_markers, group.by = "immune_plot", dot.scale = 4.4) +
  labs(title = "Immune-subtype marker-gene dot plot", x = "Marker genes", y = "Immune subtype") +
  scale_x_discrete(labels = function(x) label_values(x, immune_marker_labels)) +
  theme_submission +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7, face = "italic"),
    axis.text.y = element_text(size = 8),
    legend.position = "right"
  )

fig3 <- (p3a / p3b) +
  plot_layout(heights = c(1.15, 1.00)) +
  plot_annotation(tag_levels = "a", tag_suffix = "  ") &
  theme_panel_tags
save_publication_figure(fig3, "Figure4_immune_reclustering_final", width = 11.5, height = 8.15)

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
major_counts$cell_type_plot <- factor(
  label_values(as.character(major_counts$cell_type), celltype_label_map),
  levels = celltype_order[celltype_order %in% label_values(sort(unique(as.character(major_counts$cell_type))), celltype_label_map)]
)
write.csv(major_counts, file.path(figure_source_out, "figure5_major_celltype_composition.csv"), row.names = FALSE)

# Descriptive library-level composition for immune subtypes only.
immune_md <- cell_md[cell_md$immune_subtype != "" & !is.na(cell_md$immune_subtype), ]
immune_counts <- as.data.frame(table(immune_md$library_id, immune_md$diet_group, immune_md$immune_subtype))
colnames(immune_counts) <- c("library_id", "diet_group", "immune_subtype", "n")
immune_counts <- immune_counts[immune_counts$n > 0, ]
immune_counts$prop <- ave(immune_counts$n, immune_counts$library_id, FUN = function(x) x / sum(x))
immune_counts$immune_subtype_plot <- factor(
  label_values(as.character(immune_counts$immune_subtype), immune_label_map),
  levels = immune_order[immune_order %in% label_values(sort(unique(as.character(immune_counts$immune_subtype))), immune_label_map)]
)
write.csv(immune_counts, file.path(figure_source_out, "figure5_immune_subtype_composition.csv"), row.names = FALSE)

p4a <- ggplot(major_counts, aes(x = library_id, y = prop, fill = cell_type_plot)) +
  geom_col(width = 0.82) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = celltype_palette[levels(major_counts$cell_type_plot)]) +
  labs(title = "Major intestinal lineage composition per library", x = "Library", y = "Relative abundance", fill = "Cell type") +
  theme_submission +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "right"
  )

p4b <- ggplot(immune_counts, aes(x = library_id, y = prop, fill = immune_subtype_plot)) +
  geom_col(width = 0.82) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = immune_palette[levels(immune_counts$immune_subtype_plot)]) +
  labs(title = "Immune-subtype composition per library", x = "Library", y = "Relative abundance", fill = "Immune subtype") +
  theme_submission +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "right"
  )

fig4 <- (p4a / p4b) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "a", tag_suffix = "  ") &
  theme_panel_tags
save_publication_figure(fig4, "Figure5_composition_final", width = 10.0, height = 10.0)

message("Wrote final figures to ", out_dir)
