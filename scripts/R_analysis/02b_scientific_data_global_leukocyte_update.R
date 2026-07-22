# Scientific Data 2026-07-09 global/leukocyte annotation update
#
# This script records the final Scientific Data annotation logic:
# - Figure 3 and Table 2 use broad global cell types.
# - Leukocyte-lineage nuclei are collapsed to "Leukocytes" at the global level.
# - Detailed immune states are introduced only after leukocyte-focused reclustering
#   in Figure 4 and Table 3.
#
# The script intentionally uses base R only so the table/source-data update can be
# checked without loading large Seurat objects. Upstream reconstruction is
# documented in 00_reconstruct_global_atlas_from_10x.R and
# 01_reconstruct_immune_final_object.R.

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", script_args[grepl(file_arg, script_args)])
script_dir <- if (length(script_path) > 0) dirname(normalizePath(script_path[1], mustWork = FALSE)) else getwd()
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)

path_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

write_source_csv <- function(x, filename, out_dir) {
  ensure_dir(out_dir)
  write.csv(x, file.path(out_dir, filename), row.names = FALSE, quote = TRUE, na = "")
}

out_dir <- path_env("GROUPER_SOURCE_DATA_OUT", file.path(repo_root, "source_data", "figures"))
audit_dir <- path_env("GROUPER_ANNOTATION_AUDIT_OUT", file.path(repo_root, "audit", "annotation"))
ensure_dir(out_dir)
ensure_dir(audit_dir)

table2_global_annotations <- data.frame(
  `Cell type` = c(
    "Enterocytes",
    "Lysosome-rich enterocyte-related cells (LREs)",
    "Goblet cells",
    "Putative tuft-like epithelial cells",
    "Best4+ epithelial cells",
    "Enteroendocrine cells",
    "Neuronal cells",
    "Leukocytes",
    "Fibroblasts",
    "Endothelial cells",
    "Smooth muscle cells",
    "Pancreatic acinar-like cells"
  ),
  `Marker genes` = c(
    "fabp2, cd36, apoa1, SI, CHIA, MGAM",
    "slc10a2, fabp6, CUBN, LRP2",
    "muc2, spdef, FER1L6, ST3GAL1",
    "pou2f3, avil, Pik3ap1, Adgrg2.1",
    "best4, cftr, slc20a1a",
    "neurod1, scgn, ISL1, SPATA17",
    "syt1, elavl3, phox2a, eno2",
    "ptprc, lcp1, BCL11B, SATB1, FCER1G, XCR1, Axl",
    "col1a1, col1a2, dcn, PLTP",
    "pecam1, cdh5, kdrl, flt4",
    "tagln, CNTNAP5, Pld5.1",
    "cel, cela1.1, cpa2, ctrb1"
  ),
  `Cluster(s)` = c(
    "C0, C1",
    "C2, C6, C9",
    "C4, C23",
    "C7",
    "C5, C24, C26",
    "C19",
    "C17, C18, C20, C25",
    "C3, C8, C12, C15, C22",
    "C13",
    "C10, C21",
    "C11, C16",
    "C14"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

table3_leukocyte_annotations <- data.frame(
  `Cell type` = c(
    "T cell (CCR7+)",
    "Activated T (RORA+)",
    "Activated lymphoid (CCL20hi)",
    "NK-like cytotoxic",
    "B cells",
    "cDC1-like (XCR1+)",
    "Monocytes/macrophages",
    "MoDC-like (CD209d+)",
    "Granulocyte-like",
    "Cycling (G2/M)"
  ),
  `Marker genes` = c(
    "CCR7, BCL11B, LCK",
    "roraa, FKBP5, ddit4.1",
    "Ccl20, CCR6, RORA",
    "Prf1, Gzmb.1, TYROBP",
    "CD79A, EBF1, Pax5",
    "XCR1, ZNF366, DIP2C",
    "csf1r1, Axl, CMKLR1",
    "Cd209d.1, FN1, AOC3",
    "EPX, Ncf4, CYBB.1",
    "ncaph2, ANLN, birc5.2, TOP2A"
  ),
  `Analysis level` = rep("Leukocyte reclustering", 10),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

global_marker_sets <- list(
  "Enterocytes" = c("fabp2", "cd36", "SI"),
  "LREs" = c("slc10a2", "CUBN", "LRP2"),
  "Goblet cells" = c("muc2", "spdef", "FER1L6"),
  "Tuft-like cells" = c("pou2f3", "avil", "Pik3ap1"),
  "Best4+ cells" = c("best4", "cftr", "slc20a1a"),
  "Enteroendocrine cells" = c("neurod1", "scgn", "ISL1"),
  "Neuronal cells" = c("syt1", "elavl3", "phox2a"),
  "Leukocytes" = c("ptprc", "lcp1", "BCL11B", "SATB1", "FCER1G", "XCR1", "Axl"),
  "Fibroblasts" = c("col1a1", "col1a2", "dcn"),
  "Endothelial cells" = c("pecam1", "cdh5", "kdrl"),
  "Smooth muscle cells" = c("tagln", "CNTNAP5", "Pld5.1"),
  "Acinar-like cells" = c("cel", "Cela1.1", "cpa2")
)

global_cluster_map <- c(
  "Enterocytes" = "C0, C1",
  "LREs" = "C2, C6, C9",
  "Goblet cells" = "C4, C23",
  "Tuft-like cells" = "C7",
  "Best4+ cells" = "C5, C24, C26",
  "Enteroendocrine cells" = "C19",
  "Neuronal cells" = "C17, C18, C20, C25",
  "Leukocytes" = "C3, C8, C12, C15, C22",
  "Fibroblasts" = "C13",
  "Endothelial cells" = "C10, C21",
  "Smooth muscle cells" = "C11, C16",
  "Acinar-like cells" = "C14"
)

gene_panel <- unique(unlist(global_marker_sets, use.names = FALSE))
fig3_marker_panel <- do.call(
  rbind,
  lapply(names(global_marker_sets), function(cell_type) {
    data.frame(
      figure = "Figure3c",
      cell_type = cell_type,
      marker_gene = gene_panel,
      assigned_clusters = unname(global_cluster_map[cell_type]),
      selected_support_marker = ifelse(gene_panel %in% global_marker_sets[[cell_type]], "True", "False"),
      stringsAsFactors = FALSE
    )
  })
)

write_source_csv(table2_global_annotations, "Table2_global_cell_type_annotations.csv", out_dir)
write_source_csv(table3_leukocyte_annotations, "Table3_leukocyte_reclustering_annotations.csv", out_dir)
write_source_csv(fig3_marker_panel, "Fig3_global_marker_dotplot_gene_panel.csv", out_dir)

message("Wrote Scientific Data 2026-07-09 annotation update files to: ", out_dir)
message("Figure 5 composition CSVs are deposited as checked source data in the same folder:")
message("  - Fig5a_global_celltype_composition_leukocytes_merged.csv")
message("  - Fig5b_immune_subtype_composition.csv")
