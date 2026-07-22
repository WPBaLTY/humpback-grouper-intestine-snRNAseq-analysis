options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", script_args, value = TRUE)[1])
script_dir <- if (length(script_file) == 0 || is.na(script_file)) {
  getwd()
} else {
  dirname(normalizePath(script_file, winslash = "/", mustWork = FALSE))
}
repo_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
lib_dir <- if (length(args) >= 1 && nzchar(args[1])) {
  normalizePath(args[1], winslash = "/", mustWork = FALSE)
} else {
  file.path(repo_root, ".r-library-doublet")
}

dir.create(lib_dir, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib_dir, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"))

install_if_missing <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package, lib = lib_dir)
  }
}

install_if_missing("remotes")
install_if_missing("BiocManager")

cran_versions <- c(
  SeuratObject = "5.4.0",
  Seurat = "5.5.1",
  dplyr = "1.2.1"
)
for (package in names(cran_versions)) {
  installed <- tryCatch(
    as.character(utils::packageVersion(package, lib.loc = lib_dir)),
    error = function(e) NA_character_
  )
  if (!identical(installed, cran_versions[[package]])) {
    remotes::install_version(
      package,
      version = cran_versions[[package]],
      lib = lib_dir,
      dependencies = TRUE,
      upgrade = "never"
    )
  }
}

BiocManager::install(version = "3.22", ask = FALSE, update = FALSE)
BiocManager::install(
  c("SingleCellExperiment", "scDblFinder", "BiocParallel"),
  lib = lib_dir,
  ask = FALSE,
  update = FALSE
)

expected <- c(
  Seurat = "5.5.1",
  SeuratObject = "5.4.0",
  dplyr = "1.2.1",
  SingleCellExperiment = "1.34.0",
  scDblFinder = "1.26.7",
  BiocParallel = "1.46.0"
)
installed <- installed.packages(lib.loc = lib_dir)
missing <- setdiff(names(expected), rownames(installed))
if (length(missing) > 0) {
  stop("Packages missing after installation: ", paste(missing, collapse = ", "))
}
observed <- installed[names(expected), "Version"]
mismatch <- observed != expected
if (any(mismatch)) {
  stop(
    "Version validation failed: ",
    paste(names(expected)[mismatch], observed[mismatch], "!=", expected[mismatch], collapse = "; ")
  )
}

message("Validated doublet-assessment library: ", lib_dir)
