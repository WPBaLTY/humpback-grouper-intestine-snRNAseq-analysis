options(stringsAsFactors = FALSE)

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", script_args, value = TRUE)[1])
script_dir <- if (length(script_file) == 0 || is.na(script_file)) {
  getwd()
} else {
  dirname(normalizePath(script_file, winslash = "/", mustWork = FALSE))
}
repo_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)

cran_repo <- Sys.getenv("GROUPER_CRAN_REPO", unset = "https://cloud.r-project.org")
target_library <- Sys.getenv(
  "GROUPER_R_LIBRARY",
  unset = file.path(repo_root, ".r-library")
)
dir.create(target_library, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = cran_repo, lib = target_library)
}
.libPaths(c(target_library, .libPaths()))

required_versions <- c(
  SeuratObject = "5.2.0",
  harmony = "1.2.3",
  Seurat = "5.3.0"
)

for (package in names(required_versions)) {
  required_version <- required_versions[[package]]
  installed_version <- if (requireNamespace(package, quietly = TRUE)) {
    as.character(packageVersion(package))
  } else {
    NA_character_
  }
  if (is.na(installed_version) || installed_version != required_version) {
    remotes::install_version(
      package,
      version = required_version,
      repos = cran_repo,
      lib = target_library,
      dependencies = NA,
      upgrade = "never"
    )
  }
}

installed_matrix <- installed.packages(lib.loc = target_library)
missing_from_target <- setdiff(names(required_versions), rownames(installed_matrix))
if (length(missing_from_target) > 0) {
  stop("Packages missing from exact target library: ", paste(missing_from_target, collapse = ", "))
}
installed_versions <- installed_matrix[names(required_versions), "Version"]
if (!identical(unname(installed_versions), unname(required_versions))) {
  stop("Exact package version installation did not validate")
}

message("Exact reconstruction library: ", normalizePath(target_library, winslash = "/"))
message(paste(names(installed_versions), installed_versions, collapse = "; "))
