# Optional setup helper.
#
# Run this script only if required packages are missing:
#   Rscript 00_setup/00_install_required_packages.R
#
# The original manuscript analysis used R 4.4.0 with Seurat 5.3.0 and Harmony
# 0.1. The script installs current CRAN binaries sufficient to rerun the
# repository analysis scripts from the deposited RDS objects.

options(stringsAsFactors = FALSE)

cran_repo <- Sys.getenv("CRAN_REPO", "https://cloud.r-project.org")
options(repos = c(CRAN = cran_repo))

required_packages <- c("Seurat", "harmony", "ggplot2", "patchwork", "scales")
available <- vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
missing <- names(available)[!available]

cat("R version: ", R.version.string, "\n", sep = "")
cat("CRAN mirror: ", cran_repo, "\n", sep = "")

if (length(missing) > 0) {
  cat("Installing missing packages: ", paste(missing, collapse = ", "), "\n", sep = "")
  install.packages(missing, dependencies = TRUE)
} else {
  cat("All required packages are already installed.\n")
}

cat("\nInstalled package versions:\n")
for (pkg in required_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(pkg, as.character(utils::packageVersion(pkg)), "\n")
  } else {
    cat(pkg, "NOT INSTALLED\n")
  }
}
