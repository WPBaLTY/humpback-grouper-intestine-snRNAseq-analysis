# Write local R session information for the analysis environment.
#
# By default this writes to 04_environment/R_sessionInfo.txt. The all-in-one
# driver overrides GROUPER_SESSION_INFO_FILE so session information is stored
# inside that run's output directory instead of modifying the code repository.

out_dir <- file.path(getwd(), "04_environment")
session_file <- Sys.getenv("GROUPER_SESSION_INFO_FILE", file.path(out_dir, "R_sessionInfo.txt"))
dir.create(dirname(session_file), recursive = TRUE, showWarnings = FALSE)

sink(session_file)
print(sessionInfo())
sink()

message("Wrote R session information to: ", session_file)
