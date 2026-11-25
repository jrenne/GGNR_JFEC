# Small provenance records prevent mixing estimates, data, or model versions.
ggnr_cache_provenance <- function(estimate_file = Sys.getenv(
    "GGNR_START_FILE", "estimates/baseline_r2.rds")) {
  files <- c(estimate_file, "data/input/paper_data.rds",
             list.files("code/model", full.names = TRUE, pattern = "\\.(R|cpp|h)$"),
             "code/data/load_data.R", "code/inference/parameterization.R")
  list(files = files, md5 = unname(tools::md5sum(files)))
}

ggnr_check_provenance <- function(result, label) {
  expected <- ggnr_cache_provenance()
  if (is.null(result$provenance) || !identical(result$provenance, expected)) {
    stop(label, " does not match the selected estimate, data, and model code. ",
         "Recompute that exercise before rebuilding figures or tables.")
  }
  invisible(TRUE)
}
