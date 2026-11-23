# Resume the 500 + 500 exercise. Use a new GGNR_MC_DIR for a fresh experiment;
# archived samples are never silently recomputed or excluded.
source("code/validation/parameter_coverage_helpers.R")
Rcpp::sourceCpp("code/validation/small_parameter_model.cpp")
output_directory <- Sys.getenv("GGNR_MC_DIR", "outputs/diagnostics/parameter_coverage")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
replications <- as.integer(Sys.getenv("GGNR_MC_REPS", "500"))
case <- Sys.getenv("GGNR_MC_CASE", "all")
cores <- as.integer(Sys.getenv("GGNR_MC_CORES", "4"))
stopifnot(is.finite(cores), cores >= 1)
design <- coverage_design(replications, case)
archive <- file.path(output_directory, "replications.csv")
results <- if (file.exists(archive)) read.csv(archive) else NULL
if (!is.null(results)) {
  row_keys <- paste(coverage_key(results), results$parameter, sep = ":")
  if (anyDuplicated(row_keys)) stop("Duplicate cached Monte Carlo rows.")
  counts <- table(coverage_key(results))
  if (any(counts != 5L)) stop("Incomplete cached replications must be investigated.")
}
pending <- design[!coverage_key(design) %in% coverage_key(results), , drop = FALSE]
started <- Sys.time()
if (nrow(pending)) for (first in seq.int(1L, nrow(pending), by = 40L)) {
  batch <- pending[first:min(first + 39L, nrow(pending)), , drop = FALSE]
  worker <- function(i) {
    job <- batch[i, ]
    result <- run_parameter_replication(job$target, job$rep)
    failed <- result$convergence[1] != 0L
    if (failed) {
      write.csv(result, file.path(output_directory,
        paste0("initial_", job$case, "_", job$rep, ".csv")), row.names = FALSE)
      # Repeat exactly the same seed, generated data, and initial parameter vector.
      result <- run_parameter_replication(job$target, job$rep, 1600L, 200L)
      write.csv(result, file.path(output_directory,
        paste0("retry_", job$case, "_", job$rep, ".csv")), row.names = FALSE)
    }
    result$case <- job$case
    result$initial_optimizer_failure <- failed
    result$optimizer_retry <- failed
    result
  }
  rows <- if (.Platform$OS.type == "unix" && cores > 1) {
    parallel::mclapply(seq_len(nrow(batch)), worker, mc.cores = cores)
  } else lapply(seq_len(nrow(batch)), worker)
  fresh <- do.call(rbind, rows)
  results <- if (is.null(results)) fresh else rbind(results[names(fresh)], fresh)
  write.csv(results, archive, row.names = FALSE)
  message(sprintf("Completed %d replications; %.1f seconds elapsed.",
    length(unique(coverage_key(results))),
    as.numeric(difftime(Sys.time(), started, units = "secs"))))
}
selected <- results[coverage_key(results) %in% coverage_key(design), , drop = FALSE]
summary <- summarize_parameter_coverage(selected)
summary_file <- if (case == "all" && replications == 500L) "summary.csv" else
  paste0("summary_", case, "_", replications, ".csv")
write.csv(summary, file.path(output_directory, summary_file), row.names = FALSE)
print(summary[c("case", "parameter", "replications", "valid_intervals", "coverage_hac")])
