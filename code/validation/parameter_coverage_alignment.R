# Audit the published 500 + 500 cache; optionally replay selected samples.
source("code/validation/parameter_coverage_helpers.R")
directory <- Sys.getenv("GGNR_MC_DIR", "outputs/diagnostics/parameter_coverage")
x <- read.csv(file.path(directory, "replications.csv"))
design <- coverage_design()
row_key <- paste(coverage_key(x), x$parameter, sep = ":")
stopifnot(nrow(x) == 5000L, !anyDuplicated(row_key),
          setequal(unique(coverage_key(x)), coverage_key(design)),
          all(table(coverage_key(x)) == 5L), all(x$seed == 20310100L + x$rep),
          all(x$target_binding == ifelse(x$case == "reference_15", 0.15, 0.50)),
          all(x$valid), all(x$convergence == 0L),
          all(x$coverage_hac == (abs(x$estimate - x$true) <= 1.96 * x$hac_se)),
          max(abs(x$length_hac - 2 * 1.96 * x$hac_se)) < 1e-10)
summary <- summarize_parameter_coverage(x)
saved <- read.csv(file.path(directory, "summary.csv"))
key <- function(value) paste(value$case, value$parameter)
saved <- saved[match(key(summary), key(saved)), names(summary)]
stopifnot(isTRUE(all.equal(summary, saved, check.attributes = FALSE, tolerance = 1e-12)))
if (tolower(Sys.getenv("GGNR_MC_REPLAY", "false")) %in% c("true", "1")) {
  Rcpp::sourceCpp("code/validation/small_parameter_model.cpp")
  for (pair in list(c("reference_15", "1"), c("higher_zlb_50", "151"),
                    c("higher_zlb_50", "550"))) {
    case <- pair[1]; replication <- as.integer(pair[2])
    fresh <- run_parameter_replication(if (case == "reference_15") 0.15 else 0.50, replication)
    saved <- x[x$case == case & x$rep == replication, ]
    stopifnot(max(abs(fresh$estimate - saved$estimate)) < 1e-9,
              max(abs(fresh$hac_se - saved$hac_se)) < 1e-8)
  }
  message("Three independent sample replays match the archived estimates and standard errors.")
}
message("Verified: 500 reference + 500 sensitivity; no omitted sample; all coverage denominators and diagnostics.")
