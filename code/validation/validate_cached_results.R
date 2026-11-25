# Validate the identity and internal consistency of archived R2 exercises.
source("code/validation/cache_provenance.R")
inference <- readRDS(file.path(Sys.getenv("GGNR_INFERENCE_DIR",
  "outputs/diagnostics/inference"), "hac_sandwich_covariance.rds"))
filtering <- readRDS(file.path(Sys.getenv("GGNR_FILTER_OUTPUT", "outputs/diagnostics/filtering"),
                               "filtering_accuracy.rds"))
uncertainty_directory <- Sys.getenv("GGNR_UNCERTAINTY_DIR", "outputs/diagnostics/uncertainty")
uncertainty <- readRDS(file.path(uncertainty_directory, "nonlinear_uncertainty.rds"))
ggnr_check_provenance(inference, "Empirical inference")
ggnr_check_provenance(filtering, "Filtering simulation")
ggnr_check_provenance(uncertainty, "Parameter-and-filtering uncertainty")
strip_scalar_names <- function(x) {
  if (is.list(x)) return(lapply(x, strip_scalar_names))
  if (length(x) == 1L) return(unname(x))
  x
}
estimate <- readRDS(Sys.getenv("GGNR_START_FILE", "estimates/baseline_r2.rds"))
stopifnot(isTRUE(all.equal(inference$theta, estimate$working_parameters, tolerance = 1e-12)),
          isTRUE(all.equal(uncertainty$theta, inference$theta, tolerance = 1e-12)),
          isTRUE(all.equal(strip_scalar_names(filtering$parameters), strip_scalar_names(estimate$parameters), tolerance = 1e-10)),
          isTRUE(all.equal(uncertainty$parameter_covariance, inference$hac_working,
                           tolerance = 1e-8)),
          abs(uncertainty$diagnostics$reconstructed_loglik - estimate$final_loglik) < 1e-7)
for (item in c("state_bands", "premium_bands")) {
  csv <- read.csv(file.path(uncertainty_directory, paste0(item, ".csv")))
  saved <- uncertainty[[item]]
  csv$date <- as.Date(csv$date); saved$date <- as.Date(saved$date)
  stopifnot(isTRUE(all.equal(csv, saved, check.attributes = FALSE, tolerance = 1e-10)))
}
message("Cached empirical exercises match the selected estimate and covariance.")
