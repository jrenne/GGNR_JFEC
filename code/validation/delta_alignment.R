# Independent full-covariance algebra and archived numerical step check.
source("code/validation/cache_provenance.R")
directory <- Sys.getenv("GGNR_UNCERTAINTY_DIR", "outputs/diagnostics/uncertainty")
z <- readRDS(file.path(directory, "delta_uncertainty.rds"))
v <- z$parameter_covariance
e <- eigen((v + t(v)) / 2, symmetric = TRUE)
stopifnot(min(e$values) > 0, length(z$theta) == 30L)
root <- sweep(e$vectors, 2, sqrt(e$values), "*")
errors <- numeric(ncol(z$point))
for (j in seq_len(ncol(z$point))) {
  independent <- rowSums((z$D[, j, ] %*% root)^2)
  errors[j] <- max(abs(independent - z$parameter_variance[, j]))
}
stopifnot(max(errors) < 1e-8,
          max(abs(z$total_sd^2 - z$filter_variance - z$parameter_variance)) < 1e-10)
for (item in c("state_bands", "premium_bands")) {
  bands <- z[[item]]
  for (name in sub("_estimate$", "", names(bands)[grepl("_estimate$", names(bands))])) {
    stopifnot(max(abs(bands[[paste0(name, "_estimate")]] - z$point[, name])) < 1e-12,
      max(abs(bands[[paste0(name, "_p975")]] - z$point[, name] -
        qnorm(.975) * z$total_sd[, name])) < 1e-10,
      max(abs(bands[[paste0(name, "_p025")]] - z$point[, name] +
        qnorm(.975) * z$total_sd[, name])) < 1e-10)
  }
}
half_file <- file.path(directory, "half_step.rds")
half_change <- rep(NA_real_, ncol(z$point))
if (file.exists(half_file)) {
  h <- readRDS(half_file)
  ggnr_check_provenance(h, "Uncertainty half-step check")
  stopifnot(isTRUE(all.equal(h$theta, z$theta)),
            isTRUE(all.equal(h$parameter_covariance, v)),
            h$step_scale == z$step_scale / 2)
  half_change <- apply(abs(h$total_sd / z$total_sd - 1), 2, max)
  selected <- c("r_star", "pi_star", "nominal_10y", "real_10y", "inflation_10y")
  stopifnot(max(half_change[selected]) < 2e-5)
}
write.csv(data.frame(series = colnames(z$point),
  max_spectral_reconstruction_error = errors,
  max_relative_half_step_se_change = half_change),
  file.path(directory, "delta_alignment.csv"), row.names = FALSE)
message("Local bands: full covariance, variance decomposition and numerical step checks agree.")
