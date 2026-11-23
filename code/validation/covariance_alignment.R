# Independent algebra and sample-size normalization checks for local inference.
source("code/validation/cache_provenance.R")
inference <- readRDS(file.path(Sys.getenv("GGNR_INFERENCE_DIR",
  "outputs/diagnostics/inference"), "hac_sandwich_covariance.rds"))
ggnr_check_provenance(inference, "Hessian inference")
relative_error <- function(x, y) max(abs(x - y)) / max(1, max(abs(y)))
scores <- inference$scores
centered <- sweep(scores, 2, colMeans(scores), "-")
n <- nrow(centered)
distance <- abs(outer(seq_len(n), seq_len(n), "-"))
weights <- pmax(1 - distance / (inference$lag_count + 1), 0)
S <- crossprod(centered, weights %*% centered)
B <- solve(inference$H)
V <- B %*% S %*% t(B)
mean_normalization <- solve(inference$H / n) %*% (S / n) %*%
  t(solve(inference$H / n)) / n
J <- inference$jacobian
checks <- data.frame(
  check = c("independent Toeplitz HAC", "direct versus scaled signed inverse",
            "independent covariance", "sum versus mean normalization",
            "economic-coordinate covariance"),
  relative_error = c(relative_error(S, inference$hac_meat),
    relative_error(B, inference$bread), relative_error(V, inference$hac_working),
    relative_error(mean_normalization, V),
    relative_error(J %*% V %*% t(J), inference$hac_natural))
)
stopifnot(max(checks$relative_error) < 1e-8,
          qr(J)$rank == length(inference$theta),
          all(diag(V) > 0), all(diag(inference$hac_natural) > 0))
white <- B %*% crossprod(scores) %*% t(B)
stopifnot(relative_error(white, inference$white_working) < 1e-8)
output_directory <- Sys.getenv("GGNR_INFERENCE_DIR", "outputs/diagnostics/inference")
write.csv(checks, file.path(output_directory, "covariance_audit.csv"), row.names = FALSE)
print(checks)
message("Signed-Hessian sandwich, centered HAC, units and normalization verified.")

# Reproduce the documented bandwidth checks and the archived half-step check.
reference_se <- sqrt(diag(inference$hac_natural))
sensitivity <- list()
for (lags in c(6L, 24L)) {
  weights <- pmax(1 - distance / (lags + 1), 0)
  meat <- crossprod(centered, weights %*% centered)
  se <- sqrt(diag(J %*% B %*% meat %*% t(B) %*% t(J)))
  sensitivity[[paste0("HAC_", lags)]] <- se / reference_se - 1
}
halfstep_file <- file.path(output_directory, "halfstep_check.rds")
if (file.exists(halfstep_file)) {
  half <- readRDS(halfstep_file)
  if (identical(half$theta, inference$theta)) {
    half_scores <- sweep(half$scores, 2, colMeans(half$scores), "-")
    weights <- pmax(1 - distance / (inference$lag_count + 1), 0)
    meat <- crossprod(half_scores, weights %*% half_scores)
    inverse <- solve(half$H)
    se <- sqrt(diag(J %*% inverse %*% meat %*% t(inverse) %*% t(J)))
    sensitivity$halfstep <- se / reference_se - 1
  }
}
sensitivity <- do.call(rbind, lapply(names(sensitivity), function(name) {
  value <- abs(sensitivity[[name]])
  data.frame(check = name, median_absolute_relative_change = median(value),
             maximum_absolute_relative_change = max(value))
}))
write.csv(sensitivity, file.path(output_directory, "inference_sensitivity.csv"), row.names = FALSE)
print(sensitivity)
