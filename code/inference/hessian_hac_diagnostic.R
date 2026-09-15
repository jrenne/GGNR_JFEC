# Diagnostic comparison of Hessian- and OPG-bread HAC covariance estimators.
#
# This script deliberately sources the score calculation so that the likelihood,
# parameter transformations, numerical score convention, and natural-parameter
# Jacobian are identical to those used for the published inference.

source("code/inference/sandwich_covariance.R")

hessian_step_scale <- as.numeric(Sys.getenv("GGNR_HESSIAN_STEP_SCALE", "5e-4"))
core_count <- as.integer(Sys.getenv("GGNR_HESSIAN_CORES", "4"))
if (!is.finite(hessian_step_scale) || hessian_step_scale <= 0) {
  stop("GGNR_HESSIAN_STEP_SCALE must be positive.")
}
if (!is.finite(core_count) || core_count < 1L) core_count <- 1L

objective <- function(value) sum(evaluate(value))
base_value <- objective(theta)
step <- hessian_step_scale * pmax(1, abs(theta))
parameter_pairs <- which(upper.tri(matrix(FALSE, parameter_count,
                                          parameter_count)), arr.ind = TRUE)

diagonal_value <- function(j) {
  plus <- minus <- theta
  plus[j] <- plus[j] + step[j]
  minus[j] <- minus[j] - step[j]
  -(objective(plus) - 2 * base_value + objective(minus)) / step[j]^2
}

off_diagonal_value <- function(k) {
  i <- parameter_pairs[k, 1]
  j <- parameter_pairs[k, 2]
  plus_plus <- plus_minus <- minus_plus <- minus_minus <- theta
  plus_plus[i] <- plus_plus[i] + step[i]
  plus_plus[j] <- plus_plus[j] + step[j]
  plus_minus[i] <- plus_minus[i] + step[i]
  plus_minus[j] <- plus_minus[j] - step[j]
  minus_plus[i] <- minus_plus[i] - step[i]
  minus_plus[j] <- minus_plus[j] + step[j]
  minus_minus[i] <- minus_minus[i] - step[i]
  minus_minus[j] <- minus_minus[j] - step[j]
  -(objective(plus_plus) - objective(plus_minus) - objective(minus_plus) +
      objective(minus_minus)) / (4 * step[i] * step[j])
}

cat(sprintf("Computing %d diagonal and %d off-diagonal Hessian elements\n",
            parameter_count, nrow(parameter_pairs)))
diagonal <- unlist(parallel::mclapply(
  seq_len(parameter_count), diagonal_value, mc.cores = core_count
))
off_diagonal <- unlist(parallel::mclapply(
  seq_len(nrow(parameter_pairs)), off_diagonal_value, mc.cores = core_count
))

hessian <- matrix(0, parameter_count, parameter_count,
                  dimnames = list(names(theta), names(theta)))
diag(hessian) <- diagonal
for (k in seq_len(nrow(parameter_pairs))) {
  i <- parameter_pairs[k, 1]
  j <- parameter_pairs[k, 2]
  hessian[i, j] <- hessian[j, i] <- off_diagonal[k]
}
hessian <- (hessian + t(hessian)) / 2

centered_scores <- sweep(scores, 2, colMeans(scores), FUN = "-")
hac_meat <- crossprod(centered_scores)
lag_count <- 12L
for (lag in seq_len(lag_count)) {
  weight <- 1 - lag / (lag_count + 1)
  lead <- centered_scores[(lag + 1L):observation_count, , drop = FALSE]
  lagged <- centered_scores[seq_len(observation_count - lag), , drop = FALSE]
  lag_product <- crossprod(lead, lagged)
  hac_meat <- hac_meat + weight * (lag_product + t(lag_product))
}
hac_meat <- (hac_meat + t(hac_meat)) / 2

hessian_decomposition <- eigen(hessian, symmetric = TRUE)
positive_threshold <- max(abs(hessian_decomposition$values)) * 1e-10
positive <- hessian_decomposition$values > positive_threshold
hessian_inverse <- hessian_decomposition$vectors[, positive, drop = FALSE] %*%
  diag(1 / hessian_decomposition$values[positive], nrow = sum(positive)) %*%
  t(hessian_decomposition$vectors[, positive, drop = FALSE])
dimnames(hessian_inverse) <- dimnames(hessian)

hessian_hac_working <- hessian_inverse %*% hac_meat %*% hessian_inverse
hessian_hac_working <- (hessian_hac_working + t(hessian_hac_working)) / 2
hessian_hac_natural <- jacobian %*% hessian_hac_working %*% t(jacobian)
hessian_hac_natural <- (hessian_hac_natural + t(hessian_hac_natural)) / 2

opg_hac_file <- file.path(
  "outputs", "diagnostics", "inference", "hac_sandwich_covariance.rds"
)
opg_hac <- readRDS(opg_hac_file)
comparison <- data.frame(
  parameter = names(natural_estimate),
  estimate = as.numeric(natural_estimate),
  hessian_hac_se = sqrt(pmax(0, diag(hessian_hac_natural))),
  opg_hac_se = sqrt(pmax(0, diag(opg_hac$hac_natural)))
)
comparison$se_ratio_hessian_to_opg <-
  comparison$hessian_hac_se / comparison$opg_hac_se

diagnostic_directory <- file.path(
  "outputs", "diagnostics", "inference", "hessian_diagnostic"
)
dir.create(diagnostic_directory, recursive = TRUE, showWarnings = FALSE)
step_tag <- gsub("[^0-9A-Za-z]+", "_", format(hessian_step_scale,
                                                scientific = TRUE))
write.csv(comparison, file.path(diagnostic_directory,
                                paste0("hessian_opg_hac_comparison_",
                                       step_tag, ".csv")),
          row.names = FALSE)
saveRDS(list(
  step_scale = hessian_step_scale,
  step = step,
  hessian = hessian,
  eigenvalues = hessian_decomposition$values,
  positive_rank = sum(positive),
  positive_threshold = positive_threshold,
  hessian_inverse = hessian_inverse,
  hac_meat = hac_meat,
  hessian_hac_working = hessian_hac_working,
  hessian_hac_natural = hessian_hac_natural,
  comparison = comparison
), file.path(diagnostic_directory,
             paste0("hessian_hac_diagnostic_", step_tag, ".rds")))

cat(sprintf("Hessian positive rank: %d/%d; eigenvalue range: %.6g to %.6g\n",
            sum(positive), parameter_count,
            min(hessian_decomposition$values),
            max(hessian_decomposition$values)))
print(comparison)
