# HAC sandwich covariance based on monthly likelihood-score contributions.
# The outer-product-of-gradients (OPG) matrix supplies a stable bread;
# a Bartlett/Newey-West estimate of the long-run score covariance supplies the meat.

input_file <- Sys.getenv(
  "GGNR_COVARIANCE_FILE",
  "outputs/diagnostics/inference/sandwich_covariance.rds"
)
output_directory <- Sys.getenv(
  "GGNR_HAC_DIR",
  "outputs/diagnostics/inference"
)
lag_count <- as.integer(Sys.getenv("GGNR_HAC_LAGS", "12"))
if (!is.finite(lag_count) || lag_count < 0L) stop("Invalid HAC lag count.")

inference <- readRDS(input_file)
scores <- inference$scores
if (is.null(scores) || any(!is.finite(scores))) stop("Scores are unavailable.")
scores <- sweep(scores, 2, colMeans(scores), FUN = "-")
observation_count <- nrow(scores)
parameter_count <- ncol(scores)

opg <- crossprod(scores)
hac_meat <- opg
if (lag_count > 0L) {
  for (lag in seq_len(lag_count)) {
    weight <- 1 - lag / (lag_count + 1)
    lead <- scores[(lag + 1L):observation_count, , drop = FALSE]
    lagged <- scores[seq_len(observation_count - lag), , drop = FALSE]
    autocovariance <- crossprod(lead, lagged)
    hac_meat <- hac_meat + weight * (autocovariance + t(autocovariance))
  }
}
hac_meat <- (hac_meat + t(hac_meat)) / 2

positive_inverse <- function(value, tolerance = 1e-12) {
  decomposition <- eigen((value + t(value)) / 2, symmetric = TRUE)
  threshold <- max(abs(decomposition$values)) * tolerance
  keep <- decomposition$values > threshold
  inverse <- decomposition$vectors[, keep, drop = FALSE] %*%
    diag(1 / decomposition$values[keep], nrow = sum(keep)) %*%
    t(decomposition$vectors[, keep, drop = FALSE])
  list(inverse = inverse, rank = sum(keep), values = decomposition$values)
}
bread_decomposition <- positive_inverse(opg)
bread <- bread_decomposition$inverse
hac_working <- bread %*% hac_meat %*% bread
hac_working <- (hac_working + t(hac_working)) / 2
dimnames(hac_working) <- list(colnames(scores), colnames(scores))

jacobian <- inference$jacobian
if (is.null(jacobian) || ncol(jacobian) != parameter_count) {
  stop("Natural-parameter Jacobian is unavailable or incompatible.")
}
hac_natural <- jacobian %*% hac_working %*% t(jacobian)
hac_natural <- (hac_natural + t(hac_natural)) / 2
parameter_names <- rownames(jacobian)
parameter_names[parameter_names ==
  "liquidity_spread_gfc_pp.liquidity_spread_gfc"] <-
  "liquidity_spread_gfc_pp"
rownames(hac_natural) <- colnames(hac_natural) <- parameter_names

summary_table <- data.frame(
  parameter = parameter_names,
  estimate = as.numeric(inference$natural_estimate),
  opg_se = sqrt(pmax(0, diag(inference$opg_natural))),
  hac_sandwich_se = sqrt(pmax(0, diag(hac_natural)))
)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
write.csv(summary_table, file.path(output_directory, "parameter_inference.csv"),
          row.names = FALSE)
write.csv(hac_working,
          file.path(output_directory, "hac_sandwich_covariance_working.csv"))
write.csv(hac_natural,
          file.path(output_directory, "hac_sandwich_covariance_natural.csv"))
saveRDS(list(
  lag_count = lag_count, centered_scores = scores, opg = opg,
  hac_meat = hac_meat, bread = bread, hac_working = hac_working,
  jacobian = jacobian, hac_natural = hac_natural, summary = summary_table,
  bread_rank = bread_decomposition$rank,
  bread_eigenvalues = bread_decomposition$values
), file.path(output_directory, "hac_sandwich_covariance.rds"))
cat(sprintf("HAC sandwich: %d monthly lags; OPG bread rank %d/%d\n",
            lag_count, bread_decomposition$rank, parameter_count))
print(summary_table)
