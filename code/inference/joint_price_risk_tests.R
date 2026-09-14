# Joint Wald tests for the market-price-of-risk coefficients.

inference_file <- Sys.getenv(
  "GGNR_HAC_FILE",
  "outputs/diagnostics/inference/hac_sandwich_covariance.rds"
)
output_directory <- Sys.getenv("GGNR_RISK_TEST_DIR", dirname(inference_file))
inference <- readRDS(inference_file)

wald_test <- function(parameter_names, label, chi_square_reference = TRUE) {
  estimates <- inference$summary$estimate[
    match(parameter_names, inference$summary$parameter)
  ]
  covariance <- inference$hac_natural[parameter_names, parameter_names,
                                      drop = FALSE]
  decomposition <- eigen((covariance + t(covariance)) / 2, symmetric = TRUE)
  tolerance <- max(abs(decomposition$values)) * 1e-10
  keep <- decomposition$values > tolerance
  inverse <- decomposition$vectors[, keep, drop = FALSE] %*%
    diag(1 / decomposition$values[keep], nrow = sum(keep)) %*%
    t(decomposition$vectors[, keep, drop = FALSE])
  statistic <- as.numeric(crossprod(estimates, inverse %*% estimates))
  degrees_freedom <- sum(keep)
  log_p_value <- if (chi_square_reference) {
    pchisq(statistic, degrees_freedom, lower.tail = FALSE, log.p = TRUE)
  } else {
    NA_real_
  }
  data.frame(
    test = label, coefficient_count = length(parameter_names),
    covariance_rank = degrees_freedom, wald_statistic = statistic,
    degrees_freedom = if (chi_square_reference) degrees_freedom else NA_integer_,
    p_value = if (chi_square_reference) exp(log_p_value) else NA_real_,
    log10_p_value = if (chi_square_reference) log_p_value / log(10) else NA_real_,
    reference_distribution = if (chi_square_reference) {
      "chi-square"
    } else {
      paste0("diagnostic only: beta_r and beta_pi are unidentified under this null; ",
             "the usual chi-square reference is invalid")
    }
  )
}

names_all <- rownames(inference$hac_natural)
lambda_names <- c(
  grep("^lambda0_", names_all, value = TRUE),
  grep("^lambdaw_", names_all, value = TRUE)
)
results <- rbind(
  wald_test(lambda_names, "Lambda blocks", chi_square_reference = FALSE),
  wald_test(grep("^lambda0_", names_all, value = TRUE),
            "Constant lambda block"),
  wald_test(grep("^lambdaw_", names_all, value = TRUE),
            "State-dependent lambda block", chi_square_reference = FALSE)
)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
write.csv(results, file.path(output_directory, "joint_price_risk_tests.csv"),
          row.names = FALSE)
print(results)
