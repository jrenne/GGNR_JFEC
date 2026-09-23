# Numerical monthly scores and OPG covariance for the jointly estimated model.

source("code/model/iekf.R")
if (!requireNamespace("Rcpp", quietly = TRUE)) stop("Rcpp is required.")
Rcpp::sourceCpp("code/model/pricing_helpers.cpp")

estimate_file <- Sys.getenv(
  "GGNR_START_FILE",
  "estimates/baseline.rds"
)
output_directory <- Sys.getenv(
  "GGNR_COVARIANCE_DIR",
  "outputs/diagnostics/inference"
)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
estimate <- readRDS(estimate_file)
filter_method <- if (is.null(estimate$filter_method)) "EKF" else
  estimate$filter_method
iekf_trigger_probability <- if (is.null(estimate$iekf_trigger_probability))
  0 else estimate$iekf_trigger_probability
iekf_steps <- if (is.null(estimate$iekf_steps)) 2L else estimate$iekf_steps
p_template <- modifyList(reduced_hfi_no_output_defaults(), estimate$parameters)
p_template$include_core_cpi <- FALSE
p_template$sigma_inflation_survey_long <- 0.10 / 1200
p_template$sigma_tbill_survey_long <- 0.10 / 1200
p_template$sigma_real_yield <- 0.10 / 1200
p_template$sigma_real_yield_liquidity <- 0.20 / 1200
p_template$liquidity_spread_covid <- 0

data <- load_paper_data()
short_rate <- load_observed_real_short_rate_inputs(data$dates)
hfi_none <- list(
  standardized = rep(NaN, length(data$dates)), observed = rep(FALSE, length(data$dates)),
  source = "none", shock = NA_character_, event_count = 0L, month_count = 0L
)
targets <- reduced_hfi_exact_targets(data, short_rate)

contractive_from_raw <- function(raw) {
  decomposition <- svd(raw)
  singular <- decomposition$d / sqrt(1 + decomposition$d^2)
  decomposition$u %*% diag(singular, nrow = length(singular)) %*%
    t(decomposition$v)
}
free_measurement <- c(
  "sigma_inf", "sigma_ptr", "sigma_inflation_survey",
  "sigma_tbill_survey", "sigma_nominal_yield"
)
from_working <- function(theta) {
  p <- p_template
  p$rho_r_star <- plogis(theta["rho_r_star"])
  p$sigma_r_star <- exp(theta["sigma_r_star"]) / 1200
  p$rho_pi_star <- plogis(theta["rho_pi_star"])
  p$sigma_pi_star <- exp(theta["sigma_pi_star"]) / 1200
  raw <- matrix(theta[paste0("cyclical_raw_", seq_len(4))], 2, 2)
  p$cyclical_var <- contractive_from_raw(raw)
  dimnames(p$cyclical_var) <- list(c("m", "pi_gap"), c("m", "pi_gap"))
  p$sigma_pi_gap <- exp(theta["sigma_pi_gap"]) / 1200
  p$impact_pi_m_ratio <- unname(theta["impact_pi_m_ratio"])
  p$liquidity_spread_gfc <- exp(theta["liquidity_spread_gfc"]) / 1200
  for (name in free_measurement) {
    p[[name]] <- exp(theta[paste0("measurement_", name)])
  }
  p$lambda_0[reduced_hfi_no_output_shocks] <-
    theta[paste0("lambda0_", reduced_hfi_no_output_shocks)]
  p$lambda_w[reduced_hfi_no_output_shocks] <-
    theta[paste0("lambdaw_", reduced_hfi_no_output_shocks)]
  p$lambda_r_star_scale <- unname(theta["lambda_r_star_scale"])
  p$lambda_pi_star_scale <- unname(theta["lambda_pi_star_scale"])
  p$rho_w <- plogis(theta["rho_w"])
  p$sigma_w <- sqrt(1 - p$rho_w^2)
  p$sigma_wpi <- unname(theta["sigma_wpi"])
  p
}

theta <- estimate$working_parameters
if (is.null(theta)) stop("The estimate does not contain working parameters.")
parameter_count <- length(theta)
observation_count <- nrow(data$macro)
evaluate <- function(value) {
  p <- from_working(value)
  fit <- run_kf_reduced_hfi_no_output(
    data, p, hfi_none, short_rate,
    filter_method = filter_method,
    iekf_trigger_probability = iekf_trigger_probability,
    iekf_steps = iekf_steps
  )
  fit$fval_t
}

base_loglik <- evaluate(theta)
cat(sprintf("Likelihood at working-parameter reconstruction: %.6f (saved %.6f)\n",
            sum(base_loglik), estimate$final_loglik))
if (tolower(Sys.getenv("GGNR_REPEAT_ONLY", "false")) %in%
    c("true", "1", "yes")) {
  cat("Repeated likelihoods:",
      paste(vapply(seq_len(5), function(i) sum(evaluate(theta)), numeric(1)),
            collapse = ", "), "\n")
  for (name in c("rho_r_star", "rho_pi_star", "sigma_pi_star",
                 "sigma_pi_gap", "sigma_wpi", "liquidity_spread_gfc",
                 "lambda0_eps_r_star", "lambda0_eps_pi_star",
                 "lambdaw_eps_pi_star")) {
    step <- 1e-6 * max(1, abs(theta[name]))
    plus <- minus <- theta
    plus[name] <- plus[name] + step
    minus[name] <- minus[name] - step
    minus_values <- evaluate(minus)
    plus_values <- evaluate(plus)
    cat(name, sum(minus_values), sum(base_loglik), sum(plus_values),
        "largest monthly changes at",
        paste(order(abs(plus_values - minus_values), decreasing = TRUE)[1:3],
              collapse = ","), "\n")
  }
  quit(save = "no", status = 0)
}
if (tolower(Sys.getenv("GGNR_DIAGNOSTIC_ONLY", "false")) %in%
    c("true", "1", "yes")) quit(save = "no", status = 0)
score_step_scale <- as.numeric(Sys.getenv("GGNR_SCORE_STEP_SCALE", "5e-5"))
score_method <- match.arg(
  Sys.getenv("GGNR_SCORE_METHOD", "one-sided"), c("central", "one-sided")
)
score_step <- score_step_scale * pmax(1, abs(theta))
scores <- matrix(NA_real_, observation_count, parameter_count,
                 dimnames = list(NULL, names(theta)))
score_sides <- character(parameter_count)
for (j in seq_len(parameter_count)) {
  plus <- theta
  minus <- theta
  plus[j] <- plus[j] + score_step[j]
  minus[j] <- minus[j] - score_step[j]
  plus_values <- evaluate(plus)
  minus_values <- evaluate(minus)
  if (score_method == "central") {
    scores[, j] <- (plus_values - minus_values) / (2 * score_step[j])
    score_sides[j] <- "central"
  } else {
    # Select the locally smoother side of the piecewise-smooth lower-bound
    # likelihood. This remains useful with a fixed-iteration IEKF because the
    # lower-bound measurement equations themselves contain numerical kinks.
    plus_distance <- abs(sum(plus_values) - sum(base_loglik))
    minus_distance <- abs(sum(minus_values) - sum(base_loglik))
    if (plus_distance <= minus_distance) {
      scores[, j] <- (plus_values - base_loglik) / score_step[j]
      score_sides[j] <- "plus"
    } else {
      scores[, j] <- (base_loglik - minus_values) / score_step[j]
      score_sides[j] <- "minus"
    }
  }
  if (j %% 5L == 0L) cat(sprintf("Score derivatives: %d/%d\n", j, parameter_count))
}
meat <- crossprod(scores)
cat(sprintf("Maximum absolute numerical score: %.6f\n",
            max(abs(colSums(scores)))))

generalized_inverse_positive <- function(matrix_value, tolerance = 1e-12) {
  decomposition <- eigen((matrix_value + t(matrix_value)) / 2, symmetric = TRUE)
  threshold <- max(abs(decomposition$values)) * tolerance
  keep <- decomposition$values > threshold
  inverse <- decomposition$vectors[, keep, drop = FALSE] %*%
    diag(1 / decomposition$values[keep], nrow = sum(keep)) %*%
    t(decomposition$vectors[, keep, drop = FALSE])
  list(inverse = inverse, values = decomposition$values, rank = sum(keep),
       threshold = threshold)
}
opg_decomposition <- generalized_inverse_positive(meat)
opg_working <- opg_decomposition$inverse
dimnames(opg_working) <- list(names(theta), names(theta))

natural_parameters <- function(value) {
  p <- apply_reduced_hfi_exact_moments(from_working(value), targets)$p
  c(
    rho_r_star = unname(p$rho_r_star),
    sigma_r_star_pp = 1200 * unname(p$sigma_r_star),
    rho_pi_star = unname(p$rho_pi_star),
    sigma_pi_star_pp = 1200 * unname(p$sigma_pi_star),
    phi_m_m = p$cyclical_var["m", "m"],
    phi_m_pi = p$cyclical_var["m", "pi_gap"],
    phi_pi_m = p$cyclical_var["pi_gap", "m"],
    phi_pi_pi = p$cyclical_var["pi_gap", "pi_gap"],
    sigma_m_pp = 1200 * p$sigma_m,
    sigma_pi_gap_pp = 1200 * unname(p$sigma_pi_gap),
    sigma_u_pp = 1200 * p$sigma_u,
    impact_pi_m_ratio = p$impact_pi_m_ratio,
    liquidity_spread_gfc_pp = 1200 * unname(p$liquidity_spread_gfc),
    setNames(1200 * unlist(p[free_measurement], use.names = FALSE),
             paste0(free_measurement, "_pp")),
    setNames(as.numeric(p$lambda_0), paste0("lambda0_", reduced_hfi_no_output_shocks)),
    setNames(as.numeric(p$lambda_w), paste0("lambdaw_", reduced_hfi_no_output_shocks)),
    lambda_r_star_scale = p$lambda_r_star_scale,
    lambda_pi_star_scale = p$lambda_pi_star_scale,
    rho_w = p$rho_w,
    sigma_wpi = p$sigma_wpi
  )
}
natural_estimate <- natural_parameters(theta)
jacobian <- matrix(NA_real_, length(natural_estimate), parameter_count,
                   dimnames = list(names(natural_estimate), names(theta)))
for (j in seq_len(parameter_count)) {
  plus <- theta
  minus <- theta
  plus[j] <- plus[j] + score_step[j]
  minus[j] <- minus[j] - score_step[j]
  jacobian[, j] <- (natural_parameters(plus) - natural_parameters(minus)) /
    (2 * score_step[j])
}
opg_natural <- jacobian %*% opg_working %*% t(jacobian)

summary_table <- data.frame(
  parameter = names(natural_estimate), estimate = as.numeric(natural_estimate),
  opg_se = sqrt(pmax(0, diag(opg_natural)))
)
write.csv(summary_table, file.path(output_directory, "parameter_inference.csv"),
          row.names = FALSE)
write.csv(meat, file.path(output_directory, "outer_product_scores.csv"))
saveRDS(list(
    theta = theta, scores = scores, score_sides = score_sides,
  score_method = score_method, score_step_scale = score_step_scale,
  meat = meat, opg_working = opg_working, natural_estimate = natural_estimate,
  jacobian = jacobian,
  opg_natural = opg_natural, summary = summary_table,
  opg_rank = opg_decomposition$rank,
  filter_method = filter_method,
  iekf_trigger_probability = iekf_trigger_probability,
  iekf_steps = iekf_steps
), file.path(output_directory, "sandwich_covariance.rds"))
cat(sprintf("OPG retained rank: %d/%d\n",
            opg_decomposition$rank, parameter_count))
print(summary_table)
