# Joint final estimation of the six-state specification. Restrictions are
# imposed through smooth changes of variables; four unconditional moments and
# the selected long-survey/real-yield measurement errors remain exact/fixed.

source("code/model/iekf.R")
if (!requireNamespace("Rcpp", quietly = TRUE)) stop("Rcpp is required.")
Rcpp::sourceCpp("code/model/pricing.cpp")
Rcpp::sourceCpp("code/model/pricing_helpers.cpp")

input_file <- Sys.getenv(
  "GGNR_START_FILE",
  "estimates/start_values.rds"
)
output_file <- Sys.getenv(
  "GGNR_OUTPUT_FILE",
  "estimates/baseline.rds"
)
budget <- as.integer(Sys.getenv("GGNR_EVAL_BUDGET", "1500"))
optimizer_method <- match.arg(
  Sys.getenv("GGNR_OPTIMIZER", "bobyqa"),
  c("nlminb", "Nelder-Mead", "bobyqa")
)
nelder_mead_scale <- as.numeric(Sys.getenv("GGNR_NM_SCALE", "1"))
if (!is.finite(nelder_mead_scale) || nelder_mead_scale <= 0) {
  stop("GGNR_NM_SCALE must be a positive finite number.")
}
filter_method <- match.arg(
  Sys.getenv("GGNR_FILTER_METHOD", "IEKF"), c("EKF", "IEKF")
)
iekf_trigger_probability <- as.numeric(
  Sys.getenv("GGNR_IEKF_TRIGGER", "0")
)
iekf_steps <- as.integer(Sys.getenv("GGNR_IEKF_STEPS", "2"))
liquid_tips_only <- identical(
  tolower(Sys.getenv("GGNR_LIQUID_TIPS_ONLY", "false")), "true"
)
risk_state_scaling <- match.arg(
  Sys.getenv("GGNR_RISK_STATE_SCALING", "r_pi"),
  c("r", "r_pi", "w_only")
)

previous <- readRDS(input_file)
p_template <- modifyList(reduced_hfi_no_output_defaults(), previous$parameters)
p_template$include_core_cpi <- FALSE
p_template$sigma_inflation_survey_long <- 0.10 / 1200
p_template$sigma_tbill_survey_long <- 0.10 / 1200
p_template$sigma_real_yield <- 0.10 / 1200
p_template$sigma_real_yield_liquidity <- 0.20 / 1200
p_template$liquidity_spread_covid <- 0
pi_risk_start <- Sys.getenv("GGNR_PI_RISK_START", "")
p_template$lambda_pi_star_scale <- if (risk_state_scaling == "r_pi") {
  if (nzchar(pi_risk_start)) as.numeric(pi_risk_start) else
    previous$parameters$lambda_pi_star_scale
} else {
  0
}
if (!is.finite(p_template$lambda_pi_star_scale)) {
  stop("The starting pi-star risk-price scale must be finite.")
}
if (risk_state_scaling == "w_only") p_template$lambda_r_star_scale <- 0
if (liquid_tips_only) p_template$liquidity_spread_gfc <- 0

data <- load_true_release_data()
short_rate <- load_observed_real_short_rate_inputs(data$dates)
liquid_real_yield_dates <- data$dates >= as.Date("2004-01-01") &
  !(data$dates >= as.Date("2008-01-01") &
      data$dates < as.Date("2010-01-01")) &
  !(data$dates >= as.Date("2020-03-01") &
      data$dates < as.Date("2021-03-01"))
if (liquid_tips_only) {
  data$yields_r[!liquid_real_yield_dates, ] <- NaN
  cat(sprintf(
    paste0("Liquid-TIPS robustness: retaining real yields in %d of %d months; ",
           "the GFC liquidity intercept is fixed at zero.\n"),
    sum(liquid_real_yield_dates), length(liquid_real_yield_dates)
  ))
}
omit_start_text <- Sys.getenv("GGNR_OMIT_START", "")
omit_end_text <- Sys.getenv("GGNR_OMIT_END", "")
omitted_dates <- rep(FALSE, length(data$dates))
if (nzchar(omit_start_text) || nzchar(omit_end_text)) {
  if (!nzchar(omit_start_text) || !nzchar(omit_end_text)) {
    stop("Both GGNR_OMIT_START and GGNR_OMIT_END must be supplied.")
  }
  omit_start <- as.Date(omit_start_text)
  omit_end <- as.Date(omit_end_text)
  if (is.na(omit_start) || is.na(omit_end) || omit_end < omit_start) {
    stop("Invalid omitted-date interval.")
  }
  omitted_dates <- data$dates >= omit_start & data$dates <= omit_end
  for (name in c("macro", "yields_n", "yields_r", "surv_infexp",
                 "surv_gdpexp", "surv_tbexp")) {
    data[[name]][omitted_dates, ] <- NaN
  }
  data$core_inflation[omitted_dates] <- NaN
  data$any_survey_release[omitted_dates] <- FALSE
  data$any_nonmonthly_release[omitted_dates] <- FALSE
  short_rate[omitted_dates] <- NaN
  cat(sprintf("Omitting all observations from %s through %s (%d months).\n",
              omit_start, omit_end, sum(omitted_dates)))
}
hfi_none <- list(
  standardized = rep(NaN, length(data$dates)), observed = rep(FALSE, length(data$dates)),
  source = "none", shock = NA_character_, event_count = 0L, month_count = 0L
)

run_estimation_filter <- function(p) {
  run_kf_reduced_hfi_no_output(
    data, p, hfi_none, short_rate,
    filter_method = filter_method,
    iekf_trigger_probability = iekf_trigger_probability,
    iekf_steps = iekf_steps
  )
}

contractive_from_raw <- function(raw) {
  decomposition <- svd(raw)
  singular <- decomposition$d / sqrt(1 + decomposition$d^2)
  decomposition$u %*% diag(singular, nrow = length(singular)) %*%
    t(decomposition$v)
}
raw_from_contractive <- function(coefficient) {
  decomposition <- svd(coefficient)
  if (max(decomposition$d) >= 1) stop("Starting VAR is not contractive.")
  singular <- decomposition$d / sqrt(1 - decomposition$d^2)
  decomposition$u %*% diag(singular, nrow = length(singular)) %*%
    t(decomposition$v)
}

free_measurement <- c(
  "sigma_inf", "sigma_ptr", "sigma_inflation_survey",
  "sigma_tbill_survey", "sigma_nominal_yield"
)
to_working <- function(p) c(
  rho_r_star = qlogis(unname(p$rho_r_star)),
  sigma_r_star = log(1200 * unname(p$sigma_r_star)),
  rho_pi_star = qlogis(unname(p$rho_pi_star)),
  sigma_pi_star = log(1200 * unname(p$sigma_pi_star)),
  setNames(as.vector(raw_from_contractive(p$cyclical_var)),
           paste0("cyclical_raw_", seq_len(4))),
  sigma_pi_gap = log(1200 * unname(p$sigma_pi_gap)),
  impact_pi_m_ratio = unname(p$impact_pi_m_ratio),
  if (!liquid_tips_only)
    c(liquidity_spread_gfc = log(1200 * unname(p$liquidity_spread_gfc))),
  setNames(log(unname(unlist(p[free_measurement], use.names = FALSE))),
           paste0("measurement_", free_measurement)),
  setNames(as.numeric(p$lambda_0),
           paste0("lambda0_", reduced_hfi_no_output_shocks)),
  setNames(as.numeric(p$lambda_w),
           paste0("lambdaw_", reduced_hfi_no_output_shocks)),
  if (risk_state_scaling %in% c("r", "r_pi"))
    setNames(unname(p$lambda_r_star_scale), "lambda_r_star_scale"),
  if (risk_state_scaling == "r_pi")
    setNames(unname(p$lambda_pi_star_scale), "lambda_pi_star_scale"),
  rho_w = qlogis(unname(p$rho_w)),
  sigma_wpi = unname(p$sigma_wpi)
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
  p$liquidity_spread_gfc <- if (liquid_tips_only) 0 else
    exp(theta["liquidity_spread_gfc"]) / 1200
  for (name in free_measurement) {
    p[[name]] <- exp(theta[paste0("measurement_", name)])
  }
  p$lambda_0[reduced_hfi_no_output_shocks] <-
    theta[paste0("lambda0_", reduced_hfi_no_output_shocks)]
  p$lambda_w[reduced_hfi_no_output_shocks] <-
    theta[paste0("lambdaw_", reduced_hfi_no_output_shocks)]
  p$lambda_r_star_scale <- if (risk_state_scaling %in% c("r", "r_pi")) {
    unname(theta["lambda_r_star_scale"])
  } else {
    0
  }
  p$lambda_pi_star_scale <- if (risk_state_scaling == "r_pi") {
    unname(theta["lambda_pi_star_scale"])
  } else {
    0
  }
  p$rho_w <- plogis(theta["rho_w"])
  p$sigma_w <- sqrt(1 - p$rho_w^2)
  p$sigma_wpi <- unname(theta["sigma_wpi"])
  p
}

evaluations <- 0L
best_loglik <- -Inf
best_parameters <- p_template
objective <- function(theta) {
  evaluations <<- evaluations + 1L
  if (evaluations > budget) stop("Explicit evaluation budget reached.")
  if (any(!is.finite(theta)) || any(abs(theta) > 50)) return(1e12)
  p <- tryCatch(from_working(theta), error = function(error) {
    if (evaluations <= 3L) cat("Parameter transformation failed: ",
                               conditionMessage(error), "\n", sep = "")
    NULL
  })
  if (is.null(p)) return(1e12)
  fit <- tryCatch(
    run_estimation_filter(p),
    error = function(error) {
      if (evaluations <= 3L) cat("Filter evaluation failed: ",
                                 conditionMessage(error), "\n", sep = "")
      NULL
    }
  )
  if (is.null(fit) || !is.finite(fit$fval)) return(1e12)
  if (fit$fval > best_loglik) {
    best_loglik <<- fit$fval
    best_parameters <<- fit$objects$pars
  }
  if (evaluations %% 50L == 0L) {
    cat(sprintf("joint evaluation %d: log likelihood %.6f; best %.6f\n",
                evaluations, fit$fval, best_loglik))
  }
  -fit$fval
}

theta_start <- to_working(p_template)
initial_fit <- run_estimation_filter(p_template)
best_loglik <- initial_fit$fval
objective_named <- function(theta) {
  names(theta) <- names(theta_start)
  objective(theta)
}
optimizer <- tryCatch(
  if (optimizer_method == "nlminb") {
    nlminb(theta_start, objective_named,
           control = list(eval.max = budget, iter.max = budget,
                          rel.tol = 1e-8, x.tol = 1e-7))
  } else if (optimizer_method == "Nelder-Mead") {
    fit <- optim(theta_start, objective_named, method = "Nelder-Mead",
                 control = list(
                   maxit = budget, reltol = 1e-8,
                   parscale = rep(nelder_mead_scale, length(theta_start))
                 ))
    list(par = fit$par, objective = fit$value,
         convergence = fit$convergence, message = fit$message,
         iterations = NA_integer_, evaluations = fit$counts)
  } else {
    if (!requireNamespace("nloptr", quietly = TRUE)) {
      stop("The nloptr package is required for BOBYQA.")
    }
    fit <- nloptr::nloptr(
      x0 = theta_start, eval_f = objective_named,
      opts = list(
        algorithm = "NLOPT_LN_BOBYQA", maxeval = budget,
        xtol_rel = 1e-8, ftol_rel = 1e-10, print_level = 0
      )
    )
    list(
      par = setNames(fit$solution, names(theta_start)),
      objective = fit$objective,
      convergence = fit$status, message = fit$message,
      iterations = NA_integer_, evaluations = fit$iterations
    )
  },
  error = function(error) list(
    par = to_working(best_parameters), objective = -best_loglik,
    convergence = 98L, message = conditionMessage(error)
  )
)
candidate <- tryCatch(from_working(optimizer$par), error = function(error) NULL)
if (!is.null(candidate)) {
  candidate_fit <- tryCatch(
    run_estimation_filter(candidate),
    error = function(error) NULL
  )
  if (!is.null(candidate_fit) && candidate_fit$fval >= best_loglik) {
    best_loglik <- candidate_fit$fval
    best_parameters <- candidate_fit$objects$pars
  }
}
final_fit <- run_estimation_filter(best_parameters)
best_parameters <- final_fit$objects$pars
theta_final <- to_working(best_parameters)
result <- list(
  specification = reduced_hfi_no_output_id,
  stage = "joint_final", preliminary = FALSE,
  filter_method = filter_method,
  iekf_trigger_probability = iekf_trigger_probability,
  iekf_steps = iekf_steps,
  iteration_fraction = final_fit$iteration_fraction,
  input_file = input_file,
  initial_loglik = initial_fit$fval, final_loglik = final_fit$fval,
  parameters = best_parameters, working_parameters = theta_final,
  exact_targets = final_fit$objects$targets,
  exact_moments = final_fit$objects$exact_moments,
  omitted_interval = if (any(omitted_dates))
    c(start = omit_start_text, end = omit_end_text) else NULL,
  omitted_months = sum(omitted_dates),
  liquid_tips_only = liquid_tips_only,
  real_yield_months_retained = if (liquid_tips_only)
    sum(liquid_real_yield_dates) else nrow(data$yields_r),
  risk_state_scaling = risk_state_scaling,
  optimizer = optimizer, evaluations = evaluations, timestamp = Sys.time()
)
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
saveRDS(result, output_file)
cat(sprintf("Joint estimation: %.6f -> %.6f (%d evaluations)\n",
            initial_fit$fval, final_fit$fval, evaluations))
cat(sprintf("rho(r*), rho(pi*), rho(w): %.6f, %.6f, %.6f\n",
            best_parameters$rho_r_star, best_parameters$rho_pi_star,
            best_parameters$rho_w))
if (liquid_tips_only) {
  cat("Liquid-TIPS robustness: GFC and COVID real yields omitted; ",
      "GFC and COVID spreads fixed at zero.\n", sep = "")
} else {
  cat(sprintf("GFC liquidity spread: %.4f pp; COVID spread fixed at zero\n",
              1200 * best_parameters$liquidity_spread_gfc))
}
