# Shared unconstrained-to-economic transformations for inference and uncertainty.
# This file defines functions only; sourcing it neither computes nor writes results.

source("code/model/iekf.R")

make_parameterization <- function(parameters, targets) {
    p_template <- modifyList(reduced_hfi_no_output_defaults(), parameters)
    p_template$include_core_cpi <- FALSE
    p_template$sigma_inflation_survey_long <- 0.1/1200
    p_template$sigma_tbill_survey_long <- 0.1/1200
    p_template$sigma_real_yield <- 0.1/1200
    p_template$sigma_real_yield_liquidity <- 0.2/1200
    p_template$liquidity_spread_covid <- 0
    contractive_from_raw <- function(raw) {
        decomposition <- svd(raw)
        singular <- decomposition$d/sqrt(1 + decomposition$d^2)
        decomposition$u %*% diag(singular, nrow = length(singular)) %*% t(decomposition$v)
    }
    free_measurement <- c("sigma_inf", "sigma_ptr", "sigma_inflation_survey", "sigma_tbill_survey",
        "sigma_nominal_yield")
    from_working <- function(theta) {
        p <- p_template
        p$rho_r_star <- plogis(theta["rho_r_star"])
        p$sigma_r_star <- exp(theta["sigma_r_star"])/1200
        p$rho_pi_star <- plogis(theta["rho_pi_star"])
        p$sigma_pi_star <- exp(theta["sigma_pi_star"])/1200
        raw <- matrix(theta[paste0("cyclical_raw_", seq_len(4))], 2, 2)
        p$cyclical_var <- contractive_from_raw(raw)
        dimnames(p$cyclical_var) <- list(c("m", "pi_gap"), c("m", "pi_gap"))
        p$sigma_pi_gap <- exp(theta["sigma_pi_gap"])/1200
        p$impact_pi_m_ratio <- unname(theta["impact_pi_m_ratio"])
        p$liquidity_spread_gfc <- exp(theta["liquidity_spread_gfc"])/1200
        for (name in free_measurement) {
            p[[name]] <- exp(theta[paste0("measurement_", name)])
        }
        p$lambda_0[reduced_hfi_no_output_shocks] <- theta[paste0("lambda0_", reduced_hfi_no_output_shocks)]
        p$lambda_w[reduced_hfi_no_output_shocks] <- theta[paste0("lambdaw_", reduced_hfi_no_output_shocks)]
        p$lambda_r_star_scale <- unname(theta["lambda_r_star_scale"])
        p$lambda_pi_star_scale <- unname(theta["lambda_pi_star_scale"])
        p$rho_w <- plogis(theta["rho_w"])
        p$sigma_w <- sqrt(1 - p$rho_w^2)
        p$sigma_wpi <- unname(theta["sigma_wpi"])
        p
    }
    natural_parameters <- function(value) {
        p <- apply_reduced_hfi_exact_moments(from_working(value), targets)$p
        c(rho_r_star = unname(p$rho_r_star), sigma_r_star_pp = 1200 * unname(p$sigma_r_star), rho_pi_star = unname(p$rho_pi_star),
            sigma_pi_star_pp = 1200 * unname(p$sigma_pi_star), phi_m_m = p$cyclical_var["m", "m"],
            phi_m_pi = p$cyclical_var["m", "pi_gap"], phi_pi_m = p$cyclical_var["pi_gap", "m"], phi_pi_pi = p$cyclical_var["pi_gap",
                "pi_gap"], sigma_m_pp = 1200 * p$sigma_m, sigma_pi_gap_pp = 1200 * unname(p$sigma_pi_gap),
            sigma_u_pp = 1200 * p$sigma_u, impact_pi_m_ratio = p$impact_pi_m_ratio, liquidity_spread_gfc_pp = 1200 *
                unname(p$liquidity_spread_gfc), setNames(1200 * unlist(p[free_measurement], use.names = FALSE),
                paste0(free_measurement, "_pp")), setNames(as.numeric(p$lambda_0), paste0("lambda0_",
                reduced_hfi_no_output_shocks)), setNames(as.numeric(p$lambda_w), paste0("lambdaw_",
                reduced_hfi_no_output_shocks)), lambda_r_star_scale = p$lambda_r_star_scale, lambda_pi_star_scale = p$lambda_pi_star_scale,
            rho_w = unname(p$rho_w), sigma_wpi = p$sigma_wpi)
    }
    list(from_working = from_working, natural = natural_parameters)
}

make_inference_problem <- function(estimate_file = Sys.getenv("GGNR_START_FILE", "estimates/baseline_r2.rds")) {
    estimate <- readRDS(estimate_file)
    data <- load_paper_data()
    short_rate <- load_observed_real_short_rate_inputs(data$dates)
    hfi <- no_hfi_data(data$dates)
    targets <- reduced_hfi_exact_targets(data, short_rate)
    mapping <- make_parameterization(estimate$parameters, targets)
    filter_method <- if (is.null(estimate$filter_method))
        "EKF"
    else estimate$filter_method
    trigger <- if (is.null(estimate$iekf_trigger_probability))
        0
    else estimate$iekf_trigger_probability
    steps <- if (is.null(estimate$iekf_steps))
        2L
    else estimate$iekf_steps
    evaluate_fit <- function(theta, branches = NULL, record = FALSE) {
        run_kf_reduced_hfi_no_output(data, mapping$from_working(theta), hfi, short_rate, filter_method = filter_method,
            iekf_trigger_probability = trigger, iekf_steps = steps, nominal_branches = branches,
            record_branches = record)
    }
    list(estimate = estimate, theta = estimate$working_parameters, from_working = mapping$from_working,
        natural = mapping$natural, evaluate_fit = evaluate_fit, evaluate = function(theta, branches = NULL) evaluate_fit(theta,
            branches)$fval_t, targets = targets, estimate_file = estimate_file)
}

