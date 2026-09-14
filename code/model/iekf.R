# EKF for the six-state reduced-form/HFI model without output or growth.

source("code/model/model_specification.R")
source("code/data/load_data.R")

load_observed_real_short_rate_inputs <- function(dates) {
  data <- load_paper_data()
  rate <- data$effective_federal_funds_rate[match(dates, data$dates)]
  if (sum(is.finite(rate)) < 0.9 * length(dates)) stop("Incomplete FEDFUNDS coverage.")
  rate
}

build_filter_objects_reduced_hfi_no_output <- function(
    data, p = reduced_hfi_no_output_defaults(),
    hfi = no_hfi_data(data$dates), observed_short_rate = NULL) {
  if (is.null(observed_short_rate)) {
    observed_short_rate <- load_observed_real_short_rate_inputs(data$dates)
  }
  targets <- reduced_hfi_exact_targets(data, observed_short_rate)
  restricted <- apply_reduced_hfi_exact_moments(p, targets)
  p <- restricted$p
  base <- restricted$model
  L <- reduced_hfi_no_output_loadings()
  K <- nrow(base$Phi)
  T <- nrow(data$macro)
  maxmat <- max(c(data$mats_n, data$mats_r))
  hstep_max <- max(c(data$hstep_s, data$hstep_t))
  Phi_j <- array(0, dim = c(K, K, hstep_max))
  Phi_j[, , 1] <- base$Phi
  for (j in 2:hstep_max) Phi_j[, , j] <- Phi_j[, , j - 1L] %*% base$Phi
  coefs_p <- affine_coefs_r(
    base$Phi, base$Mu, base$Sigma2_X, 0, L$shadow_rate, 0, L$inflation,
    maxmat = maxmat
  )
  coefs_q <- affine_coefs_r(
    base$PhiQ, base$MuQ, base$Sigma2_X, 0, L$shadow_rate, 0, L$inflation,
    maxmat = maxmat
  )
  s_n <- sqrt(cumsum(diag(
    t(coefs_q$B_X_for) %*% base$Sigma2_X %*% coefs_q$B_X_for
  )))
  Gamma_macro <- rbind(L$inflation, L$inflation_target)
  Gamma_s <- matrix(0, ncol(data$surv_infexp), K)
  Gamma_s0 <- matrix(0, nrow(Gamma_s), 1)
  for (j in seq_len(nrow(Gamma_s))) {
    horizon <- data$hstep_s[j]
    Gamma_s[j, ] <- rowSums(
      coefs_p$B_X_pi[, seq_len(horizon), drop = FALSE]
    ) / horizon
    Gamma_s0[j, 1] <- sum(
      coefs_p$A_X_exp_pi[seq_len(horizon)] / horizon +
        0.5 * coefs_p$Conv_pi[seq_len(horizon)] / horizon^2
    )
  }
  gfc_liquidity <-
    data$dates >= as.Date("2008-01-01") & data$dates < as.Date("2010-01-01")
  covid_liquidity <-
    data$dates >= as.Date("2020-03-01") & data$dates < as.Date("2021-03-01")
  liquidity <- data$dates < as.Date("2004-01-01") |
    gfc_liquidity | covid_liquidity
  sigma_r_t <- rep(p$sigma_real_yield, T)
  sigma_r_t[liquidity] <- p$sigma_real_yield + p$sigma_real_yield_liquidity
  liquidity_spread_gfc <- if (is.null(p$liquidity_spread_gfc)) 0 else
    p$liquidity_spread_gfc
  liquidity_spread_covid <- if (is.null(p$liquidity_spread_covid)) 0 else
    p$liquidity_spread_covid
  real_yield_offset_t <- rep(0, T)
  real_yield_offset_t[gfc_liquidity] <- liquidity_spread_gfc
  real_yield_offset_t[covid_liquidity] <- liquidity_spread_covid
  hfi_state_covariance <- p$lambda_h * base$Sigma[, "eps_m", drop = FALSE]
  list(
    Mu = base$Mu, Phi = base$Phi, Sigma_X = base$Sigma,
    Sigma2_X = base$Sigma2_X, MuQ = base$MuQ, PhiQ = base$PhiQ,
    Phi_j = Phi_j, coefs_p = coefs_p, coefs_q = coefs_q, s_n = s_n,
    pars = p, loadings = L, targets = targets,
    exact_moments = restricted$moments,
    Gamma_macro = Gamma_macro, Gamma_core = matrix(L$core_inflation, nrow = 1),
    Gamma_s = Gamma_s, Gamma_s0 = Gamma_s0,
    sigma_r_t = sigma_r_t, real_yield_offset_t = real_yield_offset_t, hfi = hfi,
    hfi_state_covariance = hfi_state_covariance,
    hfi_variance = p$lambda_h^2 + p$sigma_h^2,
    dims = list(
      T = T, K = K, N_s = ncol(data$surv_infexp),
      N_t = ncol(data$surv_tbexp), N_n = ncol(data$yields_n),
      N_r = ncol(data$yields_r)
    )
  )
}

tbill_measurement_no_output <- function(x_pred, obj, horizons) {
  K <- obj$dims$K
  prediction <- rep(NaN, length(horizons))
  Gamma <- matrix(0, length(horizons), K)
  if (!length(horizons)) return(list(y_pred = prediction, Gamma = Gamma))
  max_h <- max(horizons)
  fit_all <- rep(NaN, max_h)
  Gamma_all <- matrix(NaN, max_h, K)
  x_h <- x_pred
  if (max_h > 1L) for (j in seq_len(max_h - 1L)) {
    x_h <- obj$Mu + obj$Phi %*% x_h
    if ((j - 5L) %% 12L == 0L) {
      fit <- y_fitting_r(
        x_h, obj$coefs_q$A_X_for[1:3], obj$coefs_q$B_X_for[, 1:3],
        obj$coefs_q$A_X_exp[1:3], obj$pars$r_lb, obj$s_n[1:3],
        obj$coefs_q$A_X_for_pi[1:3], obj$coefs_q$B_X_pi[, 1:3],
        obj$Sigma2_X, obj$coefs_q$B_X_cum[, 1:3],
        obj$coefs_q$B_X_cum_pi[, 1:3]
      )
      fit_all[j + 1L] <- fit$yfit_all_n[1, 3]
      Gamma_all[j + 1L, ] <- fit$JJ_n[3, ] %*% obj$Phi_j[, , j]
    }
  }
  for (j in seq_along(horizons)) {
    prediction[j] <- mean(fit_all[seq_len(horizons[j])], na.rm = TRUE)
    Gamma[j, ] <- colMeans(Gamma_all[seq_len(horizons[j]), , drop = FALSE],
                           na.rm = TRUE)
  }
  list(y_pred = prediction, Gamma = Gamma)
}

run_kf_reduced_hfi_no_output <- function(
    data = load_paper_data(), p = reduced_hfi_no_output_defaults(),
    hfi = no_hfi_data(data$dates), observed_short_rate = NULL,
    filter_method = c("EKF", "IEKF"), iekf_trigger_probability = 0,
    iekf_steps = 2L) {
  filter_method <- match.arg(filter_method)
  iekf_steps <- as.integer(iekf_steps)
  if (!is.finite(iekf_trigger_probability) ||
      iekf_trigger_probability < 0 || iekf_trigger_probability > 1) {
    stop("Invalid IEKF trigger probability.")
  }
  if (!is.finite(iekf_steps) || iekf_steps < 2L) {
    stop("The IEKF requires at least two measurement evaluations.")
  }
  obj <- build_filter_objects_reduced_hfi_no_output(
    data, p, hfi, observed_short_rate
  )
  p <- obj$pars
  K <- obj$dims$K
  T <- obj$dims$T
  # GDP growth, the output gap, and GDP10 are intentionally absent.
  # Estimates saved before this switch was introduced used core CPI.
  include_core <- if (is.null(p$include_core_cpi)) TRUE else isTRUE(p$include_core_cpi)
  y_level_all <- cbind(data$macro[, 2], data$macro[, 4])
  if (include_core) y_level_all <- cbind(y_level_all, data$core_inflation)
  y_level_all <- cbind(
    y_level_all, data$surv_infexp, data$surv_tbexp,
    data$yields_n, data$yields_r
  )
  R_macro <- c(p$sigma_inf^2, p$sigma_ptr^2)
  R_core <- p$sigma_core_cpi^2
  R_s <- rep(p$sigma_inflation_survey^2, obj$dims$N_s)
  if (obj$dims$N_s >= 2L && is.finite(p$sigma_inflation_survey_long)) {
    R_s[2] <- p$sigma_inflation_survey_long^2
  }
  R_t <- rep(p$sigma_tbill_survey^2, obj$dims$N_t)
  if (obj$dims$N_t >= 2L && is.finite(p$sigma_tbill_survey_long)) {
    R_t[2] <- p$sigma_tbill_survey_long^2
  }
  R_n <- rep(p$sigma_nominal_yield^2, obj$dims$N_n)
  P_pred <- stationary_covariance_no_output(list(
    Phi = obj$Phi, Sigma2_X = obj$Sigma2_X
  ))
  mean_state <- solve(diag(K) - obj$Phi, obj$Mu)
  x_pred <- matrix(0, T, K, dimnames = list(NULL, reduced_hfi_no_output_states))
  x_upd <- x_pred
  P_pred_all <- array(NA_real_, dim = c(K, K, T),
                      dimnames = list(reduced_hfi_no_output_states,
                                      reduced_hfi_no_output_states, NULL))
  P_upd_all <- array(NA_real_, dim = c(K, K, T),
                     dimnames = list(reduced_hfi_no_output_states,
                                     reduced_hfi_no_output_states, NULL))
  x_pred[1, ] <- mean_state
  fval_t <- numeric(T)
  iterated <- logical(T)

  measurement_system <- function(state, tt, P_prior) {
    yf <- y_fitting_r(
      state, obj$coefs_q$A_X_for, obj$coefs_q$B_X_for,
      obj$coefs_q$A_X_exp, p$r_lb, obj$s_n,
      obj$coefs_q$A_X_for_pi, obj$coefs_q$B_X_pi, obj$Sigma2_X,
      obj$coefs_q$B_X_cum, obj$coefs_q$B_X_cum_pi
    )
    y_pred_n <- yf$yfit_all_n[1, data$mats_n]
    y_pred_r <- yf$yfit_all_r[1, data$mats_r] +
      obj$real_yield_offset_t[tt]
    Gamma_n <- yf$JJ_n[data$mats_n, , drop = FALSE]
    Gamma_r <- yf$JJ_r[data$mats_r, , drop = FALSE]
    observed_tbill <- which(is.finite(data$surv_tbexp[tt, ]))
    tbill <- tbill_measurement_no_output(
      state, obj, data$hstep_t[observed_tbill]
    )
    y_pred_t <- rep(NaN, obj$dims$N_t)
    Gamma_t <- matrix(0, obj$dims$N_t, K)
    if (length(observed_tbill)) {
      y_pred_t[observed_tbill] <- tbill$y_pred
      Gamma_t[observed_tbill, ] <- tbill$Gamma
    }
    Gamma_full <- obj$Gamma_macro
    y_pred_full <- as.numeric(obj$Gamma_macro %*% state)
    if (include_core) {
      Gamma_full <- rbind(Gamma_full, obj$Gamma_core)
      y_pred_full <- c(y_pred_full, obj$Gamma_core %*% state)
    }
    Gamma_full <- rbind(Gamma_full, obj$Gamma_s, Gamma_t, Gamma_n, Gamma_r)
    y_pred_full <- c(
      y_pred_full, obj$Gamma_s0 + obj$Gamma_s %*% state,
      y_pred_t, y_pred_n, y_pred_r
    )
    observed_level <- is.finite(y_level_all[tt, ]) &
      rowSums(abs(Gamma_full)) > 0
    G <- Gamma_full[observed_level, , drop = FALSE]
    innovation <- y_level_all[tt, observed_level] -
      y_pred_full[observed_level]
    R_r <- rep(obj$sigma_r_t[tt]^2, obj$dims$N_r)
    R_diag_all <- R_macro
    if (include_core) R_diag_all <- c(R_diag_all, R_core)
    R_diag <- c(R_diag_all, R_s, R_t, R_n, R_r)[observed_level]
    Delta <- G %*% P_prior %*% t(G) +
      diag(R_diag, nrow = length(R_diag))
    cross_state_measurement <- P_prior %*% t(G)
    if (is.finite(hfi$standardized[tt])) {
      c_h <- obj$hfi_state_covariance
      cross_level_hfi <- G %*% c_h
      Delta <- rbind(cbind(Delta, cross_level_hfi),
                     c(t(cross_level_hfi), obj$hfi_variance))
      cross_state_measurement <- cbind(cross_state_measurement, c_h)
      innovation <- c(innovation, hfi$standardized[tt])
    }
    list(
      G = G, innovation = innovation, level_count = nrow(G),
      Delta = (Delta + t(Delta)) / 2,
      cross_state_measurement = cross_state_measurement
    )
  }

  for (tt in seq_len(T)) {
    if (tt > 1L) {
      x_pred[tt, ] <- obj$Mu + obj$Phi %*% x_upd[tt - 1L, ]
      P_pred <- obj$Phi %*% P_upd %*% t(obj$Phi) + obj$Sigma2_X
      P_pred <- (P_pred + t(P_pred)) / 2
    }
    P_pred_all[, , tt] <- P_pred
    system <- measurement_system(x_pred[tt, ], tt, P_pred)
    if (!length(system$innovation)) {
      fval_t[tt] <- 0
      x_upd[tt, ] <- x_pred[tt, ]
      P_upd <- P_pred
      P_upd_all[, , tt] <- P_upd
      next
    }
    chol_delta <- chol(system$Delta)
    Delta_inv <- chol2inv(chol_delta)
    fval_t[tt] <- -length(system$innovation) * 0.5 * log(2 * pi) -
      sum(log(diag(chol_delta))) -
      0.5 * as.numeric(crossprod(
        system$innovation, Delta_inv %*% system$innovation
      ))
    cross_gain <- system$cross_state_measurement %*% Delta_inv
    updated_mean <- as.numeric(
      x_pred[tt, ] + cross_gain %*% system$innovation
    )

    if (filter_method == "IEKF") {
      shadow_loading <- as.numeric(obj$loadings$shadow_rate)
      shadow_mean <- sum(shadow_loading * x_pred[tt, ])
      shadow_variance <- as.numeric(
        crossprod(shadow_loading, P_pred %*% shadow_loading)
      )
      lower_bound_probability <- pnorm(
        (p$r_lb - shadow_mean) / sqrt(shadow_variance)
      )
      if (lower_bound_probability >= iekf_trigger_probability) {
        iterated[tt] <- TRUE
        expansion_point <- updated_mean
        for (step in seq.int(2L, iekf_steps)) {
          system <- measurement_system(expansion_point, tt, P_pred)
          adjusted_innovation <- system$innovation
          level_index <- seq_len(system$level_count)
          adjusted_innovation[level_index] <-
            adjusted_innovation[level_index] + as.numeric(
              system$G %*% (expansion_point - x_pred[tt, ])
            )
          chol_delta <- chol(system$Delta)
          Delta_inv <- chol2inv(chol_delta)
          cross_gain <- system$cross_state_measurement %*% Delta_inv
          updated_mean <- as.numeric(
            x_pred[tt, ] + cross_gain %*% adjusted_innovation
          )
          expansion_point <- updated_mean
        }
      }
    }

    x_upd[tt, ] <- updated_mean
    P_upd <- P_pred - cross_gain %*% t(system$cross_state_measurement)
    P_upd <- (P_upd + t(P_upd)) / 2
    P_upd_all[, , tt] <- P_upd
  }
  list(
    fval = sum(fval_t), fval_t = fval_t,
    x_pred = x_pred, x_upd = x_upd, P_pred = P_pred_all,
    P_upd = P_upd_all, objects = obj,
    hfi_observation_count = sum(is.finite(hfi$standardized)),
    filter_method = filter_method, iterated = iterated,
    iteration_fraction = mean(iterated)
  )
}
