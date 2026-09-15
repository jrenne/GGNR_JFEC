source("code/model/var_coefficients.R")

affine_coefs_r <- function(Phi_q, Mu_q, Sigma2_X, delta_0, delta_1,
                           delta_pi_0, delta_pi_1, maxmat = 120, freq = 1) {
  K <- length(Mu_q)
  B_X_for <- matrix(0, K, maxmat)
  B_X_for[, 1] <- delta_1
  A_X_exp <- matrix(0, maxmat, 1)
  A_X_exp[1, 1] <- delta_0
  A_X_for <- matrix(0, maxmat, 1)
  A_X_for[1, 1] <- delta_0
  B_X_cum <- matrix(0, K, maxmat)

  for (j in 2:maxmat) {
    B_X_for[, j] <- t(Phi_q) %*% B_X_for[, j - 1]
    A_X_exp[j, 1] <- A_X_exp[j - 1, 1] + crossprod(B_X_for[, j - 1], Mu_q)
    B_X_cum[, j] <- B_X_cum[, j - 1] + B_X_for[, j - 1]
    A_X_for[j, 1] <- A_X_exp[j, 1] -
      0.5 * crossprod(B_X_cum[, j], Sigma2_X %*% B_X_cum[, j]) / freq
  }

  B_X_pi <- matrix(0, K, maxmat)
  B_X_pi[, 1] <- t(Phi_q) %*% delta_pi_1
  A_X_exp_pi <- matrix(0, maxmat, 1)
  A_X_exp_pi[1, 1] <- delta_pi_0 + crossprod(delta_pi_1, Mu_q)
  B_X_cum_pi <- matrix(0, K, maxmat)
  B_X_cum_pi[, 1] <- delta_pi_1
  Conv_pi <- matrix(0, maxmat, 1)
  Conv_pi[1, 1] <- crossprod(B_X_cum_pi[, 1],
                             Sigma2_X %*% B_X_cum_pi[, 1]) / freq
  A_X_for_pi <- matrix(0, maxmat, 1)
  A_X_for_pi[1, 1] <- A_X_exp_pi[1, 1] + 0.5 * Conv_pi[1, 1]

  for (j in 2:maxmat) {
    B_X_pi[, j] <- t(Phi_q) %*% B_X_pi[, j - 1]
    A_X_exp_pi[j, 1] <- A_X_exp_pi[j - 1, 1] +
      crossprod(B_X_pi[, j - 1], Mu_q)
    B_X_cum_pi[, j] <- B_X_cum_pi[, j - 1] + B_X_pi[, j - 1]
    Conv_pi[j, 1] <- crossprod(B_X_cum_pi[, j],
                               Sigma2_X %*% B_X_cum_pi[, j]) / freq
    A_X_for_pi[j, 1] <- A_X_exp_pi[j, 1] + 0.5 * Conv_pi[j, 1]
  }

  list(
    A_X_for = A_X_for, B_X_for = B_X_for, A_X_exp = A_X_exp,
    B_X_cum = B_X_cum, A_X_for_pi = A_X_for_pi, B_X_pi = B_X_pi,
    A_X_exp_pi = A_X_exp_pi, B_X_cum_pi = B_X_cum_pi, Conv_pi = Conv_pi
  )
}

y_fitting_r_pure <- function(X0, A_X_for, B_X_for, A_X_exp, r_lb, s_n,
                             A_X_for_pi, B_X_pi, Sigma2_X, B_X_cum,
                             B_X_cum_pi) {
  X0 <- as.matrix(X0)
  if (nrow(X0) != nrow(B_X_for)) {
    X0 <- t(X0)
  }
  K <- nrow(X0)
  T0 <- ncol(X0)
  Mm <- length(A_X_for)

  y_fit_short <- as.vector(A_X_exp[1] + t(X0) %*% B_X_for[, 1])
  Probs_short <- rep(1, T0)
  Probs_short[y_fit_short < r_lb] <- 0
  y_fit_short[y_fit_short < r_lb] <- r_lb

  mu <- sweep(t(X0) %*% B_X_for[, 2:Mm, drop = FALSE],
              2, as.vector(A_X_exp[2:Mm]) - r_lb, "+")
  z_n <- sweep(mu, 2, s_n[1:(Mm - 1)], "/")
  Probs_long <- pnorm(z_n)
  Probs <- cbind(Probs_short, Probs_long)
  pdf_z <- dnorm(z_n)

  # Wu-Xia nominal term: the convexity correction belongs inside g.
  # The real-rate interaction retains the unshifted probability Probs_long.
  nominal_mu <- sweep(mu, 2,
                      as.vector(A_X_for[2:Mm] - A_X_exp[2:Mm]), "+")
  nominal_z <- sweep(nominal_mu, 2, s_n[1:(Mm - 1)], "/")
  nominal_probability <- pnorm(nominal_z)
  f_fit_long <- r_lb + nominal_mu * nominal_probability +
    sweep(dnorm(nominal_z), 2, s_n[1:(Mm - 1)], "*")
  f_fit <- cbind(y_fit_short, f_fit_long)
  yfit_all_n <- t(apply(f_fit, 1, cumsum)) / matrix(seq_len(Mm), T0, Mm,
                                                    byrow = TRUE)

  pi_fit <- sweep(t(X0) %*% B_X_pi, 2, as.vector(A_X_for_pi), "+")
  BXcumSig2BXcum_pi <- c(
    0,
    colSums(B_X_cum[, 2:Mm, drop = FALSE] *
              (Sigma2_X %*% B_X_cum_pi[, 2:Mm, drop = FALSE]))
  )
  r_fit <- f_fit - pi_fit + sweep(Probs, 2, BXcumSig2BXcum_pi, "*")
  yfit_all_r <- t(apply(r_fit, 1, cumsum)) / matrix(seq_len(Mm), T0, Mm,
                                                    byrow = TRUE)

  JJ_f <- matrix(0, K * T0, Mm)
  JJ_r <- matrix(0, K * T0, Mm)
  active_short <- as.vector(A_X_exp[1] + t(X0) %*% B_X_for[, 1]) > r_lb
  for (tt in seq_len(T0)) {
    rows <- ((tt - 1) * K + 1):(tt * K)
    if (active_short[tt]) {
      JJ_f[rows, 1] <- B_X_for[, 1]
    }
    JJ_r[rows, 1] <- JJ_f[rows, 1] - B_X_pi[, 1]
    for (m in 2:Mm) {
      JJ_f[rows, m] <- nominal_probability[tt, m - 1] * B_X_for[, m]
      JJ_r[rows, m] <- JJ_f[rows, m] - B_X_pi[, m] +
        (pdf_z[tt, m - 1] / s_n[m - 1]) * B_X_for[, m] *
        BXcumSig2BXcum_pi[m]
    }
  }
  JJ_n <- t(apply(JJ_f, 1, cumsum)) / matrix(seq_len(Mm), K * T0, Mm,
                                             byrow = TRUE)
  JJ_r <- t(apply(JJ_r, 1, cumsum)) / matrix(seq_len(Mm), K * T0, Mm,
                                             byrow = TRUE)

  list(
    yfit_all_n = yfit_all_n, yfit_all_r = yfit_all_r,
    JJ_n = t(JJ_n), JJ_r = t(JJ_r), Probs = Probs
  )
}

y_fitting_r <- function(X0, A_X_for, B_X_for, A_X_exp, r_lb, s_n,
                        A_X_for_pi, B_X_pi, Sigma2_X, B_X_cum,
                        B_X_cum_pi, use_cpp = TRUE) {
  X0 <- as.matrix(X0)
  if (nrow(X0) != nrow(B_X_for)) {
    X0 <- t(X0)
  }
  if (use_cpp && exists("y_fitting_jfec_cpp", mode = "function")) {
    return(y_fitting_jfec_cpp(
      X0, as.vector(A_X_for), B_X_for, as.vector(A_X_exp), r_lb,
      as.vector(s_n), as.vector(A_X_for_pi), B_X_pi, Sigma2_X, B_X_cum,
      B_X_cum_pi
    ))
  }
  y_fitting_r_pure(X0, A_X_for, B_X_for, A_X_exp, r_lb, s_n,
                   A_X_for_pi, B_X_pi, Sigma2_X, B_X_cum, B_X_cum_pi)
}

apply_jfec_calibrations <- function(x0, macro, calibrate_growth_mean = TRUE,
                                    calibrate_gdp_measurement = TRUE,
                                    calibrate_growth_survey_measurement = TRUE) {
  # These are the non-free quantities imposed inside the likelihood.
  if (calibrate_growth_mean) x0[6] <- mean(macro[, 1], na.rm = TRUE)
  x0[12] <- mean(macro[, 2], na.rm = TRUE)
  if (calibrate_gdp_measurement) x0[31] <- x0[30] / 100
  x0[32] <- x0[30] / 100
  x0[33] <- 5 * x0[34]
  if (calibrate_growth_survey_measurement) x0[35] <- 5 * x0[34]
  x0[36] <- x0[34] / 4
  x0
}

effective_real_variance_gaussian <- function(mean_state, covariance_state,
                                             inflation_loading,
                                             lower_bound = 0) {
  mean_shadow <- mean_state[1]
  variance_shadow <- covariance_state[1, 1]
  sd_shadow <- sqrt(variance_shadow)
  if (!is.finite(sd_shadow) || sd_shadow <= 0) stop("Invalid shadow-rate variance.")
  centered_mean <- mean_shadow - lower_bound
  d <- centered_mean / sd_shadow
  probability_unconstrained <- pnorm(d)
  density <- dnorm(d)
  first_positive <- centered_mean * probability_unconstrained + sd_shadow * density
  second_positive <- (centered_mean^2 + variance_shadow) * probability_unconstrained +
    centered_mean * sd_shadow * density
  variance_effective_nominal <- second_positive - first_positive^2
  variance_inflation <- as.numeric(
    t(inflation_loading) %*% covariance_state %*% inflation_loading
  )
  covariance_shadow_inflation <- as.numeric(
    covariance_state[1, , drop = FALSE] %*% inflation_loading
  )
  covariance_effective_inflation <-
    probability_unconstrained * covariance_shadow_inflation
  variance_effective_nominal + variance_inflation -
    2 * covariance_effective_inflation
}

apply_inflation_variance_restriction <- function(
    x0, macro, yields_n = NULL, mats_n = NULL,
    observed_short_rate = NULL,
    taylor_rule = "orphanides", natural_rate_identity = "exact",
    policy_inflation = "core", match_real_rate_variance = FALSE,
    natural_rate_volatility = "submitted_total",
    transitory_lags = "ma2") {
  if (policy_inflation != "core") {
    stop("Exact iid inflation-variance matching requires a core-inflation policy rule.")
  }
  # Under the core-inflation policy rule and a1=0, a0 affects the headline
  # inflation loading but not physical dynamics.
  x0[21] <- 0
  x0[22] <- 0
  base <- var_coeffs_jfec(
    x0 = x0, taylor_rule = taylor_rule,
    natural_rate_identity = natural_rate_identity,
    policy_inflation = policy_inflation,
    natural_rate_volatility = natural_rate_volatility,
    transitory_lags = transitory_lags
  )
  K <- nrow(base$Phi)
  inflation_state <- 9L
  radius <- max(Mod(eigen(base$Phi, only.values = TRUE)$values))
  if (!is.finite(radius) || radius >= 1 - 1e-8) {
    stop("The inflation-variance restriction requires strict stationarity.")
  }
  target_inflation_variance <- var(macro[, 2], na.rm = TRUE)
  mean_state <- solve(diag(K) - base$Phi, base$Mu)

  evaluate_volatilities <- function(sigma_i, sigma_kappa) {
    Sigma_til <- base$Sigma_til
    Sigma_til[1, 1] <- sigma_i
    Sigma_til[3, 8] <- sigma_kappa
    Sigma <- solve(base$A, Sigma_til)
    covariance_state <- matrix(
      solve(diag(K * K) - kronecker(base$Phi, base$Phi),
            as.vector(Sigma %*% t(Sigma))), K, K
    )
    core_variance <- covariance_state[inflation_state, inflation_state]
    transitory_variance <- target_inflation_variance - core_variance
    if (!is.finite(transitory_variance) || transitory_variance <= 0) return(NULL)
    sigma_transitory <- sqrt(transitory_variance)
    inflation_loading <- rep(0, K)
    inflation_loading[inflation_state] <- 1
    inflation_loading[10] <- sigma_transitory
    list(
      sigma_i = sigma_i, sigma_kappa = sigma_kappa,
      sigma_transitory = sigma_transitory,
      covariance_state = covariance_state,
      inflation_loading = inflation_loading
    )
  }

  if (match_real_rate_variance) {
    if (is.null(observed_short_rate) || length(observed_short_rate) != nrow(macro)) {
      stop("An observed policy-rate series aligned with macro is required for real-rate variance matching.")
    }
    observed_real <- as.numeric(observed_short_rate) - macro[, 2]
    target_real_variance <- var(observed_real, na.rm = TRUE)
    if (natural_rate_volatility != "free") {
      stop("Real-rate variance matching requires free natural-rate volatility.")
    }
    residual <- function(log_sigma_kappa) {
      evaluated <- evaluate_volatilities(x0[4], exp(log_sigma_kappa))
      if (is.null(evaluated)) return(NA_real_)
      model_variance <- effective_real_variance_gaussian(
        mean_state, evaluated$covariance_state,
        evaluated$inflation_loading, base$pars$r_lb
      )
      if (!is.finite(model_variance) || model_variance <= 0) return(NA_real_)
      log(model_variance / target_real_variance)
    }
    grid <- seq(log(1e-9), log(0.02), length.out = 120L)
    values <- vapply(grid, residual, numeric(1))
    crossings <- which(is.finite(values[-length(values)]) & is.finite(values[-1]) &
                         values[-length(values)] * values[-1] <= 0)
    if (!length(crossings)) {
      finite_values <- values[is.finite(values)]
      if (!length(finite_values)) stop(sprintf(
        paste0("No positive natural-rate-shock volatility matches the real-rate variance: ",
               "target %.8g; no admissible stationary point was found on the volatility grid."),
        target_real_variance
      ))
      attainable <- target_real_variance * exp(finite_values)
      stop(sprintf(
        paste0("No positive natural-rate-shock volatility matches the real-rate variance: ",
               "target %.8g, attainable grid range [%.8g, %.8g]."),
        target_real_variance, min(attainable), max(attainable)
      ))
    }
    closest <- crossings[which.min(abs((grid[crossings] + grid[crossings + 1L]) / 2 -
                                         log(max(x0[16], 1e-9))))]
    root <- uniroot(residual, c(grid[closest], grid[closest + 1L]), tol = 1e-10)$root
    evaluated <- evaluate_volatilities(x0[4], exp(root))
  } else {
    # Retain the submitted 5% forecast-error-share calibration when the new
    # real-rate variance restriction is not requested.
    h_fwd <- 12L
    Phi_j <- array(0, dim = c(K, K, h_fwd - 1L))
    Phi_j[, , 1] <- base$Phi
    for (h in 2:(h_fwd - 1L)) Phi_j[, , h] <- Phi_j[, , h - 1L] %*% base$Phi
    Sigma_without_policy <- solve(
      base$A, base$Sigma_til[, 2:ncol(base$Sigma_til), drop = FALSE]
    )
    Q_without_policy <- Sigma_without_policy %*% t(Sigma_without_policy)
    propagated_variance <- 0
    for (h in seq_len(h_fwd - 1L)) {
      row_h <- matrix(Phi_j[inflation_state, , h], nrow = 1)
      propagated_variance <- propagated_variance +
        row_h %*% Q_without_policy %*% t(row_h)
    }
    target_share <- 0.05
    numerator <- target_share *
      (Q_without_policy[inflation_state, inflation_state] + propagated_variance)
    denominator <- (1 - target_share) * sum(Phi_j[inflation_state, 1, ]^2)
    sigma_i <- as.numeric(sqrt(numerator / denominator))
    if (!is.finite(sigma_i) || sigma_i <= 0) {
      stop("The policy-shock variance restriction has no positive solution.")
    }
    evaluated <- evaluate_volatilities(sigma_i, base$pars$sigma_kappa)
    if (is.null(evaluated)) {
      stop("Core inflation variance exceeds the sample headline-inflation variance.")
    }
  }
  x0[4] <- evaluated$sigma_i
  x0[16] <- evaluated$sigma_kappa
  x0[21] <- evaluated$sigma_transitory
  x0
}

build_filter_objects_jfec <- function(macro, yields_n, yields_r, surv_infexp,
                                      surv_gdpexp, surv_tbexp, hstep_s,
                                      hstep_g, hstep_t, mats_n, mats_r,
                                      dates, x0 = jfec_x0(),
                                      observed_short_rate = NULL,
                                      growth_survey_sigma = NULL,
                                      taylor_rule = c("submitted", "orphanides"),
                                      natural_rate_identity = c("submitted", "exact"),
                                      policy_inflation = c("headline", "core"),
                                      calibrate_growth_mean = TRUE,
                                      calibrate_gdp_measurement = TRUE,
                                      calibrate_growth_survey_measurement = TRUE,
                                      match_inflation_variance = FALSE,
                                      match_real_rate_variance = FALSE,
                                      natural_rate_volatility = c("submitted_total", "free"),
                                      transitory_lags = c("ma2", "iid"),
                                      core_inflation = NULL,
                                      core_inflation_sigma = NULL) {
  taylor_rule <- match.arg(taylor_rule)
  natural_rate_identity <- match.arg(natural_rate_identity)
  policy_inflation <- match.arg(policy_inflation)
  natural_rate_volatility <- match.arg(natural_rate_volatility)
  transitory_lags <- match.arg(transitory_lags)
  x0 <- apply_jfec_calibrations(
    x0, macro, calibrate_growth_mean, calibrate_gdp_measurement,
    calibrate_growth_survey_measurement
  )
  if (match_inflation_variance) {
    x0 <- apply_inflation_variance_restriction(
      x0, macro, yields_n, mats_n, observed_short_rate,
      taylor_rule, natural_rate_identity,
      policy_inflation, match_real_rate_variance, natural_rate_volatility,
      transitory_lags
    )
  }
  base <- var_coeffs_jfec(x0 = x0, taylor_rule = taylor_rule,
                          natural_rate_identity = natural_rate_identity,
                          policy_inflation = policy_inflation,
                          natural_rate_volatility = natural_rate_volatility,
                          transitory_lags = transitory_lags)
  if (!is.null(growth_survey_sigma)) {
    if (length(growth_survey_sigma) != 1L ||
        !is.finite(growth_survey_sigma) || growth_survey_sigma <= 0) {
      stop("growth_survey_sigma must be one positive finite scalar.")
    }
    base$pars$sigma_y <- growth_survey_sigma
  }
  A <- base$A
  Mu_til <- base$Mu_til
  Phi_til <- base$Phi_til
  Sigma_til <- base$Sigma_til
  Lambda_0 <- base$Lambda_0
  Lambda_1 <- base$Lambda_1
  pars <- base$pars
  if (any(!is.finite(c(A, Mu_til, Phi_til, Sigma_til)))) {
    stop("Non-finite state-space object.")
  }
  physical_radius <- max(Mod(eigen(base$Phi, only.values = TRUE)$values))
  if (!is.finite(physical_radius) || physical_radius >= 1 - 1e-8) {
    stop("Non-stationary physical dynamics.")
  }
  K <- nrow(A)
  T <- nrow(macro)
  if (!is.null(core_inflation)) {
    if (length(core_inflation) != T) stop("core_inflation must have T observations.")
    if (length(core_inflation_sigma) != 1L ||
        !is.finite(core_inflation_sigma) || core_inflation_sigma <= 0) {
      stop("core_inflation_sigma must be one positive finite scalar.")
    }
  }
  N_m <- ncol(macro)
  N_s <- ncol(surv_infexp)
  N_g <- ncol(surv_gdpexp)
  N_t <- ncol(surv_tbexp)
  N_n <- ncol(yields_n)
  N_r <- ncol(yields_r)
  freq <- 1
  maxmat <- max(c(mats_n, mats_r)) * freq

  Mu <- solve(A, Mu_til)
  Phi <- solve(A, Phi_til)
  I_Phi_inv <- solve(diag(K) - Phi)

  hstep_max <- max(c(hstep_s, hstep_g, hstep_t))
  Phi_j <- array(0, dim = c(K, K, hstep_max))
  Phi_j[, , 1] <- Phi
  for (j in 2:hstep_max) {
    Phi_j[, , j] <- Phi_j[, , j - 1] %*% Phi
  }

  if (match_real_rate_variance) {
    pars$sigma_i <- x0[4]
  } else {
    h_fwd <- 12
    theta_pi_i <- 0.05
    ii <- 9
    Sigma_X_ii <- solve(A, Sigma_til[, 2:ncol(Sigma_til), drop = FALSE])
    Sigma2_X_ii <- Sigma_X_ii %*% t(Sigma_X_ii)
    Phi_Sigma_ii_Phi <- 0
    for (h in 1:(h_fwd - 1)) {
      row_h <- matrix(Phi_j[ii, , h], nrow = 1)
      Phi_Sigma_ii_Phi <- Phi_Sigma_ii_Phi +
        row_h %*% Sigma2_X_ii %*% t(row_h)
    }
    numerator <- theta_pi_i * (Sigma2_X_ii[ii, ii] + Phi_Sigma_ii_Phi)
    denominator <- (1 - theta_pi_i) *
      sum(Phi_j[ii, 1, 1:(h_fwd - 1)]^2)
    pars$sigma_i <- as.numeric(sqrt(numerator / denominator))
    if (!is.finite(pars$sigma_i) || pars$sigma_i <= 0) {
      stop("The inflation forecast-error restriction has no positive solution.")
    }
  }
  Sigma_til[1, 1] <- pars$sigma_i
  Sigma_X <- solve(A, Sigma_til)
  Sigma2_X <- Sigma_X %*% t(Sigma_X)

  Mu_q <- Mu + Sigma_X %*% Lambda_0
  Phi_q <- Phi + Sigma_X %*% Lambda_1

  delta_1 <- c(1, rep(0, K - 1))
  delta_pi_1 <- c(rep(0, 8), 1, pars$a_u, 0)
  coefs_p <- affine_coefs_r(Phi, Mu, Sigma2_X, 0, delta_1, 0, delta_pi_1,
                            maxmat = maxmat)
  coefs_q <- affine_coefs_r(Phi_q, Mu_q, Sigma2_X, 0, delta_1, 0, delta_pi_1,
                            maxmat = maxmat)

  s2_n <- cumsum(diag(t(coefs_q$B_X_for) %*% Sigma2_X %*% coefs_q$B_X_for))
  s_n <- sqrt(s2_n)

  Gamma_m <- matrix(0, N_m, K)
  Gamma_m[1, 4:6] <- c(1, 1, -1)
  Gamma_m[2, 9:(9 + length(pars$a_u))] <- c(1, pars$a_u)
  Gamma_m[3, 5] <- 1
  Gamma_m[4, 7] <- 1
  Gamma_m0 <- matrix(0, N_m, 1)
  Gamma_core <- matrix(0, if (is.null(core_inflation)) 0L else 1L, K)
  Gamma_core0 <- matrix(0, nrow(Gamma_core), 1)
  if (nrow(Gamma_core)) Gamma_core[1, 9] <- 1

  Gamma_s <- matrix(0, N_s, K)
  Gamma_s0 <- matrix(0, N_s, 1)
  for (j in seq_len(N_s)) {
    Gamma_s[j, ] <- rowSums(coefs_p$B_X_pi[, 1:hstep_s[j], drop = FALSE]) /
      hstep_s[j]
    Gamma_s0[j, 1] <- sum(coefs_p$A_X_exp_pi[1:hstep_s[j]] / hstep_s[j] +
                             0.5 * coefs_p$Conv_pi[1:hstep_s[j]] /
                             hstep_s[j]^2)
  }

  Gamma_g <- matrix(NaN, N_g, K)
  Gamma_g0 <- matrix(0, N_g, 1)
  for (j in seq_len(N_g)) {
    Gamma_g[j, ] <- Gamma_m[1, ] %*% Phi %*% I_Phi_inv %*%
      (diag(K) - Phi_j[, , hstep_g[j]]) / hstep_g[j]
    Gamma_g0[j, 1] <- Gamma_m[1, ] %*% I_Phi_inv %*%
      (hstep_g[j] * diag(K) - Phi %*% I_Phi_inv %*%
         (diag(K) - Phi_j[, , hstep_g[j]])) %*% Mu / hstep_g[j]
  }

  ind_r <- dates < as.Date("2004-01-01") |
    (dates >= as.Date("2008-01-01") & dates < as.Date("2010-01-01")) |
    (dates >= as.Date("2020-03-01") & dates < as.Date("2021-03-01"))
  sigma_r_t <- rep(pars$sigma_r, T)
  sigma_r_t[ind_r] <- pars$sigma_r + pars$sigma_r_liq

  list(
    A = A, Mu = Mu, Phi = Phi, Sigma2_X = Sigma2_X, Phi_j = Phi_j,
    Sigma_X = Sigma_X, MuQ = Mu_q, PhiQ = Phi_q,
    coefs_p = coefs_p, coefs_q = coefs_q, s_n = s_n, pars = pars,
    Gamma_m = Gamma_m, Gamma_m0 = Gamma_m0, Gamma_s = Gamma_s,
    Gamma_core = Gamma_core, Gamma_core0 = Gamma_core0,
    core_inflation_sigma = core_inflation_sigma,
    Gamma_s0 = Gamma_s0, Gamma_g = Gamma_g, Gamma_g0 = Gamma_g0,
    sigma_r_t = sigma_r_t, x0 = x0,
    dims = list(T = T, K = K, N_m = N_m, N_s = N_s,
    N_g = N_g, N_t = N_t, N_n = N_n, N_r = N_r)
  )
}

run_kf_jfec <- function(macro, yields_n, yields_r, surv_infexp, surv_gdpexp,
                        surv_tbexp, hstep_s, hstep_g, hstep_t, mats_n,
                        mats_r, dates, x0 = jfec_x0(),
                        observed_short_rate = NULL,
                        growth_survey_sigma = NULL,
                        taylor_rule = c("submitted", "orphanides"),
                        natural_rate_identity = c("submitted", "exact"),
                        policy_inflation = c("headline", "core"),
                        calibrate_growth_mean = TRUE,
                        calibrate_gdp_measurement = TRUE,
                        calibrate_growth_survey_measurement = TRUE,
                        match_inflation_variance = FALSE,
                        match_real_rate_variance = FALSE,
                        natural_rate_volatility = c("submitted_total", "free"),
                        transitory_lags = c("ma2", "iid"),
                        stabilize = FALSE,
                        core_inflation = NULL,
                        core_inflation_sigma = NULL) {
  taylor_rule <- match.arg(taylor_rule)
  natural_rate_identity <- match.arg(natural_rate_identity)
  policy_inflation <- match.arg(policy_inflation)
  natural_rate_volatility <- match.arg(natural_rate_volatility)
  transitory_lags <- match.arg(transitory_lags)
  obj <- build_filter_objects_jfec(macro, yields_n, yields_r, surv_infexp,
                                   surv_gdpexp, surv_tbexp, hstep_s, hstep_g,
                                   hstep_t, mats_n, mats_r, dates, x0,
                                   observed_short_rate, growth_survey_sigma,
                                   taylor_rule,
                                   natural_rate_identity,
                                   policy_inflation,
                                   calibrate_growth_mean,
                                   calibrate_gdp_measurement,
                                   calibrate_growth_survey_measurement,
                                   match_inflation_variance,
                                   match_real_rate_variance,
                                   natural_rate_volatility,
                                   transitory_lags,
                                   core_inflation,
                                   core_inflation_sigma)
  T <- obj$dims$T
  K <- obj$dims$K
  N_t <- obj$dims$N_t
  N_n <- obj$dims$N_n
  N_r <- obj$dims$N_r
  N_g <- obj$dims$N_g

  core_observations <- if (is.null(core_inflation)) {
    matrix(numeric(0), T, 0L)
  } else matrix(core_inflation, ncol = 1L)
  y_obs_all <- cbind(macro, core_observations, surv_infexp, surv_gdpexp,
                     surv_tbexp, yields_n, yields_r)
  R_m <- c(obj$pars$sigma_gdp^2, obj$pars$sigma_inf^2, obj$pars$sigma_o^2,
           obj$pars$sigma_ptr^2)
  R_core <- if (is.null(core_inflation)) numeric(0) else core_inflation_sigma^2
  R_s <- rep(obj$pars$sigma_s^2, obj$dims$N_s)
  R_g <- rep(obj$pars$sigma_y^2, N_g)
  R_t <- rep(obj$pars$sigma_t^2, N_t)
  R_n <- rep(obj$pars$sigma_n^2, N_n)

  Q <- obj$Sigma2_X
  if (!stabilize && exists("compute_uncond_Var", mode = "function")) {
    P0 <- compute_uncond_Var(obj$Phi, obj$Sigma_X)
  } else {
    P0 <- matrix(solve(diag(K * K) - kronecker(obj$Phi, obj$Phi),
                       as.vector(Q)), K, K)
    P0 <- (P0 + t(P0)) / 2
  }
  P_pred <- obj$Phi %*% P0 %*% t(obj$Phi) + Q
  P_pred <- (P_pred + t(P_pred)) / 2
  x_pred <- matrix(0, K, T)
  x_pred[, 1] <- solve(diag(K) - obj$Phi, obj$Mu)
  x_pred[, 1] <- obj$Mu + obj$Phi %*% x_pred[, 1]

  x_upd <- matrix(0, T, K)
  x_std <- matrix(0, T, K)
  P_pred_all <- array(0, dim = c(K, K, T))
  P_upd_all <- array(0, dim = c(K, K, T))
  fval_t <- rep(0, T)
  Z_all <- diag(ncol(y_obs_all))

  for (tt in seq_len(T)) {
    if (tt > 1) {
      x_pred[, tt] <- obj$Mu + obj$Phi %*% x_upd[tt - 1, ]
      P_pred <- obj$Phi %*% P_upd %*% t(obj$Phi) + Q
      P_pred <- (P_pred + t(P_pred)) / 2
    }
    predicted_tolerance <- 1e-10 * max(1, max(abs(diag(P_pred))))
    if (stabilize && inherits(try(chol(P_pred + diag(predicted_tolerance, K)), silent = TRUE), "try-error")) {
      stop("Non-positive-definite predicted state covariance.")
    }
    P_pred_all[, , tt] <- P_pred

    yf <- y_fitting_r(
      x_pred[, tt], obj$coefs_q$A_X_for, obj$coefs_q$B_X_for,
      obj$coefs_q$A_X_exp, obj$pars$r_lb, obj$s_n,
      obj$coefs_q$A_X_for_pi, obj$coefs_q$B_X_pi, obj$Sigma2_X,
      obj$coefs_q$B_X_cum, obj$coefs_q$B_X_cum_pi
    )
    y_pred_n <- yf$yfit_all_n[1, mats_n]
    y_pred_r <- yf$yfit_all_r[1, mats_r]
    Gamma_n <- yf$JJ_n[mats_n, , drop = FALSE]
    Gamma_r <- yf$JJ_r[mats_r, , drop = FALSE]

    y_pred_t <- rep(NaN, N_t)
    Gamma_t <- matrix(0, N_t, K)
    if (any(!is.na(surv_tbexp[tt, ]))) {
      obs_h <- which(!is.na(surv_tbexp[tt, ]))
      max_h <- max(hstep_t[obs_h])
      yfit_tb_all <- rep(NaN, max_h)
      Gamma_tb_all <- matrix(NaN, max_h, K)
      x_pred_h <- x_pred[, tt]
      for (j in 1:(max_h - 1)) {
        x_pred_h <- obj$Mu + obj$Phi %*% x_pred_h
        if ((j - 5) %% 12 == 0) {
          yft <- y_fitting_r(
            x_pred_h, obj$coefs_q$A_X_for[1:3],
            obj$coefs_q$B_X_for[, 1:3], obj$coefs_q$A_X_exp[1:3],
            obj$pars$r_lb, obj$s_n[1:3],
            obj$coefs_q$A_X_for_pi[1:3], obj$coefs_q$B_X_pi[, 1:3],
            obj$Sigma2_X, obj$coefs_q$B_X_cum[, 1:3],
            obj$coefs_q$B_X_cum_pi[, 1:3]
          )
          yfit_tb_all[j + 1] <- yft$yfit_all_n[1, 3]
          Gamma_tb_all[j + 1, ] <- yft$JJ_n[3, ] %*% obj$Phi_j[, , j]
        }
      }
      for (j in obs_h) {
        y_pred_t[j] <- mean(yfit_tb_all[1:hstep_t[j]], na.rm = TRUE)
        Gamma_t[j, ] <- colMeans(Gamma_tb_all[1:hstep_t[j], , drop = FALSE],
                                 na.rm = TRUE)
      }
    }

    Gamma_full <- rbind(obj$Gamma_m, obj$Gamma_core, obj$Gamma_s, obj$Gamma_g, Gamma_t,
                        Gamma_n, Gamma_r)
    y_obs <- y_obs_all[tt, ]
    y_pred <- c(
      obj$Gamma_m0 + obj$Gamma_m %*% x_pred[, tt],
      obj$Gamma_core0 + obj$Gamma_core %*% x_pred[, tt],
      obj$Gamma_s0 + obj$Gamma_s %*% x_pred[, tt],
      obj$Gamma_g0 + obj$Gamma_g %*% x_pred[, tt],
      y_pred_t, y_pred_n, y_pred_r
    )
    ind_obs <- !is.na(y_obs) & colSums(abs(t(Gamma_full))) != 0
    Z_obs <- Z_all[ind_obs, , drop = FALSE]
    Gamma_obs <- Z_obs %*% Gamma_full
    R_r <- rep(obj$sigma_r_t[tt]^2, N_r)
    R <- diag(c(R_m, R_core, R_s, R_g, R_t, R_n, R_r))
    R_obs <- Z_obs %*% R %*% t(Z_obs)
    Delta <- Gamma_obs %*% P_pred %*% t(Gamma_obs) + R_obs
    Delta <- (Delta + t(Delta)) / 2
    chol_delta <- try(chol(Delta), silent = TRUE)
    if (any(!is.finite(Delta)) || inherits(chol_delta, "try-error")) {
      stop("Non-positive-definite forecast-error covariance.")
    }
    Delta_inv <- chol2inv(chol_delta)
    innov <- y_obs[ind_obs] - y_pred[ind_obs]
    log_det <- 2 * sum(log(diag(chol_delta)))
    loglik <- -sum(ind_obs) * 0.5 * log(2 * pi) -
      0.5 * log_det - 0.5 * t(innov) %*% Delta_inv %*% innov
    fval_t[tt] <- as.numeric(loglik)
    gain <- P_pred %*% t(Gamma_obs) %*% Delta_inv
    x_upd[tt, ] <- x_pred[, tt] + gain %*% innov
    if (stabilize) {
      update_map <- diag(K) - gain %*% Gamma_obs
      P_upd <- update_map %*% P_pred %*% t(update_map) + gain %*% R_obs %*% t(gain)
    } else {
      P_upd <- P_pred - gain %*% Gamma_obs %*% P_pred
    }
    P_upd <- (P_upd + t(P_upd)) / 2
    P_upd_all[, , tt] <- P_upd
    updated_tolerance <- 1e-10 * max(1, max(abs(diag(P_upd))))
    if (stabilize && inherits(try(chol(P_upd + diag(updated_tolerance, K)), silent = TRUE), "try-error")) {
      stop("Non-positive-definite updated state covariance.")
    }
    if (abs(fval_t[tt]) > 1e12) {
      stop("Numerically implausible period likelihood contribution.")
    }
    x_std[tt, ] <- sqrt(diag(P_upd))
  }

  list(fval = sum(fval_t), fval_t = fval_t, x_pred = t(x_pred),
       x_upd = x_upd, x_std = x_std, P_pred = P_pred_all,
       P_upd = P_upd_all, objects = obj)
}

smooth_kf_jfec <- function(kf_res) {
  if (is.null(kf_res$P_pred) || is.null(kf_res$P_upd)) {
    stop("The filter result does not contain covariance histories.")
  }
  x_pred <- kf_res$x_pred
  x_upd <- kf_res$x_upd
  Phi <- kf_res$objects$Phi
  T <- nrow(x_upd)
  K <- ncol(x_upd)
  x_smooth <- x_upd
  P_smooth <- kf_res$P_upd
  symmetric_pinv <- function(M) {
    decomposition <- eigen((M + t(M)) / 2, symmetric = TRUE)
    tolerance <- max(dim(M)) * max(abs(decomposition$values)) *
      .Machine$double.eps
    inverse_values <- ifelse(
      decomposition$values > tolerance, 1 / decomposition$values, 0
    )
    decomposition$vectors %*% (inverse_values * t(decomposition$vectors))
  }
  if (T >= 2L) {
    for (tt in (T - 1L):1L) {
      gain <- kf_res$P_upd[, , tt] %*% t(Phi) %*%
        symmetric_pinv(kf_res$P_pred[, , tt + 1L])
      x_smooth[tt, ] <- x_upd[tt, ] + gain %*%
        (x_smooth[tt + 1L, ] - x_pred[tt + 1L, ])
      P_smooth[, , tt] <- kf_res$P_upd[, , tt] + gain %*%
        (P_smooth[, , tt + 1L] - kf_res$P_pred[, , tt + 1L]) %*%
        t(gain)
      P_smooth[, , tt] <- (P_smooth[, , tt] + t(P_smooth[, , tt])) / 2
    }
  }
  list(
    x_smooth = x_smooth,
    x_std_smooth = t(vapply(
      seq_len(T), function(tt) sqrt(pmax(diag(P_smooth[, , tt]), 0)),
      numeric(K)
    )),
    P_smooth = P_smooth
  )
}

postfit_jfec <- function(kf_res, surv_tbexp, hstep_t, mats_n, mats_r) {
  obj <- kf_res$objects
  x_upd <- kf_res$x_upd
  T <- nrow(x_upd)
  K <- ncol(x_upd)
  N_t <- ncol(surv_tbexp)

  macro_fit <- sweep(x_upd %*% t(obj$Gamma_m), 2, obj$Gamma_m0, "+")
  surv_infexp_fit <- sweep(x_upd %*% t(obj$Gamma_s), 2, obj$Gamma_s0, "+")
  surv_gdpexp_fit <- sweep(x_upd %*% t(obj$Gamma_g), 2, obj$Gamma_g0, "+")

  surv_tbexp_fit <- matrix(NaN, T, N_t)
  for (tt in seq_len(T)) {
    if (any(!is.na(surv_tbexp[tt, ]))) {
      obs_h <- which(!is.na(surv_tbexp[tt, ]))
      max_h <- max(hstep_t[obs_h])
      yfit_tb_all <- if (tt == T) rep(0, max_h) else rep(NaN, max_h)
      x_pred_h <- matrix(x_upd[tt, ], ncol = 1)
      y_fitting_r(
        x_pred_h, obj$coefs_q$A_X_for[1:3],
        obj$coefs_q$B_X_for[, 1:3], obj$coefs_q$A_X_exp[1:3],
        obj$pars$r_lb, obj$s_n[1:3], obj$coefs_q$A_X_for_pi[1:3],
        obj$coefs_q$B_X_pi[, 1:3], obj$Sigma2_X,
        obj$coefs_q$B_X_cum[, 1:3], obj$coefs_q$B_X_cum_pi[, 1:3]
      )
      for (j in 1:(max_h - 1)) {
        x_pred_h <- obj$Mu + obj$Phi %*% x_pred_h
        if ((j - 5) %% 12 == 0) {
          yft <- y_fitting_r(
            x_pred_h, obj$coefs_q$A_X_for[1:3],
            obj$coefs_q$B_X_for[, 1:3], obj$coefs_q$A_X_exp[1:3],
            obj$pars$r_lb, obj$s_n[1:3], obj$coefs_q$A_X_for_pi[1:3],
            obj$coefs_q$B_X_pi[, 1:3], obj$Sigma2_X,
            obj$coefs_q$B_X_cum[, 1:3], obj$coefs_q$B_X_cum_pi[, 1:3]
          )
          yfit_tb_all[j + 1] <- yft$yfit_all_n[1, 3]
        }
      }
      for (j in obs_h) {
        surv_tbexp_fit[tt, j] <- mean(yfit_tb_all[1:hstep_t[j]], na.rm = TRUE)
      }
    }
  }

  fit_q <- y_fitting_r(
    t(x_upd), obj$coefs_q$A_X_for, obj$coefs_q$B_X_for,
    obj$coefs_q$A_X_exp, obj$pars$r_lb, obj$s_n,
    obj$coefs_q$A_X_for_pi, obj$coefs_q$B_X_pi, obj$Sigma2_X,
    obj$coefs_q$B_X_cum, obj$coefs_q$B_X_cum_pi
  )
  yfit_n <- fit_q$yfit_all_n[, mats_n, drop = FALSE]
  yfit_r <- fit_q$yfit_all_r[, mats_r, drop = FALSE]

  fit_p <- y_fitting_r(
    t(x_upd), obj$coefs_p$A_X_for, obj$coefs_p$B_X_for,
    obj$coefs_p$A_X_exp, obj$pars$r_lb,
    sqrt(cumsum(diag(t(obj$coefs_p$B_X_for) %*% obj$Sigma2_X %*% obj$coefs_p$B_X_for))),
    obj$coefs_p$A_X_for_pi, obj$coefs_p$B_X_pi, obj$Sigma2_X,
    obj$coefs_p$B_X_cum, obj$coefs_p$B_X_cum_pi
  )
  yfit_n_exp <- fit_p$yfit_all_n[, mats_n, drop = FALSE]
  yfit_r_exp <- fit_p$yfit_all_r[, mats_r, drop = FALSE]

  list(
    macro_fit = macro_fit,
    surv_infexp_fit = surv_infexp_fit,
    surv_gdpexp_fit = surv_gdpexp_fit,
    surv_tbexp_fit = surv_tbexp_fit,
    yfit_n = yfit_n,
    yfit_r = yfit_r,
    term_prem_n = yfit_n - yfit_n_exp,
    term_prem_r = yfit_r - yfit_r_exp
  )
}
