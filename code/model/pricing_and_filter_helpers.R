# Analytical affine coefficients and constrained nominal/real yield pricing.
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
