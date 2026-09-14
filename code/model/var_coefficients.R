# State-space coefficient construction used by the JFEC model.
#
# Given the parameter vector, this file builds the
# A, Mu_til, Phi_til, Sigma_til, Lambda_0, Lambda_1, and derived objects.

x_bounds_vals <- function() {
  rbind(
    c(0, 0.995),      # rho_i
    c(1, 3),          # alpha_pi
    c(0, 0.1),        # alpha_z
    c(0, 1),          # sigma_i
    c(0, 0.995),      # rho_g
    c(0, 1),          # mu_g
    c(0, 0.00056),    # sigma_g
    c(0, 1),          # alpha
    c(0, 0.995),      # rho_z
    c(0, 1),          # sigma_z
    c(0, 0.995),      # rho_st
    c(-100, 100),     # mu_st
    c(0, 1),          # sigma_st
    c(0, 1),          # mu_kappa
    c(0, 0.995),      # rho_kappa
    c(-100, 100),     # sigma_kappa
    c(-100, 100),     # theta
    c(-0.5, 0.995),   # rho_bar
    c(0, 1),          # beta
    c(0, 0.005),      # sigma_pi
    c(0, 1),          # a0
    c(0, 0.005),      # a1
    c(-100, 100),     # r_lb
    c(0, 1),          # rho_w
    c(-100, 100),     # mu_w
    c(-100, 100),     # sigma_w
    c(-0.1, 0.1),     # sigma_wg
    c(-0.1, 0.1),     # sigma_wz
    c(-0.1, 0.1),     # sigma_wpi
    c(0, 0.05),       # sigma_o
    c(0, 0.05),       # sigma_gdp
    c(0, 0.05),       # sigma_inf
    c(0, 0.05),       # sigma_ptr
    c(0, 0.05),       # sigma_s
    c(0, 0.05),       # sigma_y
    c(0, 0.05),       # sigma_t
    c(0, 0.05),       # sigma_n
    c(0, 0.05),       # sigma_r
    c(0, 0.05),       # sigma_r_liq
    matrix(rep(c(-10, 10), 32), ncol = 2, byrow = TRUE)
  )
}

pars_trans_in <- function(x0, x_bounds = x_bounds_vals()) {
  log((x0 - x_bounds[, 1]) / (x_bounds[, 2] - x0))
}

pars_trans_out <- function(x1, x_bounds = x_bounds_vals()) {
  (x_bounds[, 2] - x_bounds[, 1]) / (1 + exp(-x1)) + x_bounds[, 1]
}

jfec_x0 <- function() {
  c(
    0.0104382545811124005, 1.5, 0.041666666666666699,
    0.000424503825792636974, 0.79041234946036798,
    0.00221100170331838984, 0.000507469265124473024,
    0.747219154808188013, 0.946493694514137007,
    0.00417761772895170001, 0.994999999999999996,
    0.00335487686007241012, 3.07978178192533006e-05,
    1.14215241622012997e-19, 0.994999999999999996,
    0.000236801488497918004, 1, -0.381562308802388006,
    0.00780528812560161035, 1.10528561709930992e-05,
    0.00245726191949591785, 0.00082486187531531543, 0,
    0.928721447221876995, 0, 0.370778199831249988,
    0.0933049551681942035, 0.000918623576709006969,
    0.00328229760763929001, 0.00499451540382658994,
    4.99451540382658975e-05, 4.99451540382658975e-05,
    0.00100000000000000002, 0.00020000000000000001,
    0.00100000000000000002, 5.00000000000000024e-05,
    0.000252040803131816022, 0.000225037669504208997,
    0.000781977814498796098, 0.0776018851782361019,
    0.21309572617054201, -1.08650492172258994,
    0.202307500330795004, 2.1269596752413702, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0.150750599551518011, 0.227559426649380003,
    -0.43899629971773102, 0.544441562310157945,
    3.28184089816334534, 0, 0, 0
  )
}

var_coeffs_jfec <- function(x0 = jfec_x0(), ind_sig_w = FALSE,
                            taylor_rule = c("submitted", "orphanides"),
                            natural_rate_identity = c("submitted", "exact"),
                            policy_inflation = c("headline", "core"),
                            natural_rate_volatility = c("submitted_total", "free"),
                            transitory_lags = c("ma2", "iid")) {
  taylor_rule <- match.arg(taylor_rule)
  natural_rate_identity <- match.arg(natural_rate_identity)
  policy_inflation <- match.arg(policy_inflation)
  natural_rate_volatility <- match.arg(natural_rate_volatility)
  transitory_lags <- match.arg(transitory_lags)
  v <- x0
  K <- if (transitory_lags == "iid") 11L else 13L

  pars <- list()
  pars$rho_i <- v[1]
  pars$alpha_pi <- v[2]
  pars$alpha_z <- v[3]
  pars$sigma_i <- v[4]
  pars$rho_g <- v[5]
  pars$mu_g <- v[6]
  pars$sigma_g <- v[7]
  pars$alpha <- v[8]
  pars$rho_z <- v[9]
  pars$sigma_z <- v[10]
  pars$rho_st <- v[11]
  pars$mu_st <- v[12]
  pars$sigma_st <- v[13]
  pars$mu_kappa <- v[14]
  pars$rho_kappa <- v[15]
  pars$sigma_kappa <- if (natural_rate_volatility == "submitted_total") {
    sqrt(0.00056^2 - pars$sigma_g^2)
  } else {
    v[16]
  }
  if (!is.finite(pars$sigma_kappa) || pars$sigma_kappa <= 0) {
    stop("sigma_kappa must be positive.")
  }
  v[16] <- pars$sigma_kappa
  pars$theta <- v[17]
  pars$rho_bar <- v[18]
  pars$beta <- v[19]
  pars$sigma_pi <- v[20]
  pars$a_u <- if (transitory_lags == "iid") v[21] else c(v[21], v[22], 0)
  pars$r_lb <- v[23]
  pars$rho_w <- v[24]
  pars$mu_w <- v[25]
  pars$sigma_w <- if (ind_sig_w) v[26] else sqrt(1 - pars$rho_w^2)
  v[26] <- pars$sigma_w
  pars$sigma_wg <- v[27]
  pars$sigma_wz <- v[28]
  pars$sigma_wpi <- v[29]
  pars$sigma_o <- v[30]
  pars$sigma_gdp <- v[31]
  pars$sigma_inf <- v[32]
  pars$sigma_ptr <- v[33]
  pars$sigma_s <- v[34]
  pars$sigma_y <- v[35]
  pars$sigma_t <- v[36]
  pars$sigma_n <- v[37]
  pars$sigma_r <- v[38]
  pars$sigma_r_liq <- v[39]
  pars$Lambda_0 <- v[40:47]
  pars$Lambda_z <- v[48:55]
  pars$Lambda_pi <- v[56:63]
  pars$Lambda_w <- v[64:71]

  state_names <- c(
    "s", "s_1", "kappa", "g", "z", "z_1", "pi_st", "r_st", "pi_bar",
    if (transitory_lags == "iid") "eps_u_0" else
      c("eps_u_0", "eps_u_1", "eps_u_2"),
    "w"
  )
  shock_names <- c("eps_i", "eps_g", "eps_z", "eps_st", "eps_pi", "eps_u",
                   "eps_w", "eps_kappa")

  A <- diag(K)
  Mu_til <- matrix(0, K, 1, dimnames = list(state_names, NULL))
  Phi_til <- matrix(0, K, K, dimnames = list(state_names, state_names))
  Sigma_til <- matrix(0, K, 8, dimnames = list(state_names, shock_names))
  transitory_policy_weight <- if (policy_inflation == "headline") 1 else 0

  if (taylor_rule == "submitted") {
    # s_t = (1-rho_i)s_{t-1} + rho_i[r*_t + pi*_t
    #       + alpha_pi(pi_t-pi*_t) + alpha_z z_t] + sigma_i eps_i,t
    A[1, 5] <- -pars$rho_i * pars$alpha_z
    A[1, 7] <- -pars$rho_i * (1 - pars$alpha_pi)
    A[1, 8] <- -pars$rho_i
    A[1, 9] <- -pars$rho_i * pars$alpha_pi
    A[1, 10:(9 + length(pars$a_u))] <- -transitory_policy_weight * pars$rho_i *
      pars$alpha_pi * pars$a_u
  } else {
    # R2/Orphanides alternative: rho_i is conventional rate inertia and the
    # responses to inflation and activity are not scaled by (1-rho_i).
    A[1, 5] <- -pars$alpha_z
    A[1, 7] <- pars$alpha_pi - (1 - pars$rho_i)
    A[1, 8] <- -(1 - pars$rho_i)
    A[1, 9] <- -pars$alpha_pi
    A[1, 10:(9 + length(pars$a_u))] <-
      -transitory_policy_weight * pars$alpha_pi * pars$a_u
  }
  A[8, 3] <- -1
  A[9, 7] <- -(1 - pars$rho_bar)
  dimnames(A) <- list(state_names, state_names)

  Mu_til[3, 1] <- (1 - pars$rho_kappa) * pars$mu_kappa
  Mu_til[4, 1] <- (1 - pars$rho_g) * pars$mu_g
  Mu_til[5, 1] <- pars$alpha * (1 - pars$rho_bar) *
    (1 - pars$rho_st) * pars$mu_st
  Mu_til[7, 1] <- (1 - pars$rho_st) * pars$mu_st
  # In the exact identity r*_t = kappa_t + theta*g_t, kappa_t already carries
  # its own unconditional mean through row 3 and A[8,3] = -1. The submitted
  # code adds mu_kappa once more in row 8; that is immaterial at its binding
  # estimate mu_kappa=0 but not when the restriction is relaxed.
  Mu_til[8, 1] <- pars$theta * (1 - pars$rho_g) * pars$mu_g
  if (natural_rate_identity == "submitted") {
    Mu_til[8, 1] <- Mu_til[8, 1] + pars$mu_kappa
  }
  Mu_til[K, 1] <- (1 - pars$rho_w) * pars$mu_w

  if (taylor_rule == "submitted") {
    Phi_til[1, 1] <- 1 - pars$rho_i
  } else {
    Phi_til[1, 1] <- pars$rho_i
  }
  Phi_til[2, 1] <- 1
  Phi_til[3, 3] <- pars$rho_kappa
  Phi_til[4, 4] <- pars$rho_g
  Phi_til[5, 1] <- -pars$alpha
  Phi_til[5, 5] <- pars$alpha * pars$beta + pars$rho_z
  Phi_til[5, 7] <- pars$alpha * (1 - pars$rho_bar) * pars$rho_st
  Phi_til[5, 8] <- pars$alpha
  Phi_til[5, 9] <- pars$alpha * pars$rho_bar
  Phi_til[5, 10:(9 + length(pars$a_u))] <-
    pars$alpha * c(pars$a_u[-1], 0)
  Phi_til[6, 5] <- 1
  Phi_til[7, 7] <- pars$rho_st
  Phi_til[8, 4] <- pars$theta * pars$rho_g
  Phi_til[9, 5] <- pars$beta
  Phi_til[9, 9] <- pars$rho_bar
  u_indices <- 10:(9 + length(pars$a_u))
  if (length(u_indices) > 1L) {
    Phi_til[u_indices[-1], u_indices[-length(u_indices)]] <-
      diag(length(u_indices) - 1L)
  }
  Phi_til[K, K] <- pars$rho_w

  Sigma_til[1, 1] <- pars$sigma_i
  Sigma_til[3, 8] <- pars$sigma_kappa
  Sigma_til[4, 2] <- pars$sigma_g
  Sigma_til[5, 3] <- pars$sigma_z
  Sigma_til[7, 4] <- pars$sigma_st
  Sigma_til[8, 2] <- pars$theta * pars$sigma_g
  Sigma_til[9, 5] <- pars$sigma_pi
  Sigma_til[10, 6] <- 1
  Sigma_til[K, 2] <- pars$sigma_wg
  Sigma_til[K, 3] <- pars$sigma_wz
  Sigma_til[K, 5] <- pars$sigma_wpi
  Sigma_til[K, 7] <- pars$sigma_w

  Lambda_1 <- matrix(0, 8, K)
  Lambda_1[, 5] <- pars$Lambda_z
  Lambda_1[, 9] <- pars$Lambda_pi
  Lambda_1[, K] <- pars$Lambda_w
  dimnames(Lambda_1) <- list(shock_names, state_names)
  Lambda_0 <- matrix(pars$Lambda_0, ncol = 1, dimnames = list(shock_names, NULL))

  Mu <- solve(A, Mu_til)
  Phi <- solve(A, Phi_til)
  Sigma <- solve(A, Sigma_til)
  MuQ <- Mu + Sigma %*% Lambda_0
  PhiQ <- Phi + Sigma %*% Lambda_1

  list(
    A = A, Mu_til = Mu_til, Phi_til = Phi_til, Sigma_til = Sigma_til,
    Lambda_0 = Lambda_0, Lambda_1 = Lambda_1, Mu = Mu, Phi = Phi,
    Sigma = Sigma, MuQ = MuQ, PhiQ = PhiQ, pars = pars, v = v,
    taylor_rule = taylor_rule,
    natural_rate_identity = natural_rate_identity,
    policy_inflation = policy_inflation,
    natural_rate_volatility = natural_rate_volatility,
    transitory_lags = transitory_lags
  )
}
