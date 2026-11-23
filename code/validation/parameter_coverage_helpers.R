# Two-factor Monte Carlo used in OA.3.2--OA.3.3. No bias correction.
source("code/inference/sandwich_helpers.R")

coverage_design <- function(replications = 500L, case = "all") {
  stopifnot(length(replications) == 1L, is.finite(replications),
            replications >= 1, replications == as.integer(replications))
  design <- rbind(
    data.frame(case = "reference_15", target = 0.15,
               rep = seq_len(replications)),
    data.frame(case = "higher_zlb_50", target = 0.50,
               rep = 50L + seq_len(replications))
  )
  if (!case %in% c("all", "reference_15", "higher_zlb_50")) stop("Unknown Monte Carlo design.")
  if (case != "all") design <- design[design$case == case, , drop = FALSE]
  design
}

coverage_key <- function(x) paste(x$case, x$rep, sep = ":")

run_parameter_replication <- function(binding, replication,
                                      evaluation_limit = 800L,
                                      max_iterations = 100L) {
  months <- 695L
  seed <- 20310100L + replication
  set.seed(seed)
  # Freeze these submitted-estimate measurement magnitudes. Reloading them from
  # a polished empirical estimate would silently change the Monte Carlo design.
  nominal_error <- 0.066301273257871
  inflation_error <- 2.18317492899878
  multipliers <- c(1, 0.10 / nominal_error, 1, 0.10 / nominal_error,
                   inflation_error / nominal_error)
  truth <- c(logit_rho_r = qlogis(0.97), logit_rho_pi = qlogis(0.95),
             log_sd_r = log(0.20), log_sd_pi = log(0.15),
             log_measurement = log(nominal_error))
  stationary_sd <- sqrt(0.20^2 / (1 - 0.97^2) + 0.15^2 / (1 - 0.95^2))
  means <- c(stationary_sd * qnorm(1 - binding) - 2, 2)
  rho <- plogis(truth[1:2])
  sigma <- exp(truth[3:4])
  states <- matrix(0, months, 2)
  states[1, ] <- rnorm(2, sd = sigma / sqrt(1 - rho^2))
  for (t in 2:months) states[t, ] <- rho * states[t - 1, ] + rnorm(2, sd = sigma)
  observations <- small_measure(states, truth, means) +
    sweep(matrix(rnorm(months * 5), months, 5), 2,
          exp(truth[5]) * multipliers, "*")
  initial <- truth + rnorm(5, sd = 0.10)
  evaluations <- 0L
  best <- Inf
  best_theta <- initial
  evaluate <- function(theta) small_filter(observations, theta, means,
                                          noise_multipliers = multipliers)
  objective <- function(theta) {
    evaluations <<- evaluations + 1L
    if (evaluations > evaluation_limit) stop("Evaluation limit reached.")
    if (any(!is.finite(theta)) || any(abs(theta) > 12)) return(1e12)
    values <- evaluate(theta)
    if (any(!is.finite(values))) return(1e12)
    value <- -sum(values)
    if (value < best) { best <<- value; best_theta <<- theta }
    value
  }
  started <- proc.time()["elapsed"]
  fit <- tryCatch(
    optim(initial, objective, method = "BFGS",
          control = list(maxit = max_iterations, reltol = 1e-9)),
    error = function(e) list(par = best_theta, value = best, convergence = 98L)
  )
  branches <- attr(evaluate(fit$par), "branches")
  local_evaluate <- function(theta) small_filter(observations, theta, means,
    fixed_branches = branches, noise_multipliers = multipliers)
  # The five-coordinate experiment uses the same central-difference formula as
  # the empirical exercise, with an absolute 2e-4 step in each log/logit coordinate.
  base <- local_evaluate(fit$par)
  H <- matrix(0, 5, 5)
  scores <- matrix(0, months, 5)
  step <- 2e-4
  for (i in 1:5) {
    plus <- minus <- fit$par
    plus[i] <- plus[i] + step; minus[i] <- minus[i] - step
    fp <- local_evaluate(plus); fm <- local_evaluate(minus)
    scores[, i] <- (fp - fm) / (2 * step)
    H[i, i] <- -sum(fp - 2 * base + fm) / step^2
    if (i > 1) for (j in seq_len(i - 1)) {
      pp <- pm <- mp <- mm <- fit$par
      pp[i] <- pm[i] <- fit$par[i] + step
      mp[i] <- mm[i] <- fit$par[i] - step
      pp[j] <- mp[j] <- fit$par[j] + step
      pm[j] <- mm[j] <- fit$par[j] - step
      H[i, j] <- H[j, i] <- -sum(local_evaluate(pp) - local_evaluate(pm) -
        local_evaluate(mp) + local_evaluate(mm)) / (4 * step^2)
    }
  }
  eigenvalues <- if (all(is.finite(H))) eigen(H, symmetric = TRUE)$values else NA_real_
  invertible <- all(is.finite(eigenvalues)) && min(abs(eigenvalues)) > 1e-7
  white_se <- hac_se <- rep(NA_real_, 5)
  if (invertible && all(is.finite(scores))) {
    inverse <- solve(H)
    white_variance <- diag(inverse %*% crossprod(scores) %*% t(inverse))
    hac_variance <- diag(inverse %*% bartlett_score_covariance(scores, 12) %*% t(inverse))
    if (any(white_variance < 0) || any(hac_variance < 0)) {
      invertible <- FALSE
    } else {
      white_se <- sqrt(white_variance); hac_se <- sqrt(hac_variance)
    }
  }
  valid <- invertible && fit$convergence == 0L && all(is.finite(hac_se))
  data.frame(phase = "confirmation", profile = "empirical",
    target_binding = binding, rep = replication, seed = seed, method = "IEKF",
    parameter = names(truth), true = truth, estimate = fit$par,
    se = white_se, hac_se = hac_se,
    coverage_white = abs(fit$par - truth) <= 1.96 * white_se,
    coverage_hac = abs(fit$par - truth) <= 1.96 * hac_se,
    length_hac = 2 * 1.96 * hac_se, convergence = fit$convergence,
    invertible = invertible, positive = all(eigenvalues > 1e-7), valid = valid,
    observed_binding = mean(rowSums(states) + sum(means) < 0),
    evaluations = evaluations, seconds = unname(proc.time()["elapsed"] - started))
}

summarize_parameter_coverage <- function(results) {
  groups <- split(results, list(results$case, results$parameter), drop = TRUE)
  do.call(rbind, lapply(groups, function(x) {
    valid <- x$valid
    coverage <- if (any(valid)) mean(x$coverage_hac[valid]) else NA_real_
    data.frame(case = x$case[1], parameter = x$parameter[1], replications = nrow(x),
      valid_intervals = sum(valid), optimizer_failures = sum(x$convergence != 0),
      noninvertible = sum(!x$invertible), positive_hessians = sum(x$positive),
      observed_binding = mean(x$observed_binding), bias = mean(x$estimate - x$true),
      empirical_sd = sd(x$estimate), rmse = sqrt(mean((x$estimate - x$true)^2)),
      mean_se_hac = mean(x$hac_se[valid]), coverage_white = mean(x$coverage_white[valid]),
      coverage_hac = coverage, mean_length_hac = mean(x$length_hac[valid]),
      mc_se_hac = sqrt(coverage * (1 - coverage) / sum(valid)),
      coverage_all_failures_miss = mean(valid & !is.na(x$coverage_hac) & x$coverage_hac),
      initial_optimizer_failures = sum(x$initial_optimizer_failure),
      optimizer_retries = sum(x$optimizer_retry))
  }))
}
