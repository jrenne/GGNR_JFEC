# Source from the project root. Check the Treasury-bill survey expectation
# against independent one-dimensional Gaussian quadrature and finite differences.
source("code/model/iekf.R")
estimate <- readRDS("estimates/baseline.rds")
data <- load_true_release_data()
short_rate <- load_observed_real_short_rate_inputs(data$dates)
hfi_none <- list(standardized = rep(NaN, length(data$dates)),
                 observed = rep(FALSE, length(data$dates)))
fit <- run_kf_reduced_hfi_no_output(data, estimate$parameters, hfi_none,
                                   short_rate, filter_method = "IEKF",
                                   iekf_steps = 2L)
obj <- fit$objects
stopifnot(identical(vapply(obj$tbill_terms, function(z) as.integer(z$horizon),
                         integer(1)), setNames(seq(5L, 113L, 12L),
                                               as.character(seq(5L, 113L, 12L)))))
quadrature_error <- derivative_error <- 0
for (row in unique(round(seq(1, nrow(fit$x_upd), length.out = 9)))) {
  x <- fit$x_upd[row, ]
  check <- tbill_measurement_no_output(x, obj, c(12, 120))
  independent <- vapply(obj$tbill_terms, function(term) {
    j <- term$horizon
    mu <- rep(0, length(x)); V <- matrix(0, length(x), length(x))
    for (h in seq_len(j)) {
      mu <- as.vector(obj$Mu + obj$Phi %*% mu)
      V <- obj$Phi %*% V %*% t(obj$Phi) + obj$Sigma2_X
    }
    mu <- mu + obj$Phi_j[, , j] %*% x
    mean(vapply(1:3, function(k) {
      b <- obj$coefs_q$B_X_for[, k]
      a <- obj$coefs_q$A_X_for[k] + sum(b * mu) - obj$pars$r_lb
      v <- sqrt(drop(t(b) %*% V %*% b))
      q <- if (k == 1) 0 else obj$s_n[k - 1]
      integrand <- function(z) {
        y <- a + v * z
        value <- if (q == 0) pmax(y, 0) else
          y * pnorm(y / q) + q * dnorm(y / q)
        value * dnorm(z)
      }
      obj$pars$r_lb + integrate(integrand, -Inf, Inf,
                                abs.tol = 1e-12, rel.tol = 1e-10)$value
    }, numeric(1)))
  }, numeric(1))
  quadrature_error <- max(quadrature_error,
                          abs(check$y_pred - c(independent[1], mean(independent))))
  for (k in seq_along(x)) {
    delta <- rep(0, length(x)); delta[k] <- 1e-7
    numerical <- (tbill_measurement_no_output(x + delta, obj, c(12, 120))$y_pred -
                    tbill_measurement_no_output(x - delta, obj, c(12, 120))$y_pred) / 2e-7
    derivative_error <- max(derivative_error, abs(numerical - check$Gamma[, k]))
  }
}
stopifnot(quadrature_error < 1e-9, derivative_error < 1e-7)
print(c(quadrature_error = quadrature_error, derivative_error = derivative_error))
