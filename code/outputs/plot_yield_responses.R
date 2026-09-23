# State-dependent yield-response figures for the six-state model.

source("code/model/iekf.R")

diagnostic_directory <- Sys.getenv(
  "GGNR_DIAGNOSTIC_DIR",
  "outputs/diagnostics/baseline"
)
output_directory <- Sys.getenv(
  "GGNR_FIGURE_DIR",
  "outputs/figures"
)
horizon <- as.integer(Sys.getenv("GGNR_IRF_HORIZON", "120"))
maturities <- c(24L, 120L)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

diagnostics <- readRDS(file.path(diagnostic_directory, "diagnostics.rds"))
obj <- diagnostics$fit$objects
K <- obj$dims$K
state_names <- reduced_hfi_no_output_states
shock_names <- reduced_hfi_no_output_shocks

unconditional_mean <- as.numeric(solve(diag(K) - obj$Phi, obj$Mu))
unconditional_covariance <- stationary_covariance_no_output(list(
  Phi = obj$Phi, Sigma2_X = obj$Sigma2_X
))

price_yields <- function(state) {
  fit <- y_fitting_r(
    state, obj$coefs_q$A_X_for, obj$coefs_q$B_X_for,
    obj$coefs_q$A_X_exp, obj$pars$r_lb, obj$s_n,
    obj$coefs_q$A_X_for_pi, obj$coefs_q$B_X_pi, obj$Sigma2_X,
    obj$coefs_q$B_X_cum, obj$coefs_q$B_X_cum_pi,
    use_cpp = FALSE
  )
  c(fit$yfit_all_n[1, maturities], fit$yfit_all_r[1, maturities])
}

# Second-order conditional-expectation approximation evaluated through
# covariance-aligned central differences.
expected_yields <- function(mean_state, covariance_state) {
  center <- price_yields(mean_state)
  eig <- eigen((covariance_state + t(covariance_state)) / 2, symmetric = TRUE)
  keep <- eig$values > max(eig$values) * 1e-10
  if (!any(keep)) return(center)
  directions <- sweep(eig$vectors[, keep, drop = FALSE], 2,
                      sqrt(eig$values[keep]), "*")
  correction <- numeric(length(center))
  for (j in seq_len(ncol(directions))) {
    correction <- correction + price_yields(mean_state + directions[, j]) +
      price_yields(mean_state - directions[, j]) - 2 * center
  }
  center + 0.5 * correction
}

conditional_yield_irf <- function(conditioning, values, shock_bp) {
  number_conditions <- nrow(conditioning)
  number_shocks <- ncol(obj$Sigma_X)
  cross_covariance <- cbind(
    obj$Sigma_X,
    obj$Phi %*% unconditional_covariance %*% t(conditioning)
  )
  conditioning_variance <- rbind(
    cbind(diag(number_shocks), matrix(0, number_shocks, number_conditions)),
    cbind(matrix(0, number_conditions, number_shocks),
          conditioning %*% unconditional_covariance %*% t(conditioning))
  )
  policy_impact <- sum(obj$loadings$shadow_rate *
                         obj$Sigma_X[, "eps_m"])
  shock <- rep(0, number_shocks)
  shock[match("eps_m", shock_names)] <-
    (shock_bp / 10000 / 12) / policy_impact
  centered_conditions <- values - as.numeric(conditioning %*% unconditional_mean)
  inverse_variance <- solve(conditioning_variance)
  response_with_shock <- cross_covariance %*% inverse_variance %*%
    c(shock, centered_conditions)
  response_without_shock <- cross_covariance %*% inverse_variance %*%
    c(rep(0, number_shocks), centered_conditions)
  removed_covariance <- cross_covariance %*% inverse_variance %*%
    t(cross_covariance)

  nominal <- matrix(NA_real_, length(maturities), horizon + 1L)
  real <- nominal
  for (h in 0:horizon) {
    conditional_covariance <- unconditional_covariance - removed_covariance
    conditional_covariance <- (conditional_covariance +
                                  t(conditional_covariance)) / 2
    yields_with_shock <- expected_yields(
      unconditional_mean + response_with_shock, conditional_covariance
    )
    yields_without_shock <- expected_yields(
      unconditional_mean + response_without_shock, conditional_covariance
    )
    difference <- yields_with_shock - yields_without_shock
    nominal[, h + 1L] <- difference[seq_along(maturities)]
    real[, h + 1L] <- difference[length(maturities) + seq_along(maturities)]
    response_with_shock <- obj$Phi %*% response_with_shock
    response_without_shock <- obj$Phi %*% response_without_shock
    removed_covariance <- obj$Phi %*% removed_covariance %*% t(obj$Phi)
  }
  list(nominal = nominal, real = real)
}

make_loading <- function(values) {
  out <- rep(0, K)
  names(out) <- state_names
  out[names(values)] <- values
  out
}
annual_shadow <- 12 * obj$loadings$shadow_rate
annual_r_pi <- 12 * make_loading(c(r_star = 1, pi_star = 1))
annual_r_star <- 12 * make_loading(c(r_star = 1))

experiments <- list(
  s0_rpi = list(
    shock_bp = 25,
    conditioning = rbind(shadow_rate = annual_shadow, r_plus_pi = annual_r_pi),
    values = list(c(0, 0.04), c(0, 0)),
    labels = c("s=0%, r*+pi*=4%", "s=0%, r*+pi*=0%"),
    title = "25 bp short-rate innovation at the lower bound"
  ),
  s25_rpi = list(
    shock_bp = -25,
    conditioning = rbind(shadow_rate = annual_shadow, r_plus_pi = annual_r_pi),
    values = list(c(0.0025, 0.04), c(0.0025, 0)),
    labels = c("s=0.25%, r*+pi*=4%", "s=0.25%, r*+pi*=0%"),
    title = "-25 bp short-rate innovation near the lower bound"
  ),
  s0_rstar = list(
    shock_bp = 25,
    conditioning = rbind(shadow_rate = annual_shadow, r_star = annual_r_star),
    values = list(c(0, 0.02), c(0, -0.02)),
    labels = c("s=0%, r*=2%", "s=0%, r*=-2%"),
    title = "25 bp short-rate innovation: alternative natural-rate states"
  ),
  shadow_rate = list(
    shock_bp = 25,
    conditioning = rbind(shadow_rate = annual_shadow),
    values = list(0.04, -0.02),
    labels = c("s=4%", "s=-2%"),
    title = "25 bp short-rate innovation: alternative shadow-rate states"
  )
)

draw_yield_irf <- function(results, experiment) {
  colors <- c("black", "grey45")
  line_types <- c(1, 2)
  all_values <- unlist(lapply(results, function(x) c(x$nominal, x$real))) * 1200
  limits <- range(c(0, all_values), finite = TRUE)
  padding <- 0.08 * diff(limits)
  limits <- limits + c(-padding, padding)
  par(mfrow = c(1, 2), mar = c(4.6, 5.0, 3.6, 1),
      cex.axis = 1.18, cex.lab = 1.22, cex.main = 1.18,
      font.main = 2, mgp = c(2.8, 0.85, 0), tcl = -0.35)
  for (market in c("nominal", "real")) {
    plot(0:horizon, rep(0, horizon + 1L), type = "n", ylim = limits,
         xlab = "Months after innovation", ylab = "Percentage points",
         main = paste(tools::toTitleCase(market), "yields"))
    abline(h = 0, col = "grey65", lwd = 1.2)
    grid(col = "grey88", lty = 1)
    legend_labels <- character()
    legend_colors <- integer()
    legend_types <- integer()
    for (scenario in seq_along(results)) {
      for (maturity in seq_along(maturities)) {
        lines(0:horizon, 1200 * results[[scenario]][[market]][maturity, ],
              col = colors[maturity], lty = line_types[scenario], lwd = 2.6)
        legend_labels <- c(
          legend_labels,
          paste0(maturities[maturity], "m, ", experiment$labels[scenario])
        )
        legend_colors <- c(legend_colors, maturity)
        legend_types <- c(legend_types, line_types[scenario])
      }
    }
    if (market == "real") {
      legend("topright", legend_labels, col = colors[legend_colors],
             lty = legend_types, lwd = 2.6, cex = 1.00,
             bg = adjustcolor("white", alpha.f = 0.9))
    }
  }
}

all_yield_results <- list()
for (name in names(experiments)) {
  experiment <- experiments[[name]]
  results <- lapply(experiment$values, function(values) {
    conditional_yield_irf(
      experiment$conditioning, values, experiment$shock_bp
    )
  })
  all_yield_results[[name]] <- results
  png(file.path(output_directory, paste0("irf_yields_", name, ".png")),
      width = 1800, height = 850, res = 180)
  draw_yield_irf(results, experiment)
  dev.off()
  pdf(file.path(output_directory, paste0("irf_yields_", name, ".pdf")),
      width = 10, height = 4.7, useDingbats = FALSE)
  draw_yield_irf(results, experiment)
  dev.off()
}

# Macro responses to a 25 bp innovation in the shadow rate.
policy_impulse <- obj$Sigma_X[, "eps_m"]
policy_impulse <- policy_impulse *
  ((25 / 10000 / 12) / sum(obj$loadings$shadow_rate * policy_impulse))
macro_loadings <- rbind(
  inflation = obj$loadings$inflation,
  shadow_rate = obj$loadings$shadow_rate,
  real_rate = obj$loadings$shadow_rate - obj$loadings$inflation
)
macro_irf <- matrix(NA_real_, horizon + 1L, nrow(macro_loadings),
                    dimnames = list(NULL, rownames(macro_loadings)))
response <- policy_impulse
for (h in 0:horizon) {
  macro_irf[h + 1L, ] <- 1200 * as.numeric(macro_loadings %*% response)
  response <- obj$Phi %*% response
}
write.csv(data.frame(horizon = 0:horizon, macro_irf),
          file.path(output_directory, "irf_macro_policy.csv"), row.names = FALSE)

draw_macro_irf <- function() {
  colors <- c("#D55E00", "#0072B2", "#009E73")
  titles <- c("Inflation", "Shadow short rate", "Ex-ante real short rate")
  par(mfrow = c(1, 3), mar = c(4, 3.8, 3.2, 0.8))
  for (j in seq_len(ncol(macro_irf))) {
    plot(0:horizon, macro_irf[, j], type = "l", lwd = 1.7,
         col = colors[j], xlab = "Months after innovation",
         ylab = "Percentage points", main = titles[j])
    abline(h = 0, col = "grey65", lwd = 0.8)
    grid(col = "grey88", lty = 1)
  }
}

png(file.path(output_directory, "irf_macro_policy.png"),
    width = 1800, height = 720, res = 180)
draw_macro_irf()
dev.off()
pdf(file.path(output_directory, "irf_macro_policy.pdf"),
    width = 10, height = 4, useDingbats = FALSE)
draw_macro_irf()
dev.off()

saveRDS(
  list(macro = macro_irf, yield_experiments = all_yield_results,
       experiments = experiments, maturities = maturities),
  file.path(output_directory, "irf_original_set.rds")
)
cat("Saved macro IRF and four conditional yield-IRF experiments.\n")
