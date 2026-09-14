# Probability of the shadow short rate hitting its lower bound, conditional on
# the natural real rate and the inflation target. This is the six-state R-model
# calculation reported in the paper.

source("code/model/iekf.R")

diagnostic_directory <- Sys.getenv(
  "GGNR_DIAGNOSTIC_DIR",
  "outputs/diagnostics/baseline"
)
output_directory <- Sys.getenv(
  "GGNR_FIGURE_DIR",
  "outputs/figures"
)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

diagnostics <- readRDS(file.path(diagnostic_directory, "diagnostics.rds"))
obj <- diagnostics$fit$objects
K <- obj$dims$K

unconditional_mean <- as.numeric(solve(diag(K) - obj$Phi, obj$Mu))
unconditional_covariance <- stationary_covariance_no_output(list(
  Phi = obj$Phi, Sigma2_X = obj$Sigma2_X
))

conditioning <- matrix(
  0, 2, K,
  dimnames = list(c("r_star", "pi_star"), reduced_hfi_no_output_states)
)
conditioning["r_star", "r_star"] <- 1
conditioning["pi_star", "pi_star"] <- 1
shadow_loading <- as.numeric(obj$loadings$shadow_rate)

conditioning_covariance <- conditioning %*% unconditional_covariance %*%
  t(conditioning)
shadow_conditioning_covariance <- matrix(
  shadow_loading %*% unconditional_covariance %*% t(conditioning), nrow = 1
)
conditional_slope <- shadow_conditioning_covariance %*%
  solve(conditioning_covariance)
conditional_variance <- as.numeric(
  shadow_loading %*% unconditional_covariance %*% shadow_loading -
    shadow_conditioning_covariance %*% solve(conditioning_covariance) %*%
      t(shadow_conditioning_covariance)
)
if (!is.finite(conditional_variance) || conditional_variance <= 0) {
  stop("The conditional shadow-rate variance is not positive.")
}

r_star_grid_pp <- seq(-2, 6, by = 0.05)
pi_star_grid_pp <- seq(1, 6, by = 0.05)
grid <- expand.grid(
  r_star = r_star_grid_pp / 1200,
  pi_star = pi_star_grid_pp / 1200
)
conditioning_mean <- as.numeric(conditioning %*% unconditional_mean)
shadow_mean <- sum(shadow_loading * unconditional_mean)
conditional_shadow_mean <- shadow_mean + as.numeric(
  as.matrix(grid) %*% t(conditional_slope) -
    sum(conditioning_mean * as.numeric(conditional_slope))
)
probability <- pnorm(
  (obj$pars$r_lb - conditional_shadow_mean) / sqrt(conditional_variance)
)
probability_matrix <- matrix(
  probability,
  nrow = length(r_star_grid_pp), ncol = length(pi_star_grid_pp)
)

output <- transform(
  grid,
  r_star = 1200 * r_star,
  pi_star = 1200 * pi_star,
  lower_bound_probability = probability
)
write.csv(output, file.path(output_directory, "iso_probability.csv"),
          row.names = FALSE)

levels <- c(0.001, 0.005, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30)
draw_figure <- function() {
  contour(
    r_star_grid_pp, pi_star_grid_pp, probability_matrix,
    levels = levels, drawlabels = TRUE, labcex = 1.10,
    lwd = 2.3, col = "black",
    xlab = expression(paste("Real-rate trend, ", r*"*", " (%)")),
    ylab = expression(paste("Inflation trend, ", pi*"*", " (%)")),
    axes = FALSE
  )
  axis(1)
  axis(2, las = 1)
  box()
  grid(col = "grey88", lty = 1)
  contour(
    r_star_grid_pp, pi_star_grid_pp, probability_matrix,
    levels = levels, drawlabels = TRUE, labcex = 1.10,
    lwd = 2.3, col = "black", add = TRUE
  )
}

png(file.path(output_directory, "iso_probability.png"),
    width = 1450, height = 1050, res = 180)
par(mar = c(5.0, 5.6, 1.0, 1.2), cex.axis = 1.18, cex.lab = 1.22,
    mgp = c(3.0, 0.9, 0), tcl = -0.35)
draw_figure()
dev.off()

pdf(file.path(output_directory, "iso_probability.pdf"),
    width = 8.1, height = 5.9, useDingbats = FALSE)
par(mar = c(5.0, 5.6, 1.0, 1.2), cex.axis = 1.18, cex.lab = 1.22,
    mgp = c(3.0, 0.9, 0), tcl = -0.35)
draw_figure()
dev.off()

cat(sprintf(
  "Conditional shadow-rate standard deviation: %.3f annualized pp\n",
  1200 * sqrt(conditional_variance)
))
cat(sprintf(
  "Probability range on plotted grid: %.5f to %.5f\n",
  min(probability), max(probability)
))
