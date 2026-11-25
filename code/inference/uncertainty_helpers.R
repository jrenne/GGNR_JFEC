# Nonlinear propagation of a full parameter covariance and conditional states.
ggnr_uncertainty_series <- function(fit, states = fit$x_upd) {
  state <- sweep(states, 2, c(1200,1200,1200,1200,1200,1), "*")
  cbind(state, shadow_rate = rowSums(state[,c("r_star","pi_star","m")]),
        calculate_premia(fit, states))
}

ggnr_uncertainty_draw <- function(input, problem, theta, root) {
  set.seed(input$state_seed)
  tryCatch({
    proposed <- theta + as.numeric(root %*% input$noise); names(proposed) <- names(theta)
    fit <- problem$evaluate_fit(proposed)
    if (!is.finite(fit$fval)) stop("Nonfinite filter output.")
    means <- ggnr_uncertainty_series(fit)
    state <- fit$x_upd
    minimum <- 0
    for (tt in seq_len(nrow(state))) {
      P <- fit$P_upd[,,tt]
      e <- eigen((P + t(P))/2, symmetric = TRUE)
      minimum <- min(minimum, min(e$values))
      if (min(e$values) < -1e-10) stop("Indefinite conditional state covariance.")
      # Remove only floating-point negative zeros in a conditional state covariance.
      state[tt,] <- state[tt,] + as.numeric(e$vectors %*%
        (sqrt(pmax(e$values,0)) * rnorm(ncol(state))))
    }
    paths <- ggnr_uncertainty_series(fit, state)
    if (any(!is.finite(paths))) stop("Nonfinite simulated output.")
    list(ok=TRUE, paths=paths, means=means, minimum_filter_eigenvalue=minimum)
  }, error=function(error) {
    message <- conditionMessage(error)
    domain_errors <- c("Non-stationary dynamics.",
      "The shadow-rate variance leaves no positive sigma_m squared.",
      "The inflation variance leaves no positive sigma_u squared.")
    if (!(message %in% domain_errors)) stop(error)
    list(ok=FALSE, reason=message)
  })
}

ggnr_uncertainty_quantiles <- function(paths, point, dates, columns) {
  bands <- data.frame(date=dates)
  for (j in columns) {
    q <- t(apply(paths[,j,],1,quantile,probs=c(.025,.16,.5,.84,.975),type=7))
    name <- colnames(point)[j]
    bands[[paste0(name,"_estimate")]] <- point[,j]
    for (k in 1:5) bands[[paste0(name,"_",c("p025","p16","median","p84","p975")[k])]] <- q[,k]
  }
  bands
}
