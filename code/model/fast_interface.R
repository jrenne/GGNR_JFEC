# Exact C++ implementation of the reference R filter, including both IEKF updates.
if (!requireNamespace("Rcpp", quietly = TRUE) || !requireNamespace("RcppArmadillo", quietly = TRUE))
  stop("Install Rcpp and RcppArmadillo before running the package.")
Rcpp::sourceCpp("code/model/fast_filter.cpp")
evaluate_filter_fast <- function(data, parameters, filter_method = "IEKF",
    trigger = 0, steps = 2L, branches = NULL, record = FALSE,
    observed_short_rate = data$effective_federal_funds_rate) {
  obj <- build_filter_objects_reduced_hfi_no_output(data, parameters,
    no_hfi_data(data$dates), observed_short_rate)
  Y <- cbind(data$macro[, 2], data$macro[, 4], data$surv_infexp,
    data$surv_tbexp, data$yields_n, data$yields_r)
  mu <- as.vector(solve(diag(6) - obj$Phi, obj$Mu))
  P0 <- stationary_covariance_no_output(obj)
  fit <- fast_iekf_cpp(obj, Y, data$mats_n, data$mats_r, data$hstep_t,
    mu, P0, filter_method == "IEKF", trigger, steps, branches, record)
  colnames(fit$x_pred) <- colnames(fit$x_upd) <- reduced_hfi_no_output_states
  dimnames(fit$P_pred) <- dimnames(fit$P_upd) <- list(reduced_hfi_no_output_states,
    reduced_hfi_no_output_states, NULL)
  fit$objects <- obj; fit$hfi_observation_count <- 0L
  fit$filter_method <- filter_method; fit$iteration_fraction <- mean(fit$iterated)
  fit
}
