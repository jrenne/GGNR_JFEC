# Thirty free coordinates. CPI has no extra measurement disturbance; PTR SD=20bp.
# Four unconditional moment targets determine the two means and two scales.
source("code/model/iekf.R")
source("code/model/fast_interface.R")
make_parameterization <- function(parameters, targets, liquid_tips_only = FALSE) {
  template <- modifyList(reduced_hfi_no_output_defaults(), parameters)
  template$include_core_cpi <- FALSE
  template$sigma_inf <- 0; template$sigma_ptr <- .20/1200
  template$sigma_inflation_survey_long <- template$sigma_tbill_survey_long <- .10/1200
  template$sigma_real_yield <- .10/1200; template$sigma_real_yield_liquidity <- .20/1200
  template$liquidity_spread_covid <- 0
  if (liquid_tips_only) template$liquidity_spread_gfc <- 0
  shocks <- reduced_hfi_no_output_shocks
  free_measurement <- c("sigma_inflation_survey", "sigma_tbill_survey", "sigma_nominal_yield")
  from_working <- function(theta) {
    if (any(!is.finite(theta))) stop("Nonfinite estimation coordinate.")
    p <- template
    for (n in c("rho_r_star", "rho_pi_star", "rho_w")) p[[n]] <- unname(plogis(theta[n]))
    for (n in c("sigma_r_star", "sigma_pi_star", "sigma_pi_gap")) p[[n]] <- unname(exp(theta[n])/1200)
    p$cyclical_var <- matrix(theta[c("phi_m_m", "phi_pi_m", "phi_m_pi", "phi_pi_pi")],2,2,
      dimnames=list(c("m","pi_gap"),c("m","pi_gap")))
    p$impact_pi_m_ratio <- unname(theta["impact_pi_m_ratio"])
    if (!liquid_tips_only) p$liquidity_spread_gfc <- unname(exp(theta["liquidity_spread_gfc"])/1200)
    for (n in free_measurement) p[[n]] <- unname(exp(theta[paste0("measurement_",n)]))
    for (n in setdiff(shocks,"eps_u")) {
      p$lambda_0[n] <- theta[paste0("lambda0_",n)]
      p$lambda_w[n] <- theta[paste0("lambdaw_",n)]
    }
    p$lambda_r_star_scale <- unname(theta["lambda_r_star_scale"])
    p$lambda_pi_star_scale <- unname(theta["lambda_pi_star_scale"])
    p$sigma_w <- sqrt(1-p$rho_w^2);p$sigma_wpi <- unname(theta["sigma_wpi"])
    p <- apply_reduced_hfi_exact_moments(p,targets)$p
    p$lambda_0["eps_u"] <- (theta["kappa0_u_pp"]/1200-.5*p$sigma_u^2)/p$sigma_u
    p$lambda_w["eps_u"] <- theta["kappaw_u_pp"]/(1200*p$sigma_u)
    p
  }
  natural <- function(theta) {
    p <- from_working(theta)
    c(rho_r_star=unname(p$rho_r_star),sigma_r_star_pp=1200*unname(p$sigma_r_star),
      rho_pi_star=unname(p$rho_pi_star),sigma_pi_star_pp=1200*unname(p$sigma_pi_star),
      phi_m_m=p$cyclical_var["m","m"],phi_m_pi=p$cyclical_var["m","pi_gap"],
      phi_pi_m=p$cyclical_var["pi_gap","m"],phi_pi_pi=p$cyclical_var["pi_gap","pi_gap"],
      sigma_m_pp=1200*unname(p$sigma_m),sigma_pi_gap_pp=1200*unname(p$sigma_pi_gap),
      sigma_u_pp=1200*unname(p$sigma_u),impact_pi_m_ratio=unname(p$impact_pi_m_ratio),
      liquidity_spread_gfc_pp=1200*unname(p$liquidity_spread_gfc),
      sigma_inf_pp=0,sigma_ptr_pp=.20,
      setNames(1200*unlist(p[free_measurement],use.names=FALSE),paste0(free_measurement,"_pp")),
      setNames(as.numeric(p$lambda_0),paste0("lambda0_",shocks)),
      setNames(as.numeric(p$lambda_w),paste0("lambdaw_",shocks)),
      lambda_r_star_scale=unname(p$lambda_r_star_scale),lambda_pi_star_scale=unname(p$lambda_pi_star_scale),
      rho_w=unname(p$rho_w),sigma_wpi=unname(p$sigma_wpi))
  }
  list(from_working=from_working,natural=natural)
}
make_inference_problem <- function(estimate_file=Sys.getenv("GGNR_START_FILE","estimates/baseline_r2.rds"),
    data_override=NULL, short_override=NULL, liquid_tips_only=NULL) {
  estimate<-readRDS(estimate_file)
  data<-if(is.null(data_override))load_paper_data() else data_override
  short<-if(is.null(short_override))data$effective_federal_funds_rate else short_override
  if(is.null(liquid_tips_only))liquid_tips_only<-isTRUE(estimate$liquid_tips_only)
  targets<-reduced_hfi_exact_targets(data,short)
  mapping<-make_parameterization(estimate$parameters,targets,liquid_tips_only)
  theta<-estimate$working_parameters
  if(length(theta)!=30L && !isTRUE(liquid_tips_only))stop("This archived estimate does not use the current 30-coordinate parameterization.")
  if(liquid_tips_only)theta<-theta[names(theta)!="liquidity_spread_gfc"]
  evaluate_fit<-function(theta,branches=NULL,record=FALSE) evaluate_filter_fast(data,mapping$from_working(theta),
    filter_method=if(is.null(estimate$filter_method))"IEKF" else estimate$filter_method,
    trigger=if(is.null(estimate$iekf_trigger_probability))0 else estimate$iekf_trigger_probability,
    steps=if(is.null(estimate$iekf_steps))2L else estimate$iekf_steps,
    branches=branches,record=record,observed_short_rate=short)
  list(estimate=estimate,theta=theta,from_working=mapping$from_working,natural=mapping$natural,
    evaluate_fit=evaluate_fit,evaluate=function(theta,branches=NULL)evaluate_fit(theta,branches)$fval_t,
    targets=targets,estimate_file=estimate_file,data=data,short_rate=short)
}
