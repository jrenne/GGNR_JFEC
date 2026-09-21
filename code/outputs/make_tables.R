options(stringsAsFactors = FALSE)

out <- file.path("outputs", "tables")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

inference_directory <- Sys.getenv(
  "GGNR_INFERENCE_DIR",
  file.path("outputs", "diagnostics", "inference")
)
diagnostic_directory <- Sys.getenv(
  "GGNR_DIAGNOSTIC_DIR",
  file.path("outputs", "diagnostics", "baseline")
)
inf <- read.csv(file.path(inference_directory, "parameter_inference.csv"))
fit <- read.csv(file.path(diagnostic_directory, "fit_summary.csv"))
hac_covariance <- as.matrix(read.csv(
  file.path(inference_directory, "hac_sandwich_covariance_natural.csv"),
  row.names = 1, check.names = FALSE
))
inf$hac_sandwich_se <- sqrt(pmax(
  0, diag(hac_covariance[inf$parameter, inf$parameter, drop = FALSE])
))

labels <- c(
  rho_r_star="$\\rho_{r^*}$", sigma_r_star_pp="$\\sigma_{r^*}$ (pp)",
  rho_pi_star="$\\rho_{\\pi^*}$", sigma_pi_star_pp="$\\sigma_{\\pi^*}$ (pp)",
  phi_m_m="$\\Phi_{m,m}$", phi_m_pi="$\\Phi_{m,\\widetilde\\pi}$",
  phi_pi_m="$\\Phi_{\\widetilde\\pi,m}$", phi_pi_pi="$\\Phi_{\\widetilde\\pi,\\widetilde\\pi}$",
  sigma_m_pp="$\\sigma_m$ (pp)", sigma_pi_gap_pp="$\\sigma_{\\widetilde\\pi}$ (pp)",
  sigma_u_pp="$\\sigma_u$ (pp)", impact_pi_m_ratio="$b_{\\pi m}$",
  liquidity_spread_gfc_pp="$\\ell_{GFC}$ (pp)",
  sigma_inf_pp="$\\sigma_{\\pi,me}$ (pp)", sigma_ptr_pp="$\\sigma_{PTR,me}$ (pp)",
  sigma_inflation_survey_pp="$\\sigma_{CPI1,me}$ (pp)", sigma_tbill_survey_pp="$\\sigma_{TB1,me}$ (pp)",
  sigma_nominal_yield_pp="$\\sigma_{y^N,me}$ (pp)",
  lambda0_eps_r_star="$\\lambda_{0,r^*}$", lambda0_eps_pi_star="$\\lambda_{0,\\pi^*}$",
  lambda0_eps_m="$\\lambda_{0,m}$", lambda0_eps_pi_gap="$\\lambda_{0,\\widetilde\\pi}$",
  lambda0_eps_u="$\\lambda_{0,u}$", lambda0_eps_w="$\\lambda_{0,w}$",
  lambdaw_eps_r_star="$\\lambda_{w,r^*}$", lambdaw_eps_pi_star="$\\lambda_{w,\\pi^*}$",
  lambdaw_eps_m="$\\lambda_{w,m}$", lambdaw_eps_pi_gap="$\\lambda_{w,\\widetilde\\pi}$",
  lambdaw_eps_u="$\\lambda_{w,u}$", lambdaw_eps_w="$\\lambda_{w,w}$",
  lambda_r_star_scale="$\\beta_r$", lambda_pi_star_scale="$\\beta_\\pi$",
  rho_w.rho_w="$\\rho_w$", sigma_wpi="$\\sigma_{w\\pi}$"
)

groups <- list(
  "State dynamics"=1:12,
  "Liquidity adjustment"=13,
  "Measurement-error standard deviations"=14:18,
  "Prices of risk and risk state"=19:34
)

fmt <- function(x) {
  if (abs(x) < 0.0001 && x != 0) {
    exponent <- floor(log10(abs(x)))
    return(sprintf("$%.2f\\times10^{%d}$", x / 10^exponent, exponent))
  }
  formatC(x, format = "f", digits = 4)
}
significance <- function(estimate, se) {
  p <- 2 * pnorm(-abs(estimate / se))
  stars <- if (p < 0.01) "***" else if (p < 0.05) "**" else if (p < 0.10) "*" else ""
  if (nzchar(stars)) paste0("$^{", stars, "}$") else ""
}
lines <- c(
  "\\begin{table}[!htbp]", "\\centering",
  "\\caption{Parameter estimates and score-based standard errors}",
  "\\label{tab:parameters}", "\\begin{threeparttable}", "\\small",
  "\\renewcommand{\\arraystretch}{0.92}",
  "\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}lrrr}",
  "\\toprule", "Parameter & Estimate & OPG & OPG--HAC \\\\", "\\midrule"
)
for (g in names(groups)) {
  lines <- c(lines, paste0("\\multicolumn{4}{l}{\\textit{", g, "}} \\\\"))
  for (i in groups[[g]]) lines <- c(lines, paste0(
    labels[[inf$parameter[i]]], " & ", fmt(inf$estimate[i]), significance(inf$estimate[i], inf$hac_sandwich_se[i]), " & ",
    fmt(inf$opg_se[i]), " & ", fmt(inf$hac_sandwich_se[i]), " \\\\"
  ))
  if (g != tail(names(groups),1)) lines <- c(lines, "\\addlinespace")
}
uncertainty_note <- paste0(
  "The OPG column uses the inverse outer product of monthly likelihood scores. ",
  "The OPG--HAC column uses the score-based HAC sandwich documented in ",
  "the Online Appendix section ``State-space and measurement system.''"
)
lines <- c(
  lines, "\\bottomrule", "\\end{tabular*}",
  "\\begin{tablenotes}[flushleft]", "\\footnotesize",
  paste0(
    "\\item Notes: Innovation and measurement-error standard deviations are expressed in annualized percentage points where indicated. ",
    uncertainty_note,
    " $^*$, $^{**}$, and $^{***}$ denote nominal two-sided normal-approximation significance at 10\\%, 5\\%, and 1\\%, using OPG--HAC standard errors and a zero null. These markers are not unit-root tests; for scale parameters, a zero null is on the boundary and the markers are descriptive. The four unconditional moment restrictions determine $\\mu_r$, $\\mu_\\pi$, $\\sigma_m$, and $\\sigma_u$ analytically. The GFC liquidity adjustment is a common additive intercept in real yields during 2008--2009. Long-horizon CPI and Treasury-bill survey errors are fixed at 0.10 percentage point; real-yield errors are fixed at 0.10 percentage point in normal periods and 0.30 percentage point before 2004, during 2008--2009, and from March 2020 through February 2021."
  ),
  "\\end{tablenotes}", "\\end{threeparttable}", "\\end{table}"
)
writeLines(lines, file.path(out, "table_parameters.tex"), useBytes=TRUE)

keep <- fit$observable != "headline_cpi"
f <- fit[keep, ]
pretty <- c(headline_cpi_12m="Headline CPI inflation, 12 months", ptr="PTR long-run inflation", cpi_forecast_1y="CPI forecast, 1 year", cpi_forecast_10y="CPI forecast, 10 years", tbill_forecast_1y="Treasury-bill forecast, 1 year", tbill_forecast_10y="Treasury-bill forecast, 10 years")
for (m in c(3, 12, 24, 36, 60, 84, 120)) {
  maturity <- if (m == 3) "3 months" else paste(m / 12, if (m == 12) "year" else "years")
  pretty[paste0("nominal_yield_", m, "m")] <- paste("Nominal yield,", maturity)
}
for (m in c(24, 60, 84, 120)) pretty[paste0("real_yield_", m, "m")] <- paste("Real yield,", m / 12, "years")

fl <- c("\\begin{table}[!htbp]", "\\centering", "\\caption{Fit of macroeconomic expectations and yields}", "\\label{tab:fit}", "\\begin{threeparttable}", "\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}lrrrr}", "\\toprule", "Series & Observations & RMSE (pp) & MAE (pp) & Correlation \\\\", "\\midrule")
for(i in seq_len(nrow(f))) fl <- c(fl, paste0(pretty[[f$observable[i]]], " & ", f$observations[i], " & ", sprintf("%.3f", f$rmse_annual_pp[i]), " & ", sprintf("%.3f", f$mae_annual_pp[i]), " & ", sprintf("%.3f", f$correlation[i]), " \\\\"))
fl <- c(fl, "\\bottomrule", "\\end{tabular*}", "\\begin{tablenotes}[flushleft]", "\\footnotesize", "\\item Notes: All yield maturities used in estimation are reported. RMSE and MAE compare observations with measurement functions evaluated at the updated filtered states; fitted real yields include the liquidity intercept. Correlation is the Pearson time-series correlation between each observed series and this fitted counterpart, using only dates when both are available. These are in-sample fits, not one-step-ahead forecasts. The CPI row reports trailing twelve-month inflation, as in the fit figure. Missing survey observations are left missing and are not carried forward. Inflation and interest rates are in annualized percentage points.", "\\end{tablenotes}", "\\end{threeparttable}", "\\end{table}")
writeLines(fl, file.path(out, "table_fit.tex"), useBytes=TRUE)
