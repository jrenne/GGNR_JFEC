# Exact decomposition of the survey-implied 10-year inflation-risk-premium fit.

source("code/data/load_data.R")
diagnostic_directory <- Sys.getenv(
  "GGNR_DIAGNOSTIC_DIR",
  "outputs/diagnostics/baseline"
)
premium_directory <- Sys.getenv(
  "GGNR_PREMIA_DIR",
  "outputs/diagnostics/risk_premia"
)
output_file <- Sys.getenv(
  "GGNR_IRP_DECOMP_PNG",
  file.path(premium_directory, "irp_fit_decomposition.png")
)

data <- load_paper_data()
diagnostics <- readRDS(file.path(diagnostic_directory, "diagnostics.rds"))
premia <- read.csv(file.path(premium_directory, "risk_premia.csv"))

nominal_error <- 1200 * (
  diagnostics$fitted$nominal_yields[, 7] - data$yields_n[, 7]
)
real_error <- 1200 * (
  diagnostics$fitted$real_yields[, 4] - data$yields_r[, 4]
)
cpi10_error <- 1200 * (
  diagnostics$fitted$inflation_survey[, 2] - data$surv_infexp[, 2]
)
fitted_observable_irp <- 1200 * (
  diagnostics$fitted$nominal_yields[, 7] -
    diagnostics$fitted$real_yields[, 4] -
    diagnostics$fitted$inflation_survey[, 2]
)
observed_irp <- premia$inflation_rp_10y_survey
structural_irp <- premia$inflation_rp_10y
identity_error <- fitted_observable_irp - observed_irp
decomposed_error <- nominal_error - real_error - cpi10_error
stopifnot(max(abs(identity_error - decomposed_error), na.rm = TRUE) < 1e-10)
structural_error <- structural_irp - observed_irp
definition_wedge <- structural_irp - fitted_observable_irp
stopifnot(max(abs(structural_error - identity_error - definition_wedge),
              na.rm = TRUE) < 1e-10)

observed <- is.finite(observed_irp)
decomposition <- data.frame(
  date = data$dates,
  observed_irp = observed_irp,
  fitted_observable_irp = fitted_observable_irp,
  structural_irp = structural_irp,
  structural_error = structural_error,
  definition_wedge = definition_wedge,
  irp_fit_error = identity_error,
  nominal_yield_contribution = nominal_error,
  real_yield_contribution = -real_error,
  cpi10_contribution = -cpi10_error
)
write.csv(
  decomposition,
  file.path(premium_directory, "irp_fit_decomposition.csv"),
  row.names = FALSE
)

png(output_file, width = 1800, height = 1600, res = 180)
par(mfrow = c(2, 1), mar = c(3.2, 4.2, 2.5, 1))
limits <- range(c(observed_irp[observed], fitted_observable_irp[observed],
                  structural_irp[observed]), finite = TRUE)
plot(data$dates[observed], observed_irp[observed], type = "p", pch = 4,
     col = "#222222", xlab = "", ylab = "percentage points", ylim = limits,
     main = "10-year inflation risk premium: observed and fitted")
lines(data$dates[observed], fitted_observable_irp[observed],
      col = "#D55E00", lwd = 1.5)
lines(data$dates[observed], structural_irp[observed],
      col = "#0072B2", lwd = 1.2, lty = 2)
abline(h = 0, col = "grey55", lwd = 0.7)
legend("topright", c("Survey-implied", "From fitted observables", "Structural"),
       col = c("#222222", "#D55E00", "#0072B2"),
       pch = c(4, NA, NA), lty = c(NA, 1, 2), lwd = c(NA, 1.5, 1.2),
       cex = 0.8, bg = adjustcolor("white", alpha.f = 0.85))

plot(data$dates[observed], identity_error[observed], type = "l",
     col = "#222222", lwd = 1.5, xlab = "", ylab = "percentage points",
     main = "Exact IRP fit-error decomposition")
lines(data$dates[observed], nominal_error[observed], col = "#009E73", lwd = 1)
lines(data$dates[observed], -real_error[observed], col = "#D55E00", lwd = 1)
lines(data$dates[observed], -cpi10_error[observed], col = "#0072B2", lwd = 1)
abline(h = 0, col = "grey55", lwd = 0.7)
legend("bottomleft",
       c("Total IRP error", "Nominal-yield contribution",
         "Real-yield contribution", "CPI10 contribution"),
       col = c("#222222", "#009E73", "#D55E00", "#0072B2"),
       lty = 1, lwd = c(1.5, 1, 1, 1), cex = 0.8,
       bg = adjustcolor("white", alpha.f = 0.85))
dev.off()

rmse <- function(x, mask = observed) sqrt(mean(x[mask]^2, na.rm = TRUE))
average <- function(x, mask = observed) mean(x[mask], na.rm = TRUE)
post2004 <- observed & data$dates >= as.Date("2004-01-01")
series <- list(
  structural_error, definition_wedge, identity_error,
  nominal_error, -real_error, -cpi10_error
)
summary <- data.frame(
  component = c(
    "Model premium minus survey proxy",
    "Definition and convexity wedge",
    "Fitted-observable proxy minus survey proxy", "Nominal 10y yield contribution",
    "Real 10y yield contribution", "CPI10 contribution"
  ),
  mean_all_irp_dates = vapply(series, average, numeric(1)),
  rmse_all_irp_dates = vapply(series, rmse, numeric(1)),
  mean_2004_onward = vapply(series, average, numeric(1), mask = post2004),
  rmse_2004_onward = vapply(series, rmse, numeric(1), mask = post2004)
)
write.csv(summary, file.path(premium_directory, "irp_fit_error_summary.csv"),
          row.names = FALSE)

table_directory <- Sys.getenv(
  "GGNR_TABLE_DIR", "outputs/tables"
)
dir.create(table_directory, recursive = TRUE, showWarnings = FALSE)
table_labels <- c(
  "Model premium $-$ survey proxy",
  "\\quad Definition and convexity wedge",
  "\\quad Fitted-observable proxy $-$ survey proxy",
  "\\qquad Nominal-yield fit contribution",
  "\\qquad Real-yield fit contribution",
  "\\qquad CPI-forecast fit contribution"
)
table_lines <- c(
  "\\begin{table}[!htbp]", "\\centering",
  "\\caption{Decomposition of the ten-year survey-based inflation-risk-premium discrepancy}",
  "\\label{tab:irp_decomposition}", "\\begin{threeparttable}", "\\small",
  "\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}lrrrr}",
  "\\toprule",
  " & \\multicolumn{2}{c}{All survey dates} & \\multicolumn{2}{c}{2004 onward} \\\\",
  "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5}",
  "Component & Mean & RMSE & Mean & RMSE \\\\ \\midrule"
)
for (index in seq_len(nrow(summary))) {
  table_lines <- c(table_lines, sprintf(
    "%s & %.3f & %.3f & %.3f & %.3f \\\\",
    table_labels[index], summary$mean_all_irp_dates[index],
    summary$rmse_all_irp_dates[index], summary$mean_2004_onward[index],
    summary$rmse_2004_onward[index]
  ))
}
table_lines <- c(
  table_lines, "\\bottomrule", "\\end{tabular*}",
  "\\begin{tablenotes}[flushleft]", "\\footnotesize",
  paste0(
    "\\item Notes: All entries are annualized percentage points. The survey-based expectation-hypothesis proxy is the observed ten-year nominal yield minus the observed ten-year real yield and the ten-year CPI forecast. The fitted-observable proxy replaces those three inputs with their fitted values. The model premium compares log certainty-equivalent values of cumulative inflation under the risk-neutral and physical measures. Their difference is a definition and convexity wedge generated by Jensen, covariance, and nonlinear lower-bound terms; it is not a fitting error. The fitted-observable error equals exactly the sum of the three fit contributions; their RMSEs do not add because the contributions covary. The columns contain ",
    sum(observed), " and ", sum(post2004), " survey dates, respectively."
  ),
  "\\end{tablenotes}", "\\end{threeparttable}", "\\end{table}"
)
writeLines(table_lines,
           file.path(table_directory, "table_irp_decomposition.tex"),
           useBytes = TRUE)
print(summary)
