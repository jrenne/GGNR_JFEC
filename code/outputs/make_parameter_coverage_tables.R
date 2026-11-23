# OA.3.2 reference table and OA.3.3 higher-ZLB-frequency sensitivity table.
source("code/validation/parameter_coverage_helpers.R")
directory <- Sys.getenv("GGNR_MC_DIR", "outputs/diagnostics/parameter_coverage")
results <- read.csv(file.path(directory, "replications.csv"))
summary <- summarize_parameter_coverage(results)
stopifnot(nrow(summary) == 10L, all(summary$replications == 500L))
write.csv(summary, file.path(directory, "summary.csv"), row.names = FALSE)
output_directory <- "outputs/tables"
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
parameters <- c("logit_rho_r", "logit_rho_pi", "log_sd_r", "log_sd_pi", "log_measurement")
labels <- c("$\\operatorname{logit}\\rho_r$", "$\\operatorname{logit}\\rho_\\pi$",
            "$\\log\\sigma_r$", "$\\log\\sigma_\\pi$", "$\\log\\sigma_{me}$")
for (case in c("reference_15", "higher_zlb_50")) {
  suffix <- if (case == "reference_15") "reference" else "sensitivity"
  title <- if (case == "reference_15") "Finite-sample parameter estimation: reference design" else
    "Sensitivity to a higher frequency of the zero lower bound"
  selected <- summary[summary$case == case, ]
  selected <- selected[match(parameters, selected$parameter), ]
  lines <- c("\\begin{table}[!htbp]", "\\centering", paste0("\\caption{", title, "}"),
    paste0("\\label{oa:tab_parameter_mc_", suffix, "}"), "\\begin{threeparttable}",
    "\\small", "\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}lrrrrr}", "\\toprule",
    "Parameter & Bias & Empirical SD & RMSE & Mean SE & Coverage (\\%) \\\\", "\\midrule")
  for (i in seq_along(parameters)) {
    row <- selected[i, ]
    values <- sprintf("%.4f", unlist(row[c("bias", "empirical_sd", "rmse", "mean_se_hac")]))
    lines <- c(lines, paste0(labels[i], " & ", paste(values, collapse = " & "),
                            " & ", sprintf("%.1f", 100 * row$coverage_hac), " \\\\"))
  }
  lines <- c(lines, "\\bottomrule", "\\end{tabular*}",
    "\\begin{tablenotes}[flushleft]", "\\footnotesize",
    paste0("\\item Notes: 500 replications of 695 months. All diagnostics are in logit/log estimation coordinates. ",
      "Mean SE and coverage use the local Hessian--HAC sandwich with 12 Bartlett lags. ",
      "Coverage is for nominal 95\\% Wald intervals; monotone transformation to economic coordinates preserves it. ",
      "$\\sigma_{me}$ is the common measurement-error scale, with generating value 0.0663 annualized percentage point; ",
      "relative errors are fixed as described in the text. ",
      "All 500 intervals are computable for each parameter; no sample is omitted."),
    "\\end{tablenotes}", "\\end{threeparttable}", "\\end{table}")
  writeLines(lines, file.path(output_directory, paste0("table_parameter_mc_", suffix, ".tex")))
}
