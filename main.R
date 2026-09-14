# Replication driver
#
# Open GGNR_JFEC.Rproj, edit the switches below if desired, and click Source.

run_estimation <- FALSE
run_robustness_estimation <- FALSE
run_inference <- FALSE
run_pricing_exercise <- FALSE
run_filtering_exercise <- FALSE
run_uncertainty_exercise <- FALSE

# A new baseline estimate changes every estimate-dependent downstream result.
if (run_estimation) {
  run_robustness_estimation <- TRUE
  run_inference <- TRUE
  run_filtering_exercise <- TRUE
  run_uncertainty_exercise <- TRUE
}

if (run_estimation) source("estimate_model.R")
if (run_robustness_estimation) source("estimate_robustness.R")

# Baseline filtered states and fit diagnostics are inexpensive and are always
# regenerated from the archived estimate and frozen paper data.
source("code/outputs/build_diagnostics.R")

if (run_inference) {
  source("code/inference/sandwich_covariance.R")
  source("code/inference/hac_sandwich_covariance.R")
  source("code/inference/joint_price_risk_tests.R")
}
source("code/validation/pricing_approximation.R")
if (run_filtering_exercise) source("code/validation/filtering_accuracy.R")
if (run_uncertainty_exercise) source("code/inference/hamilton_uncertainty.R")

# Every figure and table is rebuilt from the estimates and numerical results.
source("code/outputs/plot_fits_and_states.R")
source("code/outputs/plot_risk_premia.R")
source("code/outputs/plot_iso_probability.R")
source("code/outputs/plot_mundell_tobin.R")
source("code/outputs/plot_yield_responses.R")
source("code/outputs/plot_shadow_rate_comparison.R")
source("code/outputs/irp_fit_decomposition.R")

if (file.exists("estimates/omit_2020.rds")) {
  source("code/validation/covid_robustness.R")
}
if (file.exists("estimates/liquid_tips_only.rds")) {
  source("code/validation/liquid_tips_robustness.R")
}

if (file.exists("outputs/diagnostics/inference/parameter_inference.csv")) {
  source("code/outputs/make_tables.R")
}
if (file.exists("outputs/diagnostics/filtering/filtering_summary.csv")) {
  source("code/outputs/make_filtering_table.R")
}
if (file.exists("outputs/diagnostics/uncertainty/state_bands.csv")) {
  source("code/outputs/plot_uncertainty.R")
}

message("Replication outputs are available in outputs/.")
