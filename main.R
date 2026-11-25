# Open GGNR_JFEC.Rproj and Source this file. No downloads are needed.
# Default: regenerate every published table/figure from the archived R2 estimate
# and checked numerical caches. Set individual switches for fresh calculations.
run_estimation <- FALSE
run_polishing <- FALSE
run_robustness_estimation <- FALSE
run_inference <- FALSE
run_pricing_exercise <- FALSE
run_filtering_exercise <- FALSE
run_uncertainty_exercise <- FALSE
run_parameter_exercise <- FALSE
run_mundell_exercise <- FALSE
run_checks <- TRUE

replication_started <- proc.time()["elapsed"]
ggnr_run_stage <- function(file) {
  message("Running ", file)
  source(file, local = globalenv())
}
if (run_estimation) {
  ggnr_run_stage("estimate_model.R")
  Sys.setenv(GGNR_START_FILE = "estimates/reestimated.rds")
  run_inference <- run_filtering_exercise <- run_uncertainty_exercise <- TRUE
  run_robustness_estimation <- run_mundell_exercise <- TRUE
}
if (run_polishing) {
  ggnr_run_stage("polish_model.R")
  Sys.setenv(GGNR_START_FILE = Sys.getenv("GGNR_POLISH_OUTPUT", "estimates/refined.rds"))
  run_inference <- run_filtering_exercise <- run_uncertainty_exercise <- TRUE
  run_mundell_exercise <- TRUE
}
if (run_robustness_estimation) ggnr_run_stage("estimate_robustness.R")
# estimate_robustness restores the incoming environment, including the selected estimate.
if (run_inference) ggnr_run_stage("code/inference/hessian_sandwich.R")
if (run_filtering_exercise) ggnr_run_stage("code/validation/filtering_accuracy.R")
if (run_uncertainty_exercise) ggnr_run_stage("code/inference/delta_uncertainty.R")
if (run_parameter_exercise) ggnr_run_stage("code/validation/parameter_coverage.R")
ggnr_run_stage("code/validation/validate_cached_results.R")
ggnr_run_stage("code/outputs/build_diagnostics.R")
ggnr_run_stage("code/inference/joint_price_risk_tests.R")
ggnr_run_stage("code/validation/pricing_approximation.R")
# Build fitted states before the comparison figure, which reads that fresh CSV.
ggnr_run_stage("code/outputs/plot_fits_and_states.R")
ggnr_run_stage("code/outputs/plot_risk_premia.R")
ggnr_run_stage("code/outputs/plot_iso_probability.R")
options(ggnr.recompute_mundell = run_mundell_exercise)
ggnr_run_stage("code/outputs/plot_mundell_tobin.R")
ggnr_run_stage("code/outputs/plot_yield_responses.R")
ggnr_run_stage("code/outputs/plot_shadow_rate_comparison.R")
ggnr_run_stage("code/outputs/irp_fit_decomposition.R")
ggnr_run_stage("code/validation/covid_robustness.R")
ggnr_run_stage("code/validation/liquid_tips_robustness.R")
ggnr_run_stage("code/outputs/make_tables.R")
ggnr_run_stage("code/outputs/make_filtering_table.R")
ggnr_run_stage("code/outputs/make_parameter_coverage_tables.R")
ggnr_run_stage("code/outputs/plot_uncertainty.R")
if (run_checks) {
  ggnr_run_stage("code/validation/estimate_stationarity.R")
  ggnr_run_stage("code/validation/delta_alignment.R")
  ggnr_run_stage("code/validation/covariance_alignment.R")
  ggnr_run_stage("code/validation/parameter_coverage_alignment.R")
  ggnr_run_stage("code/validation/pricing_formula_alignment.R")
  ggnr_run_stage("code/validation/survey_expectation_alignment.R")
  ggnr_run_stage("code/validation/pricing_timing_alignment.R")
}
message(sprintf("Replication completed in %.1f seconds; results are in outputs/.",
                proc.time()["elapsed"] - replication_started))
