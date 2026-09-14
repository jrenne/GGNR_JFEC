# Full re-estimation of the baseline model
#
# Open this file in RStudio and click Source. The resulting estimate is saved
# as estimates/baseline.rds. The three successive rounds provide a transparent
# cross-optimizer check and each round starts from the best preceding estimate.

estimation_rounds <- data.frame(
  optimizer = c("bobyqa", "Nelder-Mead", "nlminb"),
  evaluations = c(1500L, 500L, 500L)
)

for (round in seq_len(nrow(estimation_rounds))) {
  Sys.setenv(
    GGNR_START_FILE = if (round == 1L)
      "estimates/start_values.rds" else "estimates/baseline.rds",
    GGNR_OUTPUT_FILE = "estimates/baseline.rds",
    GGNR_EVAL_BUDGET = as.character(estimation_rounds$evaluations[round]),
    GGNR_OPTIMIZER = estimation_rounds$optimizer[round]
  )
  source("code/estimation/estimate_baseline.R")
}

Sys.unsetenv(c(
  "GGNR_START_FILE", "GGNR_OUTPUT_FILE", "GGNR_EVAL_BUDGET",
  "GGNR_OPTIMIZER"
))
