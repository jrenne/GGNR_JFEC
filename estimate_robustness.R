# Full re-estimation of the two robustness specifications
#
# Open this file in RStudio and click Source. Each specification starts from
# the baseline estimate, then uses a short local refinement of its own result.

estimate_robustness <- function(output_file, omit_2020, liquid_tips_only) {
  rounds <- data.frame(
    optimizer = c("bobyqa", "nlminb"),
    evaluations = c(1000L, 500L)
  )
  for (round in seq_len(nrow(rounds))) {
    Sys.setenv(
      GGNR_START_FILE = if (round == 1L) "estimates/baseline.rds" else output_file,
      GGNR_OUTPUT_FILE = output_file,
      GGNR_EVAL_BUDGET = as.character(rounds$evaluations[round]),
      GGNR_OPTIMIZER = rounds$optimizer[round],
      GGNR_OMIT_START = if (omit_2020) "2020-01-01" else "",
      GGNR_OMIT_END = if (omit_2020) "2020-12-31" else "",
      GGNR_LIQUID_TIPS_ONLY = if (liquid_tips_only) "true" else "false"
    )
    source("code/estimation/estimate_baseline.R")
  }
}

estimate_robustness("estimates/omit_2020.rds", TRUE, FALSE)
estimate_robustness("estimates/liquid_tips_only.rds", FALSE, TRUE)

Sys.unsetenv(c(
  "GGNR_START_FILE", "GGNR_EVAL_BUDGET", "GGNR_OUTPUT_FILE",
  "GGNR_OPTIMIZER", "GGNR_OMIT_START", "GGNR_OMIT_END",
  "GGNR_LIQUID_TIPS_ONLY"
))
