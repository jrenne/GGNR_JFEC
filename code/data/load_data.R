# Load the frozen all-R database used by the paper.

load_paper_data <- function(path = "data/input/paper_data.rds") {
  readRDS(path)
}

# Retain this function name internally while the model code is being cleaned.
load_true_release_data <- load_paper_data

no_hfi_data <- function(dates) {
  list(
    standardized = rep(NaN, length(dates)),
    observed = rep(FALSE, length(dates)),
    source = "none", shock = NA_character_,
    event_count = 0L, month_count = 0L
  )
}
