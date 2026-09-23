# Load the frozen all-R database used by the paper.

load_paper_data <- function(path = "data/input/paper_data.rds") {
  readRDS(path)
}

no_hfi_data <- function(dates) {
  list(
    standardized = rep(NaN, length(dates)),
    observed = rep(FALSE, length(dates)),
    source = "none", shock = NA_character_,
    event_count = 0L, month_count = 0L
  )
}
