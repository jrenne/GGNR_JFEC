# Load the frozen all-R database used by the paper.

load_paper_data <- function(path = "data/input/paper_data.rds") {
  data <- readRDS(path)
  keep <- data$dates >= as.Date("1992-02-01")
  for (name in names(data)) {
    x <- data[[name]]
    if (is.matrix(x) && nrow(x) == length(keep)) data[[name]] <- x[keep, , drop = FALSE]
    else if (is.null(dim(x)) && length(x) == length(keep)) data[[name]] <- x[keep]
  }
  data$release_counts <- c(ptr = sum(is.finite(data$macro[, 4])),
    inflation_1y = sum(is.finite(data$surv_infexp[, 1])),
    inflation_10y = sum(is.finite(data$surv_infexp[, 2])),
    tbill_1y = sum(is.finite(data$surv_tbexp[, 1])),
    tbill_10y = sum(is.finite(data$surv_tbexp[, 2])))
  data
}

no_hfi_data <- function(dates) {
  list(
    standardized = rep(NaN, length(dates)),
    observed = rep(FALSE, length(dates)),
    source = "none", shock = NA_character_,
    event_count = 0L, month_count = 0L
  )
}
