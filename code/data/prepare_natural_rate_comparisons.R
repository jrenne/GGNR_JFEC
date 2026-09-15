# Rebuild Figure 1's comparison series from the frozen New York Fed files.
# Open the package project in RStudio and source this script if needed.
# Both series are one-sided (filtered), not two-sided (smoothed).

source_directory <- "data/external/natural_rate_sources"
lw <- readxl::read_xlsx(
  file.path(source_directory, "LW_20260827.xlsx"),
  sheet = "data", range = "A7:C268", col_names = FALSE
)
hlw <- readxl::read_xlsx(
  file.path(source_directory, "HLW_20260828.xlsx"),
  sheet = "HLW Estimates", range = "A7:K268", col_names = FALSE
)
quarter_end_month <- function(date) {
  date <- as.Date(date)
  as.Date(sprintf("%s-%02d-01", format(date, "%Y"),
                  as.integer(format(date, "%m")) + 2L))
}
benchmark <- merge(
  data.frame(date = quarter_end_month(lw[[1]]), LW = as.numeric(lw[[3]])),
  data.frame(date = quarter_end_month(hlw[[1]]), HLW = as.numeric(hlw[[11]])),
  by = "date", all = TRUE
)
attr(benchmark, "source") <- list(
  provider = "Federal Reserve Bank of New York",
  url = "https://www.newyorkfed.org/research/policy/rstar",
  retrieved = "2026-09-15",
  release = c(LW = "2026-08-27", HLW = "2026-08-28"),
  final_quarter = "2026Q2",
  status = c(LW = "one-sided", HLW = "one-sided"),
  units = "annualized percentage points",
  dating = "First day of the quarter's final month; no monthly interpolation"
)
saveRDS(benchmark, "data/external/natural_rate_comparisons.rds", version = 3)
