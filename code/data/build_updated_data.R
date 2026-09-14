# Functions used by prepare_data.R to reconstruct a current-vintage database.

required_data_packages <- "readxl"
missing_data_packages <- required_data_packages[
  !vapply(required_data_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_data_packages)) {
  stop("Please install: ", paste(missing_data_packages, collapse = ", "))
}

month_start <- function(x) as.Date(format(as.Date(x), "%Y-%m-01"))

monthly_average <- function(x, date_name = "Date") {
  x$date <- month_start(x[[date_name]])
  values <- x[vapply(x, is.numeric, logical(1))]
  aggregate(values, list(date = x$date), mean, na.rm = TRUE)
}

download_to <- function(url, directory, name) {
  path <- file.path(directory, name)
  utils::download.file(url, path, mode = "wb", quiet = TRUE)
  path
}

read_spf <- function(path, value_name, source_column) {
  x <- readxl::read_xlsx(path)
  release_month <- 2L + 3L * (as.integer(x$QUARTER) - 1L)
  date <- as.Date(sprintf("%04d-%02d-01", as.integer(x$YEAR), release_month))
  value <- suppressWarnings(as.numeric(x[[source_column]]))
  setNames(data.frame(date, value), c("date", value_name))
}

build_updated_data <- function(sample_start, sample_end) {
  sample_start <- month_start(sample_start)
  sample_end <- month_start(sample_end)
  dates <- seq(sample_start, sample_end, by = "month")
  raw_directory <- file.path("data", "updated", "downloads")
  dir.create(raw_directory, recursive = TRUE, showWarnings = FALSE)

  fred_series <- function(series_id) {
    path <- download_to(
      paste0(
        "https://fred.stlouisfed.org/graph/fredgraph.csv?id=", series_id,
        "&cosd=", format(sample_start - 40, "%Y-%m-%d"),
        "&coed=", format(sample_end, "%Y-%m-%d")
      ),
      raw_directory, paste0(series_id, ".csv")
    )
    x <- utils::read.csv(path, na.strings = ".", check.names = FALSE)
    x[[1]] <- as.Date(x[[1]])
    names(x)[1:2] <- c("Date", "value")
    monthly_average(x)[, c("date", "value")]
  }

  fedfunds <- fred_series("FEDFUNDS")
  dtb3 <- fred_series("DTB3")
  cpi <- fred_series("CPIAUCSL")
  cpi$inflation <- c(NA_real_, diff(log(cpi$value)))

  nominal_path <- download_to(
    "https://www.federalreserve.gov/data/yield-curve-tables/feds200628.csv",
    raw_directory, "feds200628.csv"
  )
  nominal_daily <- utils::read.csv(nominal_path, skip = 9, check.names = FALSE)
  nominal_daily$Date <- as.Date(nominal_daily$Date)
  nominal <- monthly_average(nominal_daily)

  real_path <- download_to(
    "https://www.federalreserve.gov/data/yield-curve-tables/feds200805.csv",
    raw_directory, "feds200805.csv"
  )
  real_daily <- utils::read.csv(real_path, skip = 18, check.names = FALSE)
  real_daily$Date <- as.Date(real_daily$Date)
  real <- monthly_average(real_daily)

  spf_root <- paste0(
    "https://www.philadelphiafed.org/-/media/frbp/assets/",
    "surveys-and-data/survey-of-professional-forecasters/data-files/files/"
  )
  spf_cpi_1 <- read_spf(
    download_to(paste0(spf_root, "mean_cpi_level.xlsx"), raw_directory,
                "mean_cpi_level.xlsx"),
    "CPI1", "CPI6"
  )
  spf_cpi_10 <- read_spf(
    download_to(paste0(spf_root, "mean_cpi10_level.xlsx"), raw_directory,
                "mean_cpi10_level.xlsx"),
    "CPI10", "CPI10"
  )
  spf_tbill_1 <- read_spf(
    download_to(paste0(spf_root, "mean_tbill_level.xlsx"), raw_directory,
                "mean_tbill_level.xlsx"),
    "BILL1", "TBILL6"
  )
  spf_tbill_10 <- read_spf(
    download_to(paste0(spf_root, "mean_bill10_level.xlsx"), raw_directory,
                "mean_bill10_level.xlsx"),
    "BILL10", "BILL10"
  )

  frbus_zip <- download_to(
    "https://www.federalreserve.gov/econres/files/data_only_package.zip",
    raw_directory, "frbus_data_only_package.zip"
  )
  frbus_directory <- file.path(raw_directory, "frbus")
  dir.create(frbus_directory, showWarnings = FALSE)
  utils::unzip(
    frbus_zip, files = "data_only_package/HISTDATA.TXT",
    exdir = frbus_directory
  )
  ptr <- utils::read.csv(file.path(
    frbus_directory, "data_only_package", "HISTDATA.TXT"
  ))
  ptr$date <- as.Date(sprintf(
    "%s-%02d-01", substr(ptr$OBS, 1, 4),
    3L * as.integer(substr(ptr$OBS, 6, 6))
  ))
  ptr <- data.frame(date = ptr$date, PTR = as.numeric(ptr$PTR))

  synthetic_path <- download_to(
    paste0("https://www.newyorkfed.org/medialibrary/media/research/",
           "blog/groen_tips/Synthetic_TIPS_Breakeven_Rates.xlsx"),
    raw_directory, "Synthetic_TIPS_Breakeven_Rates.xlsx"
  )
  synthetic <- readxl::read_xlsx(synthetic_path, skip = 13)
  synthetic <- data.frame(
    date = month_start(synthetic[[1]]),
    RR10Y = as.numeric(synthetic[["10-yr real rate, backcast"]])
  )

  align <- function(x, name) x[[name]][match(dates, x$date)]
  yields_n <- cbind(
    align(dtb3, "value"),
    nominal$SVENY01[match(dates, nominal$date)],
    nominal$SVENY02[match(dates, nominal$date)],
    nominal$SVENY03[match(dates, nominal$date)],
    nominal$SVENY05[match(dates, nominal$date)],
    nominal$SVENY07[match(dates, nominal$date)],
    nominal$SVENY10[match(dates, nominal$date)]
  ) / 1200
  yields_r <- cbind(
    real$TIPSY02[match(dates, real$date)],
    real$TIPSY05[match(dates, real$date)],
    real$TIPSY07[match(dates, real$date)],
    real$TIPSY10[match(dates, real$date)]
  ) / 1200
  synthetic_10y <- align(synthetic, "RR10Y") / 1200
  replace_10y <- !is.finite(yields_r[, 4]) & is.finite(synthetic_10y)
  yields_r[replace_10y, 4] <- synthetic_10y[replace_10y]

  inflation <- align(cpi, "inflation")
  ptr_monthly <- align(ptr, "PTR") / 1200
  macro <- cbind(
    growth = rep(NaN, length(dates)),
    inflation = inflation,
    output_gap = rep(NaN, length(dates)),
    inflation_target = ptr_monthly
  )
  surv_infexp <- cbind(
    align(spf_cpi_1, "CPI1"), align(spf_cpi_10, "CPI10")
  ) / 1200
  surv_tbexp <- cbind(
    align(spf_tbill_1, "BILL1"), align(spf_tbill_10, "BILL10")
  ) / 1200

  release_counts <- c(
    ptr = sum(is.finite(ptr_monthly)),
    inflation_1y = sum(is.finite(surv_infexp[, 1])),
    inflation_10y = sum(is.finite(surv_infexp[, 2])),
    tbill_1y = sum(is.finite(surv_tbexp[, 1])),
    tbill_10y = sum(is.finite(surv_tbexp[, 2]))
  )
  any_survey_release <- rowSums(cbind(
    is.finite(surv_infexp), is.finite(surv_tbexp)
  )) > 0

  list(
    macro = macro,
    yields_n = yields_n,
    yields_r = yields_r,
    surv_infexp = surv_infexp,
    surv_gdpexp = matrix(NaN, length(dates), 1),
    surv_tbexp = surv_tbexp,
    core_inflation = rep(NaN, length(dates)),
    mats_n = c(3L, 12L, 24L, 36L, 60L, 84L, 120L),
    mats_r = c(24L, 60L, 84L, 120L),
    hstep_s = c(12L, 120L),
    hstep_g = 120L,
    hstep_t = c(12L, 120L),
    dates = dates,
    release_counts = release_counts,
    any_survey_release = any_survey_release,
    any_nonmonthly_release = any_survey_release | is.finite(ptr_monthly),
    effective_federal_funds_rate = align(fedfunds, "value") / 1200
  )
}
