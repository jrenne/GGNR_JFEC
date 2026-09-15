# Data

The package deliberately separates exact replication from database updating.

## Frozen paper vintage

`input/paper_data.rds` contains the all-R database adopted for the September
2026 JFEC revision. It is the default input to the estimation and figure code. In
particular, survey and inflation-target observations are present only in their
actual release months; other months contain `NaN`. The sample runs from October
1968 through August 2026.

The object contains headline inflation, the FRB/US perceived inflation target,
one- and ten-year SPF inflation and Treasury-bill forecasts, nominal Treasury
yields, real Treasury yields, their maturities and forecast horizons, and the
effective federal funds rate used to impose the short-real-rate moments. Some
legacy fields are retained to keep the data interface compatible with the
estimation code, although output, GDP growth, GDP forecasts, and core CPI are
not observables in the retained specification.

## Optional update

Open `prepare_data.R` in RStudio, choose `sample_end`, and click **Source**.
The script obtains current-vintage public series and writes
`updated/updated_data.rds`. It never changes `input/paper_data.rds`.
No API key is required.

The frozen input and the default updating endpoint are August 2026, the latest
complete month shared by the central monthly series when this package was
assembled. Users can
change `sample_end` at the top of `prepare_data.R` as newer observations become
available.

Fresh downloads can differ from the paper vintage because source agencies
revise historical data. Thus, the updated file is intended for extending the
analysis, whereas the frozen file is intended for exact replication.

## Sources

- FRED: effective federal funds rate, three-month Treasury bill rate, and CPI.
- Federal Reserve Board: Gürkaynak--Sack--Wright nominal and real yield curves.
- Federal Reserve Bank of Philadelphia: Survey of Professional Forecasters
  inflation and Treasury-bill forecasts.
- Federal Reserve Board FRB/US database: perceived inflation target.
- Federal Reserve Bank of New York: synthetic pre-TIPS ten-year real rate.

`external/model_free_risk_premia.rds` contains the model-free series displayed
with the model-implied risk premia. `external/shadow_rate_comparisons.rds`
contains the Wu--Xia and Krippner series displayed in the shadow-rate
comparison. These comparison data are not used to estimate the model.

`external/natural_rate_comparisons.rds` contains the one-sided (filtered)
LW and HLW U.S. natural-rate estimates shown in Figure 1. The New York Fed
releases are August 27, 2026 (LW) and August 28, 2026 (HLW), both ending in
2026Q2. LW's input-data cutoff is August 26, 2026. These current-vintage
estimates use revised data and full-sample parameter estimates; they are not
historical real-time vintages. Quarterly values are assigned to the quarter's
final month, without monthly interpolation. The figure connects the points.
The exact source workbooks are in `external/natural_rate_sources/`.
To rebuild the comparison file, source
`code/data/prepare_natural_rate_comparisons.R` from the package project in
RStudio. The script selects column C (one-sided r-star) of LW's `data` sheet
and column K (U.S. r-star) of the all-one-sided `HLW Estimates` sheet.
Source: https://www.newyorkfed.org/research/policy/rstar.

The data and comparison series remain subject to the terms of their original
providers. The MIT license in the repository applies to the authors' code, not
to third-party data.
