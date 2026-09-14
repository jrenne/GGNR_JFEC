# Replication package for "The Shadow-Rate Model: Let's Make it Real"

**Authors:** Adam Golinski, Sophie Guilloux-Nefussi, and Jean-Paul Renne  
**Revision package:** September 2026

## Overview

This repository contains the code and data used to reproduce every numerical
table and figure in the paper and its Online Appendix. The retained empirical
specification is the six-state reduced-form model estimated with the fixed
two-step iterated extended Kalman filter (IEKF), using the all-R monthly
database from October 1968 through August 2026.

The workflow is designed around two uses:

1. **Quick replication (default).** Load the archived parameter estimates and
   cached results from lengthy simulations, then rerun all table and figure
   construction. This mode reproduces every published artifact without
   re-estimating the model or repeating the Monte Carlo exercises.
2. **Full replication.** Re-estimate the baseline and robustness models,
   recompute the covariance matrix and simulation exercises, and then recreate
   every table and figure.

Cached results are part of the replication record. Skipping a lengthy stage
never means copying a finished manuscript figure: the corresponding plotting
or table-building code is always rerun from the archived numerical output.

## Running the package in RStudio

Open `GGNR_JFEC.Rproj` in RStudio and then open `main.R`. A short block at the
beginning of `main.R` controls the lengthy computations:

```r
run_estimation <- FALSE
run_robustness_estimation <- FALSE
run_inference <- FALSE
run_pricing_exercise <- FALSE
run_filtering_exercise <- FALSE
run_uncertainty_exercise <- FALSE
```

With the default `FALSE` settings, the code loads the archived estimates and
cached results of the lengthy numerical exercises. Setting a switch to `TRUE`
recomputes that stage. Setting `run_estimation <- TRUE` also refreshes every
estimate-dependent downstream stage; the other switches can still be used
separately with the archived baseline estimate. In either case, clicking
**Source** in RStudio recreates the complete collection of manuscript and
Online Appendix tables and figures in `outputs/`.

The required contributed R packages are `Rcpp`, `RcppEigen`, `nloptr`,
`readxl`, and `ggplot2`. The estimation and simulation exercises also use R's
standard `parallel` package when several processor cores are available. The
archived results were last reproduced with R 4.5.1, Rcpp 1.1.0, RcppEigen
0.3.4.0.2, nloptr 2.2.1, readxl 1.4.5, and ggplot2 4.0.1.

## Published artifacts covered

The driver will recreate:

- Tables 1--2 and Figures 1--5 in the paper;
- the pricing-approximation table and figure and the filtering-accuracy table;
- all observable-fit figures;
- the 2020-omission figure;
- the inflation-risk-premium discrepancy table;
- the TIPS-liquidity table and pre-2004 sensitivity figure;
- the parameter-and-filtering uncertainty figures;
- the alternative shadow-rate comparison; and
- the additional shadow-rate and yield-response figures in the Online
  Appendix.

## Directory structure

```text
GGNR_JFEC/
├── README.md
├── LICENSE                        MIT license for the authors' code
├── GGNR_JFEC.Rproj                RStudio project
├── main.R                         single replication entry point
├── prepare_data.R                  optional database reconstruction
├── estimate_model.R               optional full baseline estimation
├── estimate_robustness.R          optional robustness estimations
├── code/
│   ├── model/                     state-space, pricing, and filtering code
│   ├── data/                      database-construction functions
│   ├── estimation/                parameter mappings and likelihood routines
│   ├── inference/                 OPG/HAC and uncertainty calculations
│   ├── validation/                pricing and filtering assessments
│   └── outputs/                   table and figure builders
├── data/
│   ├── input/                     frozen data vintage used in the paper
│   ├── updated/                   optional current-vintage reconstruction
│   └── external/                  published comparison series
├── estimates/                     archived baseline and robustness estimates
└── outputs/
    ├── figures/                   all paper and Online Appendix figures
    ├── tables/                    all paper and Online Appendix tables
    └── diagnostics/               intermediate numerical results and checks
```

## Data and database reconstruction

The package supports two distinct data workflows.

1. **Exact replication of the paper.** The final database actually used for
   estimation is supplied in `data/input/`. This frozen paper vintage is
   the default input to `main.R`; no downloading is required, and all reported
   results are reproduced from exactly the same observations used by the
   authors.
2. **Optional database update.** The package also provides
   `prepare_data.R`. When sourced in RStudio, it will download or import the
   required series from FRED, the Federal Reserve Board, the Federal Reserve
   Bank of Philadelphia's Survey of Professional Forecasters, the Federal
   Reserve Bank of New York, and the FRB/US data package, and then apply the
   documented transformations and release-date masks. The resulting database
   will be written to `data/updated/`; the frozen paper vintage will never be
   overwritten. No API key is required.

The updated database may not be numerically identical to the paper vintage,
because providers revise historical observations and occasionally replace
downloadable files. It permits users to extend or update the analysis, whereas
`data/input/` remains the replication record. Published comparison series used
only in figures will be stored separately in `data/external/`, with their
sources and any applicable redistribution conditions documented.

The variables, transformations, sources, and contents of the frozen files are
described in [`data/README.md`](data/README.md).

## License

The authors' code is released under the [MIT License](LICENSE). Third-party
data and comparison series remain subject to the terms of their original
providers and are not relicensed by this package.

## Package status

The retained model, all-R database construction, baseline and robustness
estimations, inference calculations, validation exercises, and output builders
are contained in this package. The original-submission package in
`GGNR_Codes/` is left unchanged. The default `main.R` workflow has been run
from the frozen data and archived numerical results to reproduce the complete
set of manuscript and Online Appendix outputs.
