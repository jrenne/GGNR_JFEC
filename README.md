# Replication: The Shadow-Rate Model: Let's Make it Real

Adam Golinski, Sophie Guilloux-Nefussi, and Jean-Paul Renne. JFEC second revision.

The application uses February 1992–August 2026 (415 months), 30 free parameters,
and a PTR measurement-error standard deviation fixed at 20 annualized basis
points. Four empirical moments determine two means and two innovation scales.
The economic transitory inflation shock is retained; there is no additional
CPI measurement disturbance or PCE–CPI adjustment.

## Run the published results

Extract the package, open `GGNR_JFEC.Rproj` in RStudio, and install:

```r
install.packages(c("Rcpp", "RcppArmadillo", "ggplot2"))
```

From the package root, open `main.R` and click **Source**, or run:

```r
source("main.R")
```

The terminal equivalent is `Rscript main.R`. The default run rebuilds all
**eight tables and 17 figures** and checks their numerical consistency. It
uses the frozen data, supplied estimates and numerical archives; it does not
download data, re-estimate parameters or generate new simulations. Every
figure is redrawn by code. Results appear under `outputs/`; the final message
is `Replication completed in ... seconds; results are in outputs/.`

See [the output index](docs/OUTPUTS.md) for the correspondence with the paper
and Online Appendix, and [REPLICATION_CHECKS.txt](REPLICATION_CHECKS.txt) for
verified runtimes and checks.

## Requirements

- R, the three packages above, and a C++ compiler configured for R. Rcpp
  compiles the model helpers on first use. On macOS use Apple's Command Line
  Tools; on Windows use the Rtools release matching R; on Linux use a working
  R development toolchain.
- PDF export needs working Cairo support in R or Poppler's `pdftocairo` on
  `PATH`. If necessary, set `GGNR_PDF_CONVERTER` to the executable's full path.
  The exported PDFs have embedded fonts.
- Optional: `nloptr` for BOBYQA estimation, and `readxl` for rebuilding source
  data or comparison series. Neither is needed for default replication.

The package was tested with R 4.5.1 on Apple Silicon macOS. Windows and Linux
were not independently tested. After installation, default replication is
offline and needs no API key, LaTeX, Python or MATLAB.

## Files

| Location | Contents |
| --- | --- |
| `main.R` | Default replication and optional calculation switches |
| `data/input/` | Frozen source data |
| `data/external/` | Frozen external comparison series |
| `estimates/baseline_r2.rds` | Published empirical estimate |
| `estimates/start_values.rds` | Same estimate as a warm start for optional estimation |
| `estimates/liquid_tips_only.rds` | Liquid-TIPS sensitivity estimate |
| `code/model/` | Model, pricing and filters |
| `code/inference/` | Parameter covariance and uncertainty propagation |
| `code/validation/` | Simulations and numerical checks |
| `code/outputs/` | Figure and table builders |
| `outputs/diagnostics/` | Numerical archives used by the default build |

The omit-2020 figure keeps baseline parameters and all four moment targets
fixed, removing only that year's measurement updates. The liquid-TIPS
exercise is a separate parameter re-estimation. Method and data details are
in [the inference notes](docs/INFERENCE_AND_DESIGN.md) and
[the data documentation](data/README.md).

## Optional fresh calculations

Use a separate package copy for fresh estimation or changing the model. Edit
switches **inside `main.R`**; setting them in the console first has no effect.
All calculation switches default to `FALSE`; `run_checks` is `TRUE`.

| Switch | Calculation |
| --- | --- |
| `run_estimation` | Estimate the empirical model and refresh dependent results |
| `run_refinement` | Refine all free coordinates and refresh dependent results |
| `run_robustness_estimation` | Optional alternative-sample re-estimations |
| `run_inference` | Recompute the Hessian and score covariance |
| `run_uncertainty_exercise` | Regenerate nonlinear state and premium bands |
| `run_filtering_exercise` | Generate 200 known-parameter filtering samples |
| `run_parameter_exercise` | Resume the two 500-sample parameter experiments |
| `run_pricing_exercise` | Regenerate high-precision pricing benchmarks |
| `run_mundell_exercise` | Regenerate conditional-correlation simulations |

Fresh estimation and refinement automatically refresh dependent empirical
inference, filtering, bands and conditional correlations. Full estimation
also re-estimates the robustness variants. These can take longer than the
short standalone uncertainty exercise. A finite optimizer exit message does
not replace convergence and Hessian checks.

### Empirical covariance and nonlinear bands

To recompute into separate directories:

```r
Sys.setenv(GGNR_INFERENCE_DIR = "outputs/diagnostics/inference_fresh",
           GGNR_HESSIAN_CORES = "4")
source("code/inference/hessian_sandwich.R")
Sys.unsetenv(c("GGNR_INFERENCE_DIR", "GGNR_HESSIAN_CORES"))

Sys.setenv(GGNR_UNCERTAINTY_DIR = "outputs/diagnostics/uncertainty_fresh",
           GGNR_UNCERTAINTY_DRAWS = "5000", GGNR_UNCERTAINTY_CORES = "4")
source("code/inference/nonlinear_uncertainty.R")
Sys.unsetenv(c("GGNR_UNCERTAINTY_DIR", "GGNR_UNCERTAINTY_DRAWS",
              "GGNR_UNCERTAINTY_CORES"))
```

The second block uses the supplied covariance. To use a fresh covariance,
keep `GGNR_INFERENCE_DIR` set to its directory until that block is complete.
`GGNR_HESSIAN_STEP` defaults to 0.001 in balanced coordinates;
`GGNR_HAC_LAGS` defaults to 12. The band simulation uses five independent
batches, 5,000 admissible joint parameter draws in total, and seed 20261126
for the first batch. The draw count must be a multiple of five and at least
100. `GGNR_UNCERTAINTY_SEED` selects a different first-batch seed.

The empirical figures display central 68% and 90% bands (16th–84th and
5th–95th percentiles). All free parameters and their cross-covariances are retained. Each accepted
draw reruns the original nonlinear filter, draws conditional states and
reprices bonds. Only draws violating model admissibility are rejected;
other numerical failures stop the calculation. The output reports proposals,
rejection reasons and batch-based Monte Carlo precision. These are pointwise,
covariance-based uncertainty bands conditional on the moment targets,
calibrations and admissibility, without an established finite-sample coverage
guarantee. They are not a bootstrap.

### Parameter Monte Carlo: 500 reference + 500 higher-ZLB samples

The supplied archive already contains all published replications. To generate
fresh samples rather than reuse them:

```r
Sys.setenv(GGNR_MC_DIR = "outputs/diagnostics/parameter_coverage_fresh",
           GGNR_MC_CASE = "all", GGNR_MC_REPS = "500", GGNR_MC_CORES = "4")
source("code/validation/parameter_coverage.R")
Sys.unsetenv(c("GGNR_MC_DIR", "GGNR_MC_CASE", "GGNR_MC_REPS", "GGNR_MC_CORES"))
```

The runner saves completed samples and resumes them when called with the same
directory and design. Select `reference_15` or `higher_zlb_50` instead of `all`
for one design. The illustrative five-parameter model is independent of the
empirical estimate; no bias correction is applied. A fresh 500-sample design
took about three minutes on four processes. To independently re-estimate
three archived samples without changing the archive:

```r
Sys.setenv(GGNR_MC_REPLAY = "true")
source("code/validation/parameter_coverage_alignment.R")
Sys.unsetenv("GGNR_MC_REPLAY")
```

### Estimation, filtering and data updates

`estimate_model.R` writes `estimates/reestimated.rds`; `refine_model.R` writes
`estimates/refined.rds`, preserving the published estimate. The refinement
budget is controlled by `GGNR_REFINE_SECONDS` and `GGNR_REFINE_MAXEVAL`.
`GGNR_START_FILE` selects a different estimate for subsequent calculations;
refresh all dependent numerical archives before building outputs with it.

For fresh known-parameter filtering, set `GGNR_FILTER_OUTPUT` to a new
directory, `GGNR_FILTER_REPS` to 200 and `GGNR_FILTER_CORES` to 4, then source
`code/validation/filtering_accuracy.R`. Unset these variables afterwards.
Optional data reconstruction uses `prepare_data.R`, which writes a new vintage
to `data/updated/` and preserves the frozen data. See `data/README.md`.

## Troubleshooting

- Run from the package root, containing `main.R`, `code/` and `data/`.
- For compilation errors, check the R-compatible C++ compiler and dependencies.
- For PDF errors, check Cairo or `Sys.which("pdftocairo")`.
- For cache mismatches, restore the supplied estimate/data/code or refresh the
  dependent calculation. Disabling `run_checks` does not disable provenance.
- Use a clean R session and unset optional `GGNR_*` variables before replaying
  the published results. Unix supports parallel execution; Windows falls back
  to serial execution for the R multiprocessing stages.

Authors' code uses the MIT [License](LICENSE). Third-party data retain their
providers' terms; sources are documented in `data/README.md`.
