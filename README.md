# Replication package: The Shadow-Rate Model: Let's Make it Real

Adam Golinski, Sophie Guilloux-Nefussi, and Jean-Paul Renne  
JFEC second revision, November 2026

## Start here

Extract the archive and open `GGNR_JFEC.Rproj` in RStudio. Install the two
required R packages once, if necessary:

```r
install.packages(c("Rcpp", "ggplot2"), repos = "https://cloud.r-project.org")
```

Then open `main.R` and click **Source**, or run this in the R console:

```r
source("main.R")
```

Alternatively, run `Rscript main.R` in a terminal from the extracted
`GGNR_JFEC` directory. All commands in this README use that working directory;
`file.exists("main.R")` should return `TRUE` in R.

**The default run regenerates the eight tables and 17 figures used in the
paper and Online Appendix, and runs the numerical consistency checks.** It
uses the frozen data, supplied R2 parameter estimate, and archived simulation
results. It does not download data, estimate parameters, or generate new
Monte Carlo samples. Figures are drawn by the code rather than copied from
finished manuscript files. The last clean default run took 18.3 seconds on
the tested machine; compilation and other hardware can change this time.

Outputs are written to `outputs/tables/`, `outputs/figures/`, and
`outputs/diagnostics/`. The final console message is `Replication completed
in ... seconds; results are in outputs/.` See [the output index](docs/OUTPUTS.md)
for the correspondence with every paper table and figure.

## Software requirements

- R, with `Rcpp` and `ggplot2`, and a C++ compiler configured for R. Model
  helpers are compiled through Rcpp when a script first needs them.
- On macOS, install Apple's Command Line Tools if needed; see
  [R's macOS tools instructions](https://mac.r-project.org/tools/).
  On Windows, install the [Rtools version matching your R release](https://cran.r-project.org/bin/windows/Rtools/).
  On Linux, use a working R development toolchain with a C++ compiler and `make`.
- PDF export requires a working Cairo device in R or Poppler's `pdftocairo`
  executable on `PATH`. If it is installed elsewhere, set
  `GGNR_PDF_CONVERTER` to its full executable path before running the scripts.
  This provides embedded fonts without any machine-specific paths in the code.
- Optional packages: `nloptr` for empirical estimation or local refinement;
  `readxl` for database reconstruction or rebuilding the LW/HLW comparison
  series. Install them with `install.packages(c("nloptr", "readxl"))` if needed.
  Internet access is needed for package installation and optional data updates.

The tested environment was R 4.5.1 on Apple Silicon macOS, with Rcpp 1.1.0,
ggplot2 4.0.1, nloptr 2.2.1, and readxl 1.4.5. These are recorded versions,
not enforced version pins. Windows and Linux were not independently tested.
No LaTeX installation is needed to generate the `.tex` tables or PDF figures;
no Python or MATLAB is required. After setup, default replication is offline
and requires no API key.

## What is supplied

| Location | Contents |
| --- | --- |
| `data/input/paper_data.rds` | Frozen paper data, October 1968–August 2026 |
| `data/external/` | Archived comparison series and source workbooks |
| `estimates/baseline_r2.rds` | Published R2 empirical parameter estimate |
| `estimates/baseline.rds` | Preserved first-revision empirical estimate |
| `estimates/start_values.rds` | Starting values for fresh estimation |
| `estimates/omit_2020.rds`, `estimates/liquid_tips_only.rds` | Supplied robustness estimates |
| `code/model/` | Pricing, state-space model, and filters, including Rcpp helpers |
| `code/inference/`, `code/validation/` | Uncertainty calculations, experiments, and audits |
| `code/outputs/` | Table and figure builders |
| `outputs/diagnostics/` | Archived numerical exercises and generated diagnostics |

The default scripts use `baseline_r2.rds`. The supplied robustness variants
are separate estimates; their filtered paths and comparisons are rebuilt
against the R2 baseline. The individual simulations, seeds, optimizer retries,
and summaries are included in the parameter-coverage archive.

## Optional fresh calculations

Work in a separate copy of the package when changing estimates or the model.
The driver rebuilds outputs in place; fresh empirical exercises replace their
default diagnostic files unless a different output directory is selected.
The examples below keep fresh exercise results in separate directories.
They do not change the published tables until you explicitly select them for
an output build.

### Driver switches

Edit the switches **inside `main.R`** before sourcing it. Assigning a switch
in the console first has no effect because `main.R` sets its own values.
All switches below default to `FALSE`; `run_checks` defaults to `TRUE`.

| Switch | Calculation requested |
| --- | --- |
| `run_estimation` | Fresh empirical estimation and robustness estimation; refresh dependent empirical exercises |
| `run_polishing` | Local refinement of all 32 empirical coordinates; refresh dependent empirical exercises |
| `run_robustness_estimation` | Re-estimate the omit-2020 and liquid-TIPS variants |
| `run_inference` | Recompute the empirical Hessian and score-based covariance |
| `run_filtering_exercise` | Generate 200 known-parameter samples and compare EKF, IEKF, and UKF |
| `run_uncertainty_exercise` | Recompute local delta-method state and premium bands |
| `run_parameter_exercise` | Resume the parameter-coverage archive, targeting 500 samples per design |
| `run_pricing_exercise` | Recompute the high-precision Monte Carlo pricing benchmarks |
| `run_mundell_exercise` | Recompute simulated conditional Mundell–Tobin correlations |

`run_estimation` and `run_polishing` automatically refresh inference, filtering,
uncertainty bands, and conditional correlations. Only `run_estimation` also
automatically re-estimates the robustness variants. Default cache checks
compare fingerprints of the selected estimate, frozen data, and model code
and stop on a mismatch. The small-model parameter experiment has its own
fixed design and does not depend on the empirical estimate.

### Parameter coverage: 500 reference + 500 higher-ZLB samples

Setting `run_parameter_exercise = TRUE` reuses completed samples in the
supplied archive. **To generate new samples, select a new output directory:**

```r
Sys.setenv(GGNR_MC_DIR = "outputs/diagnostics/parameter_coverage_fresh",
           GGNR_MC_CASE = "all", GGNR_MC_REPS = "500", GGNR_MC_CORES = "4")
source("code/validation/parameter_coverage.R")
Sys.unsetenv(c("GGNR_MC_DIR", "GGNR_MC_CASE", "GGNR_MC_REPS", "GGNR_MC_CORES"))
```

This generates 500 samples in the reference design and 500 in the appendix
sensitivity with a higher frequency of the zero lower bound. Select
`reference_15` or `higher_zlb_50` instead of `all` to run just one design.
The runner saves completed batches and resumes them when called again with
the same directory and settings. Restore environment variables as shown
before running the published table builders, which expect both full designs.

Pricing, filtering, and generation use C++ through Rcpp; five-coordinate
estimation uses R's `optim`. A fresh 500-sample design took about three minutes
on four processes (about six minutes for both designs). See
[the design and inference notes](docs/INFERENCE_AND_DESIGN.md) for the data
generator, intervals, optimizer retries, and interpretation of coverage.

### Empirical Hessian and uncertainty bands

Recompute the Hessian and score covariance at the supplied estimate:

```r
Sys.setenv(GGNR_INFERENCE_DIR = "outputs/diagnostics/inference_fresh",
           GGNR_HESSIAN_CORES = "4")
source("code/inference/hessian_sandwich.R")
Sys.unsetenv(c("GGNR_INFERENCE_DIR", "GGNR_HESSIAN_CORES"))
```

This took about 98 seconds on four processes. `GGNR_HESSIAN_STEP` selects the
relative derivative step (default `1e-4`); `GGNR_HAC_LAGS` selects the score
covariance lag length (default 12). Recompute bands using the supplied
covariance, including the automatic half-step derivative check:

```r
Sys.setenv(GGNR_UNCERTAINTY_DIR = "outputs/diagnostics/uncertainty_fresh",
           GGNR_DELTA_CORES = "4")
source("code/inference/delta_uncertainty.R")
Sys.unsetenv(c("GGNR_UNCERTAINTY_DIR", "GGNR_DELTA_CORES"))
```

The band calculation took about 13 seconds including its check. To use the
fresh covariance instead, set `GGNR_INFERENCE_DIR` to its directory before
this block and unset it afterwards. `GGNR_DELTA_STEP` defaults to `1e-5`.
These are local numerical inference calculations; the Hessian's curvature
and the derivative conventions are documented in
[the technical notes](docs/INFERENCE_AND_DESIGN.md).

### Known-parameter filtering comparison

```r
Sys.setenv(GGNR_FILTER_OUTPUT = "outputs/diagnostics/filtering_fresh",
           GGNR_FILTER_REPS = "200", GGNR_FILTER_CORES = "4")
source("code/validation/filtering_accuracy.R")
Sys.unsetenv(c("GGNR_FILTER_OUTPUT", "GGNR_FILTER_REPS", "GGNR_FILTER_CORES"))
```

This exercise compares EKF, IEKF, and UKF with known latent states under the
empirical six-state specification; parameters are not re-estimated. It took
about 52 seconds on four processes. All timings are illustrative. On Windows,
these parallel exercises fall back to serial execution and can take longer.
The full pricing benchmark and conditional-correlation simulations are
separate optional calculations; the default build reuses their archives.

### Fresh empirical estimation or new data

`estimate_model.R` writes `estimates/reestimated.rds`; `polish_model.R` writes
`estimates/polished.rds`. Both preserve `baseline.rds` and `baseline_r2.rds`.
The latter script uses the original likelihood criterion, with a default
180-second budget configurable through `GGNR_POLISH_SECONDS` and
`GGNR_POLISH_MAXEVAL`. Full estimation can take longer and may reach a
different local solution. `GGNR_START_FILE` selects an alternative estimate
for downstream empirical scripts. Changing it requires recomputing the
dependent empirical caches before rebuilding outputs.

`estimate_robustness.R` replaces `omit_2020.rds` and `liquid_tips_only.rds` in
the working copy. Use a separate package copy for a full fresh estimation.
See [the data documentation](data/README.md) for sources and optional updating.
`prepare_data.R` writes a new vintage to `data/updated/updated_data.rds`;
it leaves the frozen input and default analysis unchanged. An analysis using
that new vintage requires explicitly changing the data input and re-estimating.

## Checks and troubleshooting

The default driver checks cache identity, covariance algebra and normalization,
delta-method bands, Monte Carlo statistics and denominators, pricing timing,
R/C++ agreement, yield derivatives, and survey integration. The recorded
validation is in [REPLICATION_CHECKS.txt](REPLICATION_CHECKS.txt). To additionally
re-estimate three archived parameter samples and compare their stored results:

```r
Sys.setenv(GGNR_MC_REPLAY = "true")
source("code/validation/parameter_coverage_alignment.R")
Sys.unsetenv("GGNR_MC_REPLAY")
```

- **File not found:** open the extracted project and run from its root,
  containing `main.R`, `data/`, `code/`, and `estimates/`.
- **Package or compilation error:** install the packages and the R-compatible
  compiler described above; Rcpp needs that compiler even for default replication.
- **PDF/font error:** check the Cairo device or `Sys.which("pdftocairo")`;
  alternatively set `GGNR_PDF_CONVERTER` to a working executable.
- **Cache mismatch:** restore the supplied estimate/data/model, or refresh the
  dependent exercises for the intended specification. Disabling `run_checks`
  does not disable the provenance checks.
- **Unexpected exercise or output selection:** start a clean R session and
  check for `GGNR_*` environment variables left from optional runs.

## License

Authors' code uses the MIT [License](LICENSE). Third-party data retain their
providers' terms; their sources and construction are listed in
[data/README.md](data/README.md).
