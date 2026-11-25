# Output index

Run `main.R` from the package root to regenerate the items below with the
supplied estimates and numerical archives. Paths in the tables are relative
to `outputs/tables/` or `outputs/figures/`, as indicated. OA denotes the
Online Appendix. Other intermediate outputs may also be written.

## Tables: `outputs/tables/`

| Paper item | Contents | File |
| --- | --- | --- |
| Table 1 | Parameter estimates and standard errors | `table_parameters.tex` |
| Table 2 | Fit of expectations and yields | `table_fit.tex` |
| Table OA.1 | Pricing-approximation accuracy | `table_pricing_accuracy.tex` |
| Table OA.2 | Known-parameter filtering accuracy | `table_filtering_accuracy.tex` |
| Table OA.3 | Parameter coverage: 500 reference samples | `table_parameter_mc_reference.tex` |
| Table OA.4 | Parameter coverage: 500 higher-ZLB samples | `table_parameter_mc_sensitivity.tex` |
| Table OA.5 | Inflation-risk-premium discrepancy decomposition | `table_irp_decomposition.tex` |
| Table OA.6 | Liquid-TIPS robustness | `table_liquid_tips_robustness.tex` |

## Figures: `outputs/figures/`

| Paper item | Contents | File |
| --- | --- | --- |
| Figure 1 | Filtered states, uncertainty, and benchmarks | `fig_states.pdf` |
| Figure 2 | Nominal and real term premia and inflation risk premia | `risk_premia.pdf` |
| Figure 3 | Yield responses below and above the lower bound | `irf_yields_shadow_rate.pdf` |
| Figure 4 | Yield responses at the same initial shadow rate | `irf_yields_s0_rpi.pdf` |
| Figure 5 | Conditional probability of a negative shadow rate | `iso_probability.pdf` |
| Figure 6 | Conditional Mundell–Tobin correlations | `mundell_tobin_effect.pdf` |
| Figure OA.1 | State-by-state pricing-approximation errors | `pricing_approximation_errors.pdf` |
| Figure OA.2 | Inflation and survey fit | `fig_macro_fit.pdf` |
| Figure OA.3 | Nominal yield fit | `fig_nominal_yield_fit.pdf` |
| Figure OA.4 | Real yield fit | `fig_real_yield_fit.pdf` |
| Figure OA.5 | Real-rate trend omitting 2020 observations, with parameters fixed | `rstar_omit2020.pdf` |
| Figure OA.6 | Historical real-yield sensitivity | `pre2004_real_yield_sensitivity.pdf` |
| Figure OA.7 | Persistent-state uncertainty | `state_uncertainty_bands.pdf` |
| Figure OA.8 | Ten-year premium uncertainty | `premium_uncertainty_bands.pdf` |
| Figure OA.9 | Alternative shadow-rate benchmarks | `shadow_rate_benchmarks.pdf` |
| Figure OA.10 | Observed and shadow nominal short rates | `fig_shadow_rate.pdf` |
| Figure OA.11 | Yield responses to a negative innovation | `irf_yields_s25_rpi.pdf` |

## Numerical archives: `outputs/diagnostics/`

| Directory | Contents |
| --- | --- |
| `baseline/` | Empirical fit, filtered states, and model diagnostics |
| `inference/` | Full Hessian, scores, HAC covariance, parameter standard errors, and sensitivity audits |
| `uncertainty/` | Local delta-method bands, derivatives, and half-step comparison |
| `filtering/` | Known-state filtering experiment and summary |
| `parameter_coverage/` | Individual parameter Monte Carlo results, retries, and summaries for both 500-sample designs |
| `pricing/` | Monte Carlo pricing benchmarks and approximation-error summaries |
| `risk_premia/` | Model premia and comparison diagnostics |
| `omit_2020/`, `liquid_tips/` | Robustness paths and comparisons |
| `mundell_tobin/` | Archived simulated conditional correlations |

The default build redraws figures and rewrites tables. It reuses the archived
estimates and expensive numerical exercises. For fresh calculations, see
[the README](../README.md); new output directories keep experimental results
separate from the published archive. Methodological details are in
[INFERENCE_AND_DESIGN.md](INFERENCE_AND_DESIGN.md).
