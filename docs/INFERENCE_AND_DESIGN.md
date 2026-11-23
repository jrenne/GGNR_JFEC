# Inference and simulation design

These notes document the calculations behind the supplied numerical results.
See [the main README](../README.md) for dependencies, execution examples, and
output-directory options, and [the output index](OUTPUTS.md) for paper items.

## Empirical Hessian and local uncertainty bands

`code/inference/hessian_sandwich.R` separately computes the negative Hessian
H and monthly score contributions. Both standard-error columns use the full
signed inverse of H. The White column uses the contemporaneous score outer
product; Hessian--HAC uses a Bartlett covariance of centered scores with 12
monthly lags. All matrices use sum normalization: V = H^-1 S H^-T, without
an additional factor of sample size. The full transformation Jacobian maps
32 free coordinates to 34 reported economic quantities, including innovation
scales determined by the four fixed empirical moment targets.

For numerical derivatives, the branches of the known current short-rate
max term are held at their values at the evaluation point. Estimation and
ordinary filtering retain the original max. Future Gaussian pricing
terms are unchanged. The local Hessian has three negative eigenvalues; their
signed inverse is retained without clipping or deleting directions. These
are local numerical standard errors, not a proof of regular asymptotic
validity at regime boundaries. The archived half-step check and the
independent algebra audit document numerical stability. The default derivative
step is 1e-4 times max(1, abs(theta)); override `GGNR_HESSIAN_STEP` to recompute
at another step, and `GGNR_INFERENCE_DIR` to keep a separate result.

The primary state and premium bands use a first-order delta calculation,
combining conditional filtering variance with D V D'. Here D is the total
parameter derivative of each filtered state or premium, obtained by rerunning
the filter in all 32 coordinates. For premia it includes direct pricing effects.
The same complete Hessian--HAC covariance, including every cross-covariance,
is retained. No parameters are frozen or removed from the covariance.
The four empirical moment targets are fixed. Local branch conventions are
used only for derivatives. Normal pointwise bands are centered on the existing
point estimate; finite-sample coverage is not established for these bands.
The local approximation does not eliminate weak identification or the
nonregular-curvature qualification above.

`code/inference/delta_uncertainty.R` recomputes these bands in about 7 seconds
on four processes (approximately 13 seconds including the automatic half-step
check). `GGNR_DELTA_CORES` and `GGNR_DELTA_STEP` select the core
count and relative central-difference step (default 1e-5 times max(1, abs(theta))).
Set `GGNR_UNCERTAINTY_DIR` to a separate directory for a fresh result.
The archived half-step result is independently compared in
`code/validation/delta_alignment.R`, which also reconstructs the full
parameter variance through a spectral root and checks normal-band endpoints.
Half-step total standard errors differ by less than 0.002% for the five main
state and ten-year premium series.

The earlier nonlinear Gaussian-draw experiment remains available in
`code/inference/hamilton_uncertainty.R` and
`outputs/diagnostics/hamilton_reference/` as a diagnostic, not as the published
bands. It contains 1,000 accepted draws from 1,124 proposals. Such global
draws can leave the local neighborhood, including changing risk-neutral
dynamics from stable to explosive, and produce extreme tails even when
finite-horizon prices remain computable. No new stationarity restriction,
tail trimming, or draw rejection has been introduced to narrow the primary
bands. The delta calculation is explicitly local rather than a global
probability assessment of these nonlinear tails.

## Small-model parameter-coverage experiment

The independent AR(1) factors have persistences 0.97 and 0.95 and innovation
scales 0.20 and 0.15 annualized percentage point. Samples have 695 months,
known means, and identical physical and risk-neutral dynamics. The five
estimated coordinates are two logit persistences and logs of two innovation
scales and a measurement-error scale. Observations are nominal and real
yields at 12 and 60 months and inflation. The measurement profile is frozen
at submitted-estimate magnitudes: 0.0663 pp for nominal yields, 0.10 pp for
real yields, and 2.183 pp for inflation. The common scale is re-estimated.

The reference targets a 15% stationary probability of negative shadow rates;
the sensitivity changes only the real-rate mean to target 50%. The inflation
mean is fixed at 2 annualized percent. Realized frequencies are 15.0% and
49.8%. Coverage is for nominal 95% Hessian--HAC Wald intervals in logit/log
coordinates. The reported reference rates are 89.0/93.8/89.8/92.0/94.0%; the
sensitivity rates are 82.6/93.4/70.0/89.8/94.0%. No bias correction is applied.
The generator and estimator share the pricing approximation, so this
experiment does not measure pricing misspecification or establish inference
accuracy for all 32 empirical parameters.

The runner retains all samples, including failures and optimizer retries. An
initial limit of 800 likelihood evaluations is increased to 1,600 if required,
using the same simulated panel and initial parameter vector. Three fits in
the higher-ZLB design required this retry. All 1,000 published fits have
converged estimates and computable intervals. The archive contains five rows
per replication; the summary reports the full denominator for each parameter.

The exercise is resumable, with completed batches written to
`replications.csv`. Reusing an existing directory with unchanged settings
resumes its samples; use a new directory to run a fresh experiment. A replay
check can re-estimate three archived samples without replacing the archive.

## Filtering versus parameter coverage

The known-parameter filtering exercise in
`code/validation/filtering_accuracy.R` uses the empirical six-state model and
its observation-availability masks. It compares EKF, IEKF, and UKF on 200
simulated samples with known states and fixed parameters. This differs from
the two-factor parameter experiment, which re-estimates five coordinates in
each of 500 reference and 500 higher-ZLB samples. Both use the same pricing
approximation in generation and estimation; the separate pricing benchmark
addresses approximation error.
