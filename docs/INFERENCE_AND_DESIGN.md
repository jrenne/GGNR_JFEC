# Inference and simulation design

These notes document the calculations behind the supplied numerical results.
See [the main README](../README.md) for dependencies, execution examples, and
output-directory options, and [the output index](OUTPUTS.md) for paper items.

## Empirical Hessian and local uncertainty bands

`code/inference/hessian_sandwich.R` separately computes the negative Hessian
H and monthly score contributions. Both standard-error columns use the full
inverse of the positive-definite negative Hessian H. The White column uses the contemporaneous score outer
product; Hessian--HAC uses a Bartlett covariance of centered scores with 12
monthly lags. All matrices use sum normalization: V = H^-1 S H^-T, without
an additional factor of sample size. The full transformation Jacobian maps
30 free coordinates to 34 reported economic quantities, including innovation
scales determined by the four fixed empirical moment targets.

For numerical derivatives, the current nominal-rate formula is kept unchanged
at the evaluation point: below zero it stays zero, while above zero it follows
the perturbed shadow rate. All states and prices are recomputed. Optimization
retains the original max specification. The negative Hessian is evaluated in a
nonsingular linear rescaling of the 30 free coordinates and transformed back
by the full chain rule. It is positive definite at the supplied estimate. No
eigenvalues are clipped or directions removed. The default balanced-coordinate
step is 0.001, with an automatic check at half that step; override
`GGNR_HESSIAN_STEP` and `GGNR_INFERENCE_DIR` for a separate calculation.

The two calibrated measurement-error quantities are also tabulated, with zero
conditional variance: additional CPI measurement error is zero and PTR SD is
0.20 annualized percentage point. Four empirical moment targets determine the
two unconditional means and the innovation scales of m and u. Their estimation
uncertainty and uncertainty about the calibrations are not included in these
conditional standard errors. The full Jacobian retains the dependence of the
moment-derived scales on every free parameter.

State and premium bands use nonlinear propagation of the full Hessian--HAC
covariance. Five independently seeded batches provide 5,000 admissible joint
Gaussian parameter draws in the estimation coordinates. Each draw reruns the
original nonlinear IEKF, draws the full state from its conditional Gaussian
filtering distribution at each date, and recomputes prices and premia.
The pointwise 68% and 95% bands are empirical quantiles; the solid curve is
the unchanged plug-in estimate and is not a simulated median.

The distribution is restricted to the model's admissible domain: covariance
stationarity and positive innovation variances consistent with the four fixed
moment targets. Rejected proposals and their reasons are recorded. Numerical
errors other than these explicit domain violations stop the calculation.
The Gaussian parameter covariance is not diagonalized by deleting directions:
all positive eigenvalues and all cross-covariances are retained. Only tiny
floating-point negative zeros in conditional *state* covariances are removed.

The complete simulation is implemented in `code/inference/nonlinear_uncertainty.R`.
The default build uses its compact archive and replays three stored draws
independently in `code/validation/uncertainty_alignment.R`. Fresh generation
records batch-based Monte Carlo precision as well as admissibility counts.
The calculation preserves the empirical estimate and parameter standard errors;
it is not a bootstrap or a nominal finite-sample coverage guarantee.

## Small-model parameter-coverage experiment

The independent AR(1) factors have persistences 0.97 and 0.95 and innovation
scales 0.20 and 0.15 annualized percentage point. Samples have 695 months,
known means, and identical physical and risk-neutral dynamics. The five
estimated coordinates are two logit persistences and logs of two innovation
scales and a measurement-error scale. Observations are nominal and real
yields at 12 and 60 months and inflation. The measurement profile is fixed
in this illustrative design: 0.0663 pp for nominal yields, 0.10 pp for
real yields, and 2.183 pp for inflation. The common scale is re-estimated.

The reference targets a 15% stationary probability of negative shadow rates;
the sensitivity changes only the real-rate mean to target 50%. The inflation
mean is fixed at 2 annualized percent. Realized frequencies are 15.0% and
49.8%. Coverage is for nominal 95% Hessian--HAC Wald intervals in logit/log
coordinates. The reported reference rates are 89.0/93.8/89.8/92.0/94.0%; the
sensitivity rates are 82.6/93.4/70.0/89.8/94.0%. No bias correction is applied.
The generator and estimator share the pricing approximation, so this
experiment does not measure pricing misspecification or establish inference
accuracy for all 30 empirical parameters.

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

## Omitting 2020

The published omit-2020 check holds all baseline parameters and all four moment
targets fixed. It removes the 2020 measurement updates while continuing the
same transition law. It isolates filtering sensitivity to that year, rather
than claiming robustness under parameter re-estimation. Optional alternative-
sample re-estimation remains available through estimate_robustness.R and is
saved separately; its optimizer message is not a convergence certificate.
