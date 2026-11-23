baseline.rds: immutable empirical estimate from the submitted first revision.
baseline_r2.rds: published 32-coordinate empirical estimate used by the second revision.
start_values.rds: archived starting values for optional fresh estimation.
omit_2020.rds and liquid_tips_only.rds: previously computed robustness estimates;
R2 comparisons recompute their paths and differences against baseline_r2.rds.
The default driver reads baseline_r2.rds. Fresh estimation writes reestimated.rds.
Local refinement with polish_model.R writes estimates/polished.rds by default.
Fresh robustness estimation replaces omit_2020.rds and liquid_tips_only.rds
in the working copy; use a separate package copy for a full fresh estimation.
