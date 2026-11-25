baseline_r2.rds: published empirical estimate, February 1992–August 2026;
30 free coordinates, PTR measurement-error SD fixed at 20 annualized bp.
start_values.rds: the same estimate, supplied as a warm start for optional estimation.
liquid_tips_only.rds: separate liquid-TIPS sensitivity estimate.

Default replication preserves these files. estimate_model.R writes reestimated.rds;
refine_model.R writes refined.rds. The published omit-2020 exercise holds the
baseline parameters and four moment targets fixed and omits only the 2020 updates.
Optional robustness re-estimation writes its separate estimates in the working
copy; use a separate copy when performing new estimation or updating data.
