# Source from the project root after inference has been computed.
# Ordinary OPG and OPG-HAC must share the same uncentered OPG bread.
ordinary <- readRDS("outputs/diagnostics/inference/sandwich_covariance.rds")
hac <- readRDS("outputs/diagnostics/inference/hac_sandwich_covariance.rds")
G <- crossprod(ordinary$scores)
centered <- sweep(ordinary$scores, 2, colMeans(ordinary$scores), "-")
S <- crossprod(centered)
for (lag in seq_len(hac$lag_count)) {
  lead <- centered[(lag + 1L):nrow(centered), , drop = FALSE]
  lagged <- centered[seq_len(nrow(centered) - lag), , drop = FALSE]
  cross <- crossprod(lead, lagged)
  S <- S + (1 - lag / (hac$lag_count + 1)) * (cross + t(cross))
}
stopifnot(isTRUE(all.equal(G, ordinary$meat)),
          isTRUE(all.equal(G, hac$opg)),
          isTRUE(all.equal(unname(hac$bread), unname(ordinary$opg_working))),
          isTRUE(all.equal(S, hac$hac_meat)),
          isTRUE(all.equal(unname(hac$hac_working),
                           unname(hac$bread %*% S %*% hac$bread))))
message("Shared raw-score OPG bread and centered Newey-West meat verified.")
