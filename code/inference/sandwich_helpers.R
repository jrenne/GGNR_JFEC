# Matrix calculations shared by empirical and small-model inference.
# H and score covariances use sums, not sample means.
bartlett_score_covariance <- function(scores, lags = 12L) {
  if (any(!is.finite(scores)) || lags < 0L || lags >= nrow(scores)) {
    stop("Invalid scores or HAC bandwidth.")
  }
  centered <- sweep(scores, 2, colMeans(scores), "-")
  value <- crossprod(centered)
  if (lags > 0L) for (lag in seq_len(lags)) {
    cross <- crossprod(centered[(lag + 1L):nrow(centered), , drop = FALSE],
                       centered[seq_len(nrow(centered) - lag), , drop = FALSE])
    value <- value + (1 - lag / (lags + 1)) * (cross + t(cross))
  }
  (value + t(value)) / 2
}

signed_hessian_inverse <- function(H) {
  # Diagonal scaling improves inversion accuracy without changing any sign or
  # removing a direction. The full inverse is used even when H is indefinite.
  scale <- sqrt(abs(diag(H)))
  if (any(!is.finite(H)) || any(!is.finite(scale)) || any(scale == 0)) {
    stop("The Hessian is unavailable or cannot be scaled.")
  }
  scaled <- H / outer(scale, scale)
  decomposition <- eigen((scaled + t(scaled)) / 2, symmetric = TRUE)
  if (min(abs(decomposition$values)) <= 1e-12 * max(abs(decomposition$values))) {
    stop("The local Hessian is numerically singular; no directions are discarded.")
  }
  inverse <- (decomposition$vectors %*%
    diag(1 / decomposition$values, nrow = nrow(H)) %*%
    t(decomposition$vectors)) / outer(scale, scale)
  dimnames(inverse) <- dimnames(H)
  if (max(abs(H %*% inverse - diag(nrow(H)))) > 1e-5) {
    stop("The signed Hessian inverse failed its residual check.")
  }
  inverse
}

numerical_local_derivatives <- function(theta, evaluate, step_scale = 1e-4,
                                        cores = 4L) {
  base <- evaluate(theta)
  p <- length(theta)
  step <- step_scale * pmax(1, abs(theta))
  map <- function(indices, fun) {
    if (.Platform$OS.type == "unix" && cores > 1L) {
      parallel::mclapply(indices, fun, mc.cores = cores)
    } else lapply(indices, fun)
  }
  diagonal <- map(seq_len(p), function(i) {
    plus <- minus <- theta
    plus[i] <- plus[i] + step[i]
    minus[i] <- minus[i] - step[i]
    fp <- evaluate(plus)
    fm <- evaluate(minus)
    list(curvature = -sum(fp - 2 * base + fm) / step[i]^2,
         score = (fp - fm) / (2 * step[i]))
  })
  H <- matrix(0, p, p, dimnames = list(names(theta), names(theta)))
  diag(H) <- vapply(diagonal, `[[`, numeric(1), "curvature")
  scores <- do.call(cbind, lapply(diagonal, `[[`, "score"))
  colnames(scores) <- names(theta)
  pairs <- which(upper.tri(H), arr.ind = TRUE)
  H[pairs] <- unlist(map(seq_len(nrow(pairs)), function(k) {
    i <- pairs[k, 1]; j <- pairs[k, 2]
    pp <- pm <- mp <- mm <- theta
    pp[i] <- pm[i] <- theta[i] + step[i]
    mp[i] <- mm[i] <- theta[i] - step[i]
    pp[j] <- mp[j] <- theta[j] + step[j]
    pm[j] <- mm[j] <- theta[j] - step[j]
    -sum(evaluate(pp) - evaluate(pm) - evaluate(mp) + evaluate(mm)) /
      (4 * step[i] * step[j])
  }))
  H[lower.tri(H)] <- t(H)[lower.tri(H)]
  if (any(!is.finite(H)) || any(!is.finite(scores))) {
    stop("Incomplete derivatives; existing inference files are not replaced.")
  }
  list(H = H, scores = scores, step = step)
}

transformation_jacobian <- function(theta, natural, step_scale = 1e-5) {
  value <- natural(theta)
  J <- matrix(0, length(value), length(theta),
              dimnames = list(names(value), names(theta)))
  for (j in seq_along(theta)) {
    h <- step_scale * max(1, abs(theta[j]))
    plus <- minus <- theta
    plus[j] <- plus[j] + h
    minus[j] <- minus[j] - h
    J[, j] <- (natural(plus) - natural(minus)) / (2 * h)
  }
  J
}

assemble_hessian_sandwich <- function(theta, H, scores, J, natural, lags = 12L) {
  B <- signed_hessian_inverse(H)
  S <- bartlett_score_covariance(scores, lags)
  sandwich <- function(meat) {
    V <- B %*% meat %*% t(B)
    (V + t(V)) / 2
  }
  white_working <- sandwich(crossprod(scores))
  hac_working <- sandwich(S)
  white_natural <- J %*% white_working %*% t(J)
  hac_natural <- J %*% hac_working %*% t(J)
  if (any(diag(white_natural) < -1e-12) || any(diag(hac_natural) < -1e-12)) {
    stop("Nonpositive variance; inference is not silently clipped.")
  }
  list(theta = theta, H = H, scores = scores, jacobian = J, J = J,
       bread = B, hac_meat = S, lag_count = lags,
       white_working = white_working, white_natural = white_natural,
       hac_working = hac_working, hac_natural = hac_natural,
       natural_estimate = natural,
       summary = data.frame(parameter = names(natural), estimate = unname(natural),
         qmle_white_se = sqrt(pmax(diag(white_natural), 0)),
         hac_sandwich_se = sqrt(pmax(diag(hac_natural), 0))),
       method = "Hessian sandwich; four empirical targets and measurement calibrations held fixed",
       active_branches_held_fixed = TRUE,
       empirical_moment_targets_held_fixed = TRUE)
}

write_hessian_inference <- function(result, directory) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  # Serialize the complete object only after every derivative and variance check.
  saveRDS(result, file.path(directory, "hac_sandwich_covariance.rds"))
  write.csv(result$summary, file.path(directory, "parameter_inference.csv"), row.names = FALSE)
  write.csv(result$hac_working, file.path(directory, "hac_sandwich_covariance_working.csv"))
  write.csv(result$hac_natural, file.path(directory, "hac_sandwich_covariance_natural.csv"))
  write.csv(result$H, file.path(directory, "hessian.csv"))
  write.csv(result$scores, file.path(directory, "monthly_scores.csv"), row.names = FALSE)
}
