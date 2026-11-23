# Independent evaluation of the printed, corrected timing convention.

source("code/model/pricing_and_filter_helpers.R")

set.seed(146)

rows <- list()

for (K in c(3L, 6L)) {
    P <- diag(seq(0.94, 0.98, length.out = K))
    P[1, 2] <- 0.015
    P[2, 3] <- -0.01
    mu <- rnorm(K, sd = 1e-05)
    V <- diag(seq(1e-08, 3e-08, length.out = K))
    V[1, 2] <- V[2, 1] <- 2e-09
    b <- c(1, 1, 1, rep(0, K - 3))
    c <- c(0, 1, 0.5, rep(0, K - 3))
    a <- 1e-04
    gamma0 <- 2e-04
    x <- matrix(rnorm(K * 9, sd = 0.005), K, 9)
    co <- affine_coefs_r(P, mu, V, a, b, gamma0, c)
    powers <- list(diag(K))
    for (j in 1:120) powers[[j + 1]] <- powers[[j]] %*% P
    Fn <- Fr <- matrix(0, 9, 120)
    Phi_prev <- matrix(0, K, K)
    err <- 0
    for (n in 1:120) {
        Phi_n <- Phi_prev + powers[[n]]
        bn <- t(powers[[n]]) %*% b
        cn <- t(powers[[n + 1]]) %*% c
        abar <- a + drop(t(b) %*% Phi_prev %*% mu)
        an <- abar - 0.5 * drop(t(b) %*% Phi_prev %*% V %*% t(Phi_prev) %*% b)
        cdot <- drop(t(c) %*% Phi_n %*% mu) + 0.5 * drop(t(c) %*% Phi_n %*% V %*% t(Phi_n) %*% c)
        variance <- if (n == 1)
            0
        else sum(vapply(0:(n - 2), function(j) drop(t(b) %*% powers[[j + 1]] %*% V %*% t(powers[[j +
            1]]) %*% b), numeric(1)))
        interaction <- drop(t(b) %*% Phi_prev %*% V %*% t(Phi_n) %*% c)
        z <- drop(an + t(x) %*% bn)
        zbar <- drop(abar + t(x) %*% bn)
        Fn[, n] <- if (n == 1)
            pmax(z, 0)
        else z * pnorm(z/sqrt(variance)) + sqrt(variance) * dnorm(z/sqrt(variance))
        Fr[, n] <- Fn[, n] - gamma0 - cdot - drop(t(x) %*% cn) + if (n == 1)
            0
        else pnorm(zbar/sqrt(variance)) * interaction
        err <- max(err, abs(an - co$A_X_for[n]), abs(abar - co$A_X_exp[n]), abs(bn - co$B_X_for[,
            n]), abs(cn - co$B_X_pi[, n]), abs(gamma0 + cdot - co$A_X_for_pi[n]), abs(t(b) %*% Phi_prev -
            co$B_X_cum[, n]), abs(t(c) %*% Phi_n - co$B_X_cum_pi[, n]))
        Phi_prev <- Phi_n
    }
    yn <- t(apply(Fn, 1, cumsum))/rep(1:120, each = 9)
    yr <- t(apply(Fr, 1, cumsum))/rep(1:120, each = 9)
    sd_nom <- sqrt(cumsum(diag(t(co$B_X_for) %*% V %*% co$B_X_for)))
    implemented <- y_fitting_r_pure(x, co$A_X_for, co$B_X_for, co$A_X_exp, 0, sd_nom, co$A_X_for_pi,
        co$B_X_pi, V, co$B_X_cum, co$B_X_cum_pi)
    pricing_error <- max(abs(implemented$yfit_all_n - yn), abs(implemented$yfit_all_r - yr))
    stopifnot(err < 1e-12, pricing_error < 1e-12)
    rows[[length(rows) + 1]] <- data.frame(dimension = K, coefficient_error = err, yield_error = pricing_error)
}

result <- do.call(rbind, rows)

print(result)

write.csv(result, "outputs/diagnostics/pricing/timing_alignment.csv", row.names = FALSE)

