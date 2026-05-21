# =============================================================================
# 07_mc_csardl_6_CSB.R
# CS-ARDL Monte Carlo - Cross-Section Bootstrap
#
# Method   : CSB: Cross-section bootstrap
# Reference : Goncalves and Perron (2014)
#
# Monte Carlo design (manuscript Section 4.1):
#   Grid          : 5 x 5, with N and T in {40, 50, 100, 150, 200} (25 cells)
#   Replications  : R = 500 per cell
#   Factor model  : single common factor, AR(1) coefficient rho_f = 0.6
#   Truncation lag: pT = floor(T^(1/3))
#
# Input    : none (data are simulated internally with a fixed seed)
# Output   : 07_mc_csardl_6_CSB_results.csv
#
# This script is self-contained and can be run on its own in RStudio.
# Run order: see 00_RUN_ALL.R
# =============================================================================

rm(list = ls())
cat("================================================================\n")
cat(" CS-ARDL  |  Script 6 of 6  |  METHOD: CSB\n")
cat("================================================================\n\n")

set.seed(2026)

# ============================================================
# 1.  CS ESTIMATOR
# ============================================================
EstCS <- function(Y, X, N, TT, PT) {
  ybar  <- colMeans(Y);  xbar <- colMeans(X)
  dxbar <- c(0, diff(xbar))
  nobs  <- (TT - 1) - PT
  phi_i <- b0_i <- b1_i <- gam_i <- numeric(N)
  for (i in 1:N) {
    yi <- Y[i, ];  xi <- X[i, ]
    DY <- YLAG <- X0 <- XLAG <- DX <- numeric(nobs)
    YB <- XB <- DB <- matrix(0, nobs, PT + 1)
    row <- 1
    for (tt in (PT + 2):TT) {
      DY[row] <- yi[tt] - yi[tt - 1];  YLAG[row] <- yi[tt - 1]
      X0[row] <- xi[tt];  XLAG[row] <- xi[tt - 1]
      DX[row] <- xi[tt] - xi[tt - 1]
      for (ell in 0:PT) {
        YB[row, ell + 1] <- ybar[tt - ell]
        XB[row, ell + 1] <- xbar[tt - ell]
        DB[row, ell + 1] <- dxbar[tt - ell]
      }
      row <- row + 1
    }
    Z <- cbind(YLAG, X0, XLAG, DX, 1, YB, XB, DB)
    b <- solve(crossprod(Z) + 1e-8 * diag(ncol(Z)), crossprod(Z, DY))
    phi_i[i] <- b[1]; b0_i[i] <- b[2]; b1_i[i] <- b[3]; gam_i[i] <- b[4]
  }
  list(phi = phi_i, b0 = b0_i, b1 = b1_i, gam = gam_i)
}

# ============================================================
# 2.  CSB ALGORITHM
# ============================================================
CSB <- function(Y, X, N, TT, PT, B_boot = 200) {
  res_full <- EstCS(Y, X, N, TT, PT)
  theta_i_full <- -res_full$b1 / res_full$phi
  w_full       <- res_full$phi^2 / sum(res_full$phi^2)
  theta_full   <- sum(w_full * theta_i_full)
  phi_full     <- mean(res_full$phi)

  theta_boot <- phi_boot <- numeric(B_boot)
  for (b in 1:B_boot) {
    idx <- sample(1:N, N, replace = TRUE)
    Y_b <- Y[idx, , drop = FALSE]
    X_b <- X[idx, , drop = FALSE]
    res_b <- tryCatch(EstCS(Y_b, X_b, N, TT, PT), error = function(e) NULL)
    if (!is.null(res_b)) {
      th_i <- -res_b$b1 / res_b$phi
      wb   <- res_b$phi^2 / sum(res_b$phi^2)
      theta_boot[b] <- sum(wb * th_i)
      phi_boot[b]   <- mean(res_b$phi)
    } else {
      theta_boot[b] <- theta_full
      phi_boot[b]   <- phi_full
    }
  }

  theta_csb <- 2 * theta_full - mean(theta_boot)
  phi_csb   <- 2 * phi_full   - mean(phi_boot)
  list(phi = phi_csb, theta = theta_csb)
}

# ============================================================
# 3.  DGP
# ============================================================
GenerateDGP <- function(N, TT, m, rho_f, rho_x, sig_f, sig_u, sig_v,
                         phi_true, beta0_true, beta1_true, gam_true,
                         gamma_y, gamma_x) {
  Tburn <- 50;  TT_total <- TT + Tburn
  f <- matrix(0, m, TT_total)
  for (tt in 2:TT_total)
    f[, tt] <- rho_f * f[, tt - 1] + rnorm(m) * sig_f
  Y <- X <- matrix(0, N, TT_total)
  for (i in 1:N) {
    for (tt in 2:TT_total) {
      X[i, tt] <- rho_x * X[i, tt - 1] +
                  sum(gamma_x[i, ] * f[, tt]) + rnorm(1) * sig_v
      rho_i <- 1 + phi_true[i]
      Y[i, tt] <- rho_i * Y[i, tt - 1] +
                  beta0_true[i] * X[i, tt] +
                  beta1_true[i] * X[i, tt - 1] +
                  gam_true[i] * (X[i, tt] - X[i, tt - 1]) +
                  sum(gamma_y[i, ] * f[, tt]) + rnorm(1) * sig_u
    }
  }
  list(Y = Y[, (Tburn + 1):TT_total], X = X[, (Tburn + 1):TT_total])
}

# ============================================================
# 4.  GRID & MC LOOP
# ============================================================
Nvec <- c(40, 50, 100, 150, 200)
Tvec <- c(40, 50, 100, 150, 200)
Grid <- expand.grid(N = Nvec, T = Tvec)
G <- nrow(Grid);  R <- 200
m <- 1;  rho_f <- 0.6;  rho_x <- 0.6
sig_f <- 1.0;  sig_u <- 1.0;  sig_v <- 1.0;  theta_true <- 1.0
results <- data.frame()

for (g in 1:G) {
  N <- Grid$N[g];  TT <- Grid$T[g]
  PT <- max(1, min(floor(TT^(1/3)), floor((TT - 1) / 4) - 1))
  phi_true   <- rep(c(-0.2, -0.3, -0.4, -0.5, -0.6), length.out = N)
  beta1_true <- phi_true * (-theta_true)
  beta0_true <- rep(0.5, N);  gam_true <- rep(0.3, N)
  gamma_y <- matrix(runif(N * m, 0.5, 1.5), N, m)
  gamma_x <- matrix(runif(N * m, 0.3, 0.7), N, m)
  phi_hat <- theta_hat <- numeric(R)

  cat(sprintf("Grid %2d/%d  (N=%3d, T=%3d) ... ", g, G, N, TT))
  t0 <- proc.time()

  for (r in 1:R) {
    if (r %% 20 == 0) cat(sprintf("[rep %d]", r))
    dat <- GenerateDGP(N, TT, m, rho_f, rho_x, sig_f, sig_u, sig_v,
                       phi_true, beta0_true, beta1_true, gam_true,
                       gamma_y, gamma_x)
    res <- CSB(dat$Y, dat$X, N, TT, PT)
    phi_hat[r]   <- res$phi
    theta_hat[r] <- res$theta
  }

  elapsed <- (proc.time() - t0)[3]
  phi_bar  <- mean(phi_true)
  bias_phi <- mean(phi_hat) - phi_bar
  rmse_phi <- sqrt(mean((phi_hat - phi_bar)^2))
  bias_th  <- mean(theta_hat) - theta_true
  rmse_th  <- sqrt(mean((theta_hat - theta_true)^2))
  pct_phi  <- 100 * bias_phi / abs(phi_bar)

  results <- rbind(results, data.frame(
    N = N, T = TT,
    CSB_phi_bias    = round(bias_phi, 4),
    CSB_phi_rmse    = round(rmse_phi, 4),
    CSB_phi_pctbias = round(pct_phi, 2),
    CSB_theta_bias  = round(bias_th, 4),
    CSB_theta_rmse  = round(rmse_th, 4)
  ))
  cat(sprintf(" phi bias=%.4f (%.1f%%)  [%.0fs]\n", bias_phi, pct_phi, elapsed))
}

write.csv(results, "07_mc_csardl_6_CSB_results.csv", row.names = FALSE)
cat("\n=== DONE.  Saved: 07_mc_csardl_6_CSB_results.csv ===\n")
