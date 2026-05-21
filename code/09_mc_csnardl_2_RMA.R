# =============================================================================
# 09_mc_csnardl_2_RMA.R
# CS-NARDL Monte Carlo - Recursive Mean Adjustment
#
# Method   : RMA: Recursive mean adjustment
# Reference : So and Shin (1999); Choi, Mark and Sul (2010)
#
# Monte Carlo design (manuscript Section 4.1):
#   Grid          : 5 x 5, with N and T in {40, 50, 100, 150, 200} (25 cells)
#   Replications  : R = 500 per cell
#   Factor model  : single common factor, AR(1) coefficient rho_f = 0.6
#   Truncation lag: pT = floor(T^(1/3))
#
# Input    : none (data are simulated internally with a fixed seed)
# Output   : 09_mc_csnardl_2_RMA_results.csv
#
# This script is self-contained and can be run on its own in RStudio.
# Run order: see 00_RUN_ALL.R
# =============================================================================

rm(list = ls())
cat("================================================================\n")
cat(" CS-NARDL  |  Script 2 of 6  |  METHOD: RMA\n")
cat("================================================================\n\n")

set.seed(2026)

# ============================================================
# 1.  RMA-NARDL ESTIMATOR
# ============================================================
EstRMA_NARDL <- function(Y, X, N, TT, PT) {
  nobs <- (TT - 1) - PT
  phi_i <- bp_i <- bn_i <- bp1_i <- bn1_i <- numeric(N)

  for (i in 1:N) {
    yi <- Y[i, ];  xi <- X[i, ]
    DX_raw <- c(0, diff(xi))
    X_pos <- cumsum(pmax(DX_raw, 0))
    X_neg <- cumsum(pmin(DX_raw, 0))

    DY <- YLAG <- XPOS <- XNEG <- XPLAG <- XNLAG <- DX <- numeric(nobs)
    YB <- XB <- DB <- matrix(0, nobs, PT + 1)
    row <- 1
    for (tt in (PT + 2):TT) {
      DY[row] <- yi[tt] - yi[tt-1]; YLAG[row] <- yi[tt-1]
      XPOS[row] <- X_pos[tt]; XNEG[row] <- X_neg[tt]
      XPLAG[row] <- X_pos[tt-1]; XNLAG[row] <- X_neg[tt-1]
      DX[row] <- xi[tt] - xi[tt-1]
      for (ell in 0:PT) {
        s <- tt - ell
        YB[row, ell+1] <- mean(Y[, 1:s])
        XB[row, ell+1] <- mean(X[, 1:s])
        if (s >= 2) DB[row, ell+1] <- mean(diff(colMeans(X[, 1:s])))
      }
      row <- row + 1
    }
    Z <- cbind(YLAG, XPOS, XNEG, XPLAG, XNLAG, DX, 1, YB, XB, DB)
    b <- solve(crossprod(Z) + 1e-8*diag(ncol(Z)), crossprod(Z, DY))
    phi_i[i] <- b[1]; bp_i[i] <- b[2]; bn_i[i] <- b[3]
    bp1_i[i] <- b[4]; bn1_i[i] <- b[5]
  }
  list(phi = phi_i, bp = bp_i, bn = bn_i, bp1 = bp1_i, bn1 = bn1_i)
}

# ============================================================
# 2.  DGP (NARDL)
# ============================================================
GenerateDGP_NARDL <- function(N, TT, m, rho_f, rho_x, sig_f, sig_u, sig_v,
                               phi_true, theta_pos_true, theta_neg_true,
                               gamma_y, gamma_x) {
  Tburn <- 50;  TT_total <- TT + Tburn
  f <- matrix(0, m, TT_total)
  for (tt in 2:TT_total) f[, tt] <- rho_f * f[, tt-1] + rnorm(m) * sig_f
  Y <- X <- matrix(0, N, TT_total)
  for (i in 1:N) {
    for (tt in 2:TT_total) {
      X[i, tt] <- rho_x * X[i, tt-1] + sum(gamma_x[i,] * f[, tt]) + rnorm(1)*sig_v
      DX <- X[i, tt] - X[i, tt-1]
      rho_i <- 1 + phi_true[i]
      beta1_pos <- phi_true[i] * (-theta_pos_true[i])
      beta1_neg <- phi_true[i] * (-theta_neg_true[i])
      Y[i, tt] <- rho_i * Y[i, tt-1] +
                  beta1_pos * pmax(X[i, tt-1] - X[i, 1], 0) +
                  beta1_neg * pmin(X[i, tt-1] - X[i, 1], 0) +
                  0.3 * DX + sum(gamma_y[i,] * f[, tt]) + rnorm(1)*sig_u
    }
  }
  list(Y = Y[, (Tburn+1):TT_total], X = X[, (Tburn+1):TT_total])
}

# ============================================================
# 3.  GRID & MC LOOP
# ============================================================
Nvec <- c(40,50,100,150,200); Tvec <- c(40,50,100,150,200)
Grid <- expand.grid(N=Nvec, T=Tvec); G <- nrow(Grid); R <- 500
m <- 1; rho_f <- 0.6; rho_x <- 0.6
sig_f <- 1; sig_u <- 1; sig_v <- 1
theta_pos_true_val <- 1.0; theta_neg_true_val <- 0.5
results <- data.frame()

for (g in 1:G) {
  N <- Grid$N[g]; TT <- Grid$T[g]
  PT <- max(1, min(floor(TT^(1/3)), floor((TT-1)/4)-1))
  phi_true <- rep(c(-0.2,-0.3,-0.4,-0.5,-0.6), length.out=N)
  theta_pos_true <- rep(theta_pos_true_val, N)
  theta_neg_true <- rep(theta_neg_true_val, N)
  gamma_y <- matrix(runif(N*m,0.5,1.5), N, m)
  gamma_x <- matrix(runif(N*m,0.3,0.7), N, m)
  phi_hat <- thp_hat <- thn_hat <- numeric(R)

  cat(sprintf("Grid %2d/%d  (N=%3d, T=%3d) ... ", g, G, N, TT))
  for (r in 1:R) {
    dat <- GenerateDGP_NARDL(N, TT, m, rho_f, rho_x, sig_f, sig_u, sig_v,
                              phi_true, theta_pos_true, theta_neg_true,
                              gamma_y, gamma_x)
    res <- EstRMA_NARDL(dat$Y, dat$X, N, TT, PT)
    phi_hat[r] <- mean(res$phi)
    thp_i <- -(res$bp + res$bp1) / res$phi
    thn_i <- -(res$bn + res$bn1) / res$phi
    wp <- res$phi^2 / sum(res$phi^2)
    thp_hat[r] <- sum(wp * thp_i); thn_hat[r] <- sum(wp * thn_i)
  }
  phi_bar <- mean(phi_true)
  bias_phi <- mean(phi_hat) - phi_bar
  rmse_phi <- sqrt(mean((phi_hat - phi_bar)^2))
  pct_phi <- 100 * bias_phi / abs(phi_bar)
  bias_thp <- mean(thp_hat) - theta_pos_true_val
  bias_thn <- mean(thn_hat) - theta_neg_true_val

  results <- rbind(results, data.frame(
    N=N, T=TT,
    RMA_phi_bias=round(bias_phi,4), RMA_phi_rmse=round(rmse_phi,4),
    RMA_phi_pctbias=round(pct_phi,2),
    RMA_thpos_bias=round(bias_thp,4),
    RMA_thpos_rmse=round(sqrt(mean((thp_hat-theta_pos_true_val)^2)),4),
    RMA_thneg_bias=round(bias_thn,4),
    RMA_thneg_rmse=round(sqrt(mean((thn_hat-theta_neg_true_val)^2)),4)
  ))
  cat(sprintf("phi bias=%.4f (%.1f%%)\n", bias_phi, pct_phi))
}

write.csv(results, "09_mc_csnardl_2_RMA_results.csv", row.names=FALSE)
cat("\n=== DONE.  Saved: 09_mc_csnardl_2_RMA_results.csv ===\n")
