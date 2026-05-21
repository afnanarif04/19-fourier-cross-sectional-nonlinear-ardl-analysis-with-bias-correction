# =============================================================================
# 06_mc_csardl_5_BBC.R
# CS-ARDL Monte Carlo - Bootstrap Bias Correction
#
# Method   : BBC: Bootstrap bias correction
# Reference : De Vos, Everaert and Ruyssen (2015)
#
# Monte Carlo design (manuscript Section 4.1):
#   Grid          : 5 x 5, with N and T in {40, 50, 100, 150, 200} (25 cells)
#   Replications  : R = 500 per cell
#   Factor model  : single common factor, AR(1) coefficient rho_f = 0.6
#   Truncation lag: pT = floor(T^(1/3))
#
# Input    : none (data are simulated internally with a fixed seed)
# Output   : 06_mc_csardl_5_BBC_results.csv
#
# This script is self-contained and can be run on its own in RStudio.
# Run order: see 00_RUN_ALL.R
# =============================================================================

rm(list = ls())
cat("================================================================\n")
cat(" CS-ARDL  |  Script 5 of 6  |  METHOD: BBC (v3)\n")
cat("================================================================\n\n")

set.seed(2026)

# ============================================================
# 1.  CS ESTIMATOR (full version returning residuals)
# ============================================================
EstCS_Full <- function(Y, X, N, TT, PT, lambda = 1e-6) {
  ybar  <- colMeans(Y);  xbar <- colMeans(X)
  dxbar <- c(0, diff(xbar))
  nobs  <- (TT - 1) - PT
  phi_i <- b0_i <- b1_i <- gam_i <- const_i <- numeric(N)
  U <- matrix(0, N, TT)

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
    b <- solve(crossprod(Z) + lambda * diag(ncol(Z)), crossprod(Z, DY))
    phi_i[i] <- b[1]; b0_i[i] <- b[2]; b1_i[i] <- b[3]
    gam_i[i] <- b[4]; const_i[i] <- b[5]
    resid <- DY - Z %*% b
    U[i, (PT + 2):TT] <- as.numeric(resid)
  }
  list(phi = phi_i, b0 = b0_i, b1 = b1_i, gam = gam_i,
       const = const_i, U = U)
}

# ============================================================
# 2.  BBC ALGORITHM
# ============================================================
BBC <- function(Y, X, N, TT, PT, B_boot = 100, K_iter = 3) {
  # Step 0: full-sample estimates
  full <- EstCS_Full(Y, X, N, TT, PT)
  phi_curr  <- full$phi
  b0_curr   <- full$b0
  b1_curr   <- full$b1
  gam_curr  <- full$gam
  const_curr <- full$const

  for (k in 1:K_iter) {
    phi_boot <- matrix(0, B_boot, N)

    for (b in 1:B_boot) {
      # Wild bootstrap: common time multiplier
      eta <- sample(c(-1, 1), TT, replace = TRUE)

      Y_star <- matrix(0, N, TT)
      for (i in 1:N) {
        Y_star[i, 1] <- Y[i, 1]
        rho_i <- max(0.1, min(0.95, 1 + phi_curr[i]))  # v3 bound
        for (tt in 2:TT) {
          Y_star[i, tt] <- rho_i * Y_star[i, tt - 1] +
            b0_curr[i] * X[i, tt] + b1_curr[i] * X[i, tt - 1] +
            gam_curr[i] * (X[i, tt] - X[i, tt - 1]) +
            const_curr[i] + full$U[i, tt] * eta[tt]
        }
      }

      res_b <- tryCatch(
        EstCS_Full(Y_star, X, N, TT, PT, lambda = 1e-5),
        error = function(e) NULL
      )
      if (!is.null(res_b)) phi_boot[b, ] <- res_b$phi
    }

    # Bias correction with v3 parameter bounds
    for (i in 1:N) {
      phi_new <- 2 * phi_curr[i] - mean(phi_boot[, i])
      phi_curr[i] <- max(-0.95, min(-0.05, phi_new))
    }
  }

  theta_i <- -b1_curr / phi_curr
  w       <- phi_curr^2 / sum(phi_curr^2)
  list(phi = phi_curr, b1 = b1_curr, theta_mg = sum(w * theta_i))
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
    res <- BBC(dat$Y, dat$X, N, TT, PT)
    phi_hat[r]   <- mean(res$phi)
    theta_hat[r] <- res$theta_mg
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
    BBC_phi_bias    = round(bias_phi, 4),
    BBC_phi_rmse    = round(rmse_phi, 4),
    BBC_phi_pctbias = round(pct_phi, 2),
    BBC_theta_bias  = round(bias_th, 4),
    BBC_theta_rmse  = round(rmse_th, 4)
  ))
  cat(sprintf(" phi bias=%.4f (%.1f%%)  [%.0fs]\n", bias_phi, pct_phi, elapsed))
}

write.csv(results, "06_mc_csardl_5_BBC_results.csv", row.names = FALSE)
cat("\n=== DONE.  Saved: 06_mc_csardl_5_BBC_results.csv ===\n")
