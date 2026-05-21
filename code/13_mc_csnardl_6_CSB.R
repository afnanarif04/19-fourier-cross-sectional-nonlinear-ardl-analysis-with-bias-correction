# =============================================================================
# 13_mc_csnardl_6_CSB.R
# CS-NARDL Monte Carlo - Cross-Section Bootstrap
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
# Output   : 13_mc_csnardl_6_CSB_results.csv
#
# This script is self-contained and can be run on its own in RStudio.
# Run order: see 00_RUN_ALL.R
# =============================================================================

rm(list = ls())
cat("================================================================\n")
cat(" CS-NARDL  |  Script 6 of 6  |  METHOD: CSB\n")
cat("================================================================\n\n")

set.seed(2026)

EstCS_NARDL <- function(Y, X, N, TT, PT) {
  ybar <- colMeans(Y); xbar <- colMeans(X); dxbar <- c(0,diff(xbar))
  nobs <- (TT-1)-PT
  phi_i <- bp_i <- bn_i <- bp1_i <- bn1_i <- numeric(N)
  for (i in 1:N) {
    yi <- Y[i,]; xi <- X[i,]
    DX_raw <- c(0,diff(xi)); X_pos <- cumsum(pmax(DX_raw,0)); X_neg <- cumsum(pmin(DX_raw,0))
    DY <- YLAG <- XPOS <- XNEG <- XPLAG <- XNLAG <- DX <- numeric(nobs)
    YB <- XB <- DB <- matrix(0,nobs,PT+1); row <- 1
    for (tt in (PT+2):TT) {
      DY[row] <- yi[tt]-yi[tt-1]; YLAG[row] <- yi[tt-1]
      XPOS[row] <- X_pos[tt]; XNEG[row] <- X_neg[tt]
      XPLAG[row] <- X_pos[tt-1]; XNLAG[row] <- X_neg[tt-1]
      DX[row] <- xi[tt]-xi[tt-1]
      for (ell in 0:PT) { YB[row,ell+1] <- ybar[tt-ell]; XB[row,ell+1] <- xbar[tt-ell]; DB[row,ell+1] <- dxbar[tt-ell] }
      row <- row+1
    }
    Z <- cbind(YLAG,XPOS,XNEG,XPLAG,XNLAG,DX,1,YB,XB,DB)
    b <- solve(crossprod(Z)+1e-8*diag(ncol(Z)),crossprod(Z,DY))
    phi_i[i] <- b[1]; bp_i[i] <- b[2]; bn_i[i] <- b[3]; bp1_i[i] <- b[4]; bn1_i[i] <- b[5]
  }
  list(phi=phi_i,bp=bp_i,bn=bn_i,bp1=bp1_i,bn1=bn1_i)
}

CSB_NARDL <- function(Y, X, N, TT, PT, B_boot=200) {
  full <- EstCS_NARDL(Y,X,N,TT,PT)
  thp_i <- -(full$bp+full$bp1)/full$phi; thn_i <- -(full$bn+full$bn1)/full$phi
  w <- full$phi^2/sum(full$phi^2)
  thp_full <- sum(w*thp_i); thn_full <- sum(w*thn_i); phi_full <- mean(full$phi)

  thp_boot <- thn_boot <- phi_boot <- numeric(B_boot)
  for (b in 1:B_boot) {
    idx <- sample(1:N,N,replace=TRUE)
    res_b <- tryCatch(EstCS_NARDL(Y[idx,,drop=F],X[idx,,drop=F],N,TT,PT),error=function(e) NULL)
    if (!is.null(res_b)) {
      tp <- -(res_b$bp+res_b$bp1)/res_b$phi; tn <- -(res_b$bn+res_b$bn1)/res_b$phi
      wb <- res_b$phi^2/sum(res_b$phi^2)
      thp_boot[b] <- sum(wb*tp); thn_boot[b] <- sum(wb*tn); phi_boot[b] <- mean(res_b$phi)
    } else { thp_boot[b] <- thp_full; thn_boot[b] <- thn_full; phi_boot[b] <- phi_full }
  }
  list(phi=2*phi_full-mean(phi_boot), theta_pos=2*thp_full-mean(thp_boot), theta_neg=2*thn_full-mean(thn_boot))
}

GenerateDGP_NARDL <- function(N,TT,m,rho_f,rho_x,sig_f,sig_u,sig_v,
                               phi_true,theta_pos_true,theta_neg_true,gamma_y,gamma_x) {
  Tburn <- 50; TT_total <- TT+Tburn; f <- matrix(0,m,TT_total)
  for (tt in 2:TT_total) f[,tt] <- rho_f*f[,tt-1]+rnorm(m)*sig_f
  Y <- X <- matrix(0,N,TT_total)
  for (i in 1:N) for (tt in 2:TT_total) {
    X[i,tt] <- rho_x*X[i,tt-1]+sum(gamma_x[i,]*f[,tt])+rnorm(1)*sig_v
    rho_i <- 1+phi_true[i]
    Y[i,tt] <- rho_i*Y[i,tt-1]+phi_true[i]*(-theta_pos_true[i])*pmax(X[i,tt-1]-X[i,1],0)+
      phi_true[i]*(-theta_neg_true[i])*pmin(X[i,tt-1]-X[i,1],0)+
      0.3*(X[i,tt]-X[i,tt-1])+sum(gamma_y[i,]*f[,tt])+rnorm(1)*sig_u
  }
  list(Y=Y[,(Tburn+1):TT_total],X=X[,(Tburn+1):TT_total])
}

Nvec <- c(40,50,100,150,200); Tvec <- c(40,50,100,150,200)
Grid <- expand.grid(N=Nvec,T=Tvec); G <- nrow(Grid); R <- 200
m <- 1; rho_f <- 0.6; rho_x <- 0.6; sig_f <- 1; sig_u <- 1; sig_v <- 1
theta_pos_true_val <- 1.0; theta_neg_true_val <- 0.5; results <- data.frame()

for (g in 1:G) {
  N <- Grid$N[g]; TT <- Grid$T[g]
  PT <- max(1,min(floor(TT^(1/3)),floor((TT-1)/4)-1))
  phi_true <- rep(c(-0.2,-0.3,-0.4,-0.5,-0.6),length.out=N)
  theta_pos_true <- rep(theta_pos_true_val,N); theta_neg_true <- rep(theta_neg_true_val,N)
  gamma_y <- matrix(runif(N*m,0.5,1.5),N,m); gamma_x <- matrix(runif(N*m,0.3,0.7),N,m)
  phi_hat <- thp_hat <- thn_hat <- numeric(R)
  cat(sprintf("Grid %2d/%d  (N=%3d, T=%3d) ... ",g,G,N,TT)); t0 <- proc.time()
  for (r in 1:R) {
    if (r%%20==0) cat(sprintf("[rep %d]",r))
    dat <- GenerateDGP_NARDL(N,TT,m,rho_f,rho_x,sig_f,sig_u,sig_v,
                              phi_true,theta_pos_true,theta_neg_true,gamma_y,gamma_x)
    res <- CSB_NARDL(dat$Y,dat$X,N,TT,PT)
    phi_hat[r] <- res$phi; thp_hat[r] <- res$theta_pos; thn_hat[r] <- res$theta_neg
  }
  elapsed <- (proc.time()-t0)[3]; phi_bar <- mean(phi_true)
  results <- rbind(results, data.frame(
    N=N,T=TT,
    CSB_phi_bias=round(mean(phi_hat)-phi_bar,4),
    CSB_phi_rmse=round(sqrt(mean((phi_hat-phi_bar)^2)),4),
    CSB_phi_pctbias=round(100*(mean(phi_hat)-phi_bar)/abs(phi_bar),2),
    CSB_thpos_bias=round(mean(thp_hat)-theta_pos_true_val,4),
    CSB_thpos_rmse=round(sqrt(mean((thp_hat-theta_pos_true_val)^2)),4),
    CSB_thneg_bias=round(mean(thn_hat)-theta_neg_true_val,4),
    CSB_thneg_rmse=round(sqrt(mean((thn_hat-theta_neg_true_val)^2)),4)
  ))
  cat(sprintf(" phi bias=%.4f [%.0fs]\n",mean(phi_hat)-phi_bar,elapsed))
}
write.csv(results,"13_mc_csnardl_6_CSB_results.csv",row.names=FALSE)
cat("\n=== DONE.  Saved: 13_mc_csnardl_6_CSB_results.csv ===\n")
