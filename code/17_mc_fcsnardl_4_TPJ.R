# =============================================================================
# 17_mc_fcsnardl_4_TPJ.R
# FCS-NARDL Monte Carlo - Third-Panel Jackknife
#
# Method   : TPJ: Third-panel jackknife
# Reference : Dhaene and Jochmans (2015)
#
# Monte Carlo design (manuscript Section 4.1):
#   Grid          : 5 x 5, with N and T in {40, 50, 100, 150, 200} (25 cells)
#   Replications  : R = 500 per cell
#   Factor model  : single common factor, AR(1) coefficient rho_f = 0.6
#   Truncation lag: pT = floor(T^(1/3))
#
# Input    : none (data are simulated internally with a fixed seed)
# Output   : 17_mc_fcsnardl_4_TPJ_results.csv
#
# This script is self-contained and can be run on its own in RStudio.
# Run order: see 00_RUN_ALL.R
# =============================================================================

rm(list = ls())
cat("================================================================\n")
cat(" FCS-NARDL  |  Script 4 of 6  |  METHOD: TPJ\n")
cat("================================================================\n\n")

set.seed(2026)

EstFCS_NARDL <- function(Y, X, N, TT, PT, k_fourier=1) {
  ybar <- colMeans(Y); xbar <- colMeans(X); dxbar <- c(0,diff(xbar))
  nobs <- (TT-1)-PT
  phi_i <- bp_i <- bn_i <- bp1_i <- bn1_i <- numeric(N)
  for (i in 1:N) {
    yi <- Y[i,]; xi <- X[i,]
    DX_raw <- c(0,diff(xi)); X_pos <- cumsum(pmax(DX_raw,0)); X_neg <- cumsum(pmin(DX_raw,0))
    DY <- YLAG <- XPOS <- XNEG <- XPLAG <- XNLAG <- DX <- SIN_T <- COS_T <- numeric(nobs)
    YB <- XB <- DB <- matrix(0,nobs,PT+1); row <- 1
    for (tt in (PT+2):TT) {
      DY[row] <- yi[tt]-yi[tt-1]; YLAG[row] <- yi[tt-1]
      XPOS[row] <- X_pos[tt]; XNEG[row] <- X_neg[tt]
      XPLAG[row] <- X_pos[tt-1]; XNLAG[row] <- X_neg[tt-1]
      DX[row] <- xi[tt]-xi[tt-1]
      SIN_T[row] <- sin(2*pi*k_fourier*tt/TT); COS_T[row] <- cos(2*pi*k_fourier*tt/TT)
      for (ell in 0:PT) { YB[row,ell+1] <- ybar[tt-ell]; XB[row,ell+1] <- xbar[tt-ell]; DB[row,ell+1] <- dxbar[tt-ell] }
      row <- row+1
    }
    Z <- cbind(YLAG,XPOS,XNEG,XPLAG,XNLAG,DX,SIN_T,COS_T,1,YB,XB,DB)
    b <- solve(crossprod(Z)+1e-8*diag(ncol(Z)),crossprod(Z,DY))
    phi_i[i] <- b[1]; bp_i[i] <- b[2]; bn_i[i] <- b[3]; bp1_i[i] <- b[4]; bn1_i[i] <- b[5]
  }
  list(phi=phi_i,bp=bp_i,bn=bn_i,bp1=bp1_i,bn1=bn1_i)
}

EstTPJ_FCS <- function(Y, X, N, TT, PT, k_fourier=1) {
  T3a <- floor(TT/3); T3c <- TT-2*T3a
  if (T3a < PT+2 || T3c < PT+2) return(EstFCS_NARDL(Y,X,N,TT,PT,k_fourier))
  full <- EstFCS_NARDL(Y,X,N,TT,PT,k_fourier)
  e12 <- 1:(2*T3a); T12 <- length(e12)
  e23 <- (T3a+1):TT; T23 <- length(e23)
  e13 <- c(1:T3a,(2*T3a+1):TT); T13 <- length(e13)
  r12 <- EstFCS_NARDL(Y[,e12,drop=F],X[,e12,drop=F],N,T12,PT,k_fourier)
  r23 <- EstFCS_NARDL(Y[,e23,drop=F],X[,e23,drop=F],N,T23,PT,k_fourier)
  r13 <- EstFCS_NARDL(Y[,e13,drop=F],X[,e13,drop=F],N,T13,PT,k_fourier)
  list(phi=3*full$phi-0.5*(r12$phi+r23$phi+r13$phi),
       bp=3*full$bp-0.5*(r12$bp+r23$bp+r13$bp),
       bn=3*full$bn-0.5*(r12$bn+r23$bn+r13$bn),
       bp1=3*full$bp1-0.5*(r12$bp1+r23$bp1+r13$bp1),
       bn1=3*full$bn1-0.5*(r12$bn1+r23$bn1+r13$bn1))
}

GenerateDGP_FCS <- function(N,TT,m,rho_f,rho_x,sig_f,sig_u,sig_v,
                             phi_true,theta_pos_true,theta_neg_true,
                             alpha1_true,alpha2_true,k_fourier,gamma_y,gamma_x) {
  Tburn <- 50; TT_total <- TT+Tburn; f <- matrix(0,m,TT_total)
  for (tt in 2:TT_total) f[,tt] <- rho_f*f[,tt-1]+rnorm(m)*sig_f
  Y <- X <- matrix(0,N,TT_total)
  for (i in 1:N) for (tt in 2:TT_total) {
    X[i,tt] <- rho_x*X[i,tt-1]+sum(gamma_x[i,]*f[,tt])+rnorm(1)*sig_v
    rho_i <- 1+phi_true[i]
    Y[i,tt] <- rho_i*Y[i,tt-1]+phi_true[i]*(-theta_pos_true[i])*pmax(X[i,tt-1]-X[i,1],0)+
      phi_true[i]*(-theta_neg_true[i])*pmin(X[i,tt-1]-X[i,1],0)+
      0.3*(X[i,tt]-X[i,tt-1])+alpha1_true[i]*sin(2*pi*k_fourier*tt/TT_total)+
      alpha2_true[i]*cos(2*pi*k_fourier*tt/TT_total)+sum(gamma_y[i,]*f[,tt])+rnorm(1)*sig_u
  }
  list(Y=Y[,(Tburn+1):TT_total],X=X[,(Tburn+1):TT_total])
}

Nvec <- c(40,50,100,150,200); Tvec <- c(40,50,100,150,200)
Grid <- expand.grid(N=Nvec,T=Tvec); G <- nrow(Grid); R <- 500
m <- 1; rho_f <- 0.6; rho_x <- 0.6; sig_f <- 1; sig_u <- 1; sig_v <- 1
theta_pos_true_val <- 1.0; theta_neg_true_val <- 0.5; k_fourier <- 1; results <- data.frame()

for (g in 1:G) {
  N <- Grid$N[g]; TT <- Grid$T[g]
  PT <- max(1,min(floor(TT^(1/3)),floor((TT-1)/4)-1))
  phi_true <- rep(c(-0.2,-0.3,-0.4,-0.5,-0.6),length.out=N)
  theta_pos_true <- rep(theta_pos_true_val,N); theta_neg_true <- rep(theta_neg_true_val,N)
  alpha1_true <- rep(0.3,N); alpha2_true <- rep(0.2,N)
  gamma_y <- matrix(runif(N*m,0.5,1.5),N,m); gamma_x <- matrix(runif(N*m,0.3,0.7),N,m)
  phi_hat <- thp_hat <- thn_hat <- numeric(R)
  cat(sprintf("Grid %2d/%d  (N=%3d, T=%3d) ... ",g,G,N,TT))
  for (r in 1:R) {
    dat <- GenerateDGP_FCS(N,TT,m,rho_f,rho_x,sig_f,sig_u,sig_v,phi_true,theta_pos_true,theta_neg_true,
                            alpha1_true,alpha2_true,k_fourier,gamma_y,gamma_x)
    res <- EstTPJ_FCS(dat$Y,dat$X,N,TT,PT,k_fourier)
    phi_hat[r] <- mean(res$phi)
    thp_i <- -(res$bp+res$bp1)/res$phi; thn_i <- -(res$bn+res$bn1)/res$phi
    wp <- res$phi^2/sum(res$phi^2)
    thp_hat[r] <- sum(wp*thp_i); thn_hat[r] <- sum(wp*thn_i)
  }
  phi_bar <- mean(phi_true)
  results <- rbind(results, data.frame(
    N=N,T=TT, TPJ_phi_bias=round(mean(phi_hat)-phi_bar,4),
    TPJ_phi_rmse=round(sqrt(mean((phi_hat-phi_bar)^2)),4),
    TPJ_phi_pctbias=round(100*(mean(phi_hat)-phi_bar)/abs(phi_bar),2),
    TPJ_thpos_bias=round(mean(thp_hat)-theta_pos_true_val,4),
    TPJ_thpos_rmse=round(sqrt(mean((thp_hat-theta_pos_true_val)^2)),4),
    TPJ_thneg_bias=round(mean(thn_hat)-theta_neg_true_val,4),
    TPJ_thneg_rmse=round(sqrt(mean((thn_hat-theta_neg_true_val)^2)),4)
  ))
  cat(sprintf("phi bias=%.4f\n",mean(phi_hat)-phi_bar))
}
write.csv(results,"17_mc_fcsnardl_4_TPJ_results.csv",row.names=FALSE)
cat("\n=== DONE.  Saved: 17_mc_fcsnardl_4_TPJ_results.csv ===\n")
