# =============================================================================
# 18_mc_fcsnardl_5_BBC.R
# FCS-NARDL Monte Carlo - Bootstrap Bias Correction
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
# Output   : 18_mc_fcsnardl_5_BBC_results.csv
#
# This script is self-contained and can be run on its own in RStudio.
# Run order: see 00_RUN_ALL.R
# =============================================================================

rm(list = ls())
cat("================================================================\n")
cat(" FCS-NARDL  |  Script 5 of 6  |  METHOD: BBC (v3)\n")
cat("================================================================\n\n")

set.seed(2026)

EstFCS_Full <- function(Y, X, N, TT, PT, k_fourier=1, lambda=1e-6) {
  ybar <- colMeans(Y); xbar <- colMeans(X); dxbar <- c(0,diff(xbar))
  nobs <- (TT-1)-PT
  phi_i <- bp_i <- bn_i <- bp1_i <- bn1_i <- a1_i <- a2_i <- const_i <- numeric(N)
  U <- matrix(0,N,TT)
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
    b <- solve(crossprod(Z)+lambda*diag(ncol(Z)),crossprod(Z,DY))
    phi_i[i] <- b[1]; bp_i[i] <- b[2]; bn_i[i] <- b[3]
    bp1_i[i] <- b[4]; bn1_i[i] <- b[5]; a1_i[i] <- b[7]; a2_i[i] <- b[8]; const_i[i] <- b[9]
    resid <- DY - Z %*% b
    U[i,(PT+2):TT] <- as.numeric(resid)
  }
  list(phi=phi_i,bp=bp_i,bn=bn_i,bp1=bp1_i,bn1=bn1_i,a1=a1_i,a2=a2_i,const=const_i,U=U)
}

BBC_FCS <- function(Y, X, N, TT, PT, k_fourier=1, B_boot=100, K_iter=3) {
  full <- EstFCS_Full(Y,X,N,TT,PT,k_fourier)
  phi_curr <- full$phi; bp_curr <- full$bp; bn_curr <- full$bn
  bp1_curr <- full$bp1; bn1_curr <- full$bn1
  a1_curr <- full$a1; a2_curr <- full$a2; const_curr <- full$const

  for (k in 1:K_iter) {
    phi_boot <- matrix(0,B_boot,N)
    for (b in 1:B_boot) {
      eta <- sample(c(-1,1),TT,replace=TRUE)
      Y_star <- matrix(0,N,TT)
      for (i in 1:N) {
        Y_star[i,1] <- Y[i,1]
        rho_i <- max(0.1,min(0.95,1+phi_curr[i]))
        DX_raw <- c(0,diff(X[i,])); xp <- cumsum(pmax(DX_raw,0)); xn <- cumsum(pmin(DX_raw,0))
        for (tt in 2:TT) {
          Y_star[i,tt] <- rho_i*Y_star[i,tt-1]+
            bp_curr[i]*xp[tt]+bn_curr[i]*xn[tt]+bp1_curr[i]*xp[max(1,tt-1)]+bn1_curr[i]*xn[max(1,tt-1)]+
            0.3*(X[i,tt]-X[i,tt-1])+a1_curr[i]*sin(2*pi*k_fourier*tt/TT)+
            a2_curr[i]*cos(2*pi*k_fourier*tt/TT)+const_curr[i]+full$U[i,tt]*eta[tt]
        }
      }
      res_b <- tryCatch(EstFCS_Full(Y_star,X,N,TT,PT,k_fourier,lambda=1e-5),error=function(e) NULL)
      if (!is.null(res_b)) phi_boot[b,] <- res_b$phi
    }
    for (i in 1:N) {
      phi_new <- 2*phi_curr[i]-mean(phi_boot[,i])
      phi_curr[i] <- max(-0.95,min(-0.05,phi_new))
    }
  }
  thp_i <- -(bp_curr+bp1_curr)/phi_curr; thn_i <- -(bn_curr+bn1_curr)/phi_curr
  w <- phi_curr^2/sum(phi_curr^2)
  list(phi=phi_curr,theta_pos=sum(w*thp_i),theta_neg=sum(w*thn_i))
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
Grid <- expand.grid(N=Nvec,T=Tvec); G <- nrow(Grid); R <- 200
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
  cat(sprintf("Grid %2d/%d  (N=%3d, T=%3d) ... ",g,G,N,TT)); t0 <- proc.time()
  for (r in 1:R) {
    if (r%%20==0) cat(sprintf("[rep %d]",r))
    dat <- GenerateDGP_FCS(N,TT,m,rho_f,rho_x,sig_f,sig_u,sig_v,phi_true,theta_pos_true,theta_neg_true,
                            alpha1_true,alpha2_true,k_fourier,gamma_y,gamma_x)
    res <- BBC_FCS(dat$Y,dat$X,N,TT,PT,k_fourier)
    phi_hat[r] <- mean(res$phi); thp_hat[r] <- res$theta_pos; thn_hat[r] <- res$theta_neg
  }
  elapsed <- (proc.time()-t0)[3]; phi_bar <- mean(phi_true)
  results <- rbind(results, data.frame(
    N=N,T=TT,BBC_phi_bias=round(mean(phi_hat)-phi_bar,4),
    BBC_phi_rmse=round(sqrt(mean((phi_hat-phi_bar)^2)),4),
    BBC_phi_pctbias=round(100*(mean(phi_hat)-phi_bar)/abs(phi_bar),2),
    BBC_thpos_bias=round(mean(thp_hat)-theta_pos_true_val,4),
    BBC_thpos_rmse=round(sqrt(mean((thp_hat-theta_pos_true_val)^2)),4),
    BBC_thneg_bias=round(mean(thn_hat)-theta_neg_true_val,4),
    BBC_thneg_rmse=round(sqrt(mean((thn_hat-theta_neg_true_val)^2)),4)
  ))
  cat(sprintf(" phi bias=%.4f [%.0fs]\n",mean(phi_hat)-phi_bar,elapsed))
}
write.csv(results,"18_mc_fcsnardl_5_BBC_results.csv",row.names=FALSE)
cat("\n=== DONE.  Saved: 18_mc_fcsnardl_5_BBC_results.csv ===\n")
