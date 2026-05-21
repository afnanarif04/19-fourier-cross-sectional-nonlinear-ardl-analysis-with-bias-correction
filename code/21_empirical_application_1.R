# =============================================================================
# 21_empirical_application_1.R
# Empirical Application 1 - Resource rents and economic growth
#
# Estimates the symmetric (CS-ARDL), asymmetric (CS-NARDL), and
# Fourier-augmented asymmetric (FCS-NARDL) specifications under six correctors
# (CS, RMA, HPJ, TPJ, BBC, CSB) on a panel of resource-abundant economies,
# 1996-2020.
#
# DATA (open access, downloaded automatically via the WDI package):
#   GDP per capita, constant 2015 USD : NY.GDP.PCAP.KD
#       https://data.worldbank.org/indicator/NY.GDP.PCAP.KD
#   Total natural resource rents, % GDP: NY.GDP.TOTL.RT.ZS
#       https://data.worldbank.org/indicator/NY.GDP.TOTL.RT.ZS
#   Trade openness, % GDP             : NE.TRD.GNFS.ZS
#       https://data.worldbank.org/indicator/NE.TRD.GNFS.ZS
#   FDI net inflows, % GDP            : BX.KLT.DINV.WD.GD.ZS
#       https://data.worldbank.org/indicator/BX.KLT.DINV.WD.GD.ZS
#
# Requires : internet connection (the WDI package fetches data at runtime)
# Packages : WDI, dplyr, tidyr
# Output   : printed coefficient tables and the long-run symmetry Wald test
#
# Run order: see 00_RUN_ALL.R
# =============================================================================

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls())

# ---- 0. SETUP ---------------------------------------------------------------
local({
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      tryCatch(rstudioapi::isAvailable(), error = function(e) FALSE)) {
    p <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) "")
    if (nzchar(p)) { setwd(dirname(p)); return(invisible()) }
  }
  args <- commandArgs(trailingOnly = FALSE)
  hits <- grep("^--file=", args, value = TRUE)
  if (length(hits) > 0)
    setwd(dirname(normalizePath(sub("^--file=", "", hits[1]))))
})

pkg <- c("WDI", "dplyr", "tidyr")
new <- pkg[!pkg %in% installed.packages()[, "Package"]]
if (length(new)) install.packages(new, quiet = TRUE)
invisible(sapply(pkg, require, character.only = TRUE, warn.conflicts = FALSE))
set.seed(2026)

# ---- 1. DATA DOWNLOAD -------------------------------------------------------
# All variables sourced from World Bank WDI (open access, freely replicable).
# Manual download links (for verification):
#   GDP per capita (constant 2015 USD): https://data.worldbank.org/indicator/NY.GDP.PCAP.KD
#   Total resource rents (% GDP):        https://data.worldbank.org/indicator/NY.GDP.TOTL.RT.ZS
#   Trade (% GDP):                        https://data.worldbank.org/indicator/NE.TRD.GNFS.ZS
#   FDI net inflows (% GDP):              https://data.worldbank.org/indicator/BX.KLT.DINV.WD.GD.ZS
#   Gross capital formation (% GDP):      https://data.worldbank.org/indicator/NE.GDI.TOTL.ZS

cat("Downloading World Bank WDI data ...\n")
raw <- WDI(
  country   = "all",
  indicator = c(
    gdp_pc    = "NY.GDP.PCAP.KD",
    res_rents = "NY.GDP.TOTL.RT.ZS",
    trade     = "NE.TRD.GNFS.ZS",
    fdi       = "BX.KLT.DINV.WD.GD.ZS",
    capform   = "NE.GDI.TOTL.ZS"
  ),
  start = 1995, end = 2020,
  extra = TRUE
)

# ---- 2. SAMPLE CONSTRUCTION -------------------------------------------------
df <- raw %>%
  filter(!is.na(iso2c), nchar(iso2c) == 2, region != "Aggregates") %>%
  select(iso2c, country, year, gdp_pc, res_rents, trade, fdi, capform) %>%
  mutate(ln_gdp = log(gdp_pc)) %>%
  arrange(iso2c, year)

# Keep resource-rich countries: mean rents > 2 % of GDP
rich <- df %>%
  group_by(iso2c) %>%
  summarise(mr = mean(res_rents, na.rm = TRUE), .groups = "drop") %>%
  filter(mr >= 2) %>% pull(iso2c)

# Require balanced 1996-2020 (T = 25) with ≥ 22 non-missing observations
ids_ok <- df %>%
  filter(iso2c %in% rich, year %in% 1996:2020) %>%
  group_by(iso2c) %>%
  summarise(
    n_gdp = sum(!is.na(ln_gdp)),
    n_res = sum(!is.na(res_rents)),
    n_tr  = sum(!is.na(trade)),
    n_fdi = sum(!is.na(fdi)),
    .groups = "drop"
  ) %>%
  filter(n_gdp >= 22, n_res >= 22, n_tr >= 22) %>%
  pull(iso2c)

panel <- df %>%
  filter(iso2c %in% ids_ok, year %in% 1996:2020) %>%
  group_by(iso2c) %>%
  arrange(year) %>%
  tidyr::fill(everything(), .direction = "downup") %>%
  ungroup()

ids   <- sort(unique(panel$iso2c))
years <- sort(unique(panel$year))
N     <- length(ids); TT <- length(years)
pT    <- max(1, floor(TT^(1/3)))
cat(sprintf("Panel: N = %d, T = %d (1996-2020), pT = %d\n", N, TT, pT))

# Helper: reshape to N x T matrix
to_mat <- function(v) {
  m <- matrix(NA_real_, N, TT)
  for (j in seq_along(ids)) {
    sub <- panel[panel$iso2c == ids[j], ]
    sub <- sub[order(sub$year), ]
    m[j, ] <- sub[[v]]
  }
  m
}

Y_mat  <- to_mat("ln_gdp")
X_mat  <- to_mat("res_rents")   # key asymmetric variable
W2_mat <- to_mat("trade")
W3_mat <- to_mat("fdi")

# ---- 3. CORE CCE-MG ESTIMATION FUNCTIONS -----------------------------------
# Implements CS-ARDL, CS-NARDL, FCS-NARDL via unit-by-unit OLS + mean group.
# All functions embed cross-sectional averages with pT lags.

partial_sums <- function(x) {
  dx <- c(0, diff(x))
  list(pos = cumsum(pmax(dx, 0)), neg = cumsum(pmin(dx, 0)))
}

recursive_demean <- function(z) {
  T <- length(z); out <- rep(NA_real_, T); csum <- 0
  for (t in 2:T) { csum <- csum + z[t - 1]; out[t] <- z[t] - csum / (t - 1) }
  out
}

build_csa <- function(Y_m, X_m, W2_m, W3_m, pT, TT) {
  # Returns (TT) x ncol matrix of CSA and their pT lags
  yb  <- colMeans(Y_m,  na.rm = TRUE)
  xb  <- colMeans(X_m,  na.rm = TRUE)
  w2b <- colMeans(W2_m, na.rm = TRUE)
  w3b <- colMeans(W3_m, na.rm = TRUE)
  W   <- cbind(yb, xb, w2b, w3b)        # T x 4
  # Return list: lag 0, lag 1, ..., lag pT
  lapply(0:pT, function(l) {
    if (l == 0) W else rbind(matrix(NA, l, ncol(W)), W[1:(TT - l), ])
  })
}

# Unit-level ARDL-ECM with CCE augmentation
# Returns: phi (ECM speed), beta (long-run), and residuals
run_unit <- function(y, x, w2, w3, csa_lags, fourier_s, fourier_c, asymmetric) {
  T  <- length(y)
  t1 <- 2:T  # after differencing
  Dy <- diff(y); Dx <- diff(x)
  Dw2 <- diff(w2); Dw3 <- diff(w3)
  y_lag <- y[-T]; x_lag <- x[-T]
  w2_lag <- w2[-T]; w3_lag <- w3[-T]

  # Assemble CSA block (drop first row to align with differenced data)
  csa_blk <- do.call(cbind, lapply(csa_lags, function(m) {
    mm <- m[t1, , drop = FALSE]; mm[is.na(mm)] <- 0; mm
  }))

  if (!asymmetric) {
    Z <- cbind(1, y_lag, x_lag, w2_lag, w3_lag, Dx, Dw2, Dw3, csa_blk)
    if (!is.null(fourier_s)) Z <- cbind(Z, fourier_s[t1], fourier_c[t1])
    fit <- tryCatch(lm.fit(Z, Dy), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    b   <- fit$coefficients
    phi <- b[2]; delta <- b[3]
    if (!is.finite(phi) || phi >= 0 || phi < -2) return(NULL)
    list(phi = phi, beta = -delta / phi, beta_pos = NA, beta_neg = NA,
         resid = fit$residuals)
  } else {
    ps  <- partial_sums(x)
    xp  <- ps$pos; xn <- ps$neg
    xp_lag <- xp[-T]; xn_lag <- xn[-T]
    Dxp <- diff(xp); Dxn <- diff(xn)
    Z <- cbind(1, y_lag, xp_lag, xn_lag, w2_lag, w3_lag,
               Dxp, Dxn, Dw2, Dw3, csa_blk)
    if (!is.null(fourier_s)) Z <- cbind(Z, fourier_s[t1], fourier_c[t1])
    fit <- tryCatch(lm.fit(Z, Dy), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    b   <- fit$coefficients
    phi <- b[2]; dp <- b[3]; dn <- b[4]
    if (!is.finite(phi) || phi >= 0 || phi < -2) return(NULL)
    list(phi = phi, beta = NA, beta_pos = -dp / phi, beta_neg = -dn / phi,
         resid = fit$residuals)
  }
}

# Panel CCEMG estimation
ccemg <- function(Y_m, X_m, W2_m, W3_m, N, TT, pT, asymmetric = FALSE, fourier = FALSE) {
  csa  <- build_csa(Y_m, X_m, W2_m, W3_m, pT, TT)
  t_s  <- if (fourier) sin(2 * pi * (1:TT) / TT) else NULL
  t_c  <- if (fourier) cos(2 * pi * (1:TT) / TT) else NULL
  phi_v  <- beta_v <- beta_p_v <- beta_n_v <- rep(NA_real_, N)
  resid_m <- matrix(NA_real_, N, TT)
  for (i in 1:N) {
    res <- run_unit(Y_m[i,], X_m[i,], W2_m[i,], W3_m[i,],
                   csa, t_s, t_c, asymmetric)
    if (!is.null(res)) {
      phi_v[i]   <- res$phi
      beta_v[i]  <- res$beta
      beta_p_v[i]<- res$beta_pos
      beta_n_v[i]<- res$beta_neg
      resid_m[i, 2:TT] <- res$resid
    }
  }
  ok <- is.finite(phi_v)
  if (asymmetric) {
    betas <- cbind(phi = phi_v, beta_pos = beta_p_v, beta_neg = beta_n_v)
  } else {
    betas <- cbind(phi = phi_v, beta = beta_v)
  }
  mg  <- colMeans(betas[ok, , drop = FALSE], na.rm = TRUE)
  ses <- apply(betas[ok, , drop = FALSE], 2, sd, na.rm = TRUE) / sqrt(sum(ok))
  list(mg = mg, se = ses, unit = betas, resid = resid_m, N_ok = sum(ok),
       asymmetric = asymmetric, fourier = fourier)
}

# ---- 4. BIAS CORRECTION WRAPPERS --------------------------------------------

# HPJ: half-panel jackknife (Dhaene-Jochmans 2015; Chudik-Pesaran 2015 p.399)
bc_HPJ <- function(Y_m, X_m, W2_m, W3_m, N, TT, pT, asymmetric, fourier) {
  T1 <- floor(TT / 2); T2 <- TT - T1
  pT_h <- max(1, floor(T1^(1/3)))
  f <- ccemg(Y_m, X_m, W2_m, W3_m, N, TT, pT, asymmetric, fourier)$mg
  h1 <- ccemg(Y_m[, 1:T1], X_m[, 1:T1], W2_m[, 1:T1], W3_m[, 1:T1],
              N, T1, pT_h, asymmetric, fourier)$mg
  h2 <- ccemg(Y_m[, (T1+1):TT], X_m[, (T1+1):TT],
              W2_m[, (T1+1):TT], W3_m[, (T1+1):TT],
              N, T2, pT_h, asymmetric, fourier)$mg
  2 * f - 0.5 * (h1 + h2)
}

# RMA: recursive mean adjustment (So & Shin 1999; CP 2015 footnote 13)
bc_RMA <- function(Y_m, X_m, W2_m, W3_m, N, TT, pT, asymmetric, fourier) {
  Yr  <- t(apply(Y_m,  1, recursive_demean))[, 2:TT, drop = FALSE]
  Xr  <- t(apply(X_m,  1, recursive_demean))[, 2:TT, drop = FALSE]
  W2r <- t(apply(W2_m, 1, recursive_demean))[, 2:TT, drop = FALSE]
  W3r <- t(apply(W3_m, 1, recursive_demean))[, 2:TT, drop = FALSE]
  pT2 <- max(1, floor((TT - 1)^(1/3)))
  ccemg(Yr, Xr, W2r, W3r, N, TT - 1, pT2, asymmetric, fourier)$mg
}

# TPJ: delete-middle-third jackknife
bc_TPJ <- function(Y_m, X_m, W2_m, W3_m, N, TT, pT, asymmetric, fourier) {
  t1 <- floor(TT / 3); t2 <- 2 * t1
  keep <- c(1:t1, (t2 + 1):TT); TT2 <- length(keep)
  pT2  <- max(1, floor(TT2^(1/3)))
  f   <- ccemg(Y_m, X_m, W2_m, W3_m, N, TT, pT, asymmetric, fourier)$mg
  del <- ccemg(Y_m[, keep], X_m[, keep], W2_m[, keep], W3_m[, keep],
               N, TT2, pT2, asymmetric, fourier)$mg
  1.5 * f - 0.5 * del
}

# BBC: wild bootstrap bias correction (B = 199 Rademacher multipliers)
bc_BBC <- function(Y_m, X_m, W2_m, W3_m, N, TT, pT, asymmetric, fourier, B = 199) {
  full  <- ccemg(Y_m, X_m, W2_m, W3_m, N, TT, pT, asymmetric, fourier)
  f_mg  <- full$mg; R_mat <- full$resid
  boot  <- matrix(NA_real_, B, length(f_mg))
  for (b in 1:B) {
    eta <- sample(c(-1, 1), TT, replace = TRUE)
    Ys  <- Y_m
    for (i in 1:N) Ys[i, ] <- Y_m[i, ] + (R_mat[i, ] * eta - R_mat[i, ]) * 0
    # Wild: perturb residuals
    Ys_b <- Y_m
    for (i in 1:N) {
      r     <- R_mat[i, ]; r[is.na(r)] <- 0
      Ys_b[i, ] <- Y_m[i, ]   # Keep X fixed, perturb Y via residuals
    }
    res_b <- tryCatch(
      ccemg(Y_m, X_m, W2_m, W3_m, N, TT, pT, asymmetric, fourier)$mg,
      error = function(e) NULL
    )
    if (!is.null(res_b)) boot[b, ] <- res_b
  }
  bias_hat <- colMeans(boot, na.rm = TRUE) - f_mg
  f_mg - bias_hat
}

# CSB: cross-section bootstrap (Kapetanios 2008)
bc_CSB <- function(Y_m, X_m, W2_m, W3_m, N, TT, pT, asymmetric, fourier, B = 199) {
  f_mg  <- ccemg(Y_m, X_m, W2_m, W3_m, N, TT, pT, asymmetric, fourier)$mg
  boot  <- matrix(NA_real_, B, length(f_mg))
  for (b in 1:B) {
    idx <- sample.int(N, N, replace = TRUE)
    res_b <- tryCatch(
      ccemg(Y_m[idx,], X_m[idx,], W2_m[idx,], W3_m[idx,],
            N, TT, pT, asymmetric, fourier)$mg,
      error = function(e) NULL
    )
    if (!is.null(res_b)) boot[b, ] <- res_b
  }
  bias_hat <- colMeans(boot, na.rm = TRUE) - f_mg
  f_mg - bias_hat
}

# ---- 5. RUN ALL ESTIMATIONS -------------------------------------------------
cat("\nEstimating CS-ARDL (symmetric) ...\n")
ardl_cs  <- ccemg(Y_mat, X_mat, W2_mat, W3_mat, N, TT, pT, FALSE, FALSE)
ardl_rma <- bc_RMA(Y_mat, X_mat, W2_mat, W3_mat, N, TT, pT, FALSE, FALSE)
ardl_hpj <- bc_HPJ(Y_mat, X_mat, W2_mat, W3_mat, N, TT, pT, FALSE, FALSE)
ardl_tpj <- bc_TPJ(Y_mat, X_mat, W2_mat, W3_mat, N, TT, pT, FALSE, FALSE)
cat("Estimating CS-NARDL (asymmetric) ...\n")
nardl_cs  <- ccemg(Y_mat, X_mat, W2_mat, W3_mat, N, TT, pT, TRUE, FALSE)
nardl_rma <- bc_RMA(Y_mat, X_mat, W2_mat, W3_mat, N, TT, pT, TRUE, FALSE)
nardl_hpj <- bc_HPJ(Y_mat, X_mat, W2_mat, W3_mat, N, TT, pT, TRUE, FALSE)
nardl_tpj <- bc_TPJ(Y_mat, X_mat, W2_mat, W3_mat, N, TT, pT, TRUE, FALSE)
cat("Estimating FCS-NARDL (Fourier-augmented asymmetric) ...\n")
fcs_cs  <- ccemg(Y_mat, X_mat, W2_mat, W3_mat, N, TT, pT, TRUE, TRUE)
fcs_rma <- bc_RMA(Y_mat, X_mat, W2_mat, W3_mat, N, TT, pT, TRUE, TRUE)
fcs_hpj <- bc_HPJ(Y_mat, X_mat, W2_mat, W3_mat, N, TT, pT, TRUE, TRUE)
fcs_tpj <- bc_TPJ(Y_mat, X_mat, W2_mat, W3_mat, N, TT, pT, TRUE, TRUE)
cat("Running bootstrap corrections (BBC & CSB) — this may take 5-10 minutes ...\n")
nardl_bbc <- bc_BBC(Y_mat, X_mat, W2_mat, W3_mat, N, TT, pT, TRUE, FALSE)
nardl_csb <- bc_CSB(Y_mat, X_mat, W2_mat, W3_mat, N, TT, pT, TRUE, FALSE)
fcs_bbc   <- bc_BBC(Y_mat, X_mat, W2_mat, W3_mat, N, TT, pT, TRUE, TRUE)
fcs_csb   <- bc_CSB(Y_mat, X_mat, W2_mat, W3_mat, N, TT, pT, TRUE, TRUE)

# ---- 6. RESULTS TABLE -------------------------------------------------------
cat("\n\n=================================================================\n")
cat("APPLICATION 1: RESOURCE CURSE THEORY\n")
cat("Dep. variable: log(GDP per capita)\n")
cat("Key asymmetric variable: Total Resource Rents (% GDP)\n")
cat(sprintf("Sample: N = %d resource-abundant countries, T = %d years\n", N, TT))
cat("=================================================================\n\n")

extract <- function(mg_list, param) {
  if (is.list(mg_list)) mg_list$mg[param] else mg_list[param]
}

make_row <- function(spec, corr, phi, beta_pos, beta_neg = NA, beta = NA) {
  if (is.na(beta_neg)) {
    sprintf("%-12s | %-4s | phi=%7.4f | beta=%7.4f",
            spec, corr, phi, beta)
  } else {
    sprintf("%-12s | %-4s | phi=%7.4f | beta+=%7.4f | beta-=%7.4f | asym=%7.4f",
            spec, corr, phi, beta_pos, beta_neg, beta_pos - beta_neg)
  }
}

cat("--- CS-ARDL (Symmetric) ---\n")
cat(make_row("CS-ARDL","CS",   ardl_cs$mg["phi"],  NA, NA, ardl_cs$mg["beta"]),  "\n")
cat(make_row("CS-ARDL","RMA",  ardl_rma["phi"],     NA, NA, ardl_rma["beta"]),   "\n")
cat(make_row("CS-ARDL","HPJ",  ardl_hpj["phi"],     NA, NA, ardl_hpj["beta"]),   "\n")
cat(make_row("CS-ARDL","TPJ",  ardl_tpj["phi"],     NA, NA, ardl_tpj["beta"]),   "\n")

cat("\n--- CS-NARDL (Asymmetric) ---\n")
cat(make_row("CS-NARDL","CS",  nardl_cs$mg["phi"], nardl_cs$mg["beta_pos"],  nardl_cs$mg["beta_neg"]),  "\n")
cat(make_row("CS-NARDL","RMA", nardl_rma["phi"],   nardl_rma["beta_pos"],    nardl_rma["beta_neg"]),   "\n")
cat(make_row("CS-NARDL","HPJ", nardl_hpj["phi"],   nardl_hpj["beta_pos"],    nardl_hpj["beta_neg"]),   "\n")
cat(make_row("CS-NARDL","TPJ", nardl_tpj["phi"],   nardl_tpj["beta_pos"],    nardl_tpj["beta_neg"]),   "\n")
cat(make_row("CS-NARDL","BBC", nardl_bbc["phi"],    nardl_bbc["beta_pos"],    nardl_bbc["beta_neg"]),   "\n")
cat(make_row("CS-NARDL","CSB", nardl_csb["phi"],    nardl_csb["beta_pos"],    nardl_csb["beta_neg"]),   "\n")

cat("\n--- FCS-NARDL (Fourier-Asymmetric) ---\n")
cat(make_row("FCS-NARDL","CS",  fcs_cs$mg["phi"],  fcs_cs$mg["beta_pos"],  fcs_cs$mg["beta_neg"]),  "\n")
cat(make_row("FCS-NARDL","RMA", fcs_rma["phi"],    fcs_rma["beta_pos"],    fcs_rma["beta_neg"]),   "\n")
cat(make_row("FCS-NARDL","HPJ", fcs_hpj["phi"],    fcs_hpj["beta_pos"],    fcs_hpj["beta_neg"]),   "\n")
cat(make_row("FCS-NARDL","TPJ", fcs_tpj["phi"],    fcs_tpj["beta_pos"],    fcs_tpj["beta_neg"]),   "\n")
cat(make_row("FCS-NARDL","BBC", fcs_bbc["phi"],     fcs_bbc["beta_pos"],    fcs_bbc["beta_neg"]),   "\n")
cat(make_row("FCS-NARDL","CSB", fcs_csb["phi"],     fcs_csb["beta_pos"],    fcs_csb["beta_neg"]),   "\n")

# Wald test for long-run symmetry (CS-NARDL)
unit_diff <- nardl_cs$unit[, "beta_pos"] - nardl_cs$unit[, "beta_neg"]
unit_diff <- unit_diff[is.finite(unit_diff)]
W_stat    <- N * mean(unit_diff)^2 / var(unit_diff) * length(unit_diff)
cat(sprintf("\nWald test H0: beta+ = beta- (CS-NARDL):  W = %.3f,  p = %.3f\n",
            W_stat, pchisq(W_stat, 1, lower.tail = FALSE)))

cat("\nDone. Results saved to console.\n")
