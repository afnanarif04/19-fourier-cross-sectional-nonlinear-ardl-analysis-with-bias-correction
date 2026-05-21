# =============================================================================
# 22_empirical_application_2.R
# Empirical Application 2 - Income and CO2 emissions (environmental Kuznets curve)
#
# Estimates the symmetric (CS-ARDL), asymmetric (CS-NARDL), and
# Fourier-augmented asymmetric (FCS-NARDL) specifications under six correctors
# (CS, RMA, HPJ, TPJ, BBC, CSB) on a global country panel, 1990-2019.
#
# DATA (open access):
#   CO2 per capita and GDP : Our World in Data - CO2 dataset
#       File  : owid_co2_data.csv
#       Source: https://github.com/owid/co2-data
#       Direct: https://nyc3.digitaloceanspaces.com/owid-public/data/co2/owid-co2-data.csv
#   Energy per capita      : Our World in Data - Energy dataset
#       File  : owid_energy_data.csv
#       Source: https://github.com/owid/energy-data
#       Direct: https://nyc3.digitaloceanspaces.com/owid-public/data/energy/owid-energy-data.csv
#   Urban share, trade     : World Bank WDI (small, fetched at runtime)
#
# SETUP: download the two OWID CSV files from the links above and place them in
#        the ../data/ folder (they are not bundled because of their size).
#
# Requires : internet connection (WDI download); the two OWID CSV files present
# Packages : WDI, dplyr, tidyr
# Output   : printed coefficient tables and the long-run symmetry Wald test
#
# Run order: see 00_RUN_ALL.R
# =============================================================================

rm(list = ls())

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

# ---- 1. LOCAL OWID DATA (no internet needed) --------------------------------
# Place owid_co2_data.csv and owid_energy_data.csv in the SAME folder as this script.
# Source: https://github.com/owid/co2-data  and  https://github.com/owid/energy-data

co2_file <- "../data/owid_co2_data.csv"
en_file  <- "../data/owid_energy_data.csv"

if (!file.exists(co2_file)) stop(
  sprintf("File '%s' not found in working directory. Please put it next to this R script.", co2_file))
if (!file.exists(en_file)) stop(
  sprintf("File '%s' not found in working directory. Please put it next to this R script.", en_file))

cat("Reading OWID CO2 file ...\n")
co2 <- read.csv(co2_file, stringsAsFactors = FALSE) %>%
  filter(year %in% 1990:2019, !is.na(iso_code), nchar(iso_code) == 3) %>%
  select(iso = iso_code, country, year,
         co2_pc = co2_per_capita,
         gdp    = gdp,
         pop    = population)

cat("Reading OWID Energy file ...\n")
en <- read.csv(en_file, stringsAsFactors = FALSE) %>%
  filter(year %in% 1990:2019, !is.na(iso_code), nchar(iso_code) == 3) %>%
  select(iso = iso_code, year,
         energy_pc = energy_per_capita)   # in kWh/person

cat("Downloading urban + trade from WDI (small, fast) ...\n")
wdi_data <- tryCatch(
  WDI(country = "all",
      indicator = c(urban = "SP.URB.TOTL.IN.ZS", trade = "NE.TRD.GNFS.ZS"),
      start = 1990, end = 2019, extra = TRUE),
  error = function(e) { cat("WDI failed - rerun in a minute or use cached file.\n"); stop(e) }
)
wdi <- wdi_data %>%
  filter(!is.na(iso3c), nchar(iso3c) == 3, region != "Aggregates") %>%
  select(iso = iso3c, year, urban, trade)

# ---- 2. MERGE & PANEL CONSTRUCTION ------------------------------------------
df <- co2 %>%
  inner_join(en,  by = c("iso", "year")) %>%
  inner_join(wdi, by = c("iso", "year")) %>%
  mutate(gdp_pc = gdp / pop,
         ln_co2 = log(pmax(co2_pc, 0.001)),
         ln_gdp = log(pmax(gdp_pc, 1)),
         ln_en  = log(pmax(energy_pc, 1)),
         ln_urb = log(pmax(urban, 1)),
         ln_tr  = log(pmax(trade, 1))) %>%
  filter(is.finite(ln_co2), is.finite(ln_gdp), is.finite(ln_en),
         is.finite(ln_urb), is.finite(ln_tr)) %>%
  select(iso, year, ln_co2, ln_gdp, ln_en, ln_urb, ln_tr) %>%
  arrange(iso, year)

ids_ok <- df %>% group_by(iso) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n >= 28) %>% pull(iso)

panel <- df %>%
  filter(iso %in% ids_ok) %>%
  group_by(iso) %>% arrange(year) %>%
  tidyr::fill(everything(), .direction = "downup") %>% ungroup()

# Keep only countries with all 30 years after fill
final_iso <- panel %>% group_by(iso) %>% summarise(n=n(), .groups="drop") %>%
  filter(n == 30) %>% pull(iso)
panel <- panel %>% filter(iso %in% final_iso)

ids   <- sort(unique(panel$iso))
years <- sort(unique(panel$year))
N     <- length(ids); TT <- length(years); pT <- max(1, floor(TT^(1/3)))
cat(sprintf("\nPanel ready: N = %d, T = %d (1990-2019), pT = %d\n", N, TT, pT))

to_mat <- function(v) {
  m <- matrix(NA_real_, N, TT)
  for (j in seq_along(ids)) {
    sub <- panel[panel$iso == ids[j], ]
    m[j, ] <- sub[order(sub$year), v, drop = TRUE]
  }
  m
}
Y_mat  <- to_mat("ln_co2")
X_mat  <- to_mat("ln_gdp")    # asymmetric variable
W2_mat <- to_mat("ln_en")
W3_mat <- to_mat("ln_urb")
W4_mat <- to_mat("ln_tr")

# ---- 3. ESTIMATION FUNCTIONS ------------------------------------------------
partial_sums <- function(x){ dx<-c(0,diff(x)); list(pos=cumsum(pmax(dx,0)),neg=cumsum(pmin(dx,0))) }
recursive_demean <- function(z){ T<-length(z); out<-rep(NA_real_,T); cs<-0
  for(t in 2:T){ cs<-cs+z[t-1]; out[t]<-z[t]-cs/(t-1) }; out }
build_csa <- function(Ym,Xm,W2m,W3m,W4m,pT,TT){
  W <- cbind(colMeans(Ym,na.rm=TRUE), colMeans(Xm,na.rm=TRUE),
             colMeans(W2m,na.rm=TRUE), colMeans(W3m,na.rm=TRUE),
             colMeans(W4m,na.rm=TRUE))
  lapply(0:pT, function(l) if(l==0)W else rbind(matrix(NA,l,ncol(W)),W[1:(TT-l),]))
}
run_unit <- function(y,x,w2,w3,w4,csa,fs,fc,asym){
  T <- length(y); t1 <- 2:T
  Dy<-diff(y);Dx<-diff(x);Dw2<-diff(w2);Dw3<-diff(w3);Dw4<-diff(w4)
  yl<-y[-T];xl<-x[-T];w2l<-w2[-T];w3l<-w3[-T];w4l<-w4[-T]
  cb <- do.call(cbind, lapply(csa, function(m){ mm<-m[t1,,drop=FALSE]; mm[is.na(mm)]<-0; mm }))
  if(!asym){
    Z <- cbind(1,yl,xl,w2l,w3l,w4l,Dx,Dw2,Dw3,Dw4,cb)
    if(!is.null(fs)) Z <- cbind(Z,fs[t1],fc[t1])
    fit <- tryCatch(lm.fit(Z,Dy), error=function(e)NULL); if(is.null(fit)) return(NULL)
    b<-fit$coefficients; phi<-b[2]; delta<-b[3]
    if(!is.finite(phi)||phi>=0||phi< -2) return(NULL)
    list(phi=phi, beta=-delta/phi, bp=NA, bn=NA, resid=fit$residuals)
  } else {
    ps<-partial_sums(x); xpl<-ps$pos[-T]; xnl<-ps$neg[-T]
    Dxp<-diff(ps$pos); Dxn<-diff(ps$neg)
    Z <- cbind(1,yl,xpl,xnl,w2l,w3l,w4l,Dxp,Dxn,Dw2,Dw3,Dw4,cb)
    if(!is.null(fs)) Z <- cbind(Z,fs[t1],fc[t1])
    fit <- tryCatch(lm.fit(Z,Dy), error=function(e)NULL); if(is.null(fit)) return(NULL)
    b<-fit$coefficients; phi<-b[2]; dp<-b[3]; dn<-b[4]
    if(!is.finite(phi)||phi>=0||phi< -2) return(NULL)
    list(phi=phi, beta=NA, bp=-dp/phi, bn=-dn/phi, resid=fit$residuals)
  }
}
ccemg <- function(Ym,Xm,W2m,W3m,W4m,N,TT,pT,asym,fou){
  csa<-build_csa(Ym,Xm,W2m,W3m,W4m,pT,TT)
  ts<-if(fou)sin(2*pi*(1:TT)/TT) else NULL
  tc<-if(fou)cos(2*pi*(1:TT)/TT) else NULL
  pv<-bv<-bpv<-bnv<-rep(NA_real_,N); rm_<-matrix(NA_real_,N,TT)
  for(i in 1:N){
    res <- run_unit(Ym[i,],Xm[i,],W2m[i,],W3m[i,],W4m[i,],csa,ts,tc,asym)
    if(!is.null(res)){pv[i]<-res$phi;bv[i]<-res$beta;bpv[i]<-res$bp;bnv[i]<-res$bn
                       rm_[i,2:TT]<-res$resid}
  }
  ok<-is.finite(pv)
  bt <- if(asym) cbind(phi=pv,beta_pos=bpv,beta_neg=bnv) else cbind(phi=pv,beta=bv)
  mg <- colMeans(bt[ok,,drop=FALSE], na.rm=TRUE)
  list(mg=mg, unit=bt, resid=rm_, N_ok=sum(ok))
}
bc_HPJ <- function(Ym,Xm,W2m,W3m,W4m,N,TT,pT,a,f){
  T1<-floor(TT/2); T2<-TT-T1; pTh<-max(1,floor(T1^(1/3)))
  sl<-function(m,s,e)m[,s:e,drop=FALSE]
  full<-ccemg(Ym,Xm,W2m,W3m,W4m,N,TT,pT,a,f)$mg
  h1<-ccemg(sl(Ym,1,T1),sl(Xm,1,T1),sl(W2m,1,T1),sl(W3m,1,T1),sl(W4m,1,T1),N,T1,pTh,a,f)$mg
  h2<-ccemg(sl(Ym,T1+1,TT),sl(Xm,T1+1,TT),sl(W2m,T1+1,TT),sl(W3m,T1+1,TT),sl(W4m,T1+1,TT),N,T2,pTh,a,f)$mg
  2*full - 0.5*(h1+h2)
}
bc_RMA <- function(Ym,Xm,W2m,W3m,W4m,N,TT,pT,a,f){
  rd<-function(m)t(apply(m,1,recursive_demean))[,2:TT,drop=FALSE]
  ccemg(rd(Ym),rd(Xm),rd(W2m),rd(W3m),rd(W4m),N,TT-1,max(1,floor((TT-1)^(1/3))),a,f)$mg
}
bc_TPJ <- function(Ym,Xm,W2m,W3m,W4m,N,TT,pT,a,f){
  t1<-floor(TT/3); keep<-c(1:t1,(2*t1+1):TT); TT2<-length(keep)
  sl<-function(m)m[,keep,drop=FALSE]
  1.5*ccemg(Ym,Xm,W2m,W3m,W4m,N,TT,pT,a,f)$mg -
    0.5*ccemg(sl(Ym),sl(Xm),sl(W2m),sl(W3m),sl(W4m),N,TT2,max(1,floor(TT2^(1/3))),a,f)$mg
}
bc_CSB <- function(Ym,Xm,W2m,W3m,W4m,N,TT,pT,a,f,B=199){
  fm<-ccemg(Ym,Xm,W2m,W3m,W4m,N,TT,pT,a,f)$mg
  bt<-matrix(NA_real_,B,length(fm))
  for(b in 1:B){
    idx<-sample.int(N,N,replace=TRUE)
    r<-tryCatch(ccemg(Ym[idx,],Xm[idx,],W2m[idx,],W3m[idx,],W4m[idx,],N,TT,pT,a,f)$mg,
                error=function(e)NULL)
    if(!is.null(r)) bt[b,]<-r
  }
  fm - (colMeans(bt,na.rm=TRUE) - fm)
}

# ---- 4. RUN ALL ESTIMATIONS -------------------------------------------------
cat("\nCS-ARDL ...\n")
a_cs<-ccemg(Y_mat,X_mat,W2_mat,W3_mat,W4_mat,N,TT,pT,FALSE,FALSE)
a_rma<-bc_RMA(Y_mat,X_mat,W2_mat,W3_mat,W4_mat,N,TT,pT,FALSE,FALSE)
a_hpj<-bc_HPJ(Y_mat,X_mat,W2_mat,W3_mat,W4_mat,N,TT,pT,FALSE,FALSE)
a_tpj<-bc_TPJ(Y_mat,X_mat,W2_mat,W3_mat,W4_mat,N,TT,pT,FALSE,FALSE)
cat("CS-NARDL ...\n")
n_cs<-ccemg(Y_mat,X_mat,W2_mat,W3_mat,W4_mat,N,TT,pT,TRUE,FALSE)
n_rma<-bc_RMA(Y_mat,X_mat,W2_mat,W3_mat,W4_mat,N,TT,pT,TRUE,FALSE)
n_hpj<-bc_HPJ(Y_mat,X_mat,W2_mat,W3_mat,W4_mat,N,TT,pT,TRUE,FALSE)
n_tpj<-bc_TPJ(Y_mat,X_mat,W2_mat,W3_mat,W4_mat,N,TT,pT,TRUE,FALSE)
n_csb<-bc_CSB(Y_mat,X_mat,W2_mat,W3_mat,W4_mat,N,TT,pT,TRUE,FALSE)
cat("FCS-NARDL ...\n")
f_cs<-ccemg(Y_mat,X_mat,W2_mat,W3_mat,W4_mat,N,TT,pT,TRUE,TRUE)
f_rma<-bc_RMA(Y_mat,X_mat,W2_mat,W3_mat,W4_mat,N,TT,pT,TRUE,TRUE)
f_hpj<-bc_HPJ(Y_mat,X_mat,W2_mat,W3_mat,W4_mat,N,TT,pT,TRUE,TRUE)
f_tpj<-bc_TPJ(Y_mat,X_mat,W2_mat,W3_mat,W4_mat,N,TT,pT,TRUE,TRUE)
f_csb<-bc_CSB(Y_mat,X_mat,W2_mat,W3_mat,W4_mat,N,TT,pT,TRUE,TRUE)

# ---- 5. OUTPUT --------------------------------------------------------------
cat("\n=================================================================\n")
cat("APPLICATION 2: ENERGY-AUGMENTED EKC (FIXED, OWID DATA)\n")
cat("Dep. variable: log(CO2 per capita)  |  Asymmetric: log(GDP per capita)\n")
cat("Controls: log(energy_pc), log(urban), log(trade)\n")
cat(sprintf("Sample: N = %d, T = %d (1990-2019)\n", N, TT))
cat("=================================================================\n")
pr <- function(spec,corr,mg){
  if("beta"%in%names(mg))
    cat(sprintf("%-12s|%-4s| phi=%7.4f  beta(GDP)=%7.4f\n",spec,corr,mg["phi"],mg["beta"]))
  else
    cat(sprintf("%-12s|%-4s| phi=%7.4f  GDP+=%7.4f  GDP-=%7.4f  asymm=%7.4f\n",
                spec,corr,mg["phi"],mg["beta_pos"],mg["beta_neg"],mg["beta_pos"]-mg["beta_neg"]))
}
cat("\n-- CS-ARDL --\n");pr("CS-ARDL","CS",a_cs$mg);pr("CS-ARDL","RMA",a_rma);pr("CS-ARDL","HPJ",a_hpj);pr("CS-ARDL","TPJ",a_tpj)
cat("\n-- CS-NARDL --\n");pr("CS-NARDL","CS",n_cs$mg);pr("CS-NARDL","RMA",n_rma);pr("CS-NARDL","HPJ",n_hpj);pr("CS-NARDL","TPJ",n_tpj);pr("CS-NARDL","CSB",n_csb)
cat("\n-- FCS-NARDL --\n");pr("FCS-NARDL","CS",f_cs$mg);pr("FCS-NARDL","RMA",f_rma);pr("FCS-NARDL","HPJ",f_hpj);pr("FCS-NARDL","TPJ",f_tpj);pr("FCS-NARDL","CSB",f_csb)

ud<-n_cs$unit[,"beta_pos"]-n_cs$unit[,"beta_neg"]; ud<-ud[is.finite(ud)]
W<-length(ud)*mean(ud)^2/var(ud)
cat(sprintf("\nWald H0: GDP+=GDP- (CS-NARDL): W = %.3f, p = %.3f\n", W, pchisq(W,1,lower.tail=FALSE)))
cat("Done.\n")
