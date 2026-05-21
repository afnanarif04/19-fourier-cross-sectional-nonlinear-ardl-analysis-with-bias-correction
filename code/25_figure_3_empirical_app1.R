# =============================================================================
# 25_figure_3_empirical_app1.R
# Figure 3 - Country-level FCS-NARDL long-run elasticities, Application 1
#
# Two-panel line-and-dot figure. Each coloured line tracks one country's
# unit-level long-run elasticity across the three correctors (CS, RMA, HPJ).
# The thick black dashed line is the panel mean. Panel (a) is the positive
# partial sum (resource boom); panel (b) is the negative partial sum (bust).
#
# Input    : ../results/empirical/01_App1_estimates.rds
# Output   : ../figures/figure_3_empirical_app1.png
# Packages : base R only
#
# Run order: see 00_RUN_ALL.R (run after 21_empirical_application_1.R, or use
#            the bundled pre-computed estimates in ../results/empirical/)
# =============================================================================

rm(list = ls())

# ---- 1. Load pre-computed unit-level estimates ------------------------------
est <- readRDS("../results/empirical/01_App1_estimates.rds")

correctors <- c("CS", "RMA", "HPJ")
keys       <- c("F_CS", "F_RMA", "F_HPJ")   # FCS-NARDL only

# ---- 2. Helper: assemble a unit-by-corrector matrix for one coefficient -----
build_matrix <- function(est, keys, field) {
  vecs <- lapply(keys, function(k) as.numeric(est[[k]][[field]]))
  n    <- min(sapply(vecs, length))           # common set of converged units
  m    <- sapply(vecs, function(v) v[1:n])
  colnames(m) <- correctors
  m
}

# ---- 3. Helper: pick 10 representative units spanning the distribution ------
pick_units <- function(m, n_show = 10) {
  base   <- m[, "CS"]
  ord    <- order(base)
  lo     <- floor(0.02 * length(ord)) + 1     # drop extreme 2% tails
  hi     <- ceiling(0.98 * length(ord))
  valid  <- ord[lo:hi]
  qs     <- seq(0.05, 0.95, length.out = n_show)
  valid[round(qs * (length(valid) - 1)) + 1]
}

# ---- 4. Helper: trim y-axis to exclude extreme outliers ---------------------
ylim_trim <- function(vals, lo = 0.00, hi = 1.00) {
  r <- quantile(vals, c(lo, hi), na.rm = TRUE)
  pad <- 0.10 * diff(r); if (!is.finite(pad) || pad == 0) pad <- 0.01
  c(r[1] - pad, r[2] + pad)
}

# ---- 5. Plot ----------------------------------------------------------------
cols <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd",
          "#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf")
pchs <- c(16,15,18,17,25,1,0,5,2,6)

png("../figures/figure_3_empirical_app1.png",
    width = 11, height = 4.5, units = "in", res = 300)
par(mfrow = c(1, 2), family = "sans", mar = c(4.2, 4.4, 2.2, 1.0),
    oma = c(0, 0, 2.4, 0), mgp = c(2.6, 0.7, 0), bg = "white")

panels <- list(
  list(field = "unit_LR_pos", mean = "LR_pos_mean",
       ylab = expression("Long-run elasticity " * beta^"+"),
       main = "(a) Resource boom long-run effect, by country"),
  list(field = "unit_LR_neg", mean = "LR_neg_mean",
       ylab = expression("Long-run elasticity " * beta^"-"),
       main = "(b) Resource bust long-run effect, by country")
)

x <- 1:3
for (pn in panels) {
  m    <- build_matrix(est, keys, pn$field)
  pick <- pick_units(m, 10)
  yl   <- ylim_trim(as.vector(m[pick, ]))

  plot(NA, xlim = c(0.9, 3.1), ylim = yl, xaxt = "n",
       xlab = "Bias-correction method", ylab = pn$ylab, main = pn$main, cex.main = 0.95)
  axis(1, at = x, labels = correctors)
  grid(nx = NA, ny = NULL, col = "grey88", lty = 1)
  abline(h = 0, col = "grey55", lty = 3)

  for (j in seq_along(pick)) {
    lines(x, m[pick[j], ], col = cols[j], lwd = 1.6, type = "b",
          pch = pchs[j], cex = 1.1)
  }
  # panel mean
  pm <- sapply(keys, function(k) est[[k]][[pn$mean]])
  lines(x, pm, col = "black", lwd = 2.6, lty = 2, type = "b", pch = 4, cex = 1.4)

  leg <- paste0("Country ", seq_along(pick))
  legend("topleft", legend = c(leg, "Panel mean"),
         col = c(cols[seq_along(pick)], "black"),
         pch = c(pchs[seq_along(pick)], 4),
         lty = c(rep(1, length(pick)), 2), lwd = c(rep(1.6, length(pick)), 2.6),
         ncol = 2, bty = "n", cex = 0.62)
}

mtext("Figure 3. Resource Curse Application: country-level FCS-NARDL long-run elasticities under three correctors",
      outer = TRUE, line = 0.4, cex = 0.92, font = 2)
dev.off()
cat("Saved: ../figures/figure_3_empirical_app1.png\n")
