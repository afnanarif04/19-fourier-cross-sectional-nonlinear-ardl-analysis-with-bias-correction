# =============================================================================
# 23_figure_1_montecarlo.R
# Figure 1 - Comprehensive Monte Carlo bias-correction results (CS-ARDL focus)
#
# Four-panel figure: convergence of percent bias with T, method comparison,
# and a heatmap across the (N, T) grid.
#
# Input    : ../results/montecarlo/CS_ARDL_FINAL.csv
#            ../results/montecarlo/CS_NARDL_FINAL.csv
#            ../results/montecarlo/FCS_NARDL_FINAL.csv
# Output   : ../figures/figure_1_montecarlo.pdf
# Packages : fields, viridis
#
# Run order: see 00_RUN_ALL.R (run after 20_merge_montecarlo_results.R)
# =============================================================================

# Set the working directory to this script's /code folder before running.

# Load required libraries
library(fields)   # for image.plot (MUST LOAD THIS!)
library(viridis)  # for color schemes

# Load data
cs_ardl <- read.csv("../results/montecarlo/CS_ARDL_FINAL.csv")
cs_nardl <- read.csv("../results/montecarlo/CS_NARDL_FINAL.csv")
fcs_nardl <- read.csv("../results/montecarlo/FCS_NARDL_FINAL.csv")

# Set publication-quality graphics parameters
pdf("../figures/figure_1_montecarlo.pdf", width = 12, height = 10)

# Create 2x2 layout
par(mfrow = c(2, 2),
    family = "serif",
    mar = c(4.5, 4.5, 2.5, 1),
    mgp = c(3, 0.7, 0),
    cex.lab = 1.2,
    cex.axis = 1.0,
    cex.main = 1.3,
    las = 1)

# Extract diagonal cells
diagonal_indices <- c(1, 7, 13, 19, 25)
sample_sizes <- c(40, 50, 100, 150, 200)

# Method columns
method_cols_bias <- c(2, 8, 14, 20, 26, 32)
method_cols_rmse <- c(3, 9, 15, 21, 27, 33)
method_names <- c("CS", "RMA", "HPJ", "TPJ", "BBC", "CSB")
method_colors <- c("black", "blue", "red", "darkgreen", "orange", "purple")

################################################################################
# PANEL A: BIAS CONVERGENCE
################################################################################

plot(sample_sizes, cs_ardl[diagonal_indices, method_cols_bias[1]],
     type = "n",
     ylim = c(-20, 60),
     xlab = "Sample Size (N = T)",
     ylab = "Bias (%)",
     main = "(a) Bias Convergence Across Methods",
     log = "x",
     xaxt = "n")
axis(1, at = sample_sizes, labels = sample_sizes)
abline(h = 0, lty = 2, col = "gray50", lwd = 2)
grid(col = "gray90", lty = 1)

# Plot each method
for (i in 1:6) {
  bias_values <- cs_ardl[diagonal_indices, method_cols_bias[i]]
  lines(sample_sizes, bias_values,
        col = method_colors[i], lwd = 2.5,
        pch = 15 + i, type = "b", cex = 1.3)
}

legend("topright", legend = method_names,
       col = method_colors, lwd = 2.5, pch = 16:21,
       bty = "n", cex = 1.0, ncol = 2)

# Add annotation
text(45, 55, "phi parameter", pos = 4, cex = 1.0, font = 3)

################################################################################
# PANEL B: RMSE CONVERGENCE
################################################################################

plot(sample_sizes, cs_ardl[diagonal_indices, method_cols_rmse[1]],
     type = "n",
     ylim = c(0, 60),
     xlab = "Sample Size (N = T)",
     ylab = "RMSE (%)",
     main = "(b) RMSE Convergence Across Methods",
     log = "x",
     xaxt = "n")
axis(1, at = sample_sizes, labels = sample_sizes)
grid(col = "gray90", lty = 1)

for (i in 1:6) {
  rmse_values <- cs_ardl[diagonal_indices, method_cols_rmse[i]]
  lines(sample_sizes, rmse_values,
        col = method_colors[i], lwd = 2.5,
        pch = 15 + i, type = "b", cex = 1.3)
}

legend("topright", legend = method_names,
       col = method_colors, lwd = 2.5, pch = 16:21,
       bty = "n", cex = 1.0, ncol = 2)

################################################################################
# PANEL C: BIAS HEATMAP (CS ESTIMATOR)
################################################################################

# Create bias matrix
N_vals <- c(40, 50, 100, 150, 200)
T_vals <- c(40, 50, 100, 150, 200)
bias_matrix <- matrix(cs_ardl[, 2], nrow = 5, ncol = 5, byrow = FALSE)

# Set margins for heatmap with space for colorbar
par(mar = c(4.5, 4.5, 2.5, 5))

image(1:5, 1:5, bias_matrix,
      xlab = "N",
      ylab = "T",
      main = "(c) Bias Surface: CS Estimator",
      col = viridis(100),
      xaxt = "n", yaxt = "n")
axis(1, at = 1:5, labels = N_vals)
axis(2, at = 1:5, labels = T_vals)
contour(1:5, 1:5, bias_matrix, add = TRUE,
        col = "white", lwd = 1, labcex = 0.8)

# Add color scale (THIS NOW WORKS!)
image.plot(legend.only = TRUE,
           zlim = range(bias_matrix),
           col = viridis(100),
           legend.shrink = 0.8,
           legend.width = 1.5,
           axis.args = list(cex.axis = 0.9),
           legend.args = list(text = "Bias (%)",
                              side = 3,
                              line = 0.5,
                              cex = 1.0))

# Reset margins
par(mar = c(4.5, 4.5, 2.5, 1))

################################################################################
# PANEL D: METHOD COMPARISON AT N=T=100
################################################################################

idx_medium <- 13  # (100,100)
bias_values <- as.numeric(cs_ardl[idx_medium, method_cols_bias])

bp <- barplot(bias_values,
              names.arg = method_names,
              col = c("gray30", "steelblue", "coral", "darkgreen", "gold", "orchid"),
              border = "black",
              ylim = c(-8, 3),
              ylab = "Bias (%)",
              main = "(d) Method Comparison at N=T=100",
              las = 2,
              cex.names = 1.0)
abline(h = 0, lwd = 2)
grid(nx = NA, ny = NULL, col = "gray80", lty = 2)

# Add values on bars
text(bp, bias_values, labels = sprintf("%.1f", bias_values),
     pos = ifelse(bias_values > 0, 3, 1), cex = 0.9)

dev.off()

cat("\n")
cat(paste(rep("=", 70), collapse=""), "\n")
cat("FIGURE 1 CREATED: Figure1_CS_ARDL_Comprehensive.pdf\n")
cat(paste(rep("=", 70), collapse=""), "\n")
cat("\nFour-panel figure showing:\n")
cat("  (a) Bias convergence across 6 methods\n")
cat("  (b) RMSE convergence across 6 methods\n")
cat("  (c) Bias heatmap for CS estimator across (N,T) grid\n")
cat("  (d) Method comparison at N=T=100\n")
cat("\nRecommendation: Use this as your main Monte Carlo figure in text\n")
cat("\n")
