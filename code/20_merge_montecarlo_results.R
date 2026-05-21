# =============================================================================
# 20_merge_montecarlo_results.R
# Merge the 18 Monte Carlo result files into three consolidated tables
#
# Combines the per-method result files produced by scripts 02 through 19 into
# one final CSV per model specification. Each final CSV has 25 rows (one per
# N,T grid cell) with the bias, RMSE, and percent-bias columns of all six
# correctors side by side.
#
# Input    : 02_mc_csardl_1_CS_results.csv  ...  19_mc_fcsnardl_6_CSB_results.csv
#            (run scripts 02-19 first; or use the copies in ../results/montecarlo)
# Output   : CS_ARDL_FINAL.csv, CS_NARDL_FINAL.csv, FCS_NARDL_FINAL.csv
#
# Run order: see 00_RUN_ALL.R (run after all Monte Carlo scripts)
# =============================================================================

cat("================================================================\n")
cat(" MERGE ALL MONTE CARLO RESULTS\n")
cat("================================================================\n\n")

# ------------------------------------------------------------------
# Helper: read + fix column names + merge
# ------------------------------------------------------------------
merge_model <- function(model_name, file_list) {
  cat(sprintf("\n--- Merging %s ---\n", model_name))
  dfs <- list()

  for (f in file_list) {
    fname <- paste0(f, ".csv")
    if (!file.exists(fname)) {
      cat(sprintf("  SKIP (not found): %s\n", fname))
      next
    }
    d <- read.csv(fname, check.names = FALSE)
    # Fix R's auto-renaming: N.N -> N, T.T -> T
    names(d)[names(d) == "N.N"] <- "N"
    names(d)[names(d) == "T.T"] <- "T"
    cat(sprintf("  Read: %-40s  (%d rows, %d cols)\n",
                fname, nrow(d), ncol(d)))
    dfs[[length(dfs) + 1]] <- d
  }

  if (length(dfs) == 0) {
    cat(sprintf("  ERROR: no files found for %s\n", model_name))
    return(NULL)
  }

  # Start with N,T from first file
  final <- dfs[[1]][, c("N", "T")]

  # Append method-specific columns from each file

  for (k in seq_along(dfs)) {
    cols <- dfs[[k]][, !(names(dfs[[k]]) %in% c("N", "T")), drop = FALSE]
    final <- cbind(final, cols)
  }

  cat(sprintf("  MERGED: %d rows x %d columns\n", nrow(final), ncol(final)))
  return(final)
}

# ==================================================================
# CS-ARDL
# ==================================================================
cs_ardl <- merge_model("CS-ARDL", c(
  "02_mc_csardl_1_CS_results",
  "03_mc_csardl_2_RMA_results",
  "04_mc_csardl_3_HPJ_results",
  "05_mc_csardl_4_TPJ_results",
  "06_mc_csardl_5_BBC_results",
  "07_mc_csardl_6_CSB_results"
))
if (!is.null(cs_ardl)) {
  write.csv(cs_ardl, "CS_ARDL_FINAL.csv", row.names = FALSE)
  cat("  -> Saved: CS_ARDL_FINAL.csv\n")
}

# ==================================================================
# CS-NARDL
# ==================================================================
cs_nardl <- merge_model("CS-NARDL", c(
  "08_mc_csnardl_1_CS_results",
  "09_mc_csnardl_2_RMA_results",
  "10_mc_csnardl_3_HPJ_results",
  "11_mc_csnardl_4_TPJ_results",
  "12_mc_csnardl_5_BBC_results",
  "13_mc_csnardl_6_CSB_results"
))
if (!is.null(cs_nardl)) {
  write.csv(cs_nardl, "CS_NARDL_FINAL.csv", row.names = FALSE)
  cat("  -> Saved: CS_NARDL_FINAL.csv\n")
}

# ==================================================================
# FCS-NARDL
# ==================================================================
fcs_nardl <- merge_model("FCS-NARDL", c(
  "14_mc_fcsnardl_1_CS_results",
  "15_mc_fcsnardl_2_RMA_results",
  "16_mc_fcsnardl_3_HPJ_results",
  "17_mc_fcsnardl_4_TPJ_results",
  "18_mc_fcsnardl_5_BBC_results",
  "19_mc_fcsnardl_6_CSB_results"
))
if (!is.null(fcs_nardl)) {
  write.csv(fcs_nardl, "FCS_NARDL_FINAL.csv", row.names = FALSE)
  cat("  -> Saved: FCS_NARDL_FINAL.csv\n")
}

cat("\n================================================================\n")
cat(" DONE.  Upload these 3 CSVs for table generation:\n")
cat("   1. CS_ARDL_FINAL.csv\n")
cat("   2. CS_NARDL_FINAL.csv\n")
cat("   3. FCS_NARDL_FINAL.csv\n")
cat("================================================================\n")
