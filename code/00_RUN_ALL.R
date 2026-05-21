# =============================================================================
# 00_RUN_ALL.R
# Master script - reproduces every result, table, and figure in order
#
# This package is organised so that each script can also be run on its own in
# RStudio. To reproduce everything from scratch, set the working directory to
# the /code folder and source this file. The Monte Carlo block (scripts 02-19)
# is the slow part; budget several hours for the full 5x5 grid at R = 500.
#
# Folder layout expected:
#   code/      all R scripts (this folder)
#   data/      place the two OWID CSV files here (see README, Section "Data")
#   results/   pre-computed outputs are bundled here; re-running overwrites them
#   figures/   all figures are written here
#
# To reproduce only the figures from the bundled results, skip Parts 1-3 and
# run Part 4 directly.
# =============================================================================

# Make sure the working directory is the /code folder before sourcing.
# In RStudio:  Session > Set Working Directory > To Source File Location

# ---- Part 1. Monte Carlo simulations (slow; writes 18 result CSVs) ----------
# Each script is self-contained and uses a fixed seed (set.seed(2026)).
source("02_mc_csardl_1_CS.R")
source("03_mc_csardl_2_RMA.R")
source("04_mc_csardl_3_HPJ.R")
source("05_mc_csardl_4_TPJ.R")
source("06_mc_csardl_5_BBC.R")
source("07_mc_csardl_6_CSB.R")

source("08_mc_csnardl_1_CS.R")
source("09_mc_csnardl_2_RMA.R")
source("10_mc_csnardl_3_HPJ.R")
source("11_mc_csnardl_4_TPJ.R")
source("12_mc_csnardl_5_BBC.R")
source("13_mc_csnardl_6_CSB.R")

source("14_mc_fcsnardl_1_CS.R")
source("15_mc_fcsnardl_2_RMA.R")
source("16_mc_fcsnardl_3_HPJ.R")
source("17_mc_fcsnardl_4_TPJ.R")
source("18_mc_fcsnardl_5_BBC.R")
source("19_mc_fcsnardl_6_CSB.R")

# ---- Part 2. Merge Monte Carlo output into three final tables ---------------
source("20_merge_montecarlo_results.R")

# ---- Part 3. Empirical applications (need internet; App 2 needs OWID data) --
source("21_empirical_application_1.R")
source("22_empirical_application_2.R")

# ---- Part 4. Figures --------------------------------------------------------
source("23_figure_1_montecarlo.R")
source("24_figure_2_montecarlo.R")
source("25_figure_3_empirical_app1.R")
source("26_figure_4_empirical_app2.R")

cat("\nAll scripts completed. Figures are in ../figures/.\n")
