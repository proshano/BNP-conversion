#===============================================================================
# Main Orchestrator Script for BNP Conversion Analysis
#===============================================================================
# This script sources all component scripts in order to run the full analysis.
# Run this script from the project root directory.
#===============================================================================

# Source scripts in order
print("=== Starting BNP Conversion Analysis ===")

print("--- Loading configuration ---")
source("R/00_config.R")

print("--- Loading libraries ---")
source("R/01_libraries.R")

print("--- Loading utility functions ---")
source("R/02_utils.R")

print("--- Running data preparation ---")
source("R/03_data_prep.R")

print("--- Generating baseline characteristics table ---")
source("R/04_baseline_table.R")

print("--- Evaluating Kasahara model ---")
source("R/05_kasahara_eval.R")

print("--- Developing and evaluating new model ---")
source("R/06_new_model.R")

print("--- Running outcome association analysis ---")
source("R/07_outcomes.R")

print("--- Generating report ---")
source("R/08_report.R")

print("=== Analysis Complete ===")
