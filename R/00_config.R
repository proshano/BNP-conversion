#===============================================================================
# Configuration Settings
#===============================================================================

CONFIG <- list(
  seed = 20250409, # Reproducibility seed
  input_data_file_old = 'BNPdata_for_analysis_20230411_clean.csv',
  input_data_file_new = 'BNP Data Analysis Jan 9 2025.csv',
  use_old_dataset = FALSE,
  output_dir = "assets",
  cr_conversion_factor = 88.42, # umol/L per mg/dL for Creatinine
  binary_threshold = 300, # NT-proBNP threshold for binary metrics
  class_breaks = c(100, 200, 1500), # NT-proBNP thresholds for multi-class metrics
  rcs_knots = 3, # Number of knots for restricted cubic splines
  bootstrap_reps = 100, # Set to 1000 for protocol-level bootstrapping
  bootstrap_cores = 10, # Multicore bootstrapping
  kasahara_crcl_knots = c(56.5, 72.4, 93.7), # Knots from Kasahara paper
  log10_bnp_interval = log10(c(1, 50000)), # Interval for root finding (log10 scale)
  target_ntprobnp_thresholds = c(100, 200, 1500), # Thresholds for outcome analysis
  write_csv_outputs = FALSE # Disable CSV outputs; keep results in-memory/RDS
)

# --- Output Directory ---
CONFIG$output_dir <- normalizePath(CONFIG$output_dir, mustWork = FALSE)
dir.create(CONFIG$output_dir, showWarnings = FALSE, recursive = TRUE)
print(paste("Output directory set to:", CONFIG$output_dir))

# --- Set Seed ---
set.seed(CONFIG$seed)
print(paste("Random seed set to:", CONFIG$seed))
