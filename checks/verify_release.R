#===============================================================================
# Release Verification for the BNP Conversion Analysis
#===============================================================================

if (!file.exists("BNP Data Analysis Jan 9 2025.csv")) {
  stop(
    "Private input file BNP Data Analysis Jan 9 2025.csv is missing.",
    call. = FALSE
  )
}

data_current <- read.csv(
  "BNP Data Analysis Jan 9 2025.csv",
  fileEncoding = "UTF-8-BOM",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (nrow(data_current) != 460L) {
  stop(
    "Expected 460 source records; found ",
    nrow(data_current),
    ".",
    call. = FALSE
  )
}

model_columns <- c(
  "Age", "Sex", "Wt", "Ht", "Cr", "Hb", "BNP", "NTproBNP", "AF"
)
model_complete <- complete.cases(data_current[, model_columns])
outcome_complete <- model_complete &
  complete.cases(data_current[, c("MINS30", "vascular_death")])
composite_outcome <- with(
  data_current,
  ifelse(
    is.na(MINS30) | is.na(vascular_death),
    NA,
    as.integer(MINS30 == 1 | vascular_death == 1)
  )
)

if (sum(model_complete) != 448L) {
  stop("Expected 448 complete model records; found ", sum(model_complete), ".", call. = FALSE)
}
if (sum(outcome_complete) != 423L) {
  stop(
    "Expected 423 complete outcome records; found ",
    sum(outcome_complete),
    ".",
    call. = FALSE
  )
}
if (sum(composite_outcome[outcome_complete]) != 20L) {
  stop(
    "Expected 20 events in the outcome cohort; found ",
    sum(composite_outcome[outcome_complete]),
    ".",
    call. = FALSE
  )
}

expected_summaries <- list(
  "assets/ntprobnp_cat_summary.rds" = data.frame(
    n = c(129L, 105L, 164L, 25L),
    events = c(1L, 3L, 10L, 6L)
  ),
  "assets/pred_cat_summary.rds" = data.frame(
    n = c(131L, 109L, 165L, 18L),
    events = c(1L, 3L, 12L, 4L)
  )
)
for (path in names(expected_summaries)) {
  if (!file.exists(path)) {
    stop("Missing release output: ", path, call. = FALSE)
  }
  observed <- readRDS(path)
  expected <- expected_summaries[[path]]
  if (!identical(as.integer(observed$n), expected$n) ||
      !identical(as.integer(observed$events), expected$events)) {
    stop("Unexpected category counts in ", path, ".", call. = FALSE)
  }
}

bootstrap_outputs <- c(
  "assets/kasahara_validation_metrics_with_ci.rds",
  "assets/kasahara_recal_optimism_corrected.rds",
  "assets/base_model_optimism_corrected.rds",
  "assets/extended_model_optimism_corrected.rds"
)
for (path in bootstrap_outputs) {
  if (!file.exists(path)) {
    stop("Missing bootstrap release output: ", path, call. = FALSE)
  }
  bootstrap_result <- readRDS(path)
  if (!("n_boots" %in% names(bootstrap_result)) ||
      any(bootstrap_result$n_boots != 1000L)) {
    stop("Expected 1,000 completed bootstrap samples in ", path, ".", call. = FALSE)
  }
}

required_outputs <- c(
  "assets/analysis_report.Rmd",
  "assets/outcome_assoc_combined.png"
)
missing_outputs <- required_outputs[!file.exists(required_outputs)]
if (length(missing_outputs) > 0) {
  stop("Missing release outputs: ", paste(missing_outputs, collapse = ", "), call. = FALSE)
}

prose_files <- c("R/08_report.R", "assets/analysis_report.Rmd")
prose <- unlist(lapply(prose_files, readLines, warn = FALSE), use.names = FALSE)
if (any(grepl("overpredict|overestimat", prose, ignore.case = TRUE))) {
  stop("Outdated overprediction wording remains in a report source.", call. = FALSE)
}

outcome_code <- readLines("R/07_outcomes.R", warn = FALSE)
if (any(grepl("c\\(0, *1000, *2000, *3000, *4000\\)", outcome_code))) {
  stop("Figure 6 still contains the obsolete hard-coded BNP x-axis breaks.", call. = FALSE)
}
if (!any(grepl("ntprobnp_shared_max", outcome_code, fixed = TRUE))) {
  stop("Figure 6 does not enforce a shared NT-proBNP x-axis range.", call. = FALSE)
}
if (any(grepl("coord_cartesian(ylim = c(1, 100))", outcome_code, fixed = TRUE))) {
  stop("Figure 6 still clips confidence bands at an odds ratio of 100.", call. = FALSE)
}
if (!any(grepl("y_axis_upper", outcome_code, fixed = TRUE))) {
  stop("Figure 6 does not derive its upper odds-ratio limit from the data.", call. = FALSE)
}

message(
  "Release verification passed: 448 model records; ",
  "423 outcome records; 20 events."
)
