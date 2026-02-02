#===============================================================================
# Section 1.5: Baseline Characteristics Table (using gtsummary)
#===============================================================================
print("--- Starting Section 1.5: Generate Baseline Table ---")

# Variables to include in Table 1 (same list as before)
table1_vars <- c(
  # Demographics / Vitals
  "Age",
  "Sex",
  "Wt",
  "Ht",
  "BMI",
  # Labs
  "Cr",
  "CrCl",
  "Hb",
  "Hb_gdl",
  "BNP",
  "NTproBNP",
  # Comorbidities
  "AF",
  "CAD",
  "ACS",
  "CCSC_lll",
  "CCSC_lV",
  "high_risk_CAD",
  "Cardiac_arrest",
  "revascularization_6",
  "revascularization",
  "PVD",
  "stroke",
  "tia",
  "copd",
  "cancer",
  "chf",
  "hf_echo",
  "PE_DVT",
  "DM",
  "HTN",
  "OSA",
  # Surgery / Procedure Related
  "Major_GenS",
  "Major_NeurS",
  "Major_VascS",
  "Major_OrthS",
  "Major_UroGynS",
  "Major_ThorS",
  "Low_Surgery",
  "No_surgery",
  "Type of Anesthesia"
)

# Select available variables FROM THE RAW CLEANED DATA
available_table1_vars <- intersect(table1_vars, names(analysis_data_raw))
print(glue(
  "Columns available for Table 1 (from initial data): {paste(available_table1_vars, collapse=', ')}"
))
print(glue(
  "Total patients in initial dataset for Table 1: {nrow(analysis_data_raw)}"
))


# Prepare data for gtsummary using the initial cleaned data
data_for_table1 <- analysis_data_raw %>%
  select(StudyID, all_of(available_table1_vars)) %>%
  mutate(
    stroke_tia = case_when(
      "stroke" %in% names(.) && "tia" %in% names(.) ~ if_else(
        is.na(stroke) & is.na(tia),
        NA_real_,
        if_else(coalesce(stroke, 0) == 1 | coalesce(tia, 0) == 1, 1, 0)
      ),
      "stroke" %in% names(.) ~ as.numeric(stroke),
      "tia" %in% names(.) ~ as.numeric(tia),
      TRUE ~ NA_real_
    ),
    # Convert relevant comorbidities to factors
    across(
      any_of(c(
        "AF",
        "CAD",
        "ACS",
        "CCSC_lll",
        "CCSC_lV",
        "high_risk_CAD",
        "Cardiac_arrest",
        "revascularization_6",
        "revascularization",
        "PVD",
        "stroke_tia",
        "copd",
        "cancer",
        "chf",
        "hf_echo",
        "PE_DVT",
        "DM",
        "HTN",
        "OSA"
      )),
      ~ factor(., levels = c(0, 1), labels = c("No", "Yes"))
    ),
    Sex = factor(Sex, levels = c("F", "M"), labels = c("Female", "Male")),
    SurgeryCategory = case_when(
      if ("Major_VascS" %in% names(.)) Major_VascS == 1 ~ "Major Vascular",
      if ("Major_ThorS" %in% names(.)) Major_ThorS == 1 ~ "Major Thoracic",
      if ("Major_GenS" %in% names(.)) Major_GenS == 1 ~ "Major General",
      if ("Major_NeurS" %in% names(.)) Major_NeurS == 1 ~ "Major Neurological",
      if ("Major_OrthS" %in% names(.)) Major_OrthS == 1 ~ "Major Orthopedic",
      if ("Major_UroGynS" %in% names(.)) Major_UroGynS == 1 ~ "Major Uro/Gyn",
      if ("Low_Surgery" %in% names(.)) Low_Surgery == 1 ~ "Low Risk Surgery",
      if ("No_surgery" %in% names(.)) No_surgery == 1 ~ "Surgery cancelled",
      TRUE ~ "Other/Unknown"
    ) %>%
      factor(),
    across(any_of("Type of Anesthesia"), as.factor)
  ) %>%
  # Select final columns for the table
  select(
    Age,
    Sex,
    BMI,
    CrCl,
    Hb_gdl,
    BNP,
    NTproBNP,
    any_of(c(
      "AF",
      "CAD",
      "PVD",
      "stroke_tia",
      "copd",
      "cancer",
      "chf",
      "DM",
      "HTN",
      "OSA"
    )),
    SurgeryCategory,
    any_of("Type of Anesthesia")
  )

# Define variable labels for the table AS A NAMED LIST
var_labels_list <- list(
  Age = "Age (years)",
  BMI = "Body Mass Index (kg/m²)",
  CrCl = "Creatinine Clearance (mL/min)",
  Hb_gdl = "Hemoglobin (g/dL)",
  BNP = "BNP (pg/mL)",
  NTproBNP = "NT-proBNP (pg/mL)",
  AF = "Atrial Fibrillation",
  CAD = "Coronary Artery Disease",
  PVD = "Peripheral Vascular Disease",
  stroke_tia = "Stroke/TIA History",
  copd = "COPD",
  cancer = "Cancer History",
  chf = "Congestive Heart Failure History",
  DM = "Diabetes Mellitus",
  HTN = "Hypertension",
  OSA = "Obstructive Sleep Apnea",
  SurgeryCategory = "Surgery Category"
)
if ("Type of Anesthesia" %in% names(data_for_table1)) {
  var_labels_list <- c(
    var_labels_list,
    list("Type of Anesthesia" = "Type of Anesthesia")
  )
}

# Filter labels to only those present in the final table data
final_var_labels <- var_labels_list[
  names(var_labels_list) %in% names(data_for_table1)
]

# --- Calculate Missing Counts for Footnote ---
# Map variable names to reader-friendly labels for footnote
footnote_var_labels <- c(
  CrCl = "Creatinine Clearance",
  Hb_gdl = "Hemoglobin",
  DM = "Diabetes",
  Age = "Age",
  BMI = "Body Mass Index",
  BNP = "BNP",
  NTproBNP = "NT-proBNP",
  AF = "Atrial Fibrillation",
  CAD = "Coronary Artery Disease",
  PVD = "Peripheral Vascular Disease",
  stroke_tia = "Stroke/TIA",
  copd = "COPD",
  cancer = "Cancer",
  chf = "Heart Failure",
  HTN = "Hypertension",
  OSA = "Sleep Apnea",
  SurgeryCategory = "Surgery Category"
)

missing_counts <- sapply(data_for_table1, function(x) sum(is.na(x)))
missing_counts_text <- names(missing_counts[missing_counts > 0]) %>%
  {
    # Use friendly names if available, otherwise use original name
    friendly_names <- ifelse(. %in% names(footnote_var_labels),
                              footnote_var_labels[.], .)
    paste0(friendly_names, " (", missing_counts[missing_counts > 0], ")")
  } %>%
  paste(collapse = "; ")

missing_footnote_text <- if (nchar(missing_counts_text) > 0) {
  paste("Number missing: ", missing_counts_text)
} else {
  "No missing data for variables shown."
}

# Determine stats footnote based on presence of biomarkers
biomarkers_present <- length(intersect(
  c("BNP", "NTproBNP"),
  names(data_for_table1)
)) >
  0
stats_footnote_text <- if (biomarkers_present) {
  "Median (IQR) for continuous variables. Statistics calculated on non-missing data."
} else {
  "Median (IQR) for continuous variables. Statistics calculated on non-missing data."
}

# --- Generate Table 1 using gtsummary ---
tryCatch(
  {
    theme_gtsummary_compact() # Use a compact theme

    # 1. Create the base table object, setting missing="no"
    baseline_table_base <- data_for_table1 %>%
      tbl_summary(
        label = final_var_labels,
        statistic = list(
          all_continuous() ~ "{median} ({p25}, {p75})",
          all_categorical() ~ "{n} ({p}%)"
        ),
        digits = list(
          all_continuous() ~ 1,
          all_of(intersect(c("BNP", "NTproBNP"), names(data_for_table1))) ~ 0,
          all_categorical() ~ c(0, 1)
        ),
        type = list(all_continuous() ~ "continuous"),
        missing = "no" # <<<<<<<<<<<<<<<<<<<<<<<<<<< CHANGE: Hide missing rows
      )

    # 2. Add footnotes (Stats description AND Missing counts)
    #    Addressing the deprecation warning by using direct arguments
    baseline_table_footnote <- baseline_table_base %>%
      modify_footnote(
        all_stat_cols() ~ stats_footnote_text, # Footnote for statistics shown
        abbreviation = TRUE
      ) %>%
      modify_footnote(
        all_stat_cols() ~ missing_footnote_text # Separate footnote for missing counts
      )

    # 3. Apply remaining modifications and print/save
    baseline_table_final <- baseline_table_footnote %>%
      modify_header(
        label ~ "**Characteristic**",
        all_stat_cols() ~ "**Overall (N = {N})**"
      ) %>%
      modify_caption(
        "**Table 1: Baseline Characteristics of Study Cohort**"
      ) %>%
      bold_labels()

    print(baseline_table_final)
  },
  error = function(e) {
    print(paste("Error generating Baseline Table 1:", e$message))
    print(traceback())
  }
)

print("--- Finished Section 1.5 ---")
