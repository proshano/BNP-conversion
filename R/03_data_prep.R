#===============================================================================
# Section 1: Data Loading and Preparation
#===============================================================================
print("--- Starting Section 1: Data Loading and Preparation ---")

# --- Load Data ---
d_new_raw <- read_input_csv(
  CONFIG$input_data_file_new,
  label = "input_data_file_new"
)
if (isTRUE(CONFIG$use_old_dataset)) {
  d_old_raw <- read_input_csv(
    CONFIG$input_data_file_old,
    label = "input_data_file_old"
  )
} else {
  d_old_raw <- NULL
}

required_cols_new <- c(
  "Age",
  "Sex",
  "Wt",
  "Ht",
  "Cr",
  "Hb",
  "BNP",
  "NTproBNP",
  "AF"
)
check_required_columns(d_new_raw, required_cols_new, "input_data_file_new")
missing_outcomes <- setdiff(c("MINS30", "vascular_death"), names(d_new_raw))
if (length(missing_outcomes) > 0) {
  message(
    "Outcome columns not found in input_data_file_new: ",
    paste(missing_outcomes, collapse = ", "),
    ". Outcome analyses will be skipped if composite outcome is missing."
  )
}

# --- Clean Data Function ---
clean_data <- function(df, cr_conversion_factor) {
  if (!"StudyID" %in% names(df)) {
    df$StudyID <- seq_len(nrow(df))
  }
  if (!"MINS30" %in% names(df)) {
    df$MINS30 <- NA_real_
  }
  if (!"vascular_death" %in% names(df)) {
    df$vascular_death <- NA_real_
  }
  df %>%
    mutate(
      Sex_M = as.numeric(Sex == 'M'),
      Ht_m = Ht / 100,
      BMI = Wt / (Ht_m^2),
      Hb_gdl = as.numeric(Hb) / 10,
      Cr_umol_input = as.numeric(Cr),
      Cr_mgdl = ifelse(
        is.finite(Cr_umol_input),
        Cr_umol_input / cr_conversion_factor,
        NA_real_
      ),
      CrCl = case_when(
        !is.finite(Cr_mgdl) |
          Cr_mgdl <= 0 |
          !is.finite(Sex_M) |
          !is.finite(Age) |
          Age <= 0 |
          !is.finite(Wt) |
          Wt <= 0 ~ NA_real_,
        TRUE ~ nephro::CG(creatinine = Cr_mgdl, sex = Sex_M, age = Age, wt = Wt)
      ),
      CrCl = ifelse(is.finite(CrCl) & CrCl > 0, CrCl, NA_real_),
      across(
        c(BNP, NTproBNP),
        ~ ifelse(is.finite(.x) & .x > 0, as.numeric(.x), NA_real_)
      ),
      across(
        any_of(c(
          "AF",
          "CAD",
          "MI",
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
          "Insulin",
          "HTN",
          "OSA",
          "Major_GenS",
          "Major_NeurS",
          "Major_VascS",
          "Major_OrthS",
          "Major_UroGynS",
          "Major_ThorS",
          "Low_Surgery",
          "No_surgery"
        )),
        ~ as.numeric(.)
      ),
      composite_outcome = case_when(
        is.na(MINS30) | is.na(vascular_death) ~ NA_real_,
        MINS30 == 1 | vascular_death == 1 ~ 1,
        TRUE ~ 0
      ),
      # Log10 transformations (matching Kasahara scale)
      log10_bnp = ifelse(is.finite(BNP), log10(BNP), NA_real_),
      log10_ntprobnp = ifelse(is.finite(NTproBNP), log10(NTproBNP), NA_real_)
    ) %>%
    select(-any_of(c("Ht_m", "Cr_umol_input", "Cr_mgdl")))
}

# --- Apply the cleaning function ---
analysis_data_raw <- clean_data(d_new_raw, CONFIG$cr_conversion_factor)

# --- Calculate Old Model Predictions ---
analysis_data_raw <- analysis_data_raw %>%
  mutate(
    pred_ntprobnp_kasahara = predict_ntprobnp_kasahara(
      BNP = BNP,
      Age = Age,
      BMI = BMI,
      Hb_gdl = Hb_gdl,
      CrCl = CrCl,
      Sex_M = Sex_M,
      AF = AF
    )
  )

# --- Define Final Analysis Cohort ---
required_vars_new_model <- c(
  "log10_ntprobnp",
  "Age",
  "CrCl",
  "log10_bnp",
  "Sex_M",
  "BMI",
  "AF",
  "Hb_gdl"
)
required_vars_biomarkers <- c("NTproBNP", "BNP")
required_numeric <- unique(c(required_vars_new_model, required_vars_biomarkers))

analysis_data <- analysis_data_raw %>%
  filter(if_all(all_of(required_numeric), is.finite)) %>%
  {
    if ("StudyID" %in% names(.)) filter(., !is.na(StudyID)) else .
  }
n_analysis <- nrow(analysis_data)
print(glue(
  "Number of subjects in the model-development cohort (complete predictors/biomarkers): {n_analysis}"
))
if (n_analysis < 50) {
  stop(
    "Insufficient data in the final analysis cohort (<50 subjects). Stopping analysis."
  )
}

# Define datadist for rms functions based on the final analysis cohort
dd_final <- datadist(analysis_data)
options(datadist = 'dd_final')

print("--- Finished Section 1 ---")
