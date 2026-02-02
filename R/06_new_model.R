#===============================================================================
# Section 4: Develop, Recalibrate, and Evaluate New Model
# Uses log10 scale throughout (matching Kasahara) and PaBa recalibration (per protocol)
#===============================================================================
print("--- Starting Section 4: New Model Development & Evaluation ---")

# Use the final 'analysis_data' cohort
md <- analysis_data

#===============================================================================
# Section 4A: Nested Model Comparison - Do interactions help?
# Per protocol: "The model will be extended with an interaction term between
# creatinine clearance, age, and BNP to determine if the interaction term
# improves the precision of the NT-proBNP estimates."
#
# Compares: (1) Recalibrated Kasahara (2-param adjustment)
#           (2) Kasahara-like (re-estimated, additive)
#           (3) Protocol model (+ CrCl*Age*BNP interactions)
#===============================================================================
print("--- Section 4A: Nested Model Comparison ---")

# --- Fit Kasahara-like model (same predictors as Kasahara, re-estimated) ---
# Kasahara predictors: log10(BNP), Age, BMI, Hb, CrCl (RCS), Sex, AF (ADDITIVE)
print("Fitting Kasahara-like model (same predictors, re-estimated on this study's data)...")
fit_base <- fit_kasahara_like_model(md, knots = CONFIG$rcs_knots)
if (is.null(fit_base)) {
  stop("Kasahara-like model fitting failed. Cannot proceed with comparison.")
}
print("Kasahara-like model (base, additive) summary:")
print(fit_base)

# --- Fit Protocol model (Kasahara + CrCl*Age*BNP interactions per protocol) ---
print("Fitting protocol model (Kasahara predictors + CrCl*Age*BNP interactions)...")
fit_protocol <- fit_protocol_model(md, knots = CONFIG$rcs_knots)
if (is.null(fit_protocol)) {
  warning("Protocol model fitting failed. Skipping nested comparison.")
} else {
  print("Protocol model (with interactions) summary:")
  print(fit_protocol)
}

# --- Likelihood Ratio Test: Additive vs Interaction model ---
if (!is.null(fit_base) && !is.null(fit_protocol)) {
  print("--- Likelihood Ratio Test: Additive vs Interaction Model ---")

  # Get log-likelihoods and df
  ll_base <- logLik(fit_base)
  ll_protocol <- logLik(fit_protocol)
  lr_stat <- 2 * (as.numeric(ll_protocol) - as.numeric(ll_base))
  df_diff <- attr(ll_protocol, "df") - attr(ll_base, "df")
  p_value <- pchisq(lr_stat, df_diff, lower.tail = FALSE)

  # Calculate R² change
  r2_base <- fit_base$stats["R2"]
  r2_protocol <- fit_protocol$stats["R2"]
  r2_change <- r2_protocol - r2_base

  # Calculate RMSE change (on log10 scale)
  sigma_base <- fit_base$stats["Sigma"]
  sigma_protocol <- fit_protocol$stats["Sigma"]

  # Summary table
  nested_model_comparison <- tibble(
    comparison = "Protocol (interactions) vs Additive (Kasahara-like)",
    base_predictors = "log10(BNP), Age, CrCl, Sex, AF, BMI, Hb (additive)",
    added_terms = "CrCl * Age * log10(BNP) interactions",
    n = fit_base$stats["n"],
    df_base = attr(ll_base, "df"),
    df_protocol = attr(ll_protocol, "df"),
    df_difference = df_diff,
    lr_chi_sq = lr_stat,
    p_value = p_value,
    r2_base = r2_base,
    r2_protocol = r2_protocol,
    r2_change = r2_change,
    rmse_base_log10 = sigma_base,
    rmse_protocol_log10 = sigma_protocol,
    rmse_change_log10 = sigma_protocol - sigma_base,
    aic_base = AIC(fit_base),
    aic_protocol = AIC(fit_protocol),
    bic_base = BIC(fit_base),
    bic_protocol = BIC(fit_protocol)
  )

  print("Nested Model Comparison Results:")
  print(nested_model_comparison %>% select(
    comparison, added_terms, df_difference, lr_chi_sq, p_value,
    r2_base, r2_protocol, r2_change
  ))

  # Interpretation
  cat("\n--- Interpretation ---\n")
  cat(sprintf("Additive model R²: %.3f\n", r2_base))
  cat(sprintf("Protocol model (with interactions) R²: %.3f\n", r2_protocol))
  cat(sprintf("R² improvement from adding interactions: %.4f (%.2f%%)\n",
              r2_change, r2_change * 100))
  cat(sprintf("LRT chi-square: %.2f (df=%d), p-value: %.4f\n",
              lr_stat, df_diff, p_value))
  if (p_value < 0.05) {
    cat("Interactions are statistically significant (p < 0.05),\n")
    cat("but the improvement in R² is minimal.\n")
  } else {
    cat("Interactions do NOT significantly improve prediction (p >= 0.05).\n")
  }

  # Save results
  write_rds(
    nested_model_comparison,
    output_path('nested_model_comparison.rds')
  )
  if (isTRUE(CONFIG$write_csv_outputs)) {
    write_csv(
      nested_model_comparison,
      output_path('nested_model_comparison.csv')
    )
  }

  # --- Compare all three approaches: Recalibrated Kasahara vs Re-estimated vs Protocol ---
  print("--- Comparison: Recalibrated Kasahara vs Re-estimated vs Protocol ---")

  # Use analysis_data_raw which has recalibrated Kasahara predictions from Section 2
  md_comparison <- analysis_data_raw %>%
    filter(
      is.finite(log10_ntprobnp) &
        is.finite(pred_ntprobnp_kasahara) &
        is.finite(pred_ntprobnp_kasahara_recal_ols)
    )

  # Recalibrated Kasahara (from Section 2)
  md_comparison$pred_recal_kasahara <- md_comparison$pred_ntprobnp_kasahara_recal_ols

  # Re-estimated Kasahara-like model predictions (additive)
  md_comparison$log10_pred_base <- predict(fit_base, newdata = md_comparison)
  md_comparison$pred_base <- 10^md_comparison$log10_pred_base

  # Protocol model predictions (with interactions)
  md_comparison$log10_pred_protocol <- predict(fit_protocol, newdata = md_comparison)
  md_comparison$pred_protocol <- 10^md_comparison$log10_pred_protocol

  # Calculate metrics for each
  metrics_recal_kasahara <- compute_metrics(
    md_comparison, "NTproBNP", "pred_recal_kasahara"
  ) %>% mutate(model = "Recalibrated Kasahara (2 params)")

  metrics_base <- compute_metrics(
    md_comparison, "NTproBNP", "pred_base"
  ) %>% mutate(model = "Re-estimated (additive)")

  metrics_protocol <- compute_metrics(
    md_comparison, "NTproBNP", "pred_protocol"
  ) %>% mutate(model = "Protocol (+ interactions)")

  # Combine and compare key metrics
  all_metrics <- bind_rows(metrics_recal_kasahara, metrics_base, metrics_protocol)

  comparison_summary <- all_metrics %>%
    filter(.metric %in% c("rmse_log10", "rsq_log10", "paba_slope_log10")) %>%
    select(model, .metric, .estimate) %>%
    pivot_wider(names_from = .metric, values_from = .estimate) %>%
    arrange(desc(rsq_log10))

  print("Performance Comparison (all approaches):")
  print(comparison_summary)

  write_rds(
    comparison_summary,
    output_path('model_comparison_summary.rds')
  )
  if (isTRUE(CONFIG$write_csv_outputs)) {
    write_csv(
      comparison_summary,
      output_path('model_comparison_summary.csv')
    )
  }
}

print("--- Finished Section 4A ---")

# --- Fit Initial OLS Model (log10 scale) ---
print("Fitting initial OLS model on log10(NTproBNP)...")
dd <- datadist(md)
options(datadist = 'dd')
fit_ols <- fit_new_model_ols(md, knots = CONFIG$rcs_knots)
if (is.null(fit_ols)) {
  stop("Initial OLS model fitting failed. Cannot proceed.")
}

# --- Compare Null vs Full Model (log10 scale) ---
print("Comparing null model vs full model (log10 scale)...")
required_cols_cmp <- c(
  "log10_ntprobnp",
  "Age",
  "CrCl",
  "log10_bnp",
  "Sex_M",
  "BMI",
  "AF",
  "Hb_gdl"
)
data_model_cmp <- md[, required_cols_cmp, drop = FALSE]
dd_cmp <- tryCatch(
  datadist(data_model_cmp),
  error = function(e) {
    warning("datadist failed for null model comparison: ", e$message)
    return(NULL)
  }
)
if (!is.null(dd_cmp)) {
  options(datadist = 'dd_cmp')
  fit_ols_null <- tryCatch(
    stats::lm(log10_ntprobnp ~ 1, data = data_model_cmp, na.action = na.omit),
    error = function(e) {
      warning("Null model fitting failed: ", e$message)
      return(NULL)
    }
  )

  null_vs_full <- NULL
  if (!is.null(fit_ols_null)) {
    null_vs_full <- tryCatch(
      {
        ll_null <- logLik(fit_ols_null)
        ll_full <- logLik(fit_ols)
        lr_stat <- 2 * (as.numeric(ll_full) - as.numeric(ll_null))
        df_lr <- attr(ll_full, "df") - attr(ll_null, "df")
        tibble(
          test = "LR",
          df = df_lr,
          lr_chisq = lr_stat,
          p_value = pchisq(lr_stat, df_lr, lower.tail = FALSE),
          AIC_null = AIC(fit_ols_null),
          BIC_null = BIC(fit_ols_null),
          AIC_full = AIC(fit_ols),
          BIC_full = BIC(fit_ols),
          R2_full = fit_ols$stats["R2"]
        )
      },
      error = function(e) NULL
    )

    if (is.null(null_vs_full)) {
      null_vs_full <- tibble(
        model = c("null", "full"),
        AIC = c(AIC(fit_ols_null), AIC(fit_ols)),
        BIC = c(BIC(fit_ols_null), BIC(fit_ols)),
        R2 = c(NA_real_, fit_ols$stats["R2"])
      )
    }

    if (isTRUE(CONFIG$write_csv_outputs)) {
      write_csv(null_vs_full, output_path("null_vs_full_model_comparison.csv"))
      print("Null vs full model comparison saved.")
    }
  }
}
options(datadist = 'dd')

# --- Bias-corrected calibration curve (rms::calibrate; log10 scale) ---
print("Generating bias-corrected calibration curve (rms::calibrate)...")
calibrate_ols_boot <- tryCatch(
  calibrate(fit_ols, method = "boot", B = CONFIG$bootstrap_reps),
  error = function(e) {
    warning("Bias-corrected calibration failed: ", e$message)
    return(NULL)
  }
)

# Residual SE on log10 scale
sigma_ols <- fit_ols$stats['Sigma']
if (!is.finite(sigma_ols)) {
  sigma_ols <- 0
}

# Add initial predictions (log10 scale) to the dataset
md$log10_pred_initial <- predict(fit_ols)

# --- Recalibration using Passing-Bablok (per protocol) ---
# Recalibration equation (log10 scale):
# log10(NTproBNP) = I + S * log10_pred_initial
# pred_median_recal = 10^(I + S * log10_pred_initial)
print("Running Log10-Log10 Passing-Bablok recalibration (per protocol)...")
paba_log10_model <- run_mcreg(
  md,
  x_col = log10_pred_initial,
  y_col = log10_ntprobnp
)
I_log10_full <- 0
S_log10_full <- 1
if (!is.null(paba_log10_model)) {
  paba_log10_results_para <- paba_log10_model@para
  I_log10_full <- paba_log10_results_para["Intercept", "EST"]
  S_log10_full <- paba_log10_results_para["Slope", "EST"]
  print("Log10-Log10 PaBa Results (Full Sample):")
  print(paba_log10_results_para)
} else {
  print("Log10-Log10 Passing-Bablok failed or insufficient data (full sample).")
}

# Bootstrap PaBa to estimate CI and recalibration parameters (per protocol: 1000 iterations)
print(glue(
  "--- Bootstrapping Log10-Log10 PaBa ({CONFIG$bootstrap_reps} iterations) ---"
))
paba_log10_boot <- bootstrap_paba(
  md,
  x_col = log10_pred_initial,
  y_col = log10_ntprobnp,
  reps = CONFIG$bootstrap_reps
) %>%
  filter(is.finite(intercept) & is.finite(slope))

if (nrow(paba_log10_boot) > 10) {
  paba_log10_boot_summary <- tibble(
    term = c("Intercept", "Slope"),
    estimate_median = c(
      median(paba_log10_boot$intercept, na.rm = TRUE),
      median(paba_log10_boot$slope, na.rm = TRUE)
    ),
    estimate_mean = c(
      mean(paba_log10_boot$intercept, na.rm = TRUE),
      mean(paba_log10_boot$slope, na.rm = TRUE)
    ),
    ci_lower_95 = c(
      quantile(paba_log10_boot$intercept, probs = 0.025, na.rm = TRUE),
      quantile(paba_log10_boot$slope, probs = 0.025, na.rm = TRUE)
    ),
    ci_upper_95 = c(
      quantile(paba_log10_boot$intercept, probs = 0.975, na.rm = TRUE),
      quantile(paba_log10_boot$slope, probs = 0.975, na.rm = TRUE)
    ),
    n_boot = nrow(paba_log10_boot)
  )
  print("Bootstrapped Log10-Log10 PaBa Calibration (Median with 95% CI):")
  print(paba_log10_boot_summary)

  # Use bootstrapped median for recalibration (more robust)
  I_log10_paba <- paba_log10_boot_summary$estimate_median[
    paba_log10_boot_summary$term == "Intercept"
  ]
  S_log10_paba <- paba_log10_boot_summary$estimate_median[
    paba_log10_boot_summary$term == "Slope"
  ]
  print(glue(
    "Using bootstrapped median for PaBa recalibration: I={round(I_log10_paba, 4)}, S={round(S_log10_paba, 4)}"
  ))
} else {
  print(
    "Insufficient bootstrap PaBa results; using full-sample PaBa estimates."
  )
  I_log10_paba <- I_log10_full
  S_log10_paba <- S_log10_full
}

if (
  !is.finite(I_log10_paba) || !is.finite(S_log10_paba) || abs(S_log10_paba) < 1e-6
) {
  warning(
    "Invalid/near-zero PaBa slope. Resetting to S=1, I=0 (no recalibration)."
  )
  I_log10_paba <- 0
  S_log10_paba <- 1
}
log10_recal_params <- c(I_log10_paba, S_log10_paba)

# --- Calculate Recalibrated Predictions (Median) ---
md <- md %>%
  mutate(
    log10_pred_recal = I_log10_paba + S_log10_paba * log10_pred_initial,
    pred_ntprobnp_new_median = 10^log10_pred_recal
  ) %>%
  filter(is.finite(pred_ntprobnp_new_median))

# --- Evaluate APPARENT Performance of RECALIBRATED Model ---
print("--- Evaluating APPARENT Performance (Recalibrated New Model) ---")
if (nrow(md) > 10) {
  apparent_performance_new <- compute_metrics(
    md,
    target_col = "NTproBNP",
    predicted_col = "pred_ntprobnp_new_median"
  )
  print("Apparent Performance Metrics (Recalibrated New Model - PaBa):")
  print(apparent_performance_new)
} else {
  print("Insufficient data after processing for new model apparent evaluation.")
  apparent_performance_new <- NULL
}

# --- Optimism Correction Bootstrapping ---
print(glue(
  "--- Starting Optimism Correction Bootstrapping ({CONFIG$bootstrap_reps} iterations) ---"
))

# Function to get performance metrics for one bootstrap iteration
calculate_optimism_step <- function(split, knots) {
  bs_analysis <- analysis(split)
  bs_assessment <- assessment(split)

  # Fit model on bootstrap sample
  bs_fit <- fit_new_model_ols(bs_analysis, knots = knots)
  if (is.null(bs_fit)) {
    return(NULL)
  }

  bs_analysis <- bs_analysis %>%
    mutate(log10_pred_initial_bs = predict(bs_fit, newdata = .))

  # Fit PaBa recalibration within bootstrap sample
  bs_paba <- run_mcreg(
    bs_analysis,
    x_col = log10_pred_initial_bs,
    y_col = log10_ntprobnp
  )
  if (is.null(bs_paba)) {
    I_bs <- 0
    S_bs <- 1
  } else {
    I_bs <- bs_paba@para["Intercept", "EST"]
    S_bs <- bs_paba@para["Slope", "EST"]
  }

  # 1. Performance on Bootstrap Sample (Apparent)

  bs_analysis_pred <- bs_analysis %>%
    mutate(
      log10_pred_recal_bs = I_bs + S_bs * log10_pred_initial_bs,
      pred_median_recal_bs = 10^log10_pred_recal_bs
    ) %>%
    filter(is.finite(pred_median_recal_bs))
  if (nrow(bs_analysis_pred) < 10) {
    return(NULL)
  }
  perf_on_bootstrap <- compute_metrics(
    bs_analysis_pred,
    "NTproBNP",
    "pred_median_recal_bs"
  ) %>%
    select(.metric, .estimator, .estimate_on_bootstrap = .estimate)

  # 2. Performance on Original Sample (Test)
  orig_assessment_pred <- bs_assessment %>%
    mutate(
      log10_pred_initial_orig = predict(bs_fit, newdata = .),
      log10_pred_recal_orig = I_bs + S_bs * log10_pred_initial_orig,
      pred_median_recal_on_orig = 10^log10_pred_recal_orig
    ) %>%
    filter(is.finite(pred_median_recal_on_orig))
  if (nrow(orig_assessment_pred) < 10) {
    return(NULL)
  }
  perf_on_original <- compute_metrics(
    orig_assessment_pred,
    "NTproBNP",
    "pred_median_recal_on_orig"
  ) %>%
    select(.metric, .estimator, .estimate_on_original = .estimate)

  # Calculate Optimism
  metrics_step <- inner_join(
    perf_on_bootstrap,
    perf_on_original,
    by = c(".metric", ".estimator")
  ) %>%
    mutate(optimism = .estimate_on_bootstrap - .estimate_on_original) %>%
    filter(is.finite(optimism))
  return(metrics_step)
}

# Perform bootstrapping
boots <- bootstraps(md, times = CONFIG$bootstrap_reps, apparent = FALSE)
optimism_cores <- get_bootstrap_cores()
message(glue("Running optimism correction on {optimism_cores} core(s)..."))
optimism_estimates_list <- parallel::mclapply(
  boots$splits,
  function(split) {
    calculate_optimism_step(split, CONFIG$rcs_knots)
  },
  mc.cores = optimism_cores,
  mc.set.seed = TRUE
)
optimism_estimates <- bind_rows(keep(optimism_estimates_list, ~ !is.null(.)))
print("Finished Optimism Correction Bootstrapping.")

# Summarize Optimism and Calculate Corrected Estimates
if (nrow(optimism_estimates) > 0 && !is.null(apparent_performance_new)) {
  optimism_summary <- optimism_estimates %>%
    group_by(.metric, .estimator) %>%
    summarise(
      mean_optimism = mean(optimism, na.rm = TRUE),
      n_boots = n(),
      lower_ci_95 = quantile(.estimate_on_original, probs = 0.025, na.rm = TRUE),
      upper_ci_95 = quantile(.estimate_on_original, probs = 0.975, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    filter(is.finite(mean_optimism))

  optimism_corrected_final <- apparent_performance_new %>%
    select(
      .metric,
      .estimator,
      threshold_label,
      apparent_estimate = .estimate
    ) %>%
    inner_join(optimism_summary, by = c(".metric", ".estimator")) %>%
    mutate(optimism_corrected = apparent_estimate - mean_optimism) %>%
    select(
      .metric,
      .estimator,
      apparent_estimate,
      mean_optimism,
      optimism_corrected,
      lower_ci_95,
      upper_ci_95,
      n_boots
    ) %>%
    arrange(.metric, .estimator)

  print("Optimism-Corrected Performance Metrics (Recalibrated New Model):")
  print(optimism_corrected_final)
} else {
  print(
    "Could not calculate optimism-corrected estimates (insufficient bootstrap results or apparent performance missing)."
  )
  optimism_corrected_final <- NULL
}

print("--- Finished Section 4 ---")
