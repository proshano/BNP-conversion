#===============================================================================
# Section 4: Develop, Recalibrate, and Evaluate New Model
# Uses log10 scale throughout (matching Kasahara) and OLS recalibration
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
#           (3) Extended model (+ CrCl*Age*BNP interactions)
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

# --- Fit Extended model (Kasahara + CrCl*Age*BNP interactions per protocol) ---
print("Fitting extended model (Kasahara predictors + CrCl*Age*BNP interactions)...")
fit_protocol <- fit_protocol_model(md, knots = CONFIG$rcs_knots)
if (is.null(fit_protocol)) {
  warning("Extended model fitting failed. Skipping nested comparison.")
} else {
  print("Extended model (with interactions) summary:")
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
    comparison = "Extended (interactions) vs Additive (Kasahara-like)",
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
  cat(sprintf("Extended model (with interactions) R²: %.3f\n", r2_protocol))
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

  # --- Compare all three approaches: Recalibrated Kasahara vs Re-estimated vs Extended ---
  print("--- Comparison: Recalibrated Kasahara vs Re-estimated vs Extended ---")

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

  # Extended model predictions (with interactions)
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
  ) %>% mutate(model = "Extended (+ interactions)")

  # Combine and compare key metrics
  all_metrics <- bind_rows(metrics_recal_kasahara, metrics_base, metrics_protocol)

  comparison_summary <- all_metrics %>%
    filter(.metric %in% c("rmse_log10", "rsq_log10")) %>%
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

# --- Recalibration using OLS (used downstream for Section 5 outcomes) ---
# Recalibration equation (log10 scale):
# log10(NTproBNP) = I + S * log10_pred_initial
# pred_median_recal = 10^(I + S * log10_pred_initial)
print("Running Log10-Log10 OLS recalibration for downstream predictions...")
recal_data <- md %>%
  filter(is.finite(log10_ntprobnp) & is.finite(log10_pred_initial))

if (nrow(recal_data) >= 10) {
  recal_fit <- lm(log10_ntprobnp ~ log10_pred_initial, data = recal_data)
  I_log10_ols <- unname(coef(recal_fit)[1])
  S_log10_ols <- unname(coef(recal_fit)[2])
  print(glue(
    "Log10-Log10 OLS recalibration parameters: I={round(I_log10_ols, 4)}, S={round(S_log10_ols, 4)}"
  ))
} else {
  warning("Insufficient data for OLS recalibration; using I=0, S=1.")
  I_log10_ols <- 0
  S_log10_ols <- 1
}

if (
  !is.finite(I_log10_ols) || !is.finite(S_log10_ols) || abs(S_log10_ols) < 1e-6
) {
  warning(
    "Invalid/near-zero OLS slope. Resetting to S=1, I=0 (no recalibration)."
  )
  I_log10_ols <- 0
  S_log10_ols <- 1
}
log10_recal_params <- c(I_log10_ols, S_log10_ols)

# --- Calculate Recalibrated Predictions (Median) ---
md <- md %>%
  mutate(
    log10_pred_recal = I_log10_ols + S_log10_ols * log10_pred_initial,
    pred_ntprobnp_new_median = 10^log10_pred_recal
  ) %>%
  filter(is.finite(pred_ntprobnp_new_median))

# ===============================================================================
# Optimism Correction for Re-estimated (Additive) and Extended Models
# These are needed to demonstrate the effect of overfitting in more complex models
# ===============================================================================

if (!is.null(fit_base) && !is.null(fit_protocol)) {
  print("--- Starting Optimism Correction for Re-estimated and Extended Models ---")
  optimism_cores <- get_bootstrap_cores()

  # --- Apparent Performance for Re-estimated (Additive) Model ---
  md_with_base <- md %>%
    mutate(
      log10_pred_base = predict(fit_base, newdata = .),
      pred_base = 10^log10_pred_base
    ) %>%
    filter(is.finite(pred_base))

  apparent_performance_base <- compute_metrics(
    md_with_base,
    target_col = "NTproBNP",
    predicted_col = "pred_base"
  )
  print("Apparent Performance (Re-estimated Additive Model):")
  print(apparent_performance_base)

  # --- Apparent Performance for Extended Model ---
  md_with_protocol <- md %>%
    mutate(
      log10_pred_protocol = predict(fit_protocol, newdata = .),
      pred_protocol = 10^log10_pred_protocol
    ) %>%
    filter(is.finite(pred_protocol))

  apparent_performance_protocol <- compute_metrics(
    md_with_protocol,
    target_col = "NTproBNP",
    predicted_col = "pred_protocol"
  )
  print("Apparent Performance (Extended Model):")
  print(apparent_performance_protocol)

  # --- Optimism Correction Function for Re-estimated Model ---
  calculate_optimism_step_base <- function(split, knots) {
    bs_analysis <- analysis(split)
    bs_assessment <- assessment(split)

    # Fit additive model on bootstrap sample
    bs_fit <- fit_kasahara_like_model(bs_analysis, knots = knots)
    if (is.null(bs_fit)) {
      return(NULL)
    }

    # 1. Performance on Bootstrap Sample (Apparent)
    bs_analysis_pred <- bs_analysis %>%
      mutate(
        log10_pred_bs = predict(bs_fit, newdata = .),
        pred_bs = 10^log10_pred_bs
      ) %>%
      filter(is.finite(pred_bs))
    if (nrow(bs_analysis_pred) < 10) {
      return(NULL)
    }
    perf_on_bootstrap <- compute_metrics(
      bs_analysis_pred,
      "NTproBNP",
      "pred_bs"
    ) %>%
      select(.metric, .estimator, .estimate_on_bootstrap = .estimate)

    # 2. Performance on Original Sample (Test)
    orig_assessment_pred <- bs_assessment %>%
      mutate(
        log10_pred_orig = predict(bs_fit, newdata = .),
        pred_orig = 10^log10_pred_orig
      ) %>%
      filter(is.finite(pred_orig))
    if (nrow(orig_assessment_pred) < 10) {
      return(NULL)
    }
    perf_on_original <- compute_metrics(
      orig_assessment_pred,
      "NTproBNP",
      "pred_orig"
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

  # --- Optimism Correction Function for Extended Model ---
  calculate_optimism_step_protocol <- function(split, knots) {
    bs_analysis <- analysis(split)
    bs_assessment <- assessment(split)

    # Fit extended model on bootstrap sample
    bs_fit <- fit_protocol_model(bs_analysis, knots = knots)
    if (is.null(bs_fit)) {
      return(NULL)
    }

    # 1. Performance on Bootstrap Sample (Apparent)
    bs_analysis_pred <- bs_analysis %>%
      mutate(
        log10_pred_bs = predict(bs_fit, newdata = .),
        pred_bs = 10^log10_pred_bs
      ) %>%
      filter(is.finite(pred_bs))
    if (nrow(bs_analysis_pred) < 10) {
      return(NULL)
    }
    perf_on_bootstrap <- compute_metrics(
      bs_analysis_pred,
      "NTproBNP",
      "pred_bs"
    ) %>%
      select(.metric, .estimator, .estimate_on_bootstrap = .estimate)

    # 2. Performance on Original Sample (Test)
    orig_assessment_pred <- bs_assessment %>%
      mutate(
        log10_pred_orig = predict(bs_fit, newdata = .),
        pred_orig = 10^log10_pred_orig
      ) %>%
      filter(is.finite(pred_orig))
    if (nrow(orig_assessment_pred) < 10) {
      return(NULL)
    }
    perf_on_original <- compute_metrics(
      orig_assessment_pred,
      "NTproBNP",
      "pred_orig"
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

  # --- Run Bootstrapping for Re-estimated Model ---
  print(glue("Running optimism correction for re-estimated model ({CONFIG$bootstrap_reps} iterations)..."))
  boots_base <- bootstraps(md, times = CONFIG$bootstrap_reps, apparent = FALSE)
  optimism_estimates_base_list <- parallel::mclapply(
    boots_base$splits,
    function(split) {
      calculate_optimism_step_base(split, CONFIG$rcs_knots)
    },
    mc.cores = optimism_cores,
    mc.set.seed = TRUE
  )
  optimism_estimates_base <- bind_rows(keep(optimism_estimates_base_list, ~ !is.null(.)))

  # Summarize optimism for re-estimated model
  if (nrow(optimism_estimates_base) > 0) {
    optimism_summary_base <- optimism_estimates_base %>%
      group_by(.metric, .estimator) %>%
      summarise(
        mean_optimism = mean(optimism, na.rm = TRUE),
        n_boots = n(),
        ci_lower_95 = quantile(.estimate_on_original, probs = 0.025, na.rm = TRUE),
        ci_upper_95 = quantile(.estimate_on_original, probs = 0.975, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      filter(is.finite(mean_optimism))

    base_optimism_corrected <- apparent_performance_base %>%
      select(
        .metric,
        .estimator,
        apparent_estimate = .estimate
      ) %>%
      inner_join(optimism_summary_base, by = c(".metric", ".estimator")) %>%
      mutate(optimism_corrected = apparent_estimate - mean_optimism) %>%
      arrange(.metric, .estimator)

    print("Optimism-Corrected Performance (Re-estimated Additive Model):")
    print(base_optimism_corrected)

    write_rds(
      base_optimism_corrected,
      output_path('base_model_optimism_corrected.rds')
    )
  } else {
    print("Insufficient bootstrap results for re-estimated model optimism correction.")
    base_optimism_corrected <- NULL
  }

  # --- Run Bootstrapping for Extended Model ---
  print(glue("Running optimism correction for extended model ({CONFIG$bootstrap_reps} iterations)..."))
  boots_protocol <- bootstraps(md, times = CONFIG$bootstrap_reps, apparent = FALSE)
  optimism_estimates_protocol_list <- parallel::mclapply(
    boots_protocol$splits,
    function(split) {
      calculate_optimism_step_protocol(split, CONFIG$rcs_knots)
    },
    mc.cores = optimism_cores,
    mc.set.seed = TRUE
  )
  optimism_estimates_protocol <- bind_rows(keep(optimism_estimates_protocol_list, ~ !is.null(.)))

  # Summarize optimism for extended model
  if (nrow(optimism_estimates_protocol) > 0) {
    optimism_summary_protocol <- optimism_estimates_protocol %>%
      group_by(.metric, .estimator) %>%
      summarise(
        mean_optimism = mean(optimism, na.rm = TRUE),
        n_boots = n(),
        ci_lower_95 = quantile(.estimate_on_original, probs = 0.025, na.rm = TRUE),
        ci_upper_95 = quantile(.estimate_on_original, probs = 0.975, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      filter(is.finite(mean_optimism))

    protocol_optimism_corrected <- apparent_performance_protocol %>%
      select(
        .metric,
        .estimator,
        apparent_estimate = .estimate
      ) %>%
      inner_join(optimism_summary_protocol, by = c(".metric", ".estimator")) %>%
      mutate(optimism_corrected = apparent_estimate - mean_optimism) %>%
      arrange(.metric, .estimator)

    print("Optimism-Corrected Performance (Extended Model):")
    print(protocol_optimism_corrected)

    write_rds(
      protocol_optimism_corrected,
      output_path('extended_model_optimism_corrected.rds')
    )
  } else {
    print("Insufficient bootstrap results for extended model optimism correction.")
    protocol_optimism_corrected <- NULL
  }

  print("--- Finished Optimism Correction for Re-estimated and Extended Models ---")
}

print("--- Finished Section 4 ---")
