#===============================================================================
# Section 2: Evaluate Old (Kasahara) Model Prediction Errors
#===============================================================================
print("--- Starting Section 2: Evaluate Old (Kasahara) Model Performance ---")

# Use the raw cleaned data `analysis_data_raw` which contains 'pred_ntprobnp_kasahara'
# Filter this cohort specifically for valid actual values and Kasahara predictions
required_cols_old_eval <- c("NTproBNP", "pred_ntprobnp_kasahara")
if (!all(required_cols_old_eval %in% names(analysis_data_raw))) {
  stop(
    "Required columns for old model evaluation not found in analysis_data_raw."
  )
}

data_eval_kasahara <- analysis_data_raw %>%
  select(StudyID, NTproBNP, pred_ntprobnp_kasahara) %>%
  filter(is.finite(NTproBNP) & is.finite(pred_ntprobnp_kasahara))

n_eval_kasahara <- nrow(data_eval_kasahara)
print(glue(
  "Number of subjects for Kasahara model evaluation: {n_eval_kasahara}"
))

# --- Proceed with Standard Evaluation using compute_metrics ---
if (n_eval_kasahara > 10) {
  print("Running standard evaluation via compute_metrics...")
  # Evaluate performance (predicts MEDIAN)
  old_model_errors <- compute_metrics(
    data_eval_kasahara, # Use the correctly filtered data
    target_col = "NTproBNP",
    predicted_col = "pred_ntprobnp_kasahara",
    include_paba = FALSE
  )

  # Save results
  write_rds(
    old_model_errors,
    output_path('kasahara_model_prediction_errors.rds')
  )
  if (isTRUE(CONFIG$write_csv_outputs)) {
    write_csv(
      old_model_errors,
      output_path('kasahara_model_prediction_errors.csv')
    )
  }
  print("Kasahara (Old) model performance metrics (from compute_metrics):")
  print(old_model_errors)

  # --- Bootstrap 95% CIs for external validation metrics (original Kasahara) ---
  print(glue(
    "--- Bootstrapping 95% CIs for external validation metrics ({CONFIG$bootstrap_reps} iterations) ---"
  ))

  calculate_validation_metrics_step <- function(split) {
    bs_data <- rsample::analysis(split) %>%
      filter(
        is.finite(NTproBNP) &
          is.finite(pred_ntprobnp_kasahara) &
          NTproBNP > 0 &
          pred_ntprobnp_kasahara > 0
      )
    if (nrow(bs_data) < 20) {
      return(NULL)
    }
    compute_metrics(
      bs_data,
      target_col = "NTproBNP",
      predicted_col = "pred_ntprobnp_kasahara",
      include_paba = FALSE
    ) %>%
      select(.metric, .estimator, .estimate)
  }

  boots_validation_kasahara <- rsample::bootstraps(
    data_eval_kasahara,
    times = CONFIG$bootstrap_reps,
    apparent = FALSE
  )
  validation_cores <- get_bootstrap_cores()
  message(glue(
    "Running external validation CI bootstrap on {validation_cores} core(s)..."
  ))
  validation_boot_list <- parallel::mclapply(
    boots_validation_kasahara$splits,
    function(split) {
      calculate_validation_metrics_step(split)
    },
    mc.cores = validation_cores,
    mc.set.seed = TRUE
  )
  validation_boot <- bind_rows(purrr::keep(validation_boot_list, ~ !is.null(.)))

  if (nrow(validation_boot) > 10) {
    validation_ci <- validation_boot %>%
      group_by(.metric, .estimator) %>%
      summarise(
        ci_lower_95 = quantile(.estimate, probs = 0.025, na.rm = TRUE),
        ci_upper_95 = quantile(.estimate, probs = 0.975, na.rm = TRUE),
        n_boots = n(),
        .groups = "drop"
      )

    kasahara_validation_metrics_with_ci <- old_model_errors %>%
      rename(estimate = .estimate) %>%
      left_join(validation_ci, by = c(".metric", ".estimator")) %>%
      arrange(.metric, .estimator)

    print("Kasahara external validation metrics with 95% bootstrap CIs:")
    print(kasahara_validation_metrics_with_ci)

    write_rds(
      kasahara_validation_metrics_with_ci,
      output_path("kasahara_validation_metrics_with_ci.rds")
    )
    if (isTRUE(CONFIG$write_csv_outputs)) {
      write_csv(
        kasahara_validation_metrics_with_ci,
        output_path("kasahara_validation_metrics_with_ci.csv")
      )
    }
  } else {
    print("Insufficient bootstrap results for external validation CI estimation.")
  }

  # --- Log-Log Calibration (Kasahara): OLS primary, PaBa sensitivity ---
  # Recalibration equation (log10 scale):
  # log10(NTproBNP) = I + S*log10(pred_kasahara)
  # pred_kasahara_recal = 10^(I + S*log10(pred_kasahara))
  # Interpretation: I = calibration-in-the-large, S = calibration slope.
  loglog_data <- data_eval_kasahara %>%
    filter(
      is.finite(NTproBNP) &
        is.finite(pred_ntprobnp_kasahara) &
        NTproBNP > 0 &
        pred_ntprobnp_kasahara > 0
    )
  if (nrow(loglog_data) >= 10) {
    # Regular (OLS) log-log calibration: log10(Actual) ~ log10(Predicted)
    loglog_fit <- lm(
      log10(NTproBNP) ~ log10(pred_ntprobnp_kasahara),
      data = loglog_data
    )
    kasahara_loglog_intercept_ols <- coef(loglog_fit)[1]
    kasahara_loglog_slope_ols <- coef(loglog_fit)[2]
    loglog_ci <- tryCatch(confint(loglog_fit), error = function(e) NULL)
    intercept_ci <- if (!is.null(loglog_ci)) {
      loglog_ci["(Intercept)", ]
    } else {
      c(NA_real_, NA_real_)
    }
    slope_ci <- if (!is.null(loglog_ci)) {
      loglog_ci["log10(pred_ntprobnp_kasahara)", ]
    } else {
      c(NA_real_, NA_real_)
    }

    # Passing-Bablok sensitivity fit (reported in a terminal comparison section)
    paba_loglog_kasahara <- tryCatch(
      {
        mcreg(
          x = log10(loglog_data$pred_ntprobnp_kasahara),
          y = log10(loglog_data$NTproBNP),
          method.reg = "PaBa",
          na.rm = TRUE
        )
      },
      error = function(e) NULL
    )
    kasahara_loglog_intercept_paba <- NA_real_
    kasahara_loglog_slope_paba <- NA_real_
    paba_intercept_ci <- c(NA_real_, NA_real_)
    paba_slope_ci <- c(NA_real_, NA_real_)
    if (!is.null(paba_loglog_kasahara) && is_mcr_result(paba_loglog_kasahara)) {
      paba_para <- paba_loglog_kasahara@para
      kasahara_loglog_intercept_paba <- paba_para["Intercept", "EST"]
      kasahara_loglog_slope_paba <- paba_para["Slope", "EST"]
      paba_intercept_ci <- c(
        paba_para["Intercept", "LCI"],
        paba_para["Intercept", "UCI"]
      )
      paba_slope_ci <- c(
        paba_para["Slope", "LCI"],
        paba_para["Slope", "UCI"]
      )
    }

    kasahara_loglog_summary <- tibble(
      method = c("OLS", "PaBa"),
      intercept = c(
        kasahara_loglog_intercept_ols,
        kasahara_loglog_intercept_paba
      ),
      intercept_ci_lower = c(intercept_ci[1], paba_intercept_ci[1]),
      intercept_ci_upper = c(intercept_ci[2], paba_intercept_ci[2]),
      slope = c(kasahara_loglog_slope_ols, kasahara_loglog_slope_paba),
      slope_ci_lower = c(slope_ci[1], paba_slope_ci[1]),
      slope_ci_upper = c(slope_ci[2], paba_slope_ci[2]),
      n = nrow(loglog_data)
    )
    write_rds(
      kasahara_loglog_summary,
      output_path("kasahara_loglog_calibration_summary.rds")
    )
    if (isTRUE(CONFIG$write_csv_outputs)) {
      write_csv(
        kasahara_loglog_summary,
        output_path("kasahara_loglog_calibration_summary.csv")
      )
    }
    print("Kasahara log-log recalibration coefficients (primary OLS; PaBa sensitivity):")
    print(kasahara_loglog_summary)

    # Add recalibrated predictions to datasets (OLS primary, PaBa sensitivity)
    analysis_data_raw <- analysis_data_raw %>%
      mutate(
        pred_ntprobnp_kasahara_recal_ols = ifelse(
          is.finite(pred_ntprobnp_kasahara) & pred_ntprobnp_kasahara > 0,
          10^(kasahara_loglog_intercept_ols +
            kasahara_loglog_slope_ols * log10(pred_ntprobnp_kasahara)),
          NA_real_
        ),
        pred_ntprobnp_kasahara_recal_paba = ifelse(
          is.finite(pred_ntprobnp_kasahara) &
            pred_ntprobnp_kasahara > 0 &
            is.finite(kasahara_loglog_intercept_paba) &
            is.finite(kasahara_loglog_slope_paba),
          10^(kasahara_loglog_intercept_paba +
            kasahara_loglog_slope_paba * log10(pred_ntprobnp_kasahara)),
          NA_real_
        )
      )
    data_eval_kasahara <- data_eval_kasahara %>%
      mutate(
        pred_ntprobnp_kasahara_recal_ols = ifelse(
          is.finite(pred_ntprobnp_kasahara) & pred_ntprobnp_kasahara > 0,
          10^(kasahara_loglog_intercept_ols +
            kasahara_loglog_slope_ols * log10(pred_ntprobnp_kasahara)),
          NA_real_
        ),
        pred_ntprobnp_kasahara_recal_paba = ifelse(
          is.finite(pred_ntprobnp_kasahara) &
            pred_ntprobnp_kasahara > 0 &
            is.finite(kasahara_loglog_intercept_paba) &
            is.finite(kasahara_loglog_slope_paba),
          10^(kasahara_loglog_intercept_paba +
            kasahara_loglog_slope_paba * log10(pred_ntprobnp_kasahara)),
          NA_real_
        )
      )

    # Evaluate recalibrated Kasahara predictions (OLS)
    old_model_errors_recal_ols <- compute_metrics(
      data_eval_kasahara,
      target_col = "NTproBNP",
      predicted_col = "pred_ntprobnp_kasahara_recal_ols"
    )
    write_rds(
      old_model_errors_recal_ols,
      output_path('kasahara_model_prediction_errors_recal_ols.rds')
    )
    if (isTRUE(CONFIG$write_csv_outputs)) {
      write_csv(
        old_model_errors_recal_ols,
        output_path('kasahara_model_prediction_errors_recal_ols.csv')
      )
    }
    print("Kasahara (Recalibrated OLS) model performance metrics:")
    print(old_model_errors_recal_ols)

    # Evaluate PaBa-recalibrated predictions for end-of-report sensitivity comparison
    paba_eval_data <- data_eval_kasahara %>%
      filter(
        is.finite(NTproBNP) &
          is.finite(pred_ntprobnp_kasahara_recal_paba) &
          NTproBNP > 0 &
          pred_ntprobnp_kasahara_recal_paba > 0
      )
    if (nrow(paba_eval_data) >= 10) {
      old_model_errors_recal_paba <- compute_metrics(
        paba_eval_data,
        target_col = "NTproBNP",
        predicted_col = "pred_ntprobnp_kasahara_recal_paba"
      )
      write_rds(
        old_model_errors_recal_paba,
        output_path("kasahara_model_prediction_errors_recal_paba.rds")
      )
      if (isTRUE(CONFIG$write_csv_outputs)) {
        write_csv(
          old_model_errors_recal_paba,
          output_path("kasahara_model_prediction_errors_recal_paba.csv")
        )
      }
      print("Kasahara (Recalibrated PaBa sensitivity) model performance metrics:")
      print(old_model_errors_recal_paba)
    } else {
      warning("Insufficient valid data for PaBa recalibration sensitivity metrics.")
    }

    # --- Optimism Correction for Recalibrated Kasahara (Frank Harrell's method) ---
    # 1. Apparent performance on full sample (already computed above)
    # 2. For each bootstrap: fit on bootstrap, calc on bootstrap (θ_boot),
    #    apply to original, calc on original (θ_orig), optimism = θ_boot - θ_orig
    # 3. Optimism-corrected = Apparent - mean(Optimism)
    print(glue(
      "--- Bootstrapping Optimism Correction for Recalibrated Kasahara ({CONFIG$bootstrap_reps} iterations) ---"
    ))

    calculate_kasahara_optimism_step <- function(split) {
      bs_data <- rsample::analysis(split)
      orig_data <- rsample::assessment(split)

      # Filter to valid data
      bs_data <- bs_data %>%
        filter(
          is.finite(NTproBNP) &
            is.finite(pred_ntprobnp_kasahara) &
            NTproBNP > 0 &
            pred_ntprobnp_kasahara > 0
        )
      orig_data <- orig_data %>%
        filter(
          is.finite(NTproBNP) &
            is.finite(pred_ntprobnp_kasahara) &
            NTproBNP > 0 &
            pred_ntprobnp_kasahara > 0
        )

      if (nrow(bs_data) < 20 || nrow(orig_data) < 20) {
        return(NULL)
      }

      # Fit OLS recalibration on bootstrap sample
      bs_fit <- lm(
        log10(NTproBNP) ~ log10(pred_ntprobnp_kasahara),
        data = bs_data
      )
      I_bs <- coef(bs_fit)[1]
      S_bs <- coef(bs_fit)[2]

      if (!is.finite(I_bs) || !is.finite(S_bs)) {
        return(NULL)
      }

      # Calculate recalibrated predictions using bootstrap-derived parameters
      bs_data <- bs_data %>%
        mutate(
          pred_recal_bs = 10^(I_bs + S_bs * log10(pred_ntprobnp_kasahara))
        ) %>%
        filter(is.finite(pred_recal_bs))

      orig_data <- orig_data %>%
        mutate(
          pred_recal_bs = 10^(I_bs + S_bs * log10(pred_ntprobnp_kasahara))
        ) %>%
        filter(is.finite(pred_recal_bs))

      if (nrow(bs_data) < 10 || nrow(orig_data) < 10) {
        return(NULL)
      }

      # Performance on bootstrap sample (apparent on bootstrap)
      perf_on_bootstrap <- compute_metrics(
        bs_data,
        "NTproBNP",
        "pred_recal_bs"
      ) %>%
        select(.metric, .estimator, .estimate_on_bootstrap = .estimate)

      # Performance on original sample using bootstrap-fitted model
      perf_on_original <- compute_metrics(
        orig_data,
        "NTproBNP",
        "pred_recal_bs"
      ) %>%
        select(.metric, .estimator, .estimate_on_original = .estimate)

      # Calculate optimism = bootstrap - original
      metrics_step <- inner_join(
        perf_on_bootstrap,
        perf_on_original,
        by = c(".metric", ".estimator")
      ) %>%
        mutate(optimism = .estimate_on_bootstrap - .estimate_on_original) %>%
        filter(is.finite(optimism))

      return(metrics_step)
    }

    # Run bootstrap
    boots_kasahara <- rsample::bootstraps(
      data_eval_kasahara,
      times = CONFIG$bootstrap_reps,
      apparent = FALSE
    )
    optimism_cores <- get_bootstrap_cores()
    message(glue(
      "Running Kasahara optimism correction on {optimism_cores} core(s)..."
    ))
    kasahara_optimism_list <- parallel::mclapply(
      boots_kasahara$splits,
      function(split) {
        calculate_kasahara_optimism_step(split)
      },
      mc.cores = optimism_cores,
      mc.set.seed = TRUE
    )
    kasahara_optimism_estimates <- bind_rows(
      purrr::keep(kasahara_optimism_list, ~ !is.null(.))
    )

    # Summarize optimism and calculate corrected estimates
    if (nrow(kasahara_optimism_estimates) > 10) {
      kasahara_optimism_summary <- kasahara_optimism_estimates %>%
        group_by(.metric, .estimator) %>%
        summarise(
          mean_optimism = mean(optimism, na.rm = TRUE),
          n_boots = n(),
          ci_lower_95 = quantile(.estimate_on_original, probs = 0.025, na.rm = TRUE),
          ci_upper_95 = quantile(.estimate_on_original, probs = 0.975, na.rm = TRUE),
          .groups = 'drop'
        ) %>%
        filter(is.finite(mean_optimism))

      # Merge with apparent performance
      kasahara_recal_optimism_corrected <- old_model_errors_recal_ols %>%
        select(
          .metric,
          .estimator,
          threshold_label,
          apparent_estimate = .estimate
        ) %>%
        inner_join(kasahara_optimism_summary, by = c(".metric", ".estimator")) %>%
        mutate(optimism_corrected = apparent_estimate - mean_optimism) %>%
        select(
          .metric,
          .estimator,
          apparent_estimate,
          mean_optimism,
          optimism_corrected,
          ci_lower_95,
          ci_upper_95,
          n_boots
        ) %>%
        arrange(.metric, .estimator)

      print("Optimism-Corrected Performance Metrics (Recalibrated Kasahara):")
      print(kasahara_recal_optimism_corrected)

      # Save results
      write_rds(
        kasahara_recal_optimism_corrected,
        output_path('kasahara_recal_optimism_corrected.rds')
      )
      if (isTRUE(CONFIG$write_csv_outputs)) {
        write_csv(
          kasahara_recal_optimism_corrected,
          output_path('kasahara_recal_optimism_corrected.csv')
        )
      }
    } else {
      print("Insufficient bootstrap results for Kasahara optimism correction.")
      kasahara_recal_optimism_corrected <- NULL
    }

    # Store recalibration parameters for use elsewhere
    kasahara_recal_params <- c(kasahara_loglog_intercept_ols, kasahara_loglog_slope_ols)

  } else {
    print("Insufficient data for log-log calibration of Kasahara model.")
  }
} else {
  print(
    "Insufficient data to evaluate Kasahara model performance via compute_metrics."
  )
  old_model_errors <- NULL # Ensure it's defined for later checks
}

print("--- Finished Section 2 ---")


#===============================================================================
# Section 3: Visualize Old (Kasahara) Model Performance
#===============================================================================
print("--- Starting Section 3: Visualize Old (Kasahara) Model Performance ---")

# Use the data specifically filtered for Kasahara evaluation: data_eval_kasahara
if (exists("data_eval_kasahara") && nrow(data_eval_kasahara) > 10) {
  # Check if data exists and has rows
  # Calibration plots (Log and Linear Scale)
  cal_plot_old_log <- make_calibration_plot(
    data = data_eval_kasahara,
    actual_col = "NTproBNP",
    pred_col = "pred_ntprobnp_kasahara", # Use correct data
    model_name = "Kasahara Model",
    log_scale = TRUE
  )
  print_plot(cal_plot_old_log)
  ggsave(
    output_path('calibration_plot_kasahara_log.png'),
    plot = cal_plot_old_log,
    height = 6,
    width = 6,
    dpi = 300
  )

  # Log-Log calibration plot + actual vs predicted (Kasahara Recalibrated OLS), if available
  if ("pred_ntprobnp_kasahara_recal_ols" %in% names(data_eval_kasahara)) {
    cal_plot_old_loglog_recal_ols <- make_loglog_calibration_plot(
      data = data_eval_kasahara,
      actual_col = "NTproBNP",
      pred_col = "pred_ntprobnp_kasahara_recal_ols",
      model_name = "Kasahara Model (Recalibrated OLS)"
    )
    print_plot(cal_plot_old_loglog_recal_ols)
    ggsave(
      output_path('calibration_plot_kasahara_recal_ols_loglog.png'),
      plot = cal_plot_old_loglog_recal_ols,
      height = 6,
      width = 6,
      dpi = 300
    )

    # Bland-Altman plots for Recalibrated OLS Kasahara
    ba_kasahara_recal_ols_rel <- make_bland_altman_plot(
      data = data_eval_kasahara,
      actual_col = "NTproBNP",
      pred_col = "pred_ntprobnp_kasahara_recal_ols",
      model_name = "Kasahara Model (Recalibrated OLS)",
      type = "relative"
    )
    print_plot(ba_kasahara_recal_ols_rel)
    ggsave(
      output_path('BA_kasahara_recal_ols_relative.png'),
      plot = ba_kasahara_recal_ols_rel,
      height = 6,
      width = 6,
      dpi = 300
    )

    ba_kasahara_recal_ols_abs <- make_bland_altman_plot(
      data = data_eval_kasahara,
      actual_col = "NTproBNP",
      pred_col = "pred_ntprobnp_kasahara_recal_ols",
      model_name = "Kasahara Model (Recalibrated OLS)",
      type = "absolute"
    )
    print_plot(ba_kasahara_recal_ols_abs)
    ggsave(
      output_path('BA_kasahara_recal_ols_absolute.png'),
      plot = ba_kasahara_recal_ols_abs,
      height = 6,
      width = 6,
      dpi = 300
    )
  }

  # Bland-Altman plots (Relative and Absolute Error)
  ba_plot_old_rel <- make_bland_altman_plot(
    data = data_eval_kasahara,
    actual_col = "NTproBNP",
    pred_col = "pred_ntprobnp_kasahara", # Use correct data
    model_name = "Kasahara Model",
    type = "relative"
  )
  print_plot(ba_plot_old_rel)
  ggsave(
    output_path('BA_kasahara_relative.png'),
    plot = ba_plot_old_rel,
    height = 6,
    width = 6,
    dpi = 300
  )

  ba_plot_old_abs <- make_bland_altman_plot(
    data = data_eval_kasahara,
    actual_col = "NTproBNP",
    pred_col = "pred_ntprobnp_kasahara", # Use correct data
    model_name = "Kasahara Model",
    type = "absolute"
  )
  print_plot(ba_plot_old_abs)
  ggsave(
    output_path('BA_kasahara_absolute.png'),
    plot = ba_plot_old_abs,
    height = 6,
    width = 6,
    dpi = 300
  )
} else {
  print("Insufficient data for Kasahara model visualizations.")
}

print("--- Finished Section 3 ---")
