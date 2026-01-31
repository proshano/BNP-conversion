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
    predicted_col = "pred_ntprobnp_kasahara"
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

  # Assign the result from the compute_metrics PaBa attempt for plotting
  # (This uses the run_mcreg helper with tryCatch)
  paba_old_model_obj <- run_mcreg(
    data_eval_kasahara,
    x_col = pred_ntprobnp_kasahara,
    y_col = NTproBNP
  )

  # --- Log-Log Calibration (Kasahara): OLS + Passing-Bablok ---
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

    # Passing-Bablok log-log calibration: log10(Actual) vs log10(Predicted)
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
    if (!is.null(paba_loglog_kasahara) && is_mcr_result(paba_loglog_kasahara)) {
      kasahara_loglog_intercept_paba <- paba_loglog_kasahara@para[
        "Intercept",
        "EST"
      ]
      kasahara_loglog_slope_paba <- paba_loglog_kasahara@para["Slope", "EST"]
    }

    kasahara_loglog_summary <- tibble(
      method = c("OLS", "PaBa"),
      intercept = c(
        kasahara_loglog_intercept_ols,
        kasahara_loglog_intercept_paba
      ),
      slope = c(kasahara_loglog_slope_ols, kasahara_loglog_slope_paba),
      n = nrow(loglog_data)
    )
    if (isTRUE(CONFIG$write_csv_outputs)) {
      write_csv(
        kasahara_loglog_summary,
        output_path("kasahara_loglog_calibration_summary.csv")
      )
    }
    print("Kasahara log-log calibration coefficients (OLS and PaBa):")
    print(kasahara_loglog_summary)

    # Add recalibrated predictions to datasets (OLS and PaBa)
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

    # Evaluate recalibrated Kasahara predictions (OLS and PaBa)
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

    if (
      is.finite(kasahara_loglog_intercept_paba) &&
        is.finite(kasahara_loglog_slope_paba)
    ) {
      old_model_errors_recal_paba <- compute_metrics(
        data_eval_kasahara,
        target_col = "NTproBNP",
        predicted_col = "pred_ntprobnp_kasahara_recal_paba"
      )
      write_rds(
        old_model_errors_recal_paba,
        output_path('kasahara_model_prediction_errors_recal_paba.rds')
      )
      if (isTRUE(CONFIG$write_csv_outputs)) {
        write_csv(
          old_model_errors_recal_paba,
          output_path('kasahara_model_prediction_errors_recal_paba.csv')
        )
      }
      print("Kasahara (Recalibrated PaBa) model performance metrics:")
      print(old_model_errors_recal_paba)
    } else {
      print("Kasahara PaBa log-log recalibration failed or insufficient data.")
    }
  } else {
    print("Insufficient data for log-log calibration of Kasahara model.")
  }
} else {
  print(
    "Insufficient data to evaluate Kasahara model performance via compute_metrics."
  )
  old_model_errors <- NULL # Ensure it's defined for later checks
  paba_old_model_obj <- NULL
}

print("--- Finished Section 2 ---")


#===============================================================================
# Section 3: Visualize Old (Kasahara) Model Performance
#===============================================================================
print("--- Starting Section 3: Visualize Old (Kasahara) Model Performance ---")

# Use the data specifically filtered for Kasahara evaluation: data_eval_kasahara
if (exists("data_eval_kasahara") && nrow(data_eval_kasahara) > 10) {
  # Check if data exists and has rows
  # Actual vs Predicted (Kasahara)
  avp_kasahara <- make_actual_vs_pred_plot(
    data = data_eval_kasahara,
    actual_col = "NTproBNP",
    pred_col = "pred_ntprobnp_kasahara",
    model_name = "Kasahara Model",
    log_scale = TRUE
  )
  print_plot(avp_kasahara)
  ggsave(
    output_path('actual_vs_pred_kasahara.png'),
    plot = avp_kasahara,
    height = 6,
    width = 6,
    dpi = 300
  )

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

  cal_plot_old_lin <- make_calibration_plot(
    data = data_eval_kasahara,
    actual_col = "NTproBNP",
    pred_col = "pred_ntprobnp_kasahara", # Use correct data
    model_name = "Kasahara Model",
    log_scale = FALSE
  )
  print_plot(cal_plot_old_lin)
  ggsave(
    output_path('calibration_plot_kasahara_lin.png'),
    plot = cal_plot_old_lin,
    height = 6,
    width = 6,
    dpi = 300
  )

  # Log-Log calibration plot (Kasahara)
  cal_plot_old_loglog <- make_loglog_calibration_plot(
    data = data_eval_kasahara,
    actual_col = "NTproBNP",
    pred_col = "pred_ntprobnp_kasahara",
    model_name = "Kasahara Model"
  )
  print_plot(cal_plot_old_loglog)
  ggsave(
    output_path('calibration_plot_kasahara_loglog.png'),
    plot = cal_plot_old_loglog,
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

    avp_kasahara_recal_ols <- make_actual_vs_pred_plot(
      data = data_eval_kasahara,
      actual_col = "NTproBNP",
      pred_col = "pred_ntprobnp_kasahara_recal_ols",
      model_name = "Kasahara Model (Recalibrated OLS)",
      log_scale = TRUE
    )
    print_plot(avp_kasahara_recal_ols)
    ggsave(
      output_path('actual_vs_pred_kasahara_recal_ols.png'),
      plot = avp_kasahara_recal_ols,
      height = 6,
      width = 6,
      dpi = 300
    )
  }

  # Log-Log calibration plot + actual vs predicted (Kasahara Recalibrated PaBa), if available
  if ("pred_ntprobnp_kasahara_recal_paba" %in% names(data_eval_kasahara)) {
    cal_plot_old_loglog_recal_paba <- make_loglog_calibration_plot(
      data = data_eval_kasahara,
      actual_col = "NTproBNP",
      pred_col = "pred_ntprobnp_kasahara_recal_paba",
      model_name = "Kasahara Model (Recalibrated PaBa)"
    )
    print_plot(cal_plot_old_loglog_recal_paba)
    ggsave(
      output_path('calibration_plot_kasahara_recal_paba_loglog.png'),
      plot = cal_plot_old_loglog_recal_paba,
      height = 6,
      width = 6,
      dpi = 300
    )

    avp_kasahara_recal_paba <- make_actual_vs_pred_plot(
      data = data_eval_kasahara,
      actual_col = "NTproBNP",
      pred_col = "pred_ntprobnp_kasahara_recal_paba",
      model_name = "Kasahara Model (Recalibrated PaBa)",
      log_scale = TRUE
    )
    print_plot(avp_kasahara_recal_paba)
    ggsave(
      output_path('actual_vs_pred_kasahara_recal_paba.png'),
      plot = avp_kasahara_recal_paba,
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

  # Passing-Bablok plot (use object created in Section 2)
  plot_mcreg(
    paba_old_model_obj,
    filename = output_path('paba_plot_kasahara.png')
  )
} else {
  print("Insufficient data for Kasahara model visualizations.")
}

print("--- Finished Section 3 ---")
