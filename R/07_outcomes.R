#===============================================================================
# Section 5: Association with Clinical Outcomes (using RECALIBRATED MEDIAN Predictions)
# Uses log10 scale throughout (matching Kasahara) and OLS recalibration
#===============================================================================
print("--- Starting Section 5: Association with Clinical Outcomes ---")
# Use the 'md' dataframe which contains the final recalibrated median predictions
# --- Derive Biomarker Categories ---
print("Deriving biomarker categories...")
# Define target NT-proBNP thresholds
target_NTs <- CONFIG$target_ntprobnp_thresholds
cut_points_ntprobnp <- c(-Inf, target_NTs, Inf)
# Calculate corresponding BNP thresholds using the RECALIBRATED KASAHARA model
# (simpler 2-parameter recalibration as recommended)
print("Calculating corresponding BNP thresholds using recalibrated Kasahara model...")
median_covars_for_thresholds <- md %>%
  summarise(
    across(c(Age, BMI, CrCl, Hb_gdl), \(x) median(x, na.rm = TRUE)),
    Sex_M = 0, # Female (base case)
    AF = 0,    # No AF (base case)
    .groups = 'drop'
  )
print("Median covariate grid for threshold calculation:")
print(median_covars_for_thresholds)

# Use recalibrated Kasahara model (requires kasahara_recal_params from 05_kasahara_eval.R)
if (!exists("kasahara_recal_params") || length(kasahara_recal_params) != 2) {
  warning("Kasahara recalibration parameters not found. Using fit_ols as fallback.")
  corresponding_BNPs_new <- sapply(target_NTs, function(nt) {
    get_BNP_for_NT_new_model(
      nt,
      median_covars_for_thresholds,
      fit_ols,
      log10_recal_params,
      use_mean = FALSE
    )
  })
} else {
  corresponding_BNPs_new <- sapply(target_NTs, function(nt) {
    get_BNP_for_NT_kasahara_recal(
      nt,
      median_covars_for_thresholds,
      kasahara_recal_params[1],  # intercept
      kasahara_recal_params[2]   # slope
    )
  })
}
print("Calculated corresponding BNP thresholds (rounded to nearest integer):")
bnp_threshold_summary <- tibble(
  NTproBNP_threshold = target_NTs,
  BNP_threshold = round(corresponding_BNPs_new)
)
print(bnp_threshold_summary)
if (any(is.na(corresponding_BNPs_new))) {
  warning(
    "Could not determine all corresponding BNP thresholds using the recalibrated model."
  )
}
valid_new_bnps <- round(corresponding_BNPs_new[
  !is.na(corresponding_BNPs_new) & is.finite(corresponding_BNPs_new)
])
valid_new_bnps <- sort(unique(valid_new_bnps[valid_new_bnps > 0])) # Ensure unique and positive
cut_points_bnp_derived_new <- c(-Inf, valid_new_bnps, Inf)
print("Final derived cut points for BNP:")
print(cut_points_bnp_derived_new)
#--- Covariate-Specific BNP Thresholds (Sensitivity Analysis) ---
print("Computing sensitivity of BNP thresholds to individual covariates...")
# Base case: medians/mode as already computed
base_case <- median_covars_for_thresholds

# Define variations: one parameter at a time
# Base case uses Female (Sex_M=0) and No AF (AF=0) with median continuous values
# Table will show each covariate value as a column
sensitivity_profiles <- bind_rows(
  base_case,                             # Base case: Female, No AF, median continuous
  base_case %>% mutate(Sex_M = 1),       # Male
  base_case %>% mutate(AF = 1),          # AF = Yes
  base_case %>% mutate(Age = 50),        # Age 50
  base_case %>% mutate(CrCl = 30),       # CrCl 30
  base_case %>% mutate(Hb_gdl = 10),     # Hb 10
  base_case %>% mutate(BMI = 20)         # BMI 20
)

# Function to compute thresholds using recalibrated Kasahara
compute_threshold_row <- function(covars, recal_I, recal_S) {
  tibble(
    BNP_for_NT100 = get_BNP_for_NT_kasahara_recal(100, covars, recal_I, recal_S),
    BNP_for_NT200 = get_BNP_for_NT_kasahara_recal(200, covars, recal_I, recal_S),
    BNP_for_NT1500 = get_BNP_for_NT_kasahara_recal(1500, covars, recal_I, recal_S)
  )
}

# Compute thresholds using recalibrated Kasahara model
if (exists("kasahara_recal_params") && length(kasahara_recal_params) == 2) {
  recal_I <- kasahara_recal_params[1]
  recal_S <- kasahara_recal_params[2]

  threshold_sensitivity <- sensitivity_profiles %>%
    rowwise() %>%
    mutate(
      BNP_for_NT100 = get_BNP_for_NT_kasahara_recal(
        100, pick(Age, BMI, CrCl, Hb_gdl, Sex_M, AF), recal_I, recal_S
      ),
      BNP_for_NT200 = get_BNP_for_NT_kasahara_recal(
        200, pick(Age, BMI, CrCl, Hb_gdl, Sex_M, AF), recal_I, recal_S
      ),
      BNP_for_NT1500 = get_BNP_for_NT_kasahara_recal(
        1500, pick(Age, BMI, CrCl, Hb_gdl, Sex_M, AF), recal_I, recal_S
      )
    ) %>%
    ungroup()
} else {
  # Fallback to extended model if Kasahara recal params not available
  threshold_sensitivity <- sensitivity_profiles %>%
    rowwise() %>%
    mutate(
      BNP_for_NT100 = get_BNP_for_NT_new_model(
        100, pick(Age, BMI, CrCl, Hb_gdl, Sex_M, AF), fit_ols, log10_recal_params, FALSE
      ),
      BNP_for_NT200 = get_BNP_for_NT_new_model(
        200, pick(Age, BMI, CrCl, Hb_gdl, Sex_M, AF), fit_ols, log10_recal_params, FALSE
      ),
      BNP_for_NT1500 = get_BNP_for_NT_new_model(
        1500, pick(Age, BMI, CrCl, Hb_gdl, Sex_M, AF), fit_ols, log10_recal_params, FALSE
      )
    ) %>%
    ungroup()
}

# Round BNP thresholds and convert Sex/AF to labels
threshold_sensitivity <- threshold_sensitivity %>%
  mutate(
    across(starts_with("BNP_for"), round),
    Sex = ifelse(Sex_M == 0, "Female", ifelse(Sex_M == 1, "Male", "0.5")),
    AF = ifelse(AF == 0, "No", ifelse(AF == 1, "Yes", as.character(AF))),
    Age = round(Age, 0),
    BMI = round(BMI, 1),
    CrCl = round(CrCl, 0),
    Hb = round(Hb_gdl, 1)
  ) %>%
  select(Age, BMI, CrCl, Hb, Sex, AF, BNP_for_NT100, BNP_for_NT200, BNP_for_NT1500)

print("BNP Threshold Sensitivity to Individual Covariates:")
print(threshold_sensitivity)
if (isTRUE(CONFIG$write_csv_outputs)) {
  write_csv(
    threshold_sensitivity,
    output_path("bnp_threshold_sensitivity.csv")
  )
}
saveRDS(threshold_sensitivity, output_path("bnp_threshold_sensitivity.rds"))
# Define BNP category labels if all three thresholds are available
if (length(valid_new_bnps) == 3) {
  bnp_labels <- c(
    glue("<{valid_new_bnps[1]}"),
    glue("{valid_new_bnps[1]}-<{valid_new_bnps[2]}"),
    glue("{valid_new_bnps[2]}-<{valid_new_bnps[3]}"),
    glue(">={valid_new_bnps[3]}")
  )
} else {
  warning(
    "Derived BNP thresholds are incomplete; BNP categories will be set to NA."
  )
  bnp_labels <- NULL
}
# Prepare data for logistic regression models
logit_md <- md %>%
  mutate(
    NTproBNP_cat = make_ntprobnp_cat(NTproBNP, CONFIG$class_breaks),
    BNP_cat_derived = if (!is.null(bnp_labels)) {
      make_interval_cat(BNP, valid_new_bnps, bnp_labels)
    } else {
      NA
    },
    pred_new_median_cat = make_ntprobnp_cat(
      pred_ntprobnp_new_median,
      CONFIG$class_breaks
    )
  ) %>%
  select(
    StudyID,
    composite_outcome,
    NTproBNP,
    BNP,
    pred_ntprobnp_new_median,
    NTproBNP_cat,
    BNP_cat_derived,
    pred_new_median_cat
  ) %>%
  filter(is.finite(composite_outcome)) # Ensure outcome is not missing
if (nrow(logit_md) < 30) {
  # Arbitrary minimum for logistic regression
  warning(
    "Insufficient data (<30 subjects) for logistic regression outcome analysis. Skipping."
  )
} else {
  print(glue(
    "Number of observations for logistic regression: {nrow(logit_md)}"
  ))
  # Set factors to unordered so GLM uses category-vs-reference (not polynomial) contrasts
  logit_md <- logit_md %>%
    mutate(
      across(ends_with("_cat"), ~ factor(.x, ordered = FALSE))
    ) %>%
    mutate(
      across(ends_with("_cat"), ~ fct_relevel(.x, levels(.)[1]))
    )

  # --- Compute n/events per category for reporting ---
  ntprobnp_cat_summary <- logit_md %>%
    group_by(NTproBNP_cat) %>%
    summarise(
      n = n(),
      events = sum(composite_outcome == 1, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    rename(Category = NTproBNP_cat)

  pred_cat_summary <- logit_md %>%
    group_by(pred_new_median_cat) %>%
    summarise(
      n = n(),
      events = sum(composite_outcome == 1, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    rename(Category = pred_new_median_cat)

  print("NT-proBNP category counts (measured):")
  print(ntprobnp_cat_summary)
  print("NT-proBNP category counts (predicted):")
  print(pred_cat_summary)

  # Save for report use
  saveRDS(ntprobnp_cat_summary, output_path("ntprobnp_cat_summary.rds"))
  saveRDS(pred_cat_summary, output_path("pred_cat_summary.rds"))

  # --- Fit Logistic Regression Models ---
  # Categorical logistic regression:
  # logit(P(outcome=1)) = a0 + sum_j a_j * I(category=j), reference = lowest category.
  # Firth logistic regression uses penalized likelihood to reduce small-sample bias.
  print(
    "Fitting logistic regression models (GLM and Firth Logistic Regression)"
  )
  drop_part_column <- function(ft) {
    if (!inherits(ft, "flextable")) {
      return(ft)
    }
    # Avoid flextable::col_keys() because it may not be exported in some versions.
    keys <- if (!is.null(ft$col_keys)) {
      ft$col_keys
    } else if (!is.null(ft$body) && !is.null(ft$body$dataset)) {
      names(ft$body$dataset)
    } else {
      character(0)
    }
    if ("part" %in% keys) {
      ft <- flextable::delete_columns(ft, j = "part")
    }
    ft
  }
  # Measured NT-proBNP categories vs outcome
  model_a_glm <- glm(
    composite_outcome ~ NTproBNP_cat,
    data = logit_md,
    family = binomial
  )
  model_a_firth <- logistf(
    composite_outcome ~ NTproBNP_cat,
    data = logit_md,
    pl = FALSE
  ) # Firth for sparse data robustness
  table_a <- modelsummary(
    list('GLM' = model_a_glm, 'Firth GLM' = model_a_firth),
    estimate = "{estimate} [{conf.low}, {conf.high}]",
    statistic = "p.value",
    exponentiate = TRUE,
    gof_map = NA,
    output = "flextable",
    title = 'Association between Measured NT-proBNP Categories and Composite Outcome (Odds Ratios)'
  )
  table_a <- drop_part_column(table_a)
  print(table_a)

  # Predicted (recalibrated median, OLS) NT-proBNP categories vs outcome
  model_c_glm <- glm(
    composite_outcome ~ pred_new_median_cat,
    data = logit_md,
    family = binomial
  )
  model_c_firth <- logistf(
    composite_outcome ~ pred_new_median_cat,
    data = logit_md,
    pl = FALSE
  )
  table_c <- modelsummary(
    list('GLM' = model_c_glm, 'Firth GLM' = model_c_firth),
    estimate = "{estimate} [{conf.low}, {conf.high}]",
    statistic = "p.value",
    exponentiate = TRUE,
    gof_map = NA,
    output = "flextable",
    title = 'Association between Predicted NT-proBNP Categories (New Model) and Composite Outcome (Odds Ratios)'
  )
  table_c <- drop_part_column(table_c)
  print(table_c)

  print("Finished fitting logistic models.")

  # --- Visualize Associations (Continuous Predictors using rms::lrm and contrast) ---
  # Spline logistic regression:
  # logit(P(outcome=1)) = a0 + g(X), where g is an rcs() spline (k knots).
  print("Generating association plots (continuous predictors)...")

  # Helper functions for plotting contrasts
  get_plot_range <- function(x, probs = c(0.01, 0.99)) {
    q <- suppressWarnings(
      quantile(x, probs = probs, na.rm = TRUE, names = FALSE)
    )
    rng <- as.numeric(q)
    if (anyNA(rng) || !all(is.finite(rng)) || rng[1] >= rng[2]) {
      rng <- range(x, na.rm = TRUE)
    }
    if (anyNA(rng) || !all(is.finite(rng)) || rng[1] >= rng[2]) {
      return(NULL)
    }
    rng
  }

  get_reference_value <- function(x, prob = 0.25) {
    ref_val <- suppressWarnings(
      quantile(x, probs = prob, na.rm = TRUE, names = FALSE)
    )
    ref_val <- as.numeric(ref_val[1])
    if (!is.finite(ref_val)) {
      return(NA_real_)
    }
    ref_val
  }

  build_contrast_df <- function(
    fit,
    var_name,
    data,
    n_points = 100,
    x_range = NULL
  ) {
    if (is.null(x_range)) {
      x_range <- get_plot_range(data[[var_name]])
    }
    if (is.null(x_range)) {
      warning(glue("Could not determine range for contrast plot: {var_name}"))
      return(NULL)
    }
    seq_vals <- seq(x_range[1], x_range[2], length.out = n_points)
    if (length(seq_vals) <= 1 || anyNA(seq_vals)) {
      return(NULL)
    }

    ref_val <- get_reference_value(data[[var_name]])
    if (!is.finite(ref_val)) {
      warning(glue("Could not determine reference value for: {var_name}"))
      return(NULL)
    }

    k <- contrast(
      fit,
      setNames(list(seq_vals), var_name),
      setNames(list(ref_val), var_name)
    )
    as.data.frame(k[c('Contrast', 'Lower', 'Upper')]) %>%
      mutate(
        x = seq_vals,
        or = exp(Contrast),
        lower = exp(Lower),
        upper = exp(Upper)
      )
  }

  plot_contrast_df <- function(
    plot_df,
    x_label,
    x_axis_limits = NULL,
    y_axis_limits = NULL
  ) {
    or_breaks <- c(1, 2.5, 5, 10, 15, 30, 50, 100)
    p <- ggplot(
      plot_df,
      aes(x = x, y = or, ymin = lower, ymax = upper)
    ) +
      geom_ribbon(alpha = 0.2, fill = PLOT_STYLE$colors$ribbon) +
      geom_line(
        color = PLOT_STYLE$colors$fit,
        linewidth = PLOT_STYLE$line_width
      ) +
      geom_hline(
        yintercept = 1,
        linetype = 'dashed',
        color = PLOT_STYLE$colors$neutral
      ) +
      scale_y_log10(
        breaks = or_breaks,
        labels = scales::label_number(accuracy = 0.1)
      ) +
      scale_x_continuous(labels = lblr) + # Use lblr for biomarker scale
      labs(
        y = 'Odds Ratio',
        x = x_label
      ) +
      my_theme

    if (!is.null(x_axis_limits) || !is.null(y_axis_limits)) {
      p <- p +
        coord_cartesian(
          xlim = x_axis_limits,
          ylim = y_axis_limits
        )
    }
    p
  }

  # Ensure datadist is set for logit_md
  dd_logit <- datadist(logit_md)
  options(datadist = 'dd_logit')

  plot_labels <- list(
    NTproBNP = "NT-proBNP (measured)",
    BNP = "BNP (measured)",
    pred_ntprobnp_new_median = "NT-proBNP (predicted)"
  )
  plot_xlabels <- list(
    NTproBNP = "NTproBNP (pg/mL)",
    BNP = "BNP (pg/mL)",
    pred_ntprobnp_new_median = "Predicted NTproBNP (New Model, pg/mL)"
  )

  x_ranges <- list(
    NTproBNP = get_plot_range(logit_md$NTproBNP),
    BNP = get_plot_range(logit_md$BNP),
    pred_ntprobnp_new_median = get_plot_range(
      logit_md$pred_ntprobnp_new_median
    )
  )
  outcome_spline_ref_values <- c(
    NTproBNP = get_reference_value(logit_md$NTproBNP),
    BNP = get_reference_value(logit_md$BNP),
    pred_ntprobnp_new_median = get_reference_value(
      logit_md$pred_ntprobnp_new_median
    )
  )

  contrast_dfs <- list()
  outcome_spline_n <- list()

  # 1. Measured NTproBNP vs Outcome
  tryCatch(
    {
      fit_lrm_1 <- lrm(
        composite_outcome ~ rcs(NTproBNP, CONFIG$rcs_knots),
        data = logit_md,
        x = TRUE,
        y = TRUE
      )
      outcome_spline_n$NTproBNP <- fit_lrm_1$stats["n"]
      df1 <- build_contrast_df(
        fit_lrm_1,
        "NTproBNP",
        logit_md,
        x_range = x_ranges$NTproBNP
      )
      if (!is.null(df1)) {
        contrast_dfs$NTproBNP <- df1 %>%
          mutate(variable = plot_labels$NTproBNP)
      }
    },
    error = function(e) {
      warning("Failed to model NTproBNP association: ", e$message)
    }
  )

  # 2. Measured BNP vs Outcome
  tryCatch(
    {
      fit_lrm_2 <- lrm(
        composite_outcome ~ rcs(BNP, CONFIG$rcs_knots),
        data = logit_md,
        x = TRUE,
        y = TRUE
      )
      outcome_spline_n$BNP <- fit_lrm_2$stats["n"]
      df2 <- build_contrast_df(
        fit_lrm_2,
        "BNP",
        logit_md,
        x_range = x_ranges$BNP
      )
      if (!is.null(df2)) {
        contrast_dfs$BNP <- df2 %>%
          mutate(variable = plot_labels$BNP)
      }
    },
    error = function(e) warning("Failed to model BNP association: ", e$message)
  )

  # 3. Predicted (Recalibrated Median, OLS) NTproBNP vs Outcome
  tryCatch(
    {
      fit_lrm_3 <- lrm(
        composite_outcome ~ rcs(
          pred_ntprobnp_new_median,
          CONFIG$rcs_knots
        ),
        data = logit_md,
        x = TRUE,
        y = TRUE
      )
      outcome_spline_n$pred_ntprobnp_new_median <- fit_lrm_3$stats["n"]
      df3 <- build_contrast_df(
        fit_lrm_3,
        "pred_ntprobnp_new_median",
        logit_md,
        x_range = x_ranges$pred_ntprobnp_new_median
      )
      if (!is.null(df3)) {
        contrast_dfs$pred_ntprobnp_new_median <- df3 %>%
          mutate(variable = plot_labels$pred_ntprobnp_new_median)
      }
    },
    error = function(e) {
      warning("Failed to model predicted NTproBNP association: ", e$message)
    }
  )

  if (length(contrast_dfs) >= 2) {
    combined_df <- bind_rows(contrast_dfs)
    combined_colors <- setNames(
      c("#0072B2", "#D55E00", "#009E73"),
      unlist(plot_labels[c("NTproBNP", "BNP", "pred_ntprobnp_new_median")])
    )
    combined_colors <- combined_colors[
      names(combined_colors) %in% combined_df$variable
    ]
    combined_df <- combined_df %>%
      mutate(
        variable = factor(
          variable,
          levels = names(combined_colors)
        )
      )

    ntprobnp_panels <- unlist(
      plot_labels[c("NTproBNP", "pred_ntprobnp_new_median")]
    )
    ntprobnp_shared_max <- max(
      combined_df$x[combined_df$variable %in% ntprobnp_panels],
      na.rm = TRUE
    )
    panel_axis_data <- combined_df %>%
      group_by(variable) %>%
      dplyr::summarize(x = max(x, na.rm = TRUE), .groups = "drop") %>%
      mutate(
        x = if_else(
          variable %in% ntprobnp_panels,
          ntprobnp_shared_max,
          x
        )
      ) %>%
      bind_rows(
        tibble(
          variable = levels(combined_df$variable),
          x = 0
        )
      ) %>%
      mutate(
        variable = factor(
          variable,
          levels = levels(combined_df$variable)
        ),
        y = 1
      )

    y_axis_breaks <- c(1, 5, 10, 20, 40, 100, 200, 400, 1000)
    max_upper_ci <- max(combined_df$upper, na.rm = TRUE)
    upper_break_candidates <- y_axis_breaks[y_axis_breaks > max_upper_ci]
    y_axis_upper <- if (length(upper_break_candidates) > 0) {
      min(upper_break_candidates)
    } else {
      10^ceiling(log10(max_upper_ci))
    }

    p_combined <- ggplot(
      combined_df,
      aes(
        x = x,
        y = or,
        ymin = lower,
        ymax = upper,
        color = variable,
        fill = variable,
        group = variable
      )
    ) +
      geom_blank(
        data = panel_axis_data,
        aes(x = x, y = y),
        inherit.aes = FALSE
      ) +
      geom_ribbon(alpha = 0.15, color = NA, show.legend = FALSE) +
      geom_line(linewidth = PLOT_STYLE$line_width, show.legend = FALSE) +
      geom_hline(
        yintercept = 1,
        linetype = "dashed",
        color = PLOT_STYLE$colors$neutral
      ) +
      scale_y_log10(
        breaks = y_axis_breaks[y_axis_breaks <= y_axis_upper],
        labels = scales::label_number(accuracy = 1)
      ) +
      scale_x_continuous(
        labels = scales::label_number(accuracy = 1, big.mark = ","),
        expand = expansion(mult = c(0.03, 0.06))
      ) +
      scale_color_manual(values = combined_colors) +
      scale_fill_manual(values = combined_colors) +
      facet_wrap(~ variable, nrow = 1, scales = "free_x") +
      labs(
        x = "Biomarker concentration (pg/mL)",
        y = "Odds Ratio"
      ) +
      my_theme +
      coord_cartesian(ylim = c(1, y_axis_upper)) +
      theme(
        strip.text = element_text(size = rel(0.95)),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        panel.spacing = unit(1.5, "lines")
      )

    ggsave(
      output_path('outcome_assoc_combined.png'),
      p_combined,
      h = 4.6,
      w = 8,
      dpi = 300
    )
  }

  outcome_spline_n <- unlist(outcome_spline_n)
} # End of check for sufficient data for logistic regression

print("--- Finished Section 5 ---")
print("--- Analysis script finished successfully ---")
