#===============================================================================
# Section 0: Load Libraries, Configuration, and Utility Functions
#===============================================================================
# rm(list = ls()) # Consider restarting R session instead for a truly clean start
# --- Configuration ---
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
  bootstrap_reps = 80, # Set to 1000 for protocol-level bootstrapping
  bootstrap_cores = 8, # Multicore bootstrapping
  kasahara_crcl_knots = c(56.5, 72.4, 93.7), # Knots from Kasahara paper
  log_bnp_interval = log(c(1, 50000)), # Interval for root finding
  target_ntprobnp_thresholds = c(100, 200, 1500) # Thresholds for outcome analysis
)

# --- Output Directory ---
CONFIG$output_dir <- normalizePath(CONFIG$output_dir, mustWork = FALSE)
dir.create(CONFIG$output_dir, showWarnings = FALSE, recursive = TRUE)
print(paste("Output directory set to:", CONFIG$output_dir))

# --- Set Seed ---
set.seed(CONFIG$seed)
print(paste("Random seed set to:", CONFIG$seed))

# --- Load Libraries ---
required_packages <- c(
  "tidyverse",
  "nephro",
  "mcr",
  "santoku",
  "rms",
  "yardstick",
  "rsample",
  "marginaleffects",
  "modelsummary",
  "logistf",
  "flextable",
  "gtsummary",
  "gt",
  "glue",
  "scales",
  "Hmisc",
  "rootSolve",
  "table1"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", "),
    ". Please install them before running this script."
  )
}
message(
  "Loading required packages: ",
  paste(required_packages, collapse = ", ")
)
invisible(lapply(required_packages, library, character.only = TRUE))

# --- Output Helper ---
output_path <- function(...) {
  file.path(CONFIG$output_dir, ...)
}

get_bootstrap_cores <- function() {
  requested <- as.integer(CONFIG$bootstrap_cores)
  if (!is.finite(requested) || requested < 1) {
    requested <- 1L
  }
  available <- tryCatch(
    parallel::detectCores(logical = TRUE),
    error = function(e) NA_integer_
  )
  if (!is.finite(available) || available < 1) {
    available <- requested
  }
  cores <- max(1L, min(requested, available))
  if (.Platform$OS.type == "windows" && cores > 1) {
    warning("Multicore fork is not supported on Windows; using 1 core.")
    cores <- 1L
  }
  cores
}

check_required_columns <- function(df, required_cols, context = "data") {
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in ",
      context,
      ": ",
      paste(missing_cols, collapse = ", ")
    )
  }
}

read_input_csv <- function(path, label = path) {
  if (!file.exists(path)) {
    stop("Input file not found: ", label, " (", path, ")")
  }
  read_csv(path, show_col_types = FALSE)
}
# --- Utility Functions ---
# Metric Calculation Setup
class_metrics_set <- metric_set(accuracy, sens, spec, ppv, npv, kap) # Added kap here
continuous_metrics_set <- metric_set(rmse, rsq, mae)
# Plotting Helpers
lblr <- function(x) glue::glue('{scales::comma(x)} pg/mL')
PLOT_STYLE <- list(
  base_size = 12,
  point = list(
    alpha = 0.55,
    size = 1.6,
    shape = 21,
    fill = "grey80",
    color = "grey30",
    stroke = 0.3
  ),
  colors = list(
    identity = "#D55E00",
    fit = "#0072B2",
    mean = "#0072B2",
    limits = "#D55E00",
    ribbon = "#9ECAE1",
    neutral = "grey50"
  ),
  fills = list(
    relative = "#C6DBEF",
    absolute = "#C7E9C0"
  ),
  text = list(
    axis = "grey20",
    title = "grey10",
    subtitle = "grey30"
  ),
  grid = "grey92",
  strip_bg = "grey90",
  line_width = 0.9,
  annotation_size = 3.4
)

scatter_geom <- function(fill = PLOT_STYLE$point$fill) {
  geom_point(
    alpha = PLOT_STYLE$point$alpha,
    shape = PLOT_STYLE$point$shape,
    fill = fill,
    color = PLOT_STYLE$point$color,
    size = PLOT_STYLE$point$size,
    stroke = PLOT_STYLE$point$stroke
  )
}

my_theme <- theme_classic(base_size = PLOT_STYLE$base_size) %+replace%
  theme(
    aspect.ratio = 1,
    panel.grid.major = element_line(
      colour = PLOT_STYLE$grid,
      linewidth = 0.4
    ),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold", color = PLOT_STYLE$text$axis),
    axis.text = element_text(color = PLOT_STYLE$text$axis),
    axis.line = element_line(color = PLOT_STYLE$text$axis, linewidth = 0.4),
    plot.title = element_text(
      face = "bold",
      size = rel(1.1),
      color = PLOT_STYLE$text$title
    ),
    plot.subtitle = element_text(
      size = rel(0.9),
      color = PLOT_STYLE$text$subtitle
    ),
    plot.margin = margin(10, 12, 10, 10),
    strip.background = element_rect(
      fill = PLOT_STYLE$strip_bg,
      colour = NA
    ),
    strip.text = element_text(face = "bold", color = PLOT_STYLE$text$title),
    legend.position = "bottom",
    legend.title = element_blank()
  )
# Mode Function
Mode <- function(x, na.rm = TRUE) {
  if (na.rm) {
    x <- x[!is.na(x)]
  }
  if (length(x) == 0) {
    return(NA)
  }
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
# Helper to create categorical intervals with protocol-accurate boundaries
make_interval_cat <- function(x, breaks, labels) {
  factor(
    cut(
      x,
      breaks = c(-Inf, breaks, Inf),
      right = FALSE,
      include.lowest = TRUE,
      labels = labels
    ),
    levels = labels,
    ordered = TRUE
  )
}
# Protocol-aligned NT-proBNP categories: <100, 100-<200, 200-<1500, >=1500
make_ntprobnp_cat <- function(x, breaks = CONFIG$class_breaks) {
  labels <- c(
    glue("<{breaks[1]}"),
    glue("{breaks[1]}-<{breaks[2]}"),
    glue("{breaks[2]}-<{breaks[3]}"),
    glue(">={breaks[3]}")
  )
  make_interval_cat(x, breaks, labels)
}
# Protocol-aligned binary threshold: <threshold vs >=threshold
make_binary_cat <- function(x, threshold) {
  labels <- c(glue("<{threshold}"), glue(">={threshold}"))
  factor(
    cut(
      x,
      breaks = c(-Inf, threshold, Inf),
      right = FALSE,
      include.lowest = TRUE,
      labels = labels
    ),
    levels = labels
  )
}
# Weighted (quadratic) kappa for ordered categories
weighted_kappa_quadratic <- function(truth, estimate) {
  idx <- complete.cases(truth, estimate)
  truth <- truth[idx]
  estimate <- estimate[idx]
  if (length(truth) < 2) {
    return(NA_real_)
  }
  truth <- factor(truth, levels = levels(truth), ordered = TRUE)
  estimate <- factor(estimate, levels = levels(truth), ordered = TRUE)
  k <- length(levels(truth))
  if (k < 2) {
    return(NA_real_)
  }
  m <- table(truth, estimate)
  w <- outer(seq_len(k), seq_len(k), function(i, j) 1 - ((i - j)^2 / (k - 1)^2))
  po <- sum(w * m) / sum(m)
  rowp <- rowSums(m) / sum(m)
  colp <- colSums(m) / sum(m)
  pe <- sum(w * (rowp %o% colp))
  if (isTRUE(all.equal(1, pe))) {
    return(NA_real_)
  }
  (po - pe) / (1 - pe)
}
# Actual vs predicted scatter plot
make_actual_vs_pred_plot <- function(
  data,
  actual_col,
  pred_col,
  model_name = pred_col,
  log_scale = TRUE
) {
  actual_sym <- sym(actual_col)
  pred_sym <- sym(pred_col)
  plot_data <- data %>%
    select(actual = {{ actual_sym }}, predicted = {{ pred_sym }}) %>%
    filter(is.finite(actual) & is.finite(predicted))
  if (log_scale) {
    plot_data <- plot_data %>% filter(actual > 0 & predicted > 0)
  }
  if (nrow(plot_data) < 10) {
    warning(glue(
      "Insufficient data (<10 points) for actual vs predicted plot: {model_name}"
    ))
    return(
      ggplot() +
        labs(
          title = glue("Actual vs Predicted: {model_name}"),
          subtitle = "Insufficient Data"
        ) +
        my_theme
    )
  }
  p <- ggplot(plot_data, aes(x = predicted, y = actual)) +
    scatter_geom() +
    geom_abline(
      slope = 1,
      intercept = 0,
      color = PLOT_STYLE$colors$identity,
      linetype = "dashed",
      linewidth = PLOT_STYLE$line_width
    ) +
    labs(
      title = glue("Actual vs Predicted: {model_name}"),
      x = glue("Predicted {pred_col} (pg/mL)"),
      y = glue("Actual {actual_col} (pg/mL)")
    ) +
    my_theme
  if (log_scale) {
    p <- p + scale_x_log10(labels = lblr) + scale_y_log10(labels = lblr)
  } else {
    p <- p +
      scale_x_continuous(labels = lblr) +
      scale_y_continuous(labels = lblr)
  }
  p
}
# Helper to detect valid mcr results across return classes
is_mcr_result <- function(x) {
  isS4(x) && "para" %in% methods::slotNames(x)
}
# Master Metric Function (evaluates MEDIAN predictions by default)
# NOTE: For MEAN prediction evaluation, pass the mean prediction column as `predicted_col`
compute_metrics <- function(
  data,
  target_col,
  predicted_col,
  binary_threshold = CONFIG$binary_threshold,
  class_breaks = CONFIG$class_breaks
) {
  target_sym <- sym(target_col)
  predicted_sym <- sym(predicted_col)

  # Basic checks
  if (!all(c(target_col, predicted_col) %in% names(data))) {
    stop("Target or predicted column not found in data for compute_metrics")
  }
  d_filtered <- data %>%
    filter(is.finite({{ target_sym }}) & is.finite({{ predicted_sym }}))
  if (nrow(d_filtered) < 2) {
    warning("Less than 2 non-missing pairs for compute_metrics")
    return(tibble(
      .metric = character(),
      .estimator = character(),
      .estimate = numeric()
    ))
  }

  # Define classification breaks and levels/labels (protocol-aligned)
  class_labels <- c(
    glue("<{class_breaks[1]}"),
    glue("{class_breaks[1]}-<{class_breaks[2]}"),
    glue("{class_breaks[2]}-<{class_breaks[3]}"),
    glue(">={class_breaks[3]}")
  )

  d_processed <- d_filtered %>%
    mutate(
      binary_tgt = make_binary_cat({{ target_sym }}, binary_threshold),
      binary_pred = make_binary_cat({{ predicted_sym }}, binary_threshold),
      class_tgt = make_interval_cat(
        {{ target_sym }},
        class_breaks,
        class_labels
      ),
      class_pred = make_interval_cat(
        {{ predicted_sym }},
        class_breaks,
        class_labels
      )
    )
  # Ensure event level (>= threshold) is first for binary metrics
  binary_event_levels <- c(
    glue(">={binary_threshold}"),
    glue("<{binary_threshold}")
  )
  d_processed <- d_processed %>%
    mutate(
      binary_tgt = factor(binary_tgt, levels = binary_event_levels),
      binary_pred = factor(binary_pred, levels = binary_event_levels)
    )

  # Calculate metrics
  cont_m <- continuous_metrics_set(
    d_processed,
    truth = {{ target_sym }},
    estimate = {{ predicted_sym }}
  )
  pearson_r <- suppressWarnings(cor(
    d_processed[[target_col]],
    d_processed[[predicted_col]],
    use = "complete.obs",
    method = "pearson"
  ))
  cont_m <- bind_rows(
    cont_m,
    tibble(.metric = "pearson_r", .estimator = "pearson", .estimate = pearson_r)
  )

  bin_m <- class_metrics_set(
    d_processed,
    truth = binary_tgt,
    estimate = binary_pred
  )
  multi_m <- class_metrics_set(
    d_processed,
    truth = class_tgt,
    estimate = class_pred
  )
  w_kappa <- weighted_kappa_quadratic(
    d_processed$class_tgt,
    d_processed$class_pred
  )
  w_kappa_m <- tibble(
    .metric = "weighted_kappa_quadratic",
    .estimator = "weighted",
    .estimate = w_kappa
  )

  # OLS Calibration slope/intercept (on original scale)
  cal_fit <- lm(
    as.formula(glue("{target_col} ~ {predicted_col}")),
    data = d_processed
  )
  cal_m <- tibble(
    .metric = c('calibration_intercept', 'calibration_slope'),
    .estimator = 'OLS',
    .estimate = coef(cal_fit)
  )

  # Passing-Bablok Regression
  paba_model <- tryCatch(
    {
      mcreg(
        x = d_processed[[target_col]],
        y = d_processed[[predicted_col]],
        method.reg = "PaBa",
        na.rm = TRUE
      )
    },
    error = function(e) NULL
  )

  if (!is.null(paba_model) && is_mcr_result(paba_model)) {
    paba_results_para <- paba_model@para
    paba_m <- tibble(
      .metric = c(
        "paba_intercept",
        "paba_intercept_lower",
        "paba_intercept_upper",
        "paba_slope",
        "paba_slope_lower",
        "paba_slope_upper"
      ),
      .estimator = "PaBa",
      .estimate = c(
        paba_results_para[1, "EST"],
        paba_results_para[1, "LCI"],
        paba_results_para[1, "UCI"],
        paba_results_para[2, "EST"],
        paba_results_para[2, "LCI"],
        paba_results_para[2, "UCI"]
      )
    )
  } else {
    paba_m <- tibble(
      .metric = character(),
      .estimator = character(),
      .estimate = numeric()
    )
    warning("Passing-Bablok regression failed or insufficient data.")
  }

  bind_rows(cont_m, bin_m, multi_m, w_kappa_m, cal_m, paba_m) %>%
    mutate(model_type = predicted_col) # Add identifier based on prediction column name
}
# Bland-Altman Plot Function
make_bland_altman_plot <- function(
  data,
  actual_col,
  pred_col,
  model_name = pred_col,
  type = c("absolute", "relative")
) {
  type <- match.arg(type)
  actual_sym <- sym(actual_col)
  pred_sym <- sym(pred_col)
  plot_df <- data %>%
    filter(is.finite({{ actual_sym }}) & is.finite({{ pred_sym }})) %>%
    mutate(
      av = ({{ actual_sym }} + {{ pred_sym }}) / 2,
      diff = {{ pred_sym }} - {{ actual_sym }},
      plotted_err = if (type == "relative") {
        # Avoid division by zero or near-zero average for relative error
        ifelse(abs(av) > 1e-6, diff / av, NA)
      } else {
        diff
      }
    ) %>%
    filter(!is.na(av) & !is.na(plotted_err))
  if (nrow(plot_df) < 10) {
    # Increased minimum points for stable stats
    warning(glue(
      "Insufficient data (<10 points) for Bland-Altman plot: {model_name} ({type})"
    ))
    return(
      ggplot() +
        labs(
          title = glue("Bland-Altman Plot: {model_name} ({type})"),
          subtitle = "Insufficient Data"
        ) +
        my_theme
    )
  }
  sd_err <- sd(plot_df$plotted_err, na.rm = TRUE)
  mu_err <- mean(plot_df$plotted_err, na.rm = TRUE)
  upper_limit <- mu_err + 1.96 * sd_err
  lower_limit <- mu_err - 1.96 * sd_err
  y_label <- if (type == "relative") {
    "Relative Error (%)"
  } else {
    "Absolute Error (pg/mL)"
  }
  fill_color <- if (type == "relative") {
    PLOT_STYLE$fills$relative
  } else {
    PLOT_STYLE$fills$absolute
  }
  anno_suffix <- if (type == "relative") "%" else ""
  anno_multiplier <- if (type == "relative") 100 else 1
  p <- ggplot(plot_df, aes(x = av, y = plotted_err)) +
    scatter_geom(fill = fill_color) +
    geom_hline(
      yintercept = mu_err,
      color = PLOT_STYLE$colors$mean,
      linewidth = PLOT_STYLE$line_width
    ) +
    geom_hline(
      yintercept = c(lower_limit, upper_limit),
      linetype = "dashed",
      color = PLOT_STYLE$colors$limits,
      linewidth = PLOT_STYLE$line_width
    ) +
    annotate(
      "text",
      x = Inf,
      y = mu_err,
      label = sprintf("Mean = %.2f%s", mu_err * anno_multiplier, anno_suffix),
      hjust = 1.05,
      vjust = -0.5,
      color = PLOT_STYLE$colors$mean,
      size = PLOT_STYLE$annotation_size
    ) +
    annotate(
      "text",
      x = Inf,
      y = upper_limit,
      label = sprintf(
        "+1.96 SD = %.2f%s",
        upper_limit * anno_multiplier,
        anno_suffix
      ),
      hjust = 1.05,
      vjust = -0.5,
      color = PLOT_STYLE$colors$limits,
      size = PLOT_STYLE$annotation_size
    ) +
    annotate(
      "text",
      x = Inf,
      y = lower_limit,
      label = sprintf(
        "-1.96 SD = %.2f%s",
        lower_limit * anno_multiplier,
        anno_suffix
      ),
      hjust = 1.05,
      vjust = 1.5,
      color = PLOT_STYLE$colors$limits,
      size = PLOT_STYLE$annotation_size
    ) +
    scale_x_continuous(labels = lblr) +
    (if (type == "relative") {
      scale_y_continuous(labels = scales::label_percent(accuracy = 1))
    } else {
      scale_y_continuous(labels = lblr)
    }) +
    labs(
      x = "Mean of Actual and Predicted (pg/mL)",
      y = y_label,
      title = glue("Bland-Altman Plot: {model_name}"),
      subtitle = glue("{str_to_title(type)} Error")
    ) +
    my_theme
  return(p)
}
# Prediction Function (Kasahara Formula - OLD MODEL)
# Source: Kasahara et al., Int J Cardiol 2019; 280:184-189 (DOI: 10.1016/j.ijcard.2018.12.069)
# Predicts MEDIAN NT-proBNP
predict_ntprobnp_kasahara <- function(BNP, Age, BMI, Hb_gdl, CrCl, Sex_M, AF) {
  # Ensure numeric inputs
  BNP <- as.numeric(BNP)
  Age <- as.numeric(Age)
  BMI <- as.numeric(BMI)
  Hb_gdl <- as.numeric(Hb_gdl) # Expecting g/dL (converted in clean_data)
  CrCl <- as.numeric(CrCl) # Expecting mL/min (calculated in clean_data)
  Sex_M <- as.numeric(Sex_M) # Expecting 1 for Male, 0 for Female
  AF <- as.numeric(AF) # Expecting 1 for AF, 0 otherwise

  log10_BNP <- ifelse(BNP > 0 & is.finite(BNP), log10(BNP), NA_real_)

  # Spline knots from paper
  k1 <- CONFIG$kasahara_crcl_knots[1]
  k2 <- CONFIG$kasahara_crcl_knots[2]
  k3 <- CONFIG$kasahara_crcl_knots[3]

  # Handle potential NAs in inputs before calculation
  valid_input <- is.finite(log10_BNP) &
    is.finite(Age) &
    is.finite(BMI) &
    is.finite(Hb_gdl) &
    is.finite(CrCl) &
    is.finite(Sex_M) &
    is.finite(AF)

  pred_log10 <- ifelse(
    valid_input,
    2.05 +
      0.907 * log10_BNP -
      0.00522 * Age +
      0.00283 * BMI -
      0.00866 * Hb_gdl +
      (-0.0422 *
        CrCl +
        0.000530 * CrCl^2 -
        0.00000214 * CrCl^3 -
        0.00000278 * pmax(0, CrCl - k1)^3 +
        0.00000621 * pmax(0, CrCl - k2)^3 -
        0.00000133 * pmax(0, CrCl - k3)^3) +
      0.0164 * (1 - Sex_M) + # Adds term if Female (Sex_M=0)
      0.194 * AF, # Adds term if AF=1
    NA_real_
  )

  result <- 10^pred_log10
  result <- ifelse(is.finite(result), result, NA_real_)
  return(result)
}
# New Model Fitting Function (OLS on log scale)
fit_new_model_ols <- function(data, knots = CONFIG$rcs_knots) {
  required_cols <- c(
    "log_ntprobnp",
    "Age",
    "CrCl",
    "log_bnp",
    "Sex_M",
    "BMI",
    "AF",
    "Hb_gdl"
  )
  if (!all(required_cols %in% names(data))) {
    missing_cols <- setdiff(required_cols, names(data))
    stop(
      "Missing required columns for OLS model: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  data_model <- data[, required_cols, drop = FALSE]
  if (nrow(na.omit(data_model)) < 20) {
    # Check for sufficient data
    warning("Less than 20 complete cases for OLS model fitting.")
    return(NULL)
  }
  # Define formula using rcs() and specified knots
  # Interactions between key predictors can be important
  frml_ols <- formula(glue(
    "log_ntprobnp ~ rcs(Age, {knots}) * rcs(CrCl, {knots}) * rcs(log_bnp, {knots}) +
                    Sex_M + rcs(BMI, {knots}) + AF + rcs(Hb_gdl, {knots})"
  ))
  # Fit OLS model using rms::ols
  # Need datadist for rms functions to work properly later (like contrast)
  dd <- tryCatch(
    {
      datadist(data_model)
    },
    error = function(e) {
      warning("datadist failed for OLS model fitting: ", e$message)
      return(NULL)
    }
  )
  if (is.null(dd)) {
    return(NULL)
  }
  old_datadist <- options("datadist")
  options(datadist = 'dd')
  on.exit(options(old_datadist), add = TRUE)
  fit <- tryCatch(
    {
      ols(
        formula = frml_ols,
        data = data_model,
        x = TRUE,
        y = TRUE,
        na.action = na.omit
      )
    },
    error = function(e) {
      warning("OLS model fitting failed: ", e$message)
      return(NULL)
    }
  )
  return(fit)
}
# Root Finding Function for Thresholds (using New OLS Model)
get_BNP_for_NT_new_model <- function(
  target_NT,
  covars_grid,
  fitted_model,
  log_recal_params = NULL,
  use_mean = FALSE
) {
  if (is.null(fitted_model) || !inherits(fitted_model, "ols")) {
    warning("Invalid fitted_model provided to get_BNP_for_NT_new_model.")
    return(NA_real_)
  }
  sigma_sq_half <- (fitted_model$stats['Sigma']^2) / 2
  if (!is.finite(sigma_sq_half)) {
    sigma_sq_half <- 0
  } # Handle potential NA/Inf sigma
  objective_function <- function(log_bnp_val) {
    temp_data <- covars_grid %>% mutate(log_bnp = log_bnp_val)
    # Ensure newdata has same factor levels if Sex_M/AF are factors in covars_grid
    # temp_data$Sex_M <- factor(temp_data$Sex_M, levels = levels(fitted_model$Design$limits$Sex_M))
    # temp_data$AF <- factor(temp_data$AF, levels = levels(fitted_model$Design$limits$AF))
    predicted_log_NT_initial <- predict(fitted_model, newdata = temp_data)
    # Apply recalibration if parameters provided
    if (!is.null(log_recal_params) && length(log_recal_params) == 2) {
      I_log <- log_recal_params[1]
      S_log <- log_recal_params[2]
      if (abs(S_log) < 1e-6) {
        S_log <- sign(S_log) * 1e-6
      } # Avoid division by zero
      predicted_log_NT_final <- (predicted_log_NT_initial - I_log) / S_log
    } else {
      predicted_log_NT_final <- predicted_log_NT_initial
    }
    # Calculate NT-proBNP on original scale (median by default; mean if use_mean=TRUE)
    expected_NT_pred <- if (isTRUE(use_mean)) {
      exp(predicted_log_NT_final + sigma_sq_half)
    } else {
      exp(predicted_log_NT_final)
    }
    return(expected_NT_pred - target_NT)
  }
  # Find root (log_bnp value) where objective function is zero
  root_result <- tryCatch(
    {
      uniroot(
        objective_function,
        interval = CONFIG$log_bnp_interval,
        extendInt = "yes",
        tol = 1e-5
      )
    },
    error = function(e) {
      warning(glue(
        "Could not find root for target_NT = {target_NT}: {e$message}"
      ))
      return(NULL)
    }
  )
  if (!is.null(root_result)) {
    return(exp(root_result$root)) # Return BNP on original scale
  } else {
    return(NA_real_)
  }
}
# Calibration Plot Function
make_calibration_plot <- function(
  data,
  actual_col,
  pred_col,
  model_name = pred_col,
  log_scale = FALSE
) {
  actual_sym <- sym(actual_col)
  pred_sym <- sym(pred_col)
  plot_data <- data %>%
    select(actual = {{ actual_sym }}, predicted = {{ pred_sym }}) %>%
    filter(
      is.finite(actual) & is.finite(predicted) & actual > 0 & predicted > 0
    ) # Ensure positive for log scale
  if (nrow(plot_data) < 10) {
    warning(glue(
      "Insufficient data (<10 points) for calibration plot: {model_name}"
    ))
    return(
      ggplot() +
        labs(
          title = glue("Calibration Plot: {model_name}"),
          subtitle = "Insufficient Data"
        ) +
        my_theme
    )
  }
  # Determine plot limits
  if (log_scale) {
    lims <- range(c(plot_data$actual, plot_data$predicted), na.rm = TRUE)
    lims <- c(max(lims[1] * 0.8, 1), lims[2] * 1.2) # Ensure lower limit >= 1 for log scale, add buffer
  } else {
    lims <- range(c(plot_data$actual, plot_data$predicted), na.rm = TRUE)
    buffer <- (lims[2] - lims[1]) * 0.05
    lims <- c(lims[1] - buffer, lims[2] + buffer)
  }
  # Fit linear model for calibration assessment
  cal_fit <- lm(actual ~ predicted, data = plot_data)
  cal_slope <- coef(cal_fit)[2]
  cal_intercept <- coef(cal_fit)[1]
  cal_summary <- glue(
    "OLS Fit: Actual = {sprintf('%.2f', cal_intercept)} + {sprintf('%.2f', cal_slope)} * Predicted"
  )
  p <- ggplot(plot_data, aes(x = predicted, y = actual)) + # Swapped axes for typical calibration plot
    scatter_geom() +
    geom_abline(
      slope = 1,
      intercept = 0,
      color = PLOT_STYLE$colors$identity,
      linetype = "dashed",
      linewidth = PLOT_STYLE$line_width
    ) + # Line of identity
    geom_smooth(
      method = "lm",
      se = TRUE,
      color = PLOT_STYLE$colors$fit,
      fill = PLOT_STYLE$colors$ribbon,
      formula = y ~ x,
      linewidth = PLOT_STYLE$line_width
    ) + # Calibration line
    labs(
      title = glue("Calibration Plot: {model_name}"),
      subtitle = cal_summary, # Add OLS fit info
      x = glue("Predicted {pred_col} (pg/mL)"),
      y = glue("Actual {actual_col} (pg/mL)")
    ) +
    my_theme +
    coord_cartesian(xlim = lims, ylim = lims) # Ensure limits apply after potential expansion by geom_smooth
  if (log_scale) {
    p <- p + scale_x_log10(labels = lblr) + scale_y_log10(labels = lblr)
  } else {
    p <- p +
      scale_x_continuous(labels = lblr) +
      scale_y_continuous(labels = lblr)
  }
  return(p)
}
# Log-Log Calibration Plot Function (for log10 models)
make_loglog_calibration_plot <- function(
  data,
  actual_col,
  pred_col,
  model_name = pred_col
) {
  actual_sym <- sym(actual_col)
  pred_sym <- sym(pred_col)
  plot_data <- data %>%
    select(actual = {{ actual_sym }}, predicted = {{ pred_sym }}) %>%
    filter(
      is.finite(actual) & is.finite(predicted) & actual > 0 & predicted > 0
    )
  if (nrow(plot_data) < 10) {
    warning(glue(
      "Insufficient data (<10 points) for log-log calibration plot: {model_name}"
    ))
    return(
      ggplot() +
        labs(
          title = glue("Log-Log Calibration Plot: {model_name}"),
          subtitle = "Insufficient Data"
        ) +
        my_theme
    )
  }
  # Fit log10 calibration
  log_fit <- lm(log10(actual) ~ log10(predicted), data = plot_data)
  log_intercept <- coef(log_fit)[1]
  log_slope <- coef(log_fit)[2]
  plot_data <- plot_data %>%
    mutate(fit = 10^(log_intercept + log_slope * log10(predicted)))
  lims <- range(c(plot_data$actual, plot_data$predicted), na.rm = TRUE)
  lims <- c(max(lims[1] * 0.8, 1), lims[2] * 1.2)
  cal_summary <- glue(
    "Log10 Fit: log10(Actual) = {sprintf('%.2f', log_intercept)} + {sprintf('%.2f', log_slope)} * log10(Predicted)"
  )
  p <- ggplot(plot_data, aes(x = predicted, y = actual)) +
    scatter_geom() +
    geom_abline(
      slope = 1,
      intercept = 0,
      color = PLOT_STYLE$colors$identity,
      linetype = "dashed",
      linewidth = PLOT_STYLE$line_width
    ) +
    geom_line(
      aes(y = fit),
      color = PLOT_STYLE$colors$fit,
      linewidth = PLOT_STYLE$line_width
    ) +
    scale_x_log10(labels = lblr) +
    scale_y_log10(labels = lblr) +
    labs(
      title = glue("Log-Log Calibration Plot: {model_name}"),
      subtitle = cal_summary,
      x = glue("Predicted {pred_col} (pg/mL)"),
      y = glue("Actual {actual_col} (pg/mL)")
    ) +
    my_theme +
    coord_cartesian(xlim = lims, ylim = lims)
  return(p)
}
# Helper to run mcreg safely
run_mcreg <- function(data, x_col, y_col) {
  df_mcreg <- data %>% select(x = {{ x_col }}, y = {{ y_col }}) %>% na.omit()
  if (nrow(df_mcreg) < 10) {
    # Need sufficient points for stable PaBa
    warning(glue(
      "Skipping PaBa regression due to < 10 points for {y_col} vs {x_col}"
    ))
    return(NULL)
  }
  m <- tryCatch(
    {
      mcreg(x = df_mcreg$x, y = df_mcreg$y, method.reg = "PaBa", na.rm = TRUE)
    },
    error = function(e) {
      warning(glue(
        "PaBa regression failed for {y_col} vs {x_col}: {e$message}"
      ))
      return(NULL)
    }
  )
  if (!is.null(m) && is_mcr_result(m)) {
    return(m)
  }
  warning(glue(
    "PaBa regression returned unexpected class for {y_col} vs {x_col}."
  ))
  NULL
}
# Helper to bootstrap Passing-Bablok slope/intercept
bootstrap_paba <- function(data, x_col, y_col, reps = CONFIG$bootstrap_reps) {
  df_mcreg <- data %>% select(x = {{ x_col }}, y = {{ y_col }}) %>% na.omit()
  if (nrow(df_mcreg) < 10) {
    warning("Insufficient data (<10 points) for PaBa bootstrapping.")
    return(tibble(intercept = numeric(), slope = numeric()))
  }
  cores <- get_bootstrap_cores()
  message(glue("Running PaBa bootstrap on {cores} core(s)..."))
  boot_list <- parallel::mclapply(
    seq_len(reps),
    function(i) {
      idx <- sample.int(nrow(df_mcreg), replace = TRUE)
      d <- df_mcreg[idx, , drop = FALSE]
      m <- tryCatch(
        {
          mcreg(x = d$x, y = d$y, method.reg = "PaBa", na.rm = TRUE)
        },
        error = function(e) NULL
      )
      if (is.null(m) || !is_mcr_result(m)) {
        return(tibble(intercept = numeric(), slope = numeric()))
      }
      para <- m@para
      tibble(
        intercept = para["Intercept", "EST"],
        slope = para["Slope", "EST"]
      )
    },
    mc.cores = cores,
    mc.set.seed = TRUE
  )
  boot_results <- bind_rows(boot_list)
  if (!all(c("intercept", "slope") %in% names(boot_results))) {
    return(tibble(intercept = numeric(), slope = numeric()))
  }
  boot_results
}
# Helper to plot mcreg safely
plot_mcreg <- function(mcreg_obj, filename, title) {
  if (!is.null(mcreg_obj) && is_mcr_result(mcreg_obj)) {
    png(filename, width = 800, height = 800, res = 100)
    on.exit(dev.off(), add = TRUE)
    tryCatch(
      {
        plot(
          mcreg_obj,
          plot.type = "regression",
          ci.area = TRUE,
          main = title,
          points.col = alpha(
            PLOT_STYLE$point$color,
            PLOT_STYLE$point$alpha
          ),
          points.pch = 16,
          points.cex = 0.9
        )
      },
      error = function(e) {
        plot(1, type = "n", axes = F, xlab = "", ylab = "") # Empty plot on error
        title(main = paste("Error plotting:", title), col.main = "red")
        text(1, 1, e$message, col = "red")
      }
    )
  }
}
#===============================================================================
# Section 1: Data Loading and Preparation
#===============================================================================
print("--- Starting Section 1: Data Loading and Preparation ---")
# --- Load Data ---
# Using read_csv from readr (part of tidyverse)
# Explicit file checks for transparency
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
# --- Clean and Combine Data ---
# Goal: Create a single primary dataframe `analysis_data` for all subsequent analyses.
# We'll use the 'new' dataset as the base, assuming it's the primary cohort for outcome analysis.
# Define cleaning steps as a function for reusability if needed, but here apply directly
# --- Clean and Combine Data ---
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
      # Convert Hb from input g/L to g/dL for Kasahara formula
      Hb_gdl = as.numeric(Hb) / 10, # REINSTATED division by 10
      # Ensure input Cr is numeric (expected umol/L) and convert to mg/dL for CG
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
      # Ensure other vars are numeric/valid
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
      log_bnp = ifelse(is.finite(BNP), log(BNP), NA_real_),
      log10_bnp = ifelse(is.finite(BNP), log10(BNP), NA_real_),
      log_ntprobnp = ifelse(is.finite(NTproBNP), log(NTproBNP), NA_real_)
    ) %>%
    # Select relevant columns, drop intermediate ones if desired
    select(-any_of(c("Ht_m", "Cr_umol_input", "Cr_mgdl"))) # Remove intermediate Cr vars
}

# --- Apply the cleaning function ---
analysis_data_raw <- clean_data(d_new_raw, CONFIG$cr_conversion_factor)


# --- Calculate Old Model Predictions ---
# Apply Kasahara formula (predicts MEDIAN NTproBNP)
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
# Based on non-missing values needed for NEW model fitting and outcome analysis
required_vars_new_model <- c(
  "log_ntprobnp",
  "Age",
  "CrCl",
  "log_bnp",
  "Sex_M",
  "BMI",
  "AF",
  "Hb_gdl"
)
required_vars_biomarkers <- c("NTproBNP", "BNP") # Need actual values too
required_numeric <- unique(c(required_vars_new_model, required_vars_biomarkers))

# Filter based on required variables for new model development and outcome analysis
analysis_data <- analysis_data_raw %>%
  filter(if_all(all_of(required_numeric), is.finite)) %>% # Require predictors/biomarkers
  {
    if ("StudyID" %in% names(.)) filter(., !is.na(StudyID)) else .
  }
n_analysis <- nrow(analysis_data)
print(glue(
  "Number of subjects in the model-development cohort (complete predictors/biomarkers): {n_analysis}"
))
if (n_analysis < 50) {
  # Arbitrary minimum threshold
  stop(
    "Insufficient data in the final analysis cohort (<50 subjects). Stopping analysis."
  )
}
# Define datadist for rms functions based on the final analysis cohort
dd_final <- datadist(analysis_data)
options(datadist = 'dd_final')
# Save the final analysis dataset if needed (optional)
# write_rds(analysis_data, output_path('final_analysis_data.rds'))
# write_csv(analysis_data, output_path('final_analysis_data.csv'))
print("--- Finished Section 1 ---")

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
    # Convert relevant comorbidities to factors
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
        "OSA"
      )),
      ~ factor(., levels = c(0, 1), labels = c("No", "Yes"))
    ),
    Sex = factor(Sex, levels = c("F", "M")),
    SurgeryCategory = case_when(
      if ("Major_VascS" %in% names(.)) Major_VascS == 1 ~ "Major Vascular",
      if ("Major_ThorS" %in% names(.)) Major_ThorS == 1 ~ "Major Thoracic",
      if ("Major_GenS" %in% names(.)) Major_GenS == 1 ~ "Major General",
      if ("Major_NeurS" %in% names(.)) Major_NeurS == 1 ~ "Major Neurological",
      if ("Major_OrthS" %in% names(.)) Major_OrthS == 1 ~ "Major Orthopedic",
      if ("Major_UroGynS" %in% names(.)) Major_UroGynS == 1 ~ "Major Uro/Gyn",
      if ("Low_Surgery" %in% names(.)) Low_Surgery == 1 ~ "Low Risk Surgery",
      if ("No_surgery" %in% names(.)) No_surgery == 1 ~ "No Surgery",
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
      "MI",
      "PVD",
      "stroke",
      "tia",
      "copd",
      "cancer",
      "chf",
      "DM",
      "Insulin",
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
  MI = "Myocardial Infarction History",
  PVD = "Peripheral Vascular Disease",
  stroke = "Stroke History",
  tia = "TIA History",
  copd = "COPD",
  cancer = "Cancer History",
  chf = "Congestive Heart Failure History",
  DM = "Diabetes Mellitus",
  Insulin = "Insulin Use",
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
missing_counts <- sapply(data_for_table1, function(x) sum(is.na(x)))
missing_counts_text <- names(missing_counts[missing_counts > 0]) %>%
  {
    paste0(., " (", missing_counts[missing_counts > 0], ")")
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
  "Median (IQR) for most continuous; Mean (SD) & Median (IQR) for BNP & NTproBNP. Statistics calculated on non-missing data."
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
          all_of(intersect(c("BNP", "NTproBNP"), names(data_for_table1))) ~ c(
            "{mean} ({sd})",
            "{median} ({p25}, {p75})"
          ),
          all_categorical() ~ "{n} ({p}%)"
        ),
        digits = list(
          all_continuous() ~ 1,
          all_of(intersect(c("BNP", "NTproBNP"), names(data_for_table1))) ~ c(
            0,
            0
          ),
          all_categorical() ~ c(0, 1)
        ),
        type = list(
          all_continuous() ~ "continuous",
          all_of(intersect(
            c("BNP", "NTproBNP"),
            names(data_for_table1)
          )) ~ "continuous2"
        ),
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
        "**Table 1: Baseline Characteristics of Initial Cohort**"
      ) %>%
      bold_labels()

    print(baseline_table_final)

    # Save table as DOCX using flextable
    baseline_flextable <- as_flex_table(baseline_table_final)
    save_as_docx(
      baseline_flextable,
      path = output_path(
        "baseline_characteristics_gtsummary_initial_cohort.docx"
      )
    )
    print("Baseline Table 1 (Initial Cohort) saved successfully.")
  },
  error = function(e) {
    print(paste("Error generating Baseline Table 1:", e$message))
    print(traceback())
  }
)

print("--- Finished Section 1.5 ---")

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
  write_csv(
    old_model_errors,
    output_path('kasahara_model_prediction_errors.csv')
  )
  print("Kasahara (Old) model performance metrics (from compute_metrics):")
  print(old_model_errors)

  # Assign the result from the compute_metrics PaBa attempt for plotting
  # (This uses the run_mcreg helper with tryCatch)
  paba_old_model_obj <- run_mcreg(
    data_eval_kasahara,
    x_col = NTproBNP,
    y_col = pred_ntprobnp_kasahara
  )

  # --- Log-Log Calibration (Kasahara): OLS + Passing-Bablok ---
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
    write_csv(
      kasahara_loglog_summary,
      output_path("kasahara_loglog_calibration_summary.csv")
    )
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
    write_csv(
      old_model_errors_recal_ols,
      output_path('kasahara_model_prediction_errors_recal_ols.csv')
    )
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
      write_csv(
        old_model_errors_recal_paba,
        output_path('kasahara_model_prediction_errors_recal_paba.csv')
      )
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
  print(avp_kasahara)
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
  print(cal_plot_old_log)
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
  print(cal_plot_old_lin)
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
  print(cal_plot_old_loglog)
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
    print(cal_plot_old_loglog_recal_ols)
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
    print(avp_kasahara_recal_ols)
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
    print(cal_plot_old_loglog_recal_paba)
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
    print(avp_kasahara_recal_paba)
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
  print(ba_plot_old_rel)
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
  print(ba_plot_old_abs)
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
    filename = output_path('paba_plot_kasahara.png'),
    title = "Passing-Bablok Regression (Kasahara Model)"
  )
} else {
  print("Insufficient data for Kasahara model visualizations.")
}

print("--- Finished Section 3 ---")
#===============================================================================
# Section 4: Develop, Recalibrate, and Evaluate New Model (OLS)
#===============================================================================
print("--- Starting Section 4: New Model Development & Evaluation ---")
# Use the final 'analysis_data' cohort
md <- analysis_data
# --- Fit Initial OLS Model ---
print("Fitting initial OLS model on log(NTproBNP)...")
dd <- datadist(md) # Use the analysis cohort for datadist
options(datadist = 'dd') # Tells ols to look for an object named 'dd'
fit_ols <- fit_new_model_ols(md, knots = CONFIG$rcs_knots)
if (is.null(fit_ols)) {
  stop("Initial OLS model fitting failed. Cannot proceed.")
}

# --- Compare Null vs Full Model (log scale) ---
print("Comparing null model vs full model (log scale)...")
required_cols_cmp <- c(
  "log_ntprobnp",
  "Age",
  "CrCl",
  "log_bnp",
  "Sex_M",
  "BMI",
  "AF",
  "Hb_gdl"
)
data_model_cmp <- md[, required_cols_cmp, drop = FALSE]
dd_cmp <- tryCatch(
  {
    datadist(data_model_cmp)
  },
  error = function(e) {
    warning("datadist failed for null model comparison: ", e$message)
    return(NULL)
  }
)
if (!is.null(dd_cmp)) {
  options(datadist = 'dd_cmp')
  fit_ols_null <- tryCatch(
    {
      ols(
        log_ntprobnp ~ 1,
        data = data_model_cmp,
        x = TRUE,
        y = TRUE,
        na.action = na.omit
      )
    },
    error = function(e) {
      warning("Null model fitting failed: ", e$message)
      return(NULL)
    }
  )

  null_vs_full <- NULL
  if (!is.null(fit_ols_null)) {
    null_vs_full <- tryCatch(
      {
        anova_tbl <- anova(fit_ols_null, fit_ols)
        as.data.frame(anova_tbl)
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

    write_csv(
      null_vs_full,
      output_path("null_vs_full_model_comparison.csv")
    )
    print("Null vs full model comparison saved.")
  }
}
options(datadist = 'dd') # restore

sigma_ols <- fit_ols$stats['Sigma'] # Residual SE on log scale
if (!is.finite(sigma_ols)) {
  sigma_ols <- 0
} # Handle edge case
# Add initial predictions (log scale) to the dataset
md$log_pred_initial <- predict(fit_ols) # Predicts on the data used for fitting
# --- Recalibration using Log-Log OLS and Passing-Bablok ---
print("Running Log-Log OLS recalibration...")
ols_loglog_fit <- lm(log_ntprobnp ~ log_pred_initial, data = md)
I_log_ols <- coef(ols_loglog_fit)[1]
S_log_ols <- coef(ols_loglog_fit)[2]
print(glue(
  "OLS log-log recalibration: I={round(I_log_ols, 4)}, S={round(S_log_ols, 4)}"
))

print("Running Log-Log Passing-Bablok recalibration...")
paba_loglog_model <- run_mcreg(
  md,
  x_col = log_pred_initial,
  y_col = log_ntprobnp
)
I_log_full <- 0
S_log_full <- 1
if (!is.null(paba_loglog_model)) {
  paba_loglog_results_para <- paba_loglog_model@para
  I_log_full <- paba_loglog_results_para["Intercept", "EST"]
  S_log_full <- paba_loglog_results_para["Slope", "EST"]
  print("Log-Log PaBa Results (Full Sample):")
  print(paba_loglog_results_para)
} else {
  print("Log-Log Passing-Bablok failed or insufficient data (full sample).")
}

# Bootstrap PaBa to estimate CI and recalibration parameters
print(glue(
  "--- Bootstrapping Log-Log PaBa ({CONFIG$bootstrap_reps} iterations) ---"
))
paba_loglog_boot <- bootstrap_paba(
  md,
  x_col = log_pred_initial,
  y_col = log_ntprobnp,
  reps = CONFIG$bootstrap_reps
) %>%
  filter(is.finite(intercept) & is.finite(slope))

use_bootstrap_recal <- TRUE
if (nrow(paba_loglog_boot) > 10) {
  paba_loglog_boot_summary <- tibble(
    term = c("Intercept", "Slope"),
    estimate_median = c(
      median(paba_loglog_boot$intercept, na.rm = TRUE),
      median(paba_loglog_boot$slope, na.rm = TRUE)
    ),
    estimate_mean = c(
      mean(paba_loglog_boot$intercept, na.rm = TRUE),
      mean(paba_loglog_boot$slope, na.rm = TRUE)
    ),
    ci_lower_95 = c(
      quantile(paba_loglog_boot$intercept, probs = 0.025, na.rm = TRUE),
      quantile(paba_loglog_boot$slope, probs = 0.025, na.rm = TRUE)
    ),
    ci_upper_95 = c(
      quantile(paba_loglog_boot$intercept, probs = 0.975, na.rm = TRUE),
      quantile(paba_loglog_boot$slope, probs = 0.975, na.rm = TRUE)
    ),
    n_boot = nrow(paba_loglog_boot)
  )
  print("Bootstrapped Log-Log PaBa Calibration (Median with 95% CI):")
  print(paba_loglog_boot_summary)
  write_rds(
    paba_loglog_boot_summary,
    output_path('paba_loglog_bootstrap_summary.rds')
  )
  write_csv(
    paba_loglog_boot_summary,
    output_path('paba_loglog_bootstrap_summary.csv')
  )

  if (use_bootstrap_recal) {
    I_log_paba <- paba_loglog_boot_summary$estimate_median[
      paba_loglog_boot_summary$term == "Intercept"
    ]
    S_log_paba <- paba_loglog_boot_summary$estimate_median[
      paba_loglog_boot_summary$term == "Slope"
    ]
    print(glue(
      "Using bootstrapped median for PaBa recalibration: I={round(I_log_paba, 4)}, S={round(S_log_paba, 4)}"
    ))
  } else {
    I_log_paba <- I_log_full
    S_log_paba <- S_log_full
    print(glue(
      "Using full-sample PaBa for recalibration: I={round(I_log_paba, 4)}, S={round(S_log_paba, 4)}"
    ))
  }
} else {
  print(
    "Insufficient bootstrap PaBa results; using full-sample PaBa estimates for recalibration."
  )
  I_log_paba <- I_log_full
  S_log_paba <- S_log_full
}

if (
  !is.finite(I_log_paba) || !is.finite(S_log_paba) || abs(S_log_paba) < 1e-6
) {
  warning(
    "Invalid/near-zero PaBa slope. Resetting to S=1, I=0 (no recalibration)."
  )
  I_log_paba <- 0
  S_log_paba <- 1
}
log_recal_params <- c(I_log_paba, S_log_paba) # Used for threshold derivation

# --- Calculate Recalibrated and Final Predictions (Median) ---
md <- md %>%
  mutate(
    log_pred_recal_ols = I_log_ols + S_log_ols * log_pred_initial,
    log_pred_recal_paba = I_log_paba + S_log_paba * log_pred_initial,
    pred_ntprobnp_new_median_ols = exp(log_pred_recal_ols),
    pred_ntprobnp_new_median_paba = exp(log_pred_recal_paba)
  ) %>%
  filter(
    is.finite(pred_ntprobnp_new_median_ols) &
      is.finite(pred_ntprobnp_new_median_paba)
  )

# --- Evaluate APPARENT Performance of RECALIBRATED Model (Median Predictions) ---
print(
  "--- Evaluating APPARENT Performance (Recalibrated New Model - MEDIAN Predictions) ---"
)
if (nrow(md) > 10) {
  apparent_performance_new_median_ols <- compute_metrics(
    md,
    target_col = "NTproBNP",
    predicted_col = "pred_ntprobnp_new_median_ols"
  )
  apparent_performance_new_median_paba <- compute_metrics(
    md,
    target_col = "NTproBNP",
    predicted_col = "pred_ntprobnp_new_median_paba"
  )
  print("Apparent Performance Metrics (Recalibrated New Model - MEDIAN, OLS):")
  print(apparent_performance_new_median_ols)
  print("Apparent Performance Metrics (Recalibrated New Model - MEDIAN, PaBa):")
  print(apparent_performance_new_median_paba)
  write_rds(
    apparent_performance_new_median_ols,
    output_path('apparent_performance_new_model_median_ols.rds')
  )
  write_csv(
    apparent_performance_new_median_ols,
    output_path('apparent_performance_new_model_median_ols.csv')
  )
  write_rds(
    apparent_performance_new_median_paba,
    output_path('apparent_performance_new_model_median_paba.rds')
  )
  write_csv(
    apparent_performance_new_median_paba,
    output_path('apparent_performance_new_model_median_paba.csv')
  )

  # Visualize Apparent Performance (Median Predictions)
  cal_plot_new_median_ols_log <- make_calibration_plot(
    md,
    "NTproBNP",
    "pred_ntprobnp_new_median_ols",
    "New Model (Recal, Median - OLS)",
    TRUE
  )
  ggsave(
    output_path('calibration_plot_new_median_ols_log.png'),
    cal_plot_new_median_ols_log,
    h = 6,
    w = 6,
    dpi = 300
  )
  cal_plot_new_median_ols_lin <- make_calibration_plot(
    md,
    "NTproBNP",
    "pred_ntprobnp_new_median_ols",
    "New Model (Recal, Median - OLS)",
    FALSE
  )
  ggsave(
    output_path('calibration_plot_new_median_ols_lin.png'),
    cal_plot_new_median_ols_lin,
    h = 6,
    w = 6,
    dpi = 300
  )
  ba_plot_new_median_ols_rel <- make_bland_altman_plot(
    md,
    "NTproBNP",
    "pred_ntprobnp_new_median_ols",
    "New Model (Recal, Median - OLS)",
    "relative"
  )
  ggsave(
    output_path('BA_new_median_ols_relative.png'),
    ba_plot_new_median_ols_rel,
    h = 6,
    w = 6,
    dpi = 300
  )
  ba_plot_new_median_ols_abs <- make_bland_altman_plot(
    md,
    "NTproBNP",
    "pred_ntprobnp_new_median_ols",
    "New Model (Recal, Median - OLS)",
    "absolute"
  )
  ggsave(
    output_path('BA_new_median_ols_absolute.png'),
    ba_plot_new_median_ols_abs,
    h = 6,
    w = 6,
    dpi = 300
  )

  cal_plot_new_median_paba_log <- make_calibration_plot(
    md,
    "NTproBNP",
    "pred_ntprobnp_new_median_paba",
    "New Model (Recal, Median - PaBa)",
    TRUE
  )
  ggsave(
    output_path('calibration_plot_new_median_paba_log.png'),
    cal_plot_new_median_paba_log,
    h = 6,
    w = 6,
    dpi = 300
  )
  cal_plot_new_median_paba_lin <- make_calibration_plot(
    md,
    "NTproBNP",
    "pred_ntprobnp_new_median_paba",
    "New Model (Recal, Median - PaBa)",
    FALSE
  )
  ggsave(
    output_path('calibration_plot_new_median_paba_lin.png'),
    cal_plot_new_median_paba_lin,
    h = 6,
    w = 6,
    dpi = 300
  )
  ba_plot_new_median_paba_rel <- make_bland_altman_plot(
    md,
    "NTproBNP",
    "pred_ntprobnp_new_median_paba",
    "New Model (Recal, Median - PaBa)",
    "relative"
  )
  ggsave(
    output_path('BA_new_median_paba_relative.png'),
    ba_plot_new_median_paba_rel,
    h = 6,
    w = 6,
    dpi = 300
  )
  ba_plot_new_median_paba_abs <- make_bland_altman_plot(
    md,
    "NTproBNP",
    "pred_ntprobnp_new_median_paba",
    "New Model (Recal, Median - PaBa)",
    "absolute"
  )
  ggsave(
    output_path('BA_new_median_paba_absolute.png'),
    ba_plot_new_median_paba_abs,
    h = 6,
    w = 6,
    dpi = 300
  )

  # Actual vs Predicted plots
  avp_new_median_ols <- make_actual_vs_pred_plot(
    md,
    "NTproBNP",
    "pred_ntprobnp_new_median_ols",
    "New Model (Median - OLS)",
    TRUE
  )
  ggsave(
    output_path('actual_vs_pred_new_median_ols.png'),
    avp_new_median_ols,
    h = 6,
    w = 6,
    dpi = 300
  )
  avp_new_median_paba <- make_actual_vs_pred_plot(
    md,
    "NTproBNP",
    "pred_ntprobnp_new_median_paba",
    "New Model (Median - PaBa)",
    TRUE
  )
  ggsave(
    output_path('actual_vs_pred_new_median_paba.png'),
    avp_new_median_paba,
    h = 6,
    w = 6,
    dpi = 300
  )

  # Passing-Bablok plots (original scale) for both recalibrations
  paba_new_median_ols_obj <- run_mcreg(
    md,
    x_col = NTproBNP,
    y_col = pred_ntprobnp_new_median_ols
  )
  plot_mcreg(
    paba_new_median_ols_obj,
    filename = output_path('paba_plot_new_median_ols.png'),
    title = "Passing-Bablok Regression (New Model - Median OLS)"
  )
  paba_new_median_paba_obj <- run_mcreg(
    md,
    x_col = NTproBNP,
    y_col = pred_ntprobnp_new_median_paba
  )
  plot_mcreg(
    paba_new_median_paba_obj,
    filename = output_path('paba_plot_new_median_paba.png'),
    title = "Passing-Bablok Regression (New Model - Median PaBa)"
  )
} else {
  print("Insufficient data after processing for new model apparent evaluation.")
  apparent_performance_new_median_ols <- NULL
  apparent_performance_new_median_paba <- NULL
}
# --- Optimism Correction Bootstrapping (for RECALIBRATED MEDIAN predictions; PaBa) ---
print(glue(
  "--- Starting Optimism Correction Bootstrapping ({CONFIG$bootstrap_reps} iterations for RECALIBRATED MEDIAN PaBa model) ---"
))
# Function to get performance metrics for one bootstrap iteration
calculate_optimism_step <- function(split, recal_params, knots) {
  bs_analysis <- analysis(split)
  bs_assessment <- assessment(split) # The original data
  # Fit model on bootstrap sample
  bs_fit <- fit_new_model_ols(bs_analysis, knots = knots)
  if (is.null(bs_fit)) {
    return(NULL)
  } # Skip if fit fails
  # 1. Performance on Bootstrap Sample (Apparent performance of bs_fit)
  bs_analysis_pred <- bs_analysis %>%
    mutate(
      log_pred_initial_bs = predict(bs_fit, newdata = .),
      log_pred_recal_bs = recal_params[1] +
        recal_params[2] * log_pred_initial_bs,
      pred_median_recal_bs = exp(log_pred_recal_bs)
    ) %>%
    filter(is.finite(pred_median_recal_bs))
  if (nrow(bs_analysis_pred) < 10) {
    return(NULL)
  } # Need enough data
  perf_on_bootstrap <- compute_metrics(
    bs_analysis_pred,
    "NTproBNP",
    "pred_median_recal_bs"
  ) %>%
    select(.metric, .estimator, .estimate_on_bootstrap = .estimate)
  # 2. Performance on Original Sample (Test performance of bs_fit)
  orig_assessment_pred <- bs_assessment %>%
    mutate(
      log_pred_initial_orig = predict(bs_fit, newdata = .),
      log_pred_recal_orig = recal_params[1] +
        recal_params[2] * log_pred_initial_orig,
      pred_median_recal_on_orig = exp(log_pred_recal_orig)
    ) %>%
    filter(is.finite(pred_median_recal_on_orig))
  if (nrow(orig_assessment_pred) < 10) {
    return(NULL)
  } # Need enough data
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
    filter(is.finite(optimism)) # Filter out any non-finite optimism values
  return(metrics_step)
}
# Perform bootstrapping
boots <- bootstraps(md, times = CONFIG$bootstrap_reps, apparent = FALSE) # Use apparent=FALSE
# Run the calculation over bootstrap samples safely
# Use multicore to handle potential NULL returns from calculate_optimism_step
optimism_cores <- get_bootstrap_cores()
message(glue("Running optimism correction on {optimism_cores} core(s)..."))
optimism_estimates_list <- parallel::mclapply(
  boots$splits,
  function(split) {
    calculate_optimism_step(split, log_recal_params, CONFIG$rcs_knots)
  },
  mc.cores = optimism_cores,
  mc.set.seed = TRUE
)
# Combine results, filtering out NULLs
optimism_estimates <- bind_rows(keep(optimism_estimates_list, ~ !is.null(.)))
print("Finished Optimism Correction Bootstrapping.")
# Summarize Optimism and Calculate Corrected Estimates
if (
  nrow(optimism_estimates) > 0 && !is.null(apparent_performance_new_median_paba)
) {
  optimism_summary <- optimism_estimates %>%
    group_by(.metric, .estimator) %>%
    summarise(
      mean_optimism = mean(optimism, na.rm = TRUE),
      n_boots = n(),
      # Calculate Bootstrap CI (percentile method) for test performance
      lower_ci_95 = quantile(
        .estimate_on_original,
        probs = 0.025,
        na.rm = TRUE
      ),
      upper_ci_95 = quantile(
        .estimate_on_original,
        probs = 0.975,
        na.rm = TRUE
      ),
      .groups = 'drop'
    ) %>%
    filter(is.finite(mean_optimism))
  # Join with apparent performance
  optimism_corrected_final <- apparent_performance_new_median_paba %>%
    select(.metric, .estimator, apparent_estimate = .estimate) %>% # Rename .estimate column
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
  print(
    "Optimism-Corrected Performance Metrics (Recalibrated New Model - MEDIAN, PaBa):"
  )
  print(optimism_corrected_final)
  write_rds(
    optimism_corrected_final,
    output_path('optimism_corrected_estimates_new_model_median_paba.rds')
  )
  write_csv(
    optimism_corrected_final,
    output_path('optimism_corrected_estimates_new_model_median_paba.csv')
  )
} else {
  print(
    "Could not calculate optimism-corrected estimates (insufficient bootstrap results or apparent performance missing)."
  )
  optimism_corrected_final <- NULL
}
print("--- Finished Section 4 ---")
#===============================================================================
# Section 5: Association with Clinical Outcomes (using RECALIBRATED MEDIAN Predictions)
#===============================================================================
print("--- Starting Section 5: Association with Clinical Outcomes ---")
# Use the 'md' dataframe which contains the final recalibrated median predictions `pred_ntprobnp_new_median_paba`
# --- Derive Biomarker Categories ---
print("Deriving biomarker categories...")
# Define target NT-proBNP thresholds
target_NTs <- CONFIG$target_ntprobnp_thresholds
cut_points_ntprobnp <- c(-Inf, target_NTs, Inf)
# Calculate corresponding BNP thresholds using the RECALIBRATED MEDIAN (PaBa) model
print("Calculating corresponding BNP thresholds using median covariates...")
median_covars_for_thresholds <- md %>%
  summarise(
    across(c(Age, BMI, CrCl, Hb_gdl), \(x) median(x, na.rm = TRUE)),
    Sex_M = 0.5, # Use midpoint for sex
    AF = Mode(AF, na.rm = TRUE), # Use Mode for binary
    .groups = 'drop'
  )
print("Median covariate grid for threshold calculation:")
print(median_covars_for_thresholds)
corresponding_BNPs_new <- sapply(target_NTs, function(nt) {
  get_BNP_for_NT_new_model(
    nt,
    median_covars_for_thresholds,
    fit_ols,
    log_recal_params,
    use_mean = FALSE
  )
})
print("Calculated corresponding BNP thresholds (may contain NAs):")
print(corresponding_BNPs_new)
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
      pred_ntprobnp_new_median_paba,
      CONFIG$class_breaks
    )
  ) %>%
  select(
    StudyID,
    composite_outcome,
    NTproBNP,
    BNP,
    pred_ntprobnp_new_median_paba,
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
  # --- Fit Logistic Regression Models ---
  print(
    "Fitting logistic regression models (GLM and Firth Logistic Regression)"
  )
  # Model A: Measured NT-proBNP categories vs. Outcome
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
    title = 'Table A: Association between Measured NT-proBNP Categories and Composite Outcome (Odds Ratios)'
  )
  print(table_a)
  save_as_docx(
    table_a,
    path = output_path("outcome_table_A_ntprobnp_cat.docx")
  )

  # Model B: Derived BNP categories vs. Outcome
  model_b_glm <- glm(
    composite_outcome ~ BNP_cat_derived,
    data = logit_md,
    family = binomial
  )
  model_b_firth <- logistf(
    composite_outcome ~ BNP_cat_derived,
    data = logit_md,
    pl = FALSE
  )
  table_b <- modelsummary(
    list('GLM' = model_b_glm, 'Firth GLM' = model_b_firth),
    estimate = "{estimate} [{conf.low}, {conf.high}]",
    statistic = "p.value",
    exponentiate = TRUE,
    gof_map = NA,
    output = "flextable",
    title = 'Table B: Association between Derived BNP Categories (New Model) and Composite Outcome (Odds Ratios)'
  )
  print(table_b)
  save_as_docx(
    table_b,
    path = output_path("outcome_table_B_bnp_derived_cat.docx")
  )

  # Model C: Predicted (Recalibrated Median, PaBa) NT-proBNP categories vs. Outcome
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
    title = 'Table C: Association between Predicted NT-proBNP Categories (New Model - Median, PaBa) and Composite Outcome (Odds Ratios)'
  )
  print(table_c)
  save_as_docx(
    table_c,
    path = output_path("outcome_table_C_pred_new_median_paba_cat.docx")
  )

  print("Finished fitting logistic models.")

  # --- Visualize Associations (Continuous Predictors using rms::lrm and contrast) ---
  print("Generating association plots (continuous predictors)...")

  # Helper function for plotting contrasts
  plot_orm_contrast <- function(fit, var_name, data, n_points = 100) {
    var_sym <- sym(var_name)
    # Use quantiles for sequence generation to handle outliers
    quantiles <- quantile(data[[var_name]], probs = c(0.01, 0.99), na.rm = TRUE)
    if (anyNA(quantiles) || quantiles[1] >= quantiles[2]) {
      rng <- range(data[[var_name]], na.rm = TRUE)
      if (anyNA(rng) || !all(is.finite(rng)) || rng[1] >= rng[2]) {
        warning(glue("Could not determine range for contrast plot: {var_name}"))
        return(NULL)
      }
      seq_vals <- seq(rng[1], rng[2], length.out = n_points)
    } else {
      seq_vals <- seq(quantiles[1], quantiles[2], length.out = n_points)
    }
    if (length(seq_vals) <= 1) {
      return(NULL)
    }

    k <- contrast(
      fit,
      setNames(list(seq_vals), var_name),
      setNames(list(mean(data[[var_name]], na.rm = TRUE)), var_name)
    )
    plot_df <- as.data.frame(k[c('Contrast', 'Lower', 'Upper')]) %>%
      mutate(x = seq_vals)

    ggplot(
      plot_df,
      aes(x = x, y = exp(Contrast), ymin = exp(Lower), ymax = exp(Upper))
    ) +
      geom_ribbon(alpha = 0.2, fill = PLOT_STYLE$colors$ribbon) +
      geom_line(color = PLOT_STYLE$colors$fit, linewidth = PLOT_STYLE$line_width) +
      geom_hline(
        yintercept = 1,
        linetype = 'dashed',
        color = PLOT_STYLE$colors$neutral
      ) +
      scale_y_log10(
        n.breaks = 8,
        labels = scales::label_number(accuracy = 0.1)
      ) +
      scale_x_continuous(labels = lblr) + # Use lblr for biomarker scale
      labs(
        y = 'Odds Ratio (log scale)',
        title = glue('Association: {var_name} vs. Outcome'),
        x = glue('{var_name} (pg/mL)')
      ) +
      my_theme
  }

  # Ensure datadist is set for logit_md
  dd_logit <- datadist(logit_md)
  options(datadist = 'dd_logit')

  # 1. Measured NTproBNP vs Outcome
  tryCatch(
    {
      fit_lrm_1 <- lrm(
        composite_outcome ~ rcs(NTproBNP, CONFIG$rcs_knots),
        data = logit_md,
        x = TRUE,
        y = TRUE
      )
      p1 <- plot_orm_contrast(fit_lrm_1, "NTproBNP", logit_md)
      if (!is.null(p1)) {
        ggsave(
          output_path('outcome_assoc_ntprobnp.png'),
          p1,
          h = 5,
          w = 5,
          dpi = 300
        )
      }
    },
    error = function(e) {
      warning("Failed to plot NTproBNP association: ", e$message)
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
      p2 <- plot_orm_contrast(fit_lrm_2, "BNP", logit_md)
      if (!is.null(p2)) {
        ggsave(
          output_path('outcome_assoc_bnp.png'),
          p2,
          h = 5,
          w = 5,
          dpi = 300
        )
      }
    },
    error = function(e) warning("Failed to plot BNP association: ", e$message)
  )

  # 3. Predicted (Recalibrated Median, PaBa) NTproBNP vs Outcome
  tryCatch(
    {
      fit_lrm_3 <- lrm(
        composite_outcome ~ rcs(
          pred_ntprobnp_new_median_paba,
          CONFIG$rcs_knots
        ),
        data = logit_md,
        x = TRUE,
        y = TRUE
      )
      p3 <- plot_orm_contrast(
        fit_lrm_3,
        "pred_ntprobnp_new_median_paba",
        logit_md
      )
      if (!is.null(p3)) {
        # Adjust labels for predicted variable
        p3 <- p3 +
          labs(
            title = 'Association: Predicted NTproBNP (New Model - Median, PaBa) vs. Outcome',
            x = 'Predicted NTproBNP (New Model - Median, PaBa, pg/mL)'
          )
        ggsave(
          output_path('outcome_assoc_pred_new_median_paba.png'),
          p3,
          h = 5,
          w = 5,
          dpi = 300
        )
      }
    },
    error = function(e) {
      warning("Failed to plot predicted NTproBNP association: ", e$message)
    }
  )
} # End of check for sufficient data for logistic regression

print("--- Finished Section 5 ---")
print("--- Analysis script finished successfully ---")

#===============================================================================
# Section 6: Assemble DOCX Report (figures + tables + explanations)
#===============================================================================
print("--- Starting Section 6: Assemble DOCX Report ---")

dir.create(CONFIG$output_dir, showWarnings = FALSE, recursive = TRUE)
report_rmd_path <- output_path("analysis_report.Rmd")
report_docx_path <- output_path("analysis_report.docx")

report_lines <- c(
  "---",
  "title: \"BNP to NT-proBNP Substudy Report\"",
  "output:",
  "  word_document:",
  "    toc: true",
  "    number_sections: true",
  "params:",
  sprintf("  output_dir: \"%s\"", CONFIG$output_dir),
  "---",
  "",
  "```{r setup, include=FALSE}",
  "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE, fig.align = 'center', fig.width = 6.5, fig.height = 6.5, out.width = '100%')",
  "suppressPackageStartupMessages({",
  "  library(dplyr)",
  "  library(readr)",
  "  library(ggplot2)",
  "  library(gtsummary)",
  "  library(modelsummary)",
  "  library(knitr)",
  "  library(flextable)",
  "})",
  "dir.create(params$output_dir, showWarnings = FALSE, recursive = TRUE)",
  "get_obj <- function(name) {",
  "  if (exists(name, envir = knitr::knit_global())) {",
  "    get(name, envir = knitr::knit_global())",
  "  } else {",
  "    NULL",
  "  }",
  "}",
  "format_table <- function(df, caption, digits = 3, font_size = 9, drop_model_type = TRUE) {",
  "  if (drop_model_type && 'model_type' %in% names(df)) {",
  "    df <- dplyr::select(df, -model_type)",
  "  }",
  "  df <- df %>%",
  "    dplyr::rename_with(~ gsub('^\\\\.', '', .x)) %>%",
  "    dplyr::rename_with(~ gsub('_', ' ', .x)) %>%",
  "    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, digits)))",
  "  ft <- flextable::flextable(df)",
  "  ft <- flextable::set_caption(ft, caption = caption)",
  "  ft <- flextable::fontsize(ft, size = font_size, part = 'all')",
  "  ft <- flextable::autofit(ft)",
  "  ft",
  "}",
  "page_break <- function() knitr::asis_output('\\\\pagebreak')",
  "filter_protocol_metrics <- function(df) {",
  "  if (!all(c('.metric', '.estimator') %in% names(df))) {",
  "    return(df)",
  "  }",
  "  allowed <- c(",
  "    'pearson_r', 'rsq', 'rmse',",
  "    'weighted_kappa_quadratic',",
  "    'sens', 'spec', 'ppv', 'npv',",
  "    'paba_intercept', 'paba_intercept_lower', 'paba_intercept_upper',",
  "    'paba_slope', 'paba_slope_lower', 'paba_slope_upper'",
  "  )",
  "  df <- dplyr::filter(df, .metric %in% allowed)",
  "  df <- dplyr::filter(df, !(.metric %in% c('sens','spec','ppv','npv') & .estimator != 'binary'))",
  "  df <- dplyr::filter(df, !(.metric %in% c('weighted_kappa_quadratic') & .estimator != 'weighted'))",
  "  df <- dplyr::filter(df, !(.metric %in% c('paba_intercept','paba_intercept_lower','paba_intercept_upper','paba_slope','paba_slope_lower','paba_slope_upper') & .estimator != 'PaBa'))",
  "  df <- dplyr::filter(df, !(.metric %in% c('pearson_r') & .estimator != 'pearson'))",
  "  df <- df %>%",
  "    dplyr::mutate(.metric = factor(.metric, levels = allowed)) %>%",
  "    dplyr::arrange(.metric, .estimator)",
  "  df",
  "}",
  "```",
  "",
  "# Overview",
  "This report compiles the figures and tables produced by the analysis script and describes how each output corresponds to the published protocol.",
  "",
  "# Baseline Characteristics",
  "Baseline characteristics are summarized for the initial cohort (Table 1).",
  "",
  "```{r table-baseline}",
  "baseline_table_final <- get_obj('baseline_table_final')",
  "if (!is.null(baseline_table_final)) {",
  "  baseline_table_final <- baseline_table_final %>%",
  "    gtsummary::modify_caption('Baseline characteristics of the initial cohort.')",
  "  ft <- gtsummary::as_flex_table(baseline_table_final)",
  "  ft <- flextable::fontsize(ft, size = 8, part = 'all')",
  "  ft <- flextable::autofit(ft)",
  "  ft",
  "} else {",
  "  cat('Baseline table not available in memory.\\n')",
  "}",
  "```",
  "\\pagebreak",
  "",
  "# Primary Analysis: Kasahara Conversion Formula",
  "The Kasahara formula was applied to BNP values to predict NT-proBNP. Agreement between measured and predicted values is summarized below.",
  "",
  "```{r table-kasahara-metrics}",
  "kasahara_metrics_path <- file.path(params$output_dir, 'kasahara_model_prediction_errors.csv')",
  "if (file.exists(kasahara_metrics_path)) {",
  "  kasahara_metrics <- readr::read_csv(kasahara_metrics_path, show_col_types = FALSE)",
  "  kasahara_metrics <- filter_protocol_metrics(kasahara_metrics)",
  "  format_table(",
  "    kasahara_metrics,",
  "    'Kasahara model performance metrics (includes RMSE and R-squared).',",
  "    digits = 3,",
  "    font_size = 8",
  "  )",
  "} else {",
  "  cat('Kasahara metrics file not found.\\n')",
  "}",
  "```",
  "\\pagebreak",
  "",
  "```{r fig-kasahara-cal-log, fig.cap='Kasahara model calibration plot (log scale).'}",
  "knitr::include_graphics(file.path(params$output_dir, 'calibration_plot_kasahara_log.png'))",
  "```",
  "\\pagebreak",
  "",
  "```{r fig-kasahara-cal-lin, fig.cap='Kasahara model calibration plot (linear scale).'}",
  "knitr::include_graphics(file.path(params$output_dir, 'calibration_plot_kasahara_lin.png'))",
  "```",
  "\\pagebreak",
  "",
  "```{r fig-kasahara-ba-rel, fig.cap='Kasahara model Bland-Altman plot (relative error).'}",
  "knitr::include_graphics(file.path(params$output_dir, 'BA_kasahara_relative.png'))",
  "```",
  "\\pagebreak",
  "",
  "```{r fig-kasahara-ba-abs, fig.cap='Kasahara model Bland-Altman plot (absolute error).'}",
  "knitr::include_graphics(file.path(params$output_dir, 'BA_kasahara_absolute.png'))",
  "```",
  "\\pagebreak",
  "",
  "```{r fig-kasahara-paba, fig.cap='Passing-Bablok regression for Kasahara model (measured vs predicted).'}",
  "knitr::include_graphics(file.path(params$output_dir, 'paba_plot_kasahara.png'))",
  "```",
  "\\pagebreak",
  "",
  "# Updated Conversion Formula (Recalibrated Model)",
  "A new model was fit on log-transformed NT-proBNP with the specified predictors and interactions. Passing-Bablok regression (bootstrap) was used for recalibration.",
  "",
  "```{r table-paba-boot}",
  "paba_boot_path <- file.path(params$output_dir, 'paba_loglog_bootstrap_summary.csv')",
  "if (file.exists(paba_boot_path)) {",
  "  paba_boot <- readr::read_csv(paba_boot_path, show_col_types = FALSE)",
  "  format_table(",
  "    paba_boot,",
  "    'Bootstrapped Passing-Bablok calibration (log-log scale).',",
  "    digits = 4,",
  "    font_size = 8",
  "  )",
  "} else {",
  "  cat('Bootstrapped Passing-Bablok summary not found.\\n')",
  "}",
  "```",
  "\\pagebreak",
  "",
  "```{r table-newmodel-metrics}",
  "newmodel_metrics_path <- file.path(params$output_dir, 'apparent_performance_new_model_median_paba.csv')",
  "if (file.exists(newmodel_metrics_path)) {",
  "  newmodel_metrics <- readr::read_csv(newmodel_metrics_path, show_col_types = FALSE)",
  "  newmodel_metrics <- filter_protocol_metrics(newmodel_metrics)",
  "  format_table(",
  "    newmodel_metrics,",
  "    'Apparent performance metrics for the recalibrated model (median prediction, PaBa).',",
  "    digits = 3,",
  "    font_size = 8",
  "  )",
  "} else {",
  "  cat('New model performance file not found.\\n')",
  "}",
  "```",
  "\\pagebreak",
  "",
  "```{r table-newmodel-optimism}",
  "optimism_path <- file.path(params$output_dir, 'optimism_corrected_estimates_new_model_median_paba.csv')",
  "if (file.exists(optimism_path)) {",
  "  optimism_tbl <- readr::read_csv(optimism_path, show_col_types = FALSE)",
  "  optimism_tbl <- filter_protocol_metrics(optimism_tbl)",
  "  format_table(",
  "    optimism_tbl,",
  "    'Optimism-corrected performance metrics (bootstrap).',",
  "    digits = 3,",
  "    font_size = 8",
  "  )",
  "} else {",
  "  cat('Optimism-corrected metrics file not found.\\n')",
  "}",
  "```",
  "\\pagebreak",
  "",
  "```{r fig-new-cal-log, fig.cap='Recalibrated model calibration plot (log scale, median PaBa).'}",
  "knitr::include_graphics(file.path(params$output_dir, 'calibration_plot_new_median_paba_log.png'))",
  "```",
  "\\pagebreak",
  "",
  "```{r fig-new-cal-lin, fig.cap='Recalibrated model calibration plot (linear scale, median PaBa).'}",
  "knitr::include_graphics(file.path(params$output_dir, 'calibration_plot_new_median_paba_lin.png'))",
  "```",
  "\\pagebreak",
  "",
  "```{r fig-new-ba-rel, fig.cap='Recalibrated model Bland-Altman plot (relative error, median PaBa).'}",
  "knitr::include_graphics(file.path(params$output_dir, 'BA_new_median_paba_relative.png'))",
  "```",
  "\\pagebreak",
  "",
  "```{r fig-new-ba-abs, fig.cap='Recalibrated model Bland-Altman plot (absolute error, median PaBa).'}",
  "knitr::include_graphics(file.path(params$output_dir, 'BA_new_median_paba_absolute.png'))",
  "```",
  "\\pagebreak",
  "",
  "```{r fig-new-paba, fig.cap='Passing-Bablok regression for recalibrated model (median PaBa, measured vs predicted).'}",
  "knitr::include_graphics(file.path(params$output_dir, 'paba_plot_new_median_paba.png'))",
  "```",
  "\\pagebreak",
  "",
  "# Secondary Analysis: Outcome Associations",
  "Logistic regression models evaluate associations between biomarker categories and the 30-day composite outcome.",
  "",
  "```{r table-outcome-a}",
  "model_a_glm <- get_obj('model_a_glm')",
  "model_a_firth <- get_obj('model_a_firth')",
  "if (!is.null(model_a_glm) && !is.null(model_a_firth)) {",
  "  tbl_a <- modelsummary::modelsummary(",
  "    list('GLM' = model_a_glm, 'Firth GLM' = model_a_firth),",
  "    estimate = '{estimate} [{conf.low}, {conf.high}]',",
  "    statistic = 'p.value',",
  "    exponentiate = TRUE,",
  "    gof_map = NA,",
  "    output = 'data.frame'",
  "  )",
  "  format_table(",
  "    tbl_a,",
  "    'Association between measured NT-proBNP categories and the composite outcome (odds ratios).',",
  "    digits = 3,",
  "    font_size = 8,",
  "    drop_model_type = FALSE",
  "  )",
  "} else {",
  "  cat('Outcome model A not available.\\n')",
  "}",
  "```",
  "\\pagebreak",
  "",
  "```{r table-outcome-b}",
  "model_b_glm <- get_obj('model_b_glm')",
  "model_b_firth <- get_obj('model_b_firth')",
  "if (!is.null(model_b_glm) && !is.null(model_b_firth)) {",
  "  tbl_b <- modelsummary::modelsummary(",
  "    list('GLM' = model_b_glm, 'Firth GLM' = model_b_firth),",
  "    estimate = '{estimate} [{conf.low}, {conf.high}]',",
  "    statistic = 'p.value',",
  "    exponentiate = TRUE,",
  "    gof_map = NA,",
  "    output = 'data.frame'",
  "  )",
  "  format_table(",
  "    tbl_b,",
  "    'Association between derived BNP categories and the composite outcome (odds ratios).',",
  "    digits = 3,",
  "    font_size = 8,",
  "    drop_model_type = FALSE",
  "  )",
  "} else {",
  "  cat('Outcome model B not available.\\n')",
  "}",
  "```",
  "\\pagebreak",
  "",
  "```{r table-outcome-c}",
  "model_c_glm <- get_obj('model_c_glm')",
  "model_c_firth <- get_obj('model_c_firth')",
  "if (!is.null(model_c_glm) && !is.null(model_c_firth)) {",
  "  tbl_c <- modelsummary::modelsummary(",
  "    list('GLM' = model_c_glm, 'Firth GLM' = model_c_firth),",
  "    estimate = '{estimate} [{conf.low}, {conf.high}]',",
  "    statistic = 'p.value',",
  "    exponentiate = TRUE,",
  "    gof_map = NA,",
  "    output = 'data.frame'",
  "  )",
  "  format_table(",
  "    tbl_c,",
  "    'Association between predicted NT-proBNP categories (median PaBa) and the composite outcome (odds ratios).',",
  "    digits = 3,",
  "    font_size = 8,",
  "    drop_model_type = FALSE",
  "  )",
  "} else {",
  "  cat('Outcome model C not available.\\n')",
  "}",
  "```",
  "\\pagebreak",
  "",
  "```{r fig-outcome-ntprobnp, fig.cap='Association of measured NT-proBNP with outcome (spline model).'}",
  "knitr::include_graphics(file.path(params$output_dir, 'outcome_assoc_ntprobnp.png'))",
  "```",
  "\\pagebreak",
  "",
  "```{r fig-outcome-bnp, fig.cap='Association of measured BNP with outcome (spline model).'}",
  "knitr::include_graphics(file.path(params$output_dir, 'outcome_assoc_bnp.png'))",
  "```",
  "\\pagebreak",
  "",
  "```{r fig-outcome-pred, fig.cap='Association of predicted NT-proBNP (new model, median PaBa) with outcome (spline model).'}",
  "knitr::include_graphics(file.path(params$output_dir, 'outcome_assoc_pred_new_median_paba.png'))",
  "```",
  "\\pagebreak"
)

page_break_chunk <- c(
  "```{r, results='asis'}",
  "page_break()",
  "```"
)
report_lines <- unlist(lapply(
  report_lines,
  function(line) if (identical(line, "\\pagebreak")) page_break_chunk else line
))

writeLines(report_lines, report_rmd_path)

if (requireNamespace("rmarkdown", quietly = TRUE)) {
  tryCatch(
    {
      rmarkdown::render(
        report_rmd_path,
        output_dir = CONFIG$output_dir,
        output_file = "analysis_report.docx",
        params = list(output_dir = CONFIG$output_dir),
        envir = environment(),
        quiet = FALSE
      )
      print(glue("Report saved to: {report_docx_path}"))
    },
    error = function(e) {
      warning("Report rendering failed: ", e$message)
    }
  )
} else {
  warning("Package 'rmarkdown' is not available; skipping report generation.")
}

print("--- Finished Section 6 ---")
