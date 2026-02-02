#===============================================================================
# Utility Functions
#===============================================================================

# --- Output Helper ---
output_path <- function(...) {
  file.path(CONFIG$output_dir, ...)
}

print_plot <- function(p) {
  if (interactive()) {
    print(p)
  }
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

# --- Metric Calculation Setup (protocol-aligned) ---
binary_metrics_set <- metric_set(sens, spec, ppv, npv)
continuous_metrics_set <- metric_set(rmse, rsq)

# --- Plotting Helpers ---
# Simple comma formatting for axis tick labels (no units - units go in axis title)
lblr <- scales::label_comma()

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
    panel.grid.major = element_blank(),
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

# --- Mode Function ---
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

# --- Category Helper Functions ---
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

make_ntprobnp_cat <- function(x, breaks = CONFIG$class_breaks) {
  labels <- c(
    glue("<{breaks[1]}"),
    glue("{breaks[1]}-<{breaks[2]}"),
    glue("{breaks[2]}-<{breaks[3]}"),
    glue(">={breaks[3]}")
  )
  make_interval_cat(x, breaks, labels)
}

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

# --- Weighted Kappa ---
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

# --- MCR Helper ---
is_mcr_result <- function(x) {
  isS4(x) && "para" %in% methods::slotNames(x)
}

# --- Actual vs Predicted Plot ---
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
    return(ggplot() + my_theme + labs(x = NULL, y = NULL))
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
      x = "Predicted NT-proBNP (pg/mL)",
      y = "Measured NT-proBNP (pg/mL)"
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

# --- Compute Metrics Function ---
# Evaluates model performance with continuous metrics on log10 scale (matching Kasahara)
# and classification metrics on original scale (where clinical thresholds are defined)
compute_metrics <- function(
  data,
  target_col,
  predicted_col,
  binary_threshold = CONFIG$binary_threshold,
  class_breaks = CONFIG$class_breaks
) {
  target_sym <- sym(target_col)
  predicted_sym <- sym(predicted_col)

  if (!all(c(target_col, predicted_col) %in% names(data))) {
    stop("Target or predicted column not found in data for compute_metrics")
  }
  d_filtered <- data %>%
    filter(
      is.finite({{ target_sym }}) & is.finite({{ predicted_sym }}) &
        {{ target_sym }} > 0 & {{ predicted_sym }} > 0
    )
  if (nrow(d_filtered) < 2) {
    warning("Less than 2 non-missing pairs for compute_metrics")
    return(tibble(
      .metric = character(),
      .estimator = character(),
      .estimate = numeric()
    ))
  }

  # Create log10 transformed columns for continuous metrics
  d_filtered <- d_filtered %>%
    mutate(
      log10_target = log10({{ target_sym }}),
      log10_predicted = log10({{ predicted_sym }})
    )

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
  binary_event_levels <- c(
    glue(">={binary_threshold}"),
    glue("<{binary_threshold}")
  )
  d_processed <- d_processed %>%
    mutate(
      binary_tgt = factor(binary_tgt, levels = binary_event_levels),
      binary_pred = factor(binary_pred, levels = binary_event_levels)
    )

  # Continuous metrics on LOG10 scale (matching Kasahara)
  rmse_log10 <- sqrt(mean((d_processed$log10_target - d_processed$log10_predicted)^2))
  rsq_log10 <- cor(d_processed$log10_target, d_processed$log10_predicted)^2
  pearson_r_log10 <- cor(
    d_processed$log10_target,
    d_processed$log10_predicted,
    use = "complete.obs",
    method = "pearson"
  )
  cont_m <- tibble(
    .metric = c("rmse_log10", "rsq_log10", "pearson_r_log10"),
    .estimator = c("log10", "log10", "log10"),
    .estimate = c(rmse_log10, rsq_log10, pearson_r_log10)
  )

  # Binary metrics on ORIGINAL scale (clinical thresholds)
  bin_m <- binary_metrics_set(
    d_processed,
    truth = binary_tgt,
    estimate = binary_pred
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

  # Calibration on LOG10 scale (OLS)
  ols_fit <- tryCatch(
    stats::lm(log10_target ~ log10_predicted, data = d_processed),
    error = function(e) NULL
  )
  if (!is.null(ols_fit) && length(coef(ols_fit)) >= 2) {
    ols_m <- tibble(
      .metric = c("ols_intercept_log10", "ols_slope_log10"),
      .estimator = "OLS_log10",
      .estimate = c(unname(coef(ols_fit)[1]), unname(coef(ols_fit)[2]))
    )
  } else {
    ols_m <- tibble(
      .metric = character(),
      .estimator = character(),
      .estimate = numeric()
    )
    warning("OLS calibration failed or insufficient data.")
  }

  # Calibration on LOG10 scale (Passing-Bablok, per protocol)
  paba_model <- tryCatch(
    {
      mcreg(
        x = d_processed$log10_predicted,
        y = d_processed$log10_target,
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
        "paba_intercept_log10",
        "paba_intercept_lower_log10",
        "paba_intercept_upper_log10",
        "paba_slope_log10",
        "paba_slope_lower_log10",
        "paba_slope_upper_log10"
      ),
      .estimator = "PaBa_log10",
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

  binary_label <- glue(
    "{target_col} >= {binary_threshold} vs < {binary_threshold}"
  )
  protocol_metrics <- bind_rows(cont_m, bin_m, w_kappa_m, ols_m, paba_m) %>%
    mutate(
      threshold_label = if_else(
        .metric %in% c("sens", "spec", "ppv", "npv"),
        binary_label,
        ""
      ),
      model_type = predicted_col
    ) %>%
    filter(
      .metric %in%
        c(
          "rmse_log10",
          "rsq_log10",
          "pearson_r_log10",
          "weighted_kappa_quadratic",
          "sens",
          "spec",
          "ppv",
          "npv",
          "ols_intercept_log10",
          "ols_slope_log10",
          "paba_intercept_log10",
          "paba_intercept_lower_log10",
          "paba_intercept_upper_log10",
          "paba_slope_log10",
          "paba_slope_lower_log10",
          "paba_slope_upper_log10"
        )
    )
  protocol_metrics
}

# --- Bland-Altman Plot ---
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
        ifelse(abs(av) > 1e-6, diff / av, NA)
      } else {
        diff
      }
    ) %>%
    filter(!is.na(av) & !is.na(plotted_err))
  if (nrow(plot_df) < 10) {
    warning(glue(
      "Insufficient data (<10 points) for Bland-Altman plot: {model_name} ({type})"
    ))
    return(ggplot() + my_theme + labs(x = NULL, y = NULL))
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
      x = "Mean of Measured and Predicted (pg/mL)",
      y = y_label
    ) +
    my_theme
  return(p)
}

# --- Kasahara Prediction Function ---
predict_ntprobnp_kasahara <- function(BNP, Age, BMI, Hb_gdl, CrCl, Sex_M, AF) {
  BNP <- as.numeric(BNP)
  Age <- as.numeric(Age)
  BMI <- as.numeric(BMI)
  Hb_gdl <- as.numeric(Hb_gdl)
  CrCl <- as.numeric(CrCl)
  Sex_M <- as.numeric(Sex_M)
  AF <- as.numeric(AF)

  log10_BNP <- ifelse(BNP > 0 & is.finite(BNP), log10(BNP), NA_real_)

  k1 <- CONFIG$kasahara_crcl_knots[1]
  k2 <- CONFIG$kasahara_crcl_knots[2]
  k3 <- CONFIG$kasahara_crcl_knots[3]

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
      0.0164 * (1 - Sex_M) +
      0.194 * AF,
    NA_real_
  )

  result <- 10^pred_log10
  result <- ifelse(is.finite(result), result, NA_real_)
  return(result)
}

# --- Kasahara-like Model (same predictors as Kasahara, re-estimated on new data) ---
# Kasahara predictors: log10(BNP), Age, BMI, Hb, CrCl (RCS), Sex, AF
# This is an ADDITIVE model (no interactions) matching Kasahara's structure
fit_kasahara_like_model <- function(data, knots = CONFIG$rcs_knots) {
  required_cols <- c(
    "log10_ntprobnp",
    "Age",
    "CrCl",
    "log10_bnp",
    "Sex_M",
    "AF",
    "BMI",
    "Hb_gdl"
  )
  if (!all(required_cols %in% names(data))) {
    missing_cols <- setdiff(required_cols, names(data))
    stop(
      "Missing required columns for Kasahara-like model: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  data_model <- data[, required_cols, drop = FALSE]
  if (nrow(na.omit(data_model)) < 20) {
    warning("Less than 20 complete cases for Kasahara-like model fitting.")
    return(NULL)
  }
  # Additive model matching Kasahara structure (no interactions)
  # Kasahara uses: log10(BNP), Age, BMI, Hb, CrCl (with RCS), Sex, AF
  frml <- formula(glue(
    "log10_ntprobnp ~ rcs(log10_bnp, {knots}) + rcs(Age, {knots}) +
                      rcs(CrCl, {knots}) + Sex_M + AF +
                      rcs(BMI, {knots}) + rcs(Hb_gdl, {knots})"
  ))
  dd_base <- tryCatch(
    datadist(data_model),
    error = function(e) {
      warning("datadist failed for Kasahara-like model: ", e$message)
      return(NULL)
    }
  )
  if (is.null(dd_base)) {
    return(NULL)
  }
  # Assign to global environment for rms to find it
  assign("dd_base", dd_base, envir = .GlobalEnv)
  old_datadist <- options("datadist")
  options(datadist = "dd_base")
  on.exit({
    options(old_datadist)
    if (exists("dd_base", envir = .GlobalEnv)) {
      rm("dd_base", envir = .GlobalEnv)
    }
  }, add = TRUE)
  fit <- tryCatch(
    ols(
      formula = frml,
      data = data_model,
      x = TRUE,
      y = TRUE,
      na.action = na.omit
    ),
    error = function(e) {
      warning("Kasahara-like model fitting failed: ", e$message)
      return(NULL)
    }
  )
  return(fit)
}

# --- Protocol Model (Kasahara predictors + CrCl*Age*BNP interaction per protocol) ---
# Per protocol: "The model will be extended with an interaction term between
# creatinine clearance, age, and BNP to determine if the interaction term
# improves the precision of the NT-proBNP estimates."
fit_protocol_model <- function(data, knots = CONFIG$rcs_knots) {
  required_cols <- c(
    "log10_ntprobnp",
    "Age",
    "CrCl",
    "log10_bnp",
    "Sex_M",
    "AF",
    "BMI",
    "Hb_gdl"
  )
  if (!all(required_cols %in% names(data))) {
    missing_cols <- setdiff(required_cols, names(data))
    stop(
      "Missing required columns for protocol model: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  data_model <- data[, required_cols, drop = FALSE]
  if (nrow(na.omit(data_model)) < 20) {
    warning("Less than 20 complete cases for protocol model fitting.")
    return(NULL)
  }
  # Protocol model: Kasahara predictors + 3-way interaction (Age * CrCl * log10_bnp)
  frml <- formula(glue(
    "log10_ntprobnp ~ rcs(log10_bnp, {knots}) * rcs(Age, {knots}) * rcs(CrCl, {knots}) +
                      Sex_M + AF + rcs(BMI, {knots}) + rcs(Hb_gdl, {knots})"
  ))
  dd_prot <- tryCatch(
    datadist(data_model),
    error = function(e) {
      warning("datadist failed for protocol model: ", e$message)
      return(NULL)
    }
  )
  if (is.null(dd_prot)) {
    return(NULL)
  }
  # Assign to global environment for rms to find it
  assign("dd_prot", dd_prot, envir = .GlobalEnv)
  old_datadist <- options("datadist")
  options(datadist = "dd_prot")
  on.exit({
    options(old_datadist)
    if (exists("dd_prot", envir = .GlobalEnv)) {
      rm("dd_prot", envir = .GlobalEnv)
    }
  }, add = TRUE)
  fit <- tryCatch(
    ols(
      formula = frml,
      data = data_model,
      x = TRUE,
      y = TRUE,
      na.action = na.omit
    ),
    error = function(e) {
      warning("Protocol model fitting failed: ", e$message)
      return(NULL)
    }
  )
  return(fit)
}

# --- New Model Fitting Function (log10 scale, matching Kasahara) ---
fit_new_model_ols <- function(data, knots = CONFIG$rcs_knots) {
  required_cols <- c(
    "log10_ntprobnp",
    "Age",
    "CrCl",
    "log10_bnp",
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
    warning("Less than 20 complete cases for OLS model fitting.")
    return(NULL)
  }
  frml_ols <- formula(glue(
    "log10_ntprobnp ~ rcs(Age, {knots}) * rcs(CrCl, {knots}) * rcs(log10_bnp, {knots}) +
                      Sex_M + rcs(BMI, {knots}) + AF + rcs(Hb_gdl, {knots})"
  ))
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
        weights = NULL,
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

# --- Root Finding Function (log10 scale, matching Kasahara) ---
get_BNP_for_NT_new_model <- function(
  target_NT,
  covars_grid,
  fitted_model,
  log10_recal_params = NULL,
  use_mean = FALSE
) {
  if (is.null(fitted_model) || !inherits(fitted_model, "ols")) {
    warning("Invalid fitted_model provided to get_BNP_for_NT_new_model.")
    return(NA_real_)
  }
  # For mean prediction with lognormal: E[Y] = 10^(mu + sigma^2 * ln(10) / 2)
  # Simplified: use median (10^mu) unless use_mean = TRUE
  sigma <- fitted_model$stats['Sigma']
  if (!is.finite(sigma)) {
    sigma <- 0
  }
  objective_function <- function(log10_bnp_val) {
    temp_data <- covars_grid %>% mutate(log10_bnp = log10_bnp_val)
    predicted_log10_NT_initial <- predict(fitted_model, newdata = temp_data)
    if (!is.null(log10_recal_params) && length(log10_recal_params) == 2) {
      I_log10 <- log10_recal_params[1]
      S_log10 <- log10_recal_params[2]
      if (abs(S_log10) < 1e-6) {
        S_log10 <- sign(S_log10) * 1e-6
      }
      predicted_log10_NT_final <- I_log10 + S_log10 * predicted_log10_NT_initial
    } else {
      predicted_log10_NT_final <- predicted_log10_NT_initial
    }
    expected_NT_pred <- if (isTRUE(use_mean)) {
      # Mean of lognormal: 10^(mu + sigma^2 * ln(10) / 2)
      10^(predicted_log10_NT_final + (sigma^2) * log(10) / 2)
    } else {
      10^predicted_log10_NT_final
    }
    return(expected_NT_pred - target_NT)
  }
  root_result <- tryCatch(
    {
      uniroot(
        objective_function,
        interval = CONFIG$log10_bnp_interval,
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
    return(10^root_result$root)
  } else {
    return(NA_real_)
  }
}

# --- Root Finding Function for Recalibrated Kasahara Model ---
# Finds BNP value that produces target NT-proBNP using the recalibrated Kasahara equation
get_BNP_for_NT_kasahara_recal <- function(
  target_NT,
  covars,
  recal_intercept,
  recal_slope
) {
  # Extract covariate values
  Age <- as.numeric(covars$Age)
  BMI <- as.numeric(covars$BMI)
  CrCl <- as.numeric(covars$CrCl)
  Hb_gdl <- as.numeric(covars$Hb_gdl)
  Sex_M <- as.numeric(covars$Sex_M)
  AF <- as.numeric(covars$AF)

  # Validate inputs
  if (!all(is.finite(c(Age, BMI, CrCl, Hb_gdl, Sex_M, AF, recal_intercept, recal_slope)))) {
    warning("Invalid covariates or recalibration parameters.")
    return(NA_real_)
  }

  # Objective function: find BNP where recalibrated prediction equals target
  objective_function <- function(log10_bnp_val) {
    BNP_val <- 10^log10_bnp_val
    # Get Kasahara prediction
    kasahara_pred <- predict_ntprobnp_kasahara(
      BNP = BNP_val,
      Age = Age,
      BMI = BMI,
      Hb_gdl = Hb_gdl,
      CrCl = CrCl,
      Sex_M = Sex_M,
      AF = AF
    )
    if (!is.finite(kasahara_pred) || kasahara_pred <= 0) {
      return(Inf)
    }
    # Apply recalibration: 10^(I + S × log10(kasahara_pred))
    recal_pred <- 10^(recal_intercept + recal_slope * log10(kasahara_pred))
    return(recal_pred - target_NT)
  }

  # Use root finding
  root_result <- tryCatch(
    {
      uniroot(
        objective_function,
        interval = CONFIG$log10_bnp_interval,
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
    return(10^root_result$root)
  } else {
    return(NA_real_)
  }
}

# --- Calibration Plot ---
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
    )
  if (nrow(plot_data) < 10) {
    warning(glue(
      "Insufficient data (<10 points) for calibration plot: {model_name}"
    ))
    return(ggplot() + my_theme + labs(x = NULL, y = NULL))
  }
  if (log_scale) {
    lims <- range(c(plot_data$actual, plot_data$predicted), na.rm = TRUE)
    lims <- c(max(lims[1] * 0.8, 1), lims[2] * 1.2)
  } else {
    lims <- range(c(plot_data$actual, plot_data$predicted), na.rm = TRUE)
    buffer <- (lims[2] - lims[1]) * 0.05
    lims <- c(lims[1] - buffer, lims[2] + buffer)
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
    geom_smooth(
      method = "lm",
      se = TRUE,
      color = PLOT_STYLE$colors$fit,
      fill = PLOT_STYLE$colors$ribbon,
      formula = y ~ x,
      linewidth = PLOT_STYLE$line_width
    ) +
    labs(
      x = "Predicted NT-proBNP (pg/mL)",
      y = "Measured NT-proBNP (pg/mL)"
    ) +
    my_theme +
    coord_cartesian(xlim = lims, ylim = lims)
  if (log_scale) {
    p <- p + scale_x_log10(labels = lblr) + scale_y_log10(labels = lblr)
  } else {
    p <- p +
      scale_x_continuous(labels = lblr) +
      scale_y_continuous(labels = lblr)
  }
  return(p)
}

# --- Log-Log Calibration Plot ---
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
    return(ggplot() + my_theme + labs(x = NULL, y = NULL))
  }
  log_fit <- lm(log10(actual) ~ log10(predicted), data = plot_data)
  log_intercept <- coef(log_fit)[1]
  log_slope <- coef(log_fit)[2]
  plot_data <- plot_data %>%
    mutate(fit = 10^(log_intercept + log_slope * log10(predicted)))
  lims <- range(c(plot_data$actual, plot_data$predicted), na.rm = TRUE)
  lims <- c(max(lims[1] * 0.8, 1), lims[2] * 1.2)
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
      x = "Predicted NT-proBNP (pg/mL)",
      y = "Measured NT-proBNP (pg/mL)"
    ) +
    my_theme +
    coord_cartesian(xlim = lims, ylim = lims)
  return(p)
}

# --- Bias-Corrected Calibration Plot ---
make_bias_corrected_calibration_plot <- function(
  calibrate_obj,
  x_label = "Predicted log10(NT-proBNP)",
  y_label = "Observed log10(NT-proBNP)"
) {
  is_calibrate_obj <- inherits(calibrate_obj, "calibrate") ||
    inherits(calibrate_obj, "calibrate.default")
  if (!is_calibrate_obj) {
    warning(
      "Invalid calibrate object provided; cannot plot bias-corrected curve."
    )
    return(NULL)
  }

  cal_core <- tryCatch(unclass(calibrate_obj), error = function(e) calibrate_obj)
  if (is.matrix(cal_core)) {
    cal_df <- as.data.frame(cal_core)
  } else {
    cal_df <- as.data.frame(calibrate_obj)
    if (identical(names(cal_df), "x") && is.matrix(cal_df$x)) {
      cal_df <- as.data.frame(cal_df$x)
    }
  }
  colnames(cal_df) <- sub("^x\\.", "", colnames(cal_df))

  if (!("predy" %in% names(cal_df))) {
    pred_candidates <- c("predicted", "prediction")
    pred_match <- pred_candidates[
      tolower(pred_candidates) %in% tolower(names(cal_df))
    ]
    if (length(pred_match) > 0) {
      names(cal_df)[tolower(names(cal_df)) == tolower(pred_match[1])] <- "predy"
    }
  }
  if (!("calibrated.orig" %in% names(cal_df))) {
    orig_candidates <- c("calibrated", "calibrated.apparent")
    orig_match <- orig_candidates[
      tolower(orig_candidates) %in% tolower(names(cal_df))
    ]
    if (length(orig_match) > 0) {
      names(cal_df)[tolower(names(cal_df)) == tolower(orig_match[1])] <-
        "calibrated.orig"
    }
  }
  required_cols <- c("predy", "calibrated.corrected", "calibrated.orig")
  if (!all(required_cols %in% names(cal_df))) {
    warning("Calibrate object missing required columns; cannot plot.")
    return(NULL)
  }

  cal_df <- cal_df %>%
    filter(
      is.finite(predy) &
        is.finite(calibrated.corrected) &
        is.finite(calibrated.orig)
    )
  if (nrow(cal_df) < 5) {
    warning("Insufficient data for bias-corrected calibration plot.")
    return(NULL)
  }

  cal_long <- bind_rows(
    tibble(
      predy = cal_df$predy,
      observed = cal_df$calibrated.corrected,
      curve = "Bias-corrected"
    ),
    tibble(
      predy = cal_df$predy,
      observed = cal_df$calibrated.orig,
      curve = "Apparent"
    )
  ) %>%
    mutate(curve = factor(curve, levels = c("Bias-corrected", "Apparent")))

  p <- ggplot(cal_df, aes(x = predy)) +
    geom_abline(
      slope = 1,
      intercept = 0,
      color = PLOT_STYLE$colors$identity,
      linetype = "dashed",
      linewidth = PLOT_STYLE$line_width
    ) +
    geom_line(
      data = cal_long,
      aes(y = observed, color = curve, linetype = curve),
      linewidth = PLOT_STYLE$line_width
    ) +
    scale_color_manual(
      values = c(
        "Bias-corrected" = PLOT_STYLE$colors$fit,
        "Apparent" = PLOT_STYLE$colors$limits
      )
    ) +
    scale_linetype_manual(
      values = c("Bias-corrected" = "solid", "Apparent" = "dotted")
    ) +
    labs(x = x_label, y = y_label) +
    my_theme +
    coord_equal()

  return(p)
}

# --- MCR Helpers ---
run_mcreg <- function(data, x_col, y_col) {
  df_mcreg <- data %>% select(x = {{ x_col }}, y = {{ y_col }}) %>% na.omit()
  if (nrow(df_mcreg) < 10) {
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

plot_mcreg <- function(
  mcreg_obj,
  filename,
  x_label = "Predicted NT-proBNP (pg/mL)",
  y_label = "Measured NT-proBNP (pg/mL)"
) {
  if (is.null(mcreg_obj) || !is_mcr_result(mcreg_obj)) {
    return(NULL)
  }
  if (!("data" %in% slotNames(mcreg_obj))) {
    warning("mcreg object has no data slot; cannot plot.")
    return(NULL)
  }
  plot_data <- mcreg_obj@data
  if (!all(c("x", "y") %in% names(plot_data))) {
    warning("mcreg data missing x/y columns; cannot plot.")
    return(NULL)
  }
  plot_data <- plot_data %>%
    filter(is.finite(x) & is.finite(y))
  if (nrow(plot_data) < 5) {
    warning("Insufficient data (<5 points) for mcreg plot.")
    return(NULL)
  }

  axis_lim <- range(c(plot_data$x, plot_data$y), na.rm = TRUE, finite = TRUE)
  axis_span <- diff(axis_lim)
  if (!is.finite(axis_span) || axis_span <= 0) {
    axis_span <- max(abs(axis_lim), na.rm = TRUE)
    if (!is.finite(axis_span) || axis_span == 0) {
      axis_span <- 1
    }
  }
  axis_pad <- axis_span * 0.04
  axis_lim <- axis_lim + c(-axis_pad, axis_pad)

  para <- mcreg_obj@para
  slope <- para["Slope", "EST"]
  intercept <- para["Intercept", "EST"]
  slope_l <- para["Slope", "LCI"]
  slope_u <- para["Slope", "UCI"]
  intercept_l <- para["Intercept", "LCI"]
  intercept_u <- para["Intercept", "UCI"]

  line_df <- tibble(x = seq(axis_lim[1], axis_lim[2], length.out = 200)) %>%
    mutate(
      est = intercept + slope * x,
      lci = intercept_l + slope_l * x,
      uci = intercept_u + slope_u * x
    )

  p <- ggplot(plot_data, aes(x = x, y = y)) +
    scatter_geom() +
    geom_abline(
      slope = 1,
      intercept = 0,
      color = PLOT_STYLE$colors$identity,
      linetype = "dashed",
      linewidth = PLOT_STYLE$line_width
    )
  if (is.finite(intercept) && is.finite(slope)) {
    p <- p +
      geom_line(
        data = line_df,
        aes(y = est),
        color = PLOT_STYLE$colors$fit,
        linewidth = PLOT_STYLE$line_width
      )
  }
  if (
    is.finite(intercept_l) &&
      is.finite(slope_l) &&
      is.finite(intercept_u) &&
      is.finite(slope_u)
  ) {
    p <- p +
      geom_line(
        data = line_df,
        aes(y = lci),
        color = PLOT_STYLE$colors$limits,
        linetype = "dotted",
        linewidth = PLOT_STYLE$line_width
      ) +
      geom_line(
        data = line_df,
        aes(y = uci),
        color = PLOT_STYLE$colors$limits,
        linetype = "dotted",
        linewidth = PLOT_STYLE$line_width
      )
  }

  p <- p +
    labs(x = x_label, y = y_label) +
    my_theme +
    coord_equal(xlim = axis_lim, ylim = axis_lim) +
    scale_x_continuous(labels = lblr) +
    scale_y_continuous(labels = lblr)

  ggsave(filename, plot = p, height = 6, width = 6, dpi = 300)
  return(p)
}
