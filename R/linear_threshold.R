#' Linear Regression with Threshold Effect Analysis
#'
#' This function performs a threshold effect analysis using linear regression models,
#' including restricted cubic spline (RCS) modeling for the focal variable,
#' threshold determination, and model comparison **adjusted for confounders**.
#' Regression coefficient+95%CI visualization with explicit threshold marking.
#' Outputs include Excel results and analysis objects.
#'
#' @param formula A formula specifying the model **including both focal variable (var_name) and confounders**.
#'               Example: `outcome ~ focal_var + confounder1 + confounder2`
#' @param data Data frame containing all variables in the formula
#' @param var_name Name of the focal variable to analyze for threshold effects
#' @param rcs_nodes Number of knots for RCS (default=4) applied to focal variable
#' @return List containing:
#'   - cutoff: Optimal threshold value
#'   - mdl0: Original model with all predictors
#'   - mdl1: Segmented model (X1/X2 for focal variable + confounders)
#'   - mdl2: Threshold test model (focal variable + X2 + confounders)
#'   - plot: Visualization object of coefficient curve with threshold marker
#' @export
linear_threshold <- function(formula, data, var_name, rcs_nodes = 4) {
  # ======================
  # INITIAL SETUP AND CHECKS
  # ======================
  # Package validation
  required_packages <- c("mgcv", "lmtest", "dplyr", "broom", "openxlsx", "splines", "ggplot2")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required but not installed"))
    }
  }

  # Load libraries
  library(mgcv)
  library(lmtest)
  library(dplyr)
  library(broom)
  library(openxlsx)
  library(splines)
  library(ggplot2)

  # ======================
  # INPUT VALIDATION
  # ======================
  # Verify response variable exists in data
  response_var <- as.character(formula[[2]])  # Extract response variable name
  if (!response_var %in% names(data)) {
    stop(paste("Response variable", response_var, "not found in data"))
  }

  # Verify target variable exists in data
  if (!var_name %in% names(data)) {
    stop(paste("Variable", var_name, "not found in data"))
  }
  data[[var_name]] <- as.numeric(data[[var_name]])

  # ======================
  # GAM SMOOTHING ANALYSIS
  # ======================
  use_rcs <- TRUE
  gam_model <- tryCatch({
    # Construct RCS term for spline modeling
    rcs_term <- paste0("ns(", var_name, ", ", rcs_nodes, ")")
    new_formula <- update(formula, as.formula(paste(". ~ . +", rcs_term)))

    # Fit generalized additive model (GAM) with RCS
    gam(new_formula, data = data, family = gaussian(link = "identity"))
  }, error = function(e) {
    warning("RCS failed, using linear model: ", e$message)
    use_rcs <- FALSE
    NULL
  })

  # ======================
  # THRESHOLD DETERMINATION
  # ======================
  cut_off <- NA
  tryCatch({
    # Two-stage threshold search algorithm
    search_threshold <- function(x, formula) {
      # Ensure x is numeric
      x <- as.numeric(x)

      # Add boundary protection
      if (length(unique(x)) < 5) {
        stop("Insufficient unique values for threshold analysis")
      }

      # Phase 1: Coarse search (5%-95% quantiles, 5% step)
      tmp.ss <- seq(0.05, 0.95, 0.05)
      tp <- quantile(x, probs = tmp.ss, na.rm = TRUE)
      tmp.llk <- numeric(length(tmp.ss))

      # Calculate log-likelihood for coarse search
      for (k in seq_along(tmp.ss)) {
        tmp.X <- (x > tp[k]) * (x - tp[k])
        new_data <- cbind(data, tmp.X)
        new_formula <- update(formula, . ~ . + tmp.X)
        tmp.mdl <- glm(new_formula, data = new_data, family = gaussian)
        tmp.llk[k] <- logLik(tmp.mdl)
      }

      # Define fine search range around best coarse candidate
      tp1_idx <- which.max(tmp.llk)
      tp2.min <- max(0.05, tmp.ss[tp1_idx] - 0.04)
      tp2.max <- min(0.95, tmp.ss[tp1_idx] + 0.04)

      # Phase 2: Fine search within reduced range
      tp.pctlrange <- quantile(x, probs = c(tp2.min, tp2.max), na.rm = TRUE)
      tp.range <- unique(x[x > tp.pctlrange[1] & x < tp.pctlrange[2]])

      # Iterative range refinement
      while (length(tp.range) > 5) {
        tmp.pct3 <- quantile(tp.range, probs = c(0, 0.25, 0.5, 0.75, 1), type = 3)
        tmp.llk3 <- numeric(3)

        for (k in 2:4) {
          tmp.X <- (x > tmp.pct3[k]) * (x - tmp.pct3[k])
          new_data <- cbind(data, tmp.X)
          new_formula <- update(formula, . ~ . + tmp.X)
          tmp.mdl <- glm(new_formula, data = new_data, family = gaussian)
          tmp.llk3[k-1] <- logLik(tmp.mdl)
        }

        best_idx <- which.max(tmp.llk3)
        tp.range <- tp.range[tp.range >= tmp.pct3[best_idx] & tp.range <= tmp.pct3[best_idx+2]]
      }

      # Final threshold determination
      if (length(tp.range) > 0) {
        final.llk <- numeric(length(tp.range))
        for (k in seq_along(tp.range)) {
          tmp.X <- (x > tp.range[k]) * (x - tp.range[k])
          new_data <- cbind(data, tmp.X)
          new_formula <- update(formula, . ~ . + tmp.X)
          tmp.mdl <- glm(new_formula, data = new_data, family = gaussian)
          final.llk[k] <- logLik(tmp.mdl)
        }
        tp.val <- tp.range[which.max(final.llk)]
      } else {
        tp.val <- tp.pctlrange[1]
      }

      round(tp.val, 2)
    }

    cut_off <- search_threshold(data[[var_name]], formula)
    # Print threshold result to screen
    cat("\n===== Threshold Analysis Result =====\n")
    cat("Optimal threshold value for", var_name, ":", cut_off, "\n\n")

  }, error = function(e) {
    stop("Threshold determination failed: ", e$message)
  })

  # ======================
  # SEGMENTED VARIABLE CREATION
  # ======================
  # Create segmented variables X1 and X2
  data <- data %>%
    mutate(
      X1 = pmax(0, cut_off - .data[[var_name]]),
      X2 = pmax(0, .data[[var_name]] - cut_off)
    )

  # ======================
  # MODEL FITTING (FIXED FOR SINGLE VARIABLE CASE AND X2 DISPLAY ISSUE)
  # ======================
  mdl0 <- mdl1 <- mdl2 <- NULL

  tryCatch({
    # Extract all predictor variables from formula
    predictors <- all.vars(formula)[-1]  # Exclude response variable
    # Identify variables other than the threshold variable
    other_vars <- setdiff(predictors, var_name)

    # Model 0: Original linear regression model
    mdl0 <- glm(formula, data = data, family = gaussian)

    # Model 1: Model with segmented variables X1 and X2 (REPLACE TARGET VARIABLE COMPLETELY)
    if (length(other_vars) > 0) {
      # Create new formula without the target variable but with segmented variables
      base_formula <- as.formula(paste(response_var, "~", paste(c("X1", "X2", other_vars), collapse = " + ")))
      mdl1 <- glm(base_formula, data = data, family = gaussian)
    } else {
      # Single variable case: explicit formula without original variable
      mdl1 <- glm(
        as.formula(paste(response_var, "~ X1 + X2")),
        data = data,
        family = gaussian
      )
    }

    # Model 2: Model with X2 and original variable (KEEP ORIGINAL VARIABLE)
    if (length(other_vars) > 0) {
      # Case with additional variables: explicit inclusion of all predictors + X2
      # Reorder variables: X2 first, then original variable, then other variables
      new_formula <- as.formula(paste(response_var, "~", paste(c("X2", var_name, other_vars), collapse = " + ")))
      mdl2 <- glm(new_formula, data = data, family = gaussian)
    } else {
      # Single variable case: explicit formula with response variable
      mdl2 <- glm(
        as.formula(paste(response_var, "~ X2 +", var_name)),
        data = data,
        family = gaussian
      )
    }

  }, error = function(e) {
    stop("Model fitting failed: ", e$message)
  })

  # ======================
  # PLOT PREPARATION - beta + 95%CI CURVE
  # ======================
  # Create prediction dataframe for visualization
  pred_x <- seq(min(data[[var_name]], na.rm = TRUE),
                max(data[[var_name]], na.rm = TRUE),
                length.out = 100)

  # Create new data for predictions
  newdata <- data.frame(pred_x)
  names(newdata) <- var_name

  if (length(other_vars) > 0) {
    # Set other variables to median values
    newdata[other_vars] <- lapply(data[other_vars], function(x) median(x, na.rm = TRUE))
  }

  # Generate predictions if RCS was successful
  if (use_rcs && !is.null(gam_model)) {
    pred <- predict(gam_model, newdata = newdata, type = "terms", se.fit = TRUE)
    rcs_index <- grep(paste0("^ns\\(", var_name), colnames(pred$fit))

    pred_df <- data.frame(
      x = newdata[[var_name]],
      beta = pred$fit[, rcs_index],  # Linear predictor (beta)
      se = pred$se.fit[, rcs_index],
      low = pred$fit[, rcs_index] - 1.96*pred$se.fit[, rcs_index],
      up = pred$fit[, rcs_index] + 1.96*pred$se.fit[, rcs_index]
    )
  } else {
    # Fallback to linear model if RCS failed
    linear_model <- lm(formula, data = data)
    pred_df <- data.frame(
      x = newdata[[var_name]],
      beta = predict(linear_model, newdata = newdata),
      low = NA,
      up = NA
    )
  }

  # ======================
  # PLOTTING MODULE (beta + 95%CI CURVE)
  # ======================
  # Create beta + 95%CI visualization
  p <- ggplot(pred_df, aes(x = x)) +
    # beta curve (blue solid line)
    geom_line(aes(y = beta), color = "#2C7FBF", linewidth = 1.2) +
    # 95% CI bands (if available)
    geom_ribbon(aes(ymin = low, ymax = up),
                fill = "#3182BD", alpha = 0.2) +
    # Threshold line (red dashed)
    geom_vline(xintercept = cut_off, color = "#E63946",
               linetype = "dashed", linewidth = 1.2) +
    # Rug plot (bottom only)
    geom_rug(aes(x = .data[[var_name]]),
             data = data, sides = "b", alpha = 0.5) +
    # Annotations
    labs(x = var_name,
         y = expression(paste("beta (95% CI)")),  # 替换为ASCII字符
         title = paste("Threshold Effect Analysis for", var_name),
         subtitle = paste("Optimal cutoff:", cut_off),
         caption = "Blue line: beta coefficient; Shaded area: 95% CI; Dashed line: Threshold") +
    # Theme customization
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "#212121", linewidth = 0.5),
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "#E63946"),
      plot.caption = element_text(hjust = 0.5, face = "italic"),
      axis.title = element_text(face = "bold", size = 14),
      axis.text = element_text(size = 12)
    ) +
    # Y-axis scaling with explicit ticks
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.15)),
      breaks = scales::extended_breaks(n = 6)
    )

  # Print the plot
  print(p)

  # ======================
  # MODEL COMPARISON
  # ======================
  # Perform likelihood ratio test between models
  lrtest_result <- tryCatch({
    # Ensure models exist
    if (is.null(mdl0) || is.null(mdl1)) {
      stop("mdl0 or mdl1 not available for comparison")
    }

    # Execute likelihood ratio test
    lr_result <- lrtest(mdl0, mdl1)

    # Create results dataframe
    data.frame(
      Test = "mdl0 vs mdl1",
      Chisq = lr_result$Chisq[2],
      Df = lr_result$Df[2],
      Pr_Chisq = lr_result$`Pr(>Chisq)`[2],
      stringsAsFactors = FALSE
    )

  }, error = function(e) {
    # Return structured NA result on error
    data.frame(
      Test = "mdl0 vs mdl1",
      Chisq = NA_real_,
      Df = NA_real_,
      Pr_Chisq = NA_real_,
      stringsAsFactors = FALSE
    )
  })

  # ======================
  # RESULTS EXTRACTION
  # ======================
  # Function to extract model results with coefficients and confidence intervals
  extract_var_results <- function(model) {
    tidy_result <- tryCatch({
      tidy_model <- tidy(model, conf.int = TRUE)

      # Handle cases where coefficients might be missing
      if (nrow(tidy_model) == 0) {
        warning("Model returned empty coefficients")
        return(data.frame(
          Variable = NA_character_,
          Beta = NA_real_,
          Lower_CI = NA_real_,
          Upper_CI = NA_real_,
          P_value = NA_real_,
          stringsAsFactors = FALSE
        ))
      }

      # Rename columns for clarity
      tidy_model %>%
        rename(
          Variable = term,
          Beta = estimate,
          Lower_CI = conf.low,
          Upper_CI = conf.high,
          P_value = p.value
        ) %>%
        as.data.frame()
    }, error = function(e) {
      warning("Error in tidying model: ", e$message)
      data.frame(
        Variable = NA_character_,
        Beta = NA_real_,
        Lower_CI = NA_real_,
        Upper_CI = NA_real_,
        P_value = NA_real_,
        stringsAsFactors = FALSE
      )
    })

    return(tidy_result)
  }

  # ======================
  # RESULTS OUTPUT
  # ======================
  # Console output of results
  cat("[Model 0 Results - Original model]\n")
  print(extract_var_results(mdl0))
  cat("\n[Model 1 Results - Model with segmented variables (X1 + X2)]\n")
  print(extract_var_results(mdl1))
  cat("\n[Model 2: Model with X2 and original variable (testing for threshold effect)]\n")
  print(extract_var_results(mdl2))
  cat("\n[Likelihood Ratio Test]\n")
  print(lrtest_result)

  # Excel output of results
  wb <- createWorkbook()
  addWorksheet(wb, "Cutoff")
  addWorksheet(wb, "mdl0_Results")
  addWorksheet(wb, "mdl1_Results")
  addWorksheet(wb, "mdl2_Results")
  addWorksheet(wb, "LRTest")

  writeData(wb, "Cutoff", data.frame(Variable = var_name, Cutoff = cut_off))
  writeData(wb, "mdl0_Results", extract_var_results(mdl0))
  writeData(wb, "mdl1_Results", extract_var_results(mdl1))
  writeData(wb, "mdl2_Results", extract_var_results(mdl2))
  writeData(wb, "LRTest", lrtest_result)

  saveWorkbook(wb, "linear_threshold_results.xlsx", overwrite = TRUE)
  cat("\nResults saved to: linear_threshold_results.xlsx\n")

  # ======================
  # RETURN RESULTS
  # ======================
  # Return analysis objects as list
  list(
    cutoff = cut_off,
    gam_model = if (!is.null(gam_model)) gam_model else NULL,
    mdl0 = mdl0,
    mdl1 = mdl1,
    mdl2 = mdl2,
    mdl0_results = extract_var_results(mdl0),
    mdl1_results = extract_var_results(mdl1),
    mdl2_results = extract_var_results(mdl2),
    data = data,
    pred_data = pred_df,
    plot = p,
    lrtest_result = lrtest_result
  )
}
