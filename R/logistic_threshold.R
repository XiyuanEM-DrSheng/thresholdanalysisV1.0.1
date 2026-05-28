#' Logistic Regression with Threshold Effect Analysis
#'
#' This function performs a threshold effect analysis using logistic regression,
#' including restricted cubic spline (RCS) modeling for the focal variable,
#' threshold determination, and model comparison **adjusted for confounders**.
#' OR+95%CI visualization with explicit threshold marking.
#' Outputs include Excel results and analysis objects.
#'
#' @param formula A formula specifying the model **including both focal variable (var_name) and confounders**.
#'               Example: `outcome ~ focal_var + confounder1 + confounder2`
#' @param data Data frame containing all variables in the formula
#' @param var_name Name of the focal variable to analyze for threshold effects
#' @param rcs_nodes Number of knots for RCS (default=4) applied to focal variable
#' @return List containing:
#'   - cutoff: Optimal threshold value
#'   - mdl0: Original model with all predictors (focal + confounders)
#'   - mdl1: Segmented model (X1/X2 for focal variable + confounders)
#'   - mdl2: Threshold test model (focal variable + X2 + confounders)
#'   - plot: Visualization object of OR curve with threshold marker
#' @export
logistic_threshold <- function(formula, data, var_name, rcs_nodes = 4) {
  # ======================
  # PACKAGE & DEPENDENCY CHECK
  # ======================
  required_packages <- c("mgcv", "lmtest", "dplyr", "broom", "openxlsx",
                         "splines", "ggplot2", "purrr", "readr")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "required but not installed"))
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
  library(purrr)
  library(readr)

  # ======================
  # INPUT VALIDATION
  # ======================
  # Response variable extraction
  response_var <- as.character(formula[[2]])
  if (!response_var %in% names(data)) {
    stop(paste("Response variable", response_var, "not found in data"))
  }

  # Threshold variable validation
  if (!var_name %in% names(data)) {
    stop(paste("Variable", var_name, "not found in data"))
  }
  data[[var_name]] <- as.numeric(data[[var_name]])

  # Predictor variables identification
  predictors <- setdiff(all.vars(formula), response_var)
  other_vars <- setdiff(predictors, var_name)

  # ======================
  # RCS SMOOTHING ANALYSIS
  # ======================
  use_rcs <- TRUE
  gam_model <- tryCatch({
    # Construct RCS formula
    rcs_term <- paste0("ns(", var_name, ", ", rcs_nodes, ")")

    # Formula construction
    rcs_formula <- if (length(other_vars) > 0) {
      as.formula(paste(response_var, "~", paste(c(other_vars, rcs_term), collapse = " + ")))
    } else {
      as.formula(paste(response_var, "~", rcs_term))
    }

    # Fit GAM model
    gam_model <- gam(rcs_formula, data = data, family = binomial)

    # Generate predictions
    pred <- predict(gam_model, type = "terms", se.fit = TRUE)
    rcs_index <- grep(paste0("^ns\\(", var_name), colnames(pred$fit))

    if (length(rcs_index) == 0) stop("RCS term not found in model")

    # Add predictions to data
    data <- data %>%
      mutate(
        !!paste0(var_name, ".fit") := plogis(pred$fit[, rcs_index] + coef(gam_model)[1]),
        !!paste0(var_name, ".low") := plogis(pred$fit[, rcs_index] - 1.96*pred$se.fit[, rcs_index]),
        !!paste0(var_name, ".up") := plogis(pred$fit[, rcs_index] + 1.96*pred$se.fit[, rcs_index])
      )

    gam_model

  }, error = function(e) {
    warning("RCS analysis failed: ", e$message, "\nProceeding with linear model")
    use_rcs <- FALSE
    NULL
  })

  # ======================
  # THRESHOLD DETERMINATION
  # ======================
  cut_off <- NA
  tryCatch({
    # Prepare data for threshold search
    model_data <- data[, c(response_var, predictors), drop = FALSE]

    # Threshold search algorithm
    search_threshold <- function(x, data, response_var) {
      # Data cleaning
      data <- data %>% filter(!is.na(.data[[var_name]]))
      x <- data[[var_name]]
      predictor_vars <- setdiff(colnames(data), response_var)

      # Coarse search
      tmp.ss <- seq(0.05, 0.95, 0.05)
      tp <- quantile(x, probs = tmp.ss)
      tmp.llk <- numeric(length(tmp.ss))

      for (k in seq_along(tmp.ss)) {
        tmp.X <- ifelse(x > tp[k], x - tp[k], 0)
        temp_data <- data
        temp_data$tmp.X <- tmp.X

        # Construct formula
        new_formula <- as.formula(paste(
          response_var, "~",
          paste(c(predictor_vars, "tmp.X"), collapse = " + ")
        ))

        # Fit model
        tmp.mdl <- glm(new_formula, data = temp_data, family = binomial)
        tmp.llk[k] <- logLik(tmp.mdl)
      }

      # Fine search
      tp1_idx <- which.max(tmp.llk)
      tp2.min <- max(0.05, tmp.ss[tp1_idx] - 0.04)
      tp2.max <- min(0.95, tmp.ss[tp1_idx] + 0.04)
      tp.pctlrange <- quantile(x, probs = c(tp2.min, tp2.max))
      tp.range <- unique(x[x > tp.pctlrange[1] & x < tp.pctlrange[2]])

      # Iterative refinement
      while (length(tp.range) > 5) {
        tmp.pct3 <- quantile(tp.range, probs = c(0, 0.25, 0.5, 0.75, 1), type = 3)
        tmp.llk3 <- numeric(3)

        for (k in 2:4) {
          tmp.X <- ifelse(x > tmp.pct3[k], x - tmp.pct3[k], 0)
          temp_data <- data
          temp_data$tmp.X <- tmp.X

          # Construct formula
          new_formula <- as.formula(paste(
            response_var, "~",
            paste(c(predictor_vars, "tmp.X"), collapse = " + ")
          ))

          # Fit model
          tmp.mdl <- glm(new_formula, data = temp_data, family = binomial)
          tmp.llk3[k-1] <- logLik(tmp.mdl)
        }

        best_idx <- which.max(tmp.llk3)
        tp.range <- tp.range[tp.range >= tmp.pct3[best_idx] & tp.range <= tmp.pct3[best_idx+2]]
      }

      # Final threshold determination
      if (length(tp.range) > 0) {
        final.llk <- numeric(length(tp.range))
        for (k in seq_along(tp.range)) {
          tmp.X <- ifelse(x > tp.range[k], x - tp.range[k], 0)
          temp_data <- data
          temp_data$tmp.X <- tmp.X

          # Construct formula
          new_formula <- as.formula(paste(
            response_var, "~",
            paste(c(predictor_vars, "tmp.X"), collapse = " + ")
          ))

          # Fit model
          tmp.mdl <- glm(new_formula, data = temp_data, family = binomial)
          final.llk[k] <- logLik(tmp.mdl)
        }
        tp.val <- tp.range[which.max(final.llk)]
      } else {
        tp.val <- tp.pctlrange[1]
      }

      round(tp.val, 2)
    }

    # Execute threshold search
    cut_off <- search_threshold(model_data[[var_name]], model_data, response_var)
    cat("\n===== Threshold Analysis Result =====\n")
    cat("Optimal threshold value for", var_name, ":", cut_off, "\n\n")

  }, error = function(e) {
    stop("Threshold determination failed: ", e$message)
  })

  # ======================
  # PLOT PREPARATION
  # ======================
  # Create prediction dataframe for visualization
  pred_x <- seq(min(data[[var_name]], na.rm = TRUE),
                max(data[[var_name]], na.rm = TRUE),
                length.out = 100)

  # Create new data for predictions
  newdata <- data.frame(pred_x)
  names(newdata) <- var_name  # Set column name to the variable name

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
      fit = plogis(pred$fit[, rcs_index] + coef(gam_model)[1]),
      se = pred$se.fit[, rcs_index],
      low = plogis(pred$fit[, rcs_index] + coef(gam_model)[1] - 1.96*pred$se.fit[, rcs_index]),
      up = plogis(pred$fit[, rcs_index] + coef(gam_model)[1] + 1.96*pred$se.fit[, rcs_index])
    )
  } else {
    # Fallback to linear model if RCS failed
    linear_model <- glm(formula, data = data, family = binomial)
    pred_df <- data.frame(
      x = newdata[[var_name]],
      fit = predict(linear_model, newdata = newdata, type = "response"),
      low = NA,
      up = NA
    )
  }

  # ======================
  # PLOTTING MODULE (OR+95%CI CURVE)
  # ======================
  # Create OR+95%CI visualization
  p <- ggplot(pred_df, aes(x = x)) +
    # OR curve
    geom_line(aes(y = fit), color = "#2C7FBF", linewidth = 1.2) +
    # 95% CI bands (if available)
    geom_ribbon(aes(ymin = low, ymax = up),
                fill = "#3182BD", alpha = 0.2) +
    # Threshold line
    geom_vline(xintercept = cut_off, color = "#E63946",
               linetype = "dashed", linewidth = 1.2) +
    # Data points rug plot (bottom only)
    geom_rug(aes(x = .data[[var_name]]),
             data = data, sides = "b", alpha = 0.5) +
    # Annotations
    labs(x = var_name,
         y = "Odds Ratio (OR)",
         title = paste("Threshold Effect Analysis for", var_name),
         subtitle = paste("Optimal cutoff:", cut_off),
         caption = "Blue line: OR curve; Shaded area: 95% CI; Dashed line: Threshold") +
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
      # Set explicit tick marks for better readability
      breaks = scales::extended_breaks(n = 6)  # Show 6 ticks on y-axis
    )

  # Print the plot
  print(p)

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
  # MODEL FITTING
  # ======================
  tryCatch({
    # Model 0: Original model
    mdl0 <- glm(formula, data = data, family = binomial)

    # Model 1: Segmented variables model
    if (length(other_vars) > 0) {
      mdl1_formula <- as.formula(paste(
        response_var, "~",
        paste(c("X1", "X2", setdiff(other_vars, var_name)), collapse = " + ")
      ))
    } else {
      mdl1_formula <- as.formula(paste(response_var, "~ X1 + X2"))
    }
    mdl1 <- glm(mdl1_formula, data = data, family = binomial)

    # Model 2: Threshold effect test model
    if (length(other_vars) > 0) {
      mdl2_formula <- as.formula(paste(
        response_var, "~",
        paste(c("X2", var_name, setdiff(other_vars, var_name)), collapse = " + ")
      ))
    } else {
      mdl2_formula <- as.formula(paste(response_var, "~ X2 +", var_name))
    }
    mdl2 <- glm(mdl2_formula, data = data, family = binomial)

  }, error = function(e) {
    stop("Model fitting failed: ", e$message, "\nFormula used:\n",
         paste(response_var, "~", paste(c("X1", "X2", other_vars), collapse = " + ")))
  })

  # ======================
  # MODEL COMPARISON
  # ======================
  lrtest_result <- tryCatch({
    lrtest(mdl0, mdl1)
  }, error = function(e) {
    data.frame(
      Test = "mdl0 vs mdl1",
      Chisq = NA_real_,
      Df = NA_real_,
      Pr_Chisq = NA_real_
    )
  })

  # ======================
  # RESULTS EXTRACTION
  # ======================
  # Function to extract model results
  extract_var_results <- function(model) {
    tidy_result <- tryCatch({
      tidy_model <- tidy(model, exponentiate = TRUE, conf.int = TRUE)
      if (nrow(tidy_model) == 0) {
        warning("Model returned empty coefficients")
        return(data.frame(
          Variable = NA_character_,
          OR = NA_real_,
          Lower_CI = NA_real_,
          Upper_CI = NA_real_,
          P_value = NA_real_
        ))
      }
      rename(tidy_model,
             Variable = term,
             OR = estimate,
             Lower_CI = conf.low,
             Upper_CI = conf.high,
             P_value = p.value
      )
    }, error = function(e) {
      data.frame(
        Variable = NA_character_,
        OR = NA_real_,
        Lower_CI = NA_real_,
        Upper_CI = NA_real_,
        P_value = NA_real_
      )
    })
    return(tidy_result)
  }

  # Extract model results
  mdl0_results <- extract_var_results(mdl0)
  mdl1_results <- extract_var_results(mdl1)
  mdl2_results <- extract_var_results(mdl2)

  # ======================
  # RESULTS OUTPUT
  # ======================
  # Create workbook
  wb <- createWorkbook()
  addWorksheet(wb, "Cutoff")
  addWorksheet(wb, "mdl0_Results")
  addWorksheet(wb, "mdl1_Results")
  addWorksheet(wb, "mdl2_Results")
  addWorksheet(wb, "LRTest")

  # Write data
  writeData(wb, "Cutoff", data.frame(Variable = var_name, Cutoff = cut_off))
  writeData(wb, "mdl0_Results", mdl0_results)
  writeData(wb, "mdl1_Results", mdl1_results)
  writeData(wb, "mdl2_Results", mdl2_results)
  writeData(wb, "LRTest", lrtest_result)

  # Save workbook
  saveWorkbook(wb, "logistic_threshold_results.xlsx", overwrite = TRUE)
  cat("\nResults saved to: logistic_threshold_results.xlsx\n")

  # ======================
  # RETURN RESULTS
  # ======================
  # Return analysis objects
  list(
    cutoff = cut_off,
    mdl0 = mdl0,
    mdl1 = mdl1,
    mdl2 = mdl2,
    mdl0_results = mdl0_results,
    mdl1_results = mdl1_results,
    mdl2_results = mdl2_results,
    data = data,
    pred_data = pred_df,
    plot = p,
    lrtest_result = lrtest_result
  )
}
