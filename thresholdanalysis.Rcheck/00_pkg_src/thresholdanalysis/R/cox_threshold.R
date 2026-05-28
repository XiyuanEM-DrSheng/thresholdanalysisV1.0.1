#' Cox Proportional Hazards Model with Threshold Effect Analysis
#'
#' This function performs a threshold effect analysis using Cox proportional hazards models,
#' including restricted cubic spline (RCS) modeling for the focal variable,
#' threshold determination, and model comparison **adjusted for confounders**.
#' HR+95%CI visualization with explicit threshold marking.
#' Outputs include Excel results and analysis objects.
#'
#' @param formula A formula specifying the model **including both focal variable (var_name) and confounders**.
#'               Example: `Surv(time, status) ~ focal_var + confounder1 + confounder2`
#' @param data Data frame containing all variables in the formula
#' @param var_name Name of the focal variable to analyze for threshold effects
#' @param rcs_nodes Number of knots for RCS (default=4) applied to focal variable
#' @return List containing:
#'   - cutoff: Optimal threshold value
#'   - mdl0: Original model with all predictors (focal + confounders)
#'   - mdl1: Segmented model (X1/X2 for focal variable + confounders)
#'   - mdl2: Threshold test model (focal variable + X2 + confounders)
#'   - plot: Visualization object of HR curve with threshold marker
#' @export
cox_threshold <- function(formula, data, var_name, rcs_nodes = 4) {
  # ======================
  # INITIAL SETUP AND CHECKS
  # ======================
  required_packages <- c("rms", "survival", "mgcv", "lmtest", "dplyr", "broom", "openxlsx", "ggplot2")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required but not installed."))
    }
  }
  library(rms); library(survival); library(mgcv); library(lmtest)
  library(dplyr); library(broom); library(openxlsx); library(ggplot2)

  # ======================
  # INPUT VALIDATION
  # ======================
  if (!var_name %in% names(data)) stop(paste("Variable", var_name, "not found in data"))
  surv_obj <- formula[[2]]
  time_var <- as.character(surv_obj[[2]])
  status_var <- as.character(surv_obj[[3]])
  if (!time_var %in% names(data) || !status_var %in% names(data)) {
    stop("Survival variables not found in data")
  }
  data[[var_name]] <- as.numeric(data[[var_name]])

  # ======================
  # MODEL FORMULA PREPARATION
  # ======================
  other_terms <- setdiff(all.vars(formula)[-c(1:2)], var_name)
  rcs_term <- paste0("rcs(", var_name, ", ", rcs_nodes, ")")

  # Construct RCS formula with conditional handling
  new_formula <- if (length(other_terms) > 0) {
    reformulate(c(rcs_term, other_terms),
                response = paste0("Surv(", time_var, ", ", status_var, ")"))
  } else {
    reformulate(rcs_term,
                response = paste0("Surv(", time_var, ", ", status_var, ")"))
  }

  # ======================
  # MODEL FITTING
  # ======================
  old_dd <- options("datadist")
  dd <<- rms::datadist(data)
  options(datadist = "dd")

  cox_model <- tryCatch({
    fit <- rms::cph(new_formula, data = data, surv = TRUE, x = TRUE, y = TRUE)
    if (!any(grepl(paste0("^rcs\\("), attr(terms(fit), "term.labels")))) stop("rcs term missing")
    fit
  }, error = function(e) {
    warning("Using coxph as fallback: ", e$message)
    coxph(new_formula, data = data, x = TRUE, y = TRUE)
  })

  on.exit({
    options(datadist = old_dd$datadist)
    if (exists("dd", envir = .GlobalEnv)) rm(dd, envir = .GlobalEnv)
  })

  # ======================
  # PREDICTION PROCESSING
  # ======================
  pred <- tryCatch({
    if (inherits(cox_model, "cph")) {
      pred_df <- rms::Predict(cox_model, name = var_name, fun = exp, conf.int = 0.95)
      list(fit = pred_df$yhat,
           se.fit = (pred_df$upper - pred_df$lower)/(2*qnorm(0.975)),
           x = pred_df[[var_name]])
    } else {
      pred_data <- data.frame(seq(min(data[[var_name]]), max(data[[var_name]]), length.out = 100))
      names(pred_data) <- var_name
      if (length(other_terms) > 0) {
        for (term in other_terms) {
          pred_data[[term]] <- median(data[[term]], na.rm = TRUE)
        }
      }
      pred_result <- predict(cox_model, newdata = pred_data, type = "terms", se.fit = TRUE)
      rcs_pattern <- paste0("^rcs\\(", var_name, ".*\\)")
      var_term <- grep(rcs_pattern, attr(pred_result$fit, "dimnames")[[2]], value = TRUE)
      list(fit = pred_result$fit[,var_term],
           se.fit = pred_result$se.fit[,var_term],
           x = pred_data[[var_name]])
    }
  }, error = function(e) stop("Prediction failed: ", e$message))

  pred_df <- data.frame(
    x = pred$x,
    fit = pred$fit,
    se = pred$se.fit,
    low = pred$fit - 1.96*pred$se.fit,
    up = pred$fit + 1.96*pred$se.fit
  )
  names(pred_df) <- c(var_name, "fit", "se", "low", "up")

  # ======================
  # THRESHOLD DETERMINATION
  # ======================
  cut_off <- NA
  tryCatch({
    find_threshold <- function(data, var, time, status, formula) {
      xTMP <- data[[var]]
      tmp.ss <- seq(0.05, 0.95, 0.05)
      tp <- quantile(xTMP, probs = tmp.ss, na.rm = TRUE)
      tmp.llk <- rep(NA, length(tmp.ss))
      fml_str <- paste(deparse(formula), "+tmp.X")

      for (k in seq_along(tmp.ss)) {
        tmp.X <- (xTMP > tp[k])*(xTMP - tp[k])
        wdtmp1 <- cbind(data, tmp.X)
        fit <- tryCatch(
          coxph(formula(fml_str), data = wdtmp1),
          error = function(e) NULL
        )
        if (!is.null(fit)) tmp.llk[k] <- fit$loglik[2]
        rm(wdtmp1, tmp.X)
      }

      tp1 <- tmp.ss[which.max(tmp.llk)]
      tp2.min <- max(0.05, tp1 - 0.04)
      tp2.max <- min(0.95, tp1 + 0.04)

      tp.pctlrange <- quantile(xTMP, probs = c(tp2.min, tp2.max), na.rm = TRUE)
      tp.range <- unique(xTMP[xTMP > tp.pctlrange[1] & xTMP < tp.pctlrange[2]])

      while (length(tp.range) > 5) {
        tmp.pct3 <- quantile(tp.range, probs = c(0, 0.25, 0.5, 0.75, 1), type = 3)
        tmp.llk3 <- rep(NA, 3)

        for (k in 2:4) {
          tmp.X <- (xTMP > tmp.pct3[k])*(xTMP - tmp.pct3[k])
          wdtmp1 <- cbind(data, tmp.X)
          fit <- tryCatch(
            coxph(formula(fml_str), data = wdtmp1),
            error = function(e) NULL
          )
          if (!is.null(fit)) tmp.llk3[k-1] <- fit$loglik[2]
          rm(wdtmp1, tmp.X)
        }

        tmp.min3 <- which.max(tmp.llk3)
        tp.range <- tp.range[tp.range >= tmp.pct3[tmp.min3] & tp.range <= tmp.pct3[tmp.min3+2]]
      }

      if (length(tp.range) > 0) {
        tmp.llk <- rep(NA, length(tp.range))
        for (k in seq_along(tp.range)) {
          tmp.X <- (xTMP > tp.range[k])*(xTMP - tp.range[k])
          wdtmp1 <- cbind(data, tmp.X)
          fit <- tryCatch(
            coxph(formula(fml_str), data = wdtmp1),
            error = function(e) NULL
          )
          if (!is.null(fit)) tmp.llk[k] <- fit$loglik[2]
          rm(wdtmp1, tmp.X)
        }
        tp.val <- tp.range[which.max(tmp.llk)]
      } else {
        tp.val <- tp.pctlrange[1]
      }
      return(round(tp.val, 2))
    }

    cut_off <- find_threshold(
      data = data,
      var = var_name,
      time = time_var,
      status = status_var,
      formula = formula
    )

    # Print threshold result to screen
    cat("\n===== Threshold Analysis Result =====\n")
    cat("Optimal threshold value for", var_name, ":", cut_off, "\n\n")

    if (is.na(cut_off) || !is.numeric(cut_off)) stop("Invalid threshold")
  }, error = function(e) stop("Threshold failed: ", e$message))

  # ======================
  # SEGMENTED VARIABLE CREATION
  # ======================
  if (is.na(cut_off)) stop("Invalid threshold")
  data <- data %>% mutate(
    X1 = pmax(0, cut_off - .data[[var_name]]),
    X2 = pmax(0, .data[[var_name]] - cut_off)
  )

  # ======================
  # MODEL COMPARISON
  # ======================
  tryCatch({
    # Model 0: Original model
    mdl0 <- coxph(formula, data = data)
    response_var <- paste0("Surv(", time_var, ", ", status_var, ")")

    # Model 1: Segmented model with conditional formula
    if (length(other_terms) > 0) {
      mdl1_formula <- as.formula(paste(
        response_var, "~",
        paste(c("X1", "X2", other_terms), collapse = " + ")
      ))
    } else {
      mdl1_formula <- as.formula(paste(
        response_var, "~ X1 + X2"
      ))
    }
    mdl1 <- coxph(mdl1_formula, data = data)

    # Model 2: Threshold effect test model with conditional formula
    if (length(other_terms) > 0) {
      mdl2_formula <- as.formula(paste(
        response_var, "~",
        paste(c("X2", var_name, other_terms), collapse = " + ")
      ))
    } else {
      mdl2_formula <- as.formula(paste(
        response_var, "~ X2 +", var_name
      ))
    }
    mdl2 <- coxph(mdl2_formula, data = data)

  }, error = function(e) {
    stop("Model fitting failed: ", e$message, "\nFormula used:\n",
         paste(response_var, "~", paste(c("X1", "X2", other_terms), collapse = " + ")))
  })

  # ======================
  # RESULTS REPORTING
  # ======================
  wb <- createWorkbook()
  addWorksheet(wb, "Cutoff"); writeData(wb, "Cutoff", data.frame(Variable = var_name, Cutoff = cut_off))

  extract_results <- function(model) {
    tidy_model <- tidy(model, conf.int = TRUE, exponentiate = TRUE)
    data.frame(
      Variable = tidy_model$term,
      HR = tidy_model$estimate,
      Lower_CI = tidy_model$conf.low,
      Upper_CI = tidy_model$conf.high,
      P_value = tidy_model$p.value
    )
  }

  # Extract and print model results
  mdl0_results <- extract_results(mdl0)
  mdl1_results <- extract_results(mdl1)
  mdl2_results <- extract_results(mdl2)

  # Print results to screen
  cat("===== Model Comparison Results =====\n")
  cat("\n--- Model 0 Results - Original model ---\n")
  print(mdl0_results)

  cat("\n--- Model 1 Results - Model with segmented variables ---\n")
  print(mdl1_results)

  cat("\n--- Model with X2 and original variable (testing for threshold effect) ---\n")
  print(mdl2_results)

  # Likelihood ratio test
  lrtest_result <- lrtest(mdl0, mdl1)
  cat("\n--- Likelihood Ratio Test (Model 0 vs Model 1) ---\n")
  print(data.frame(
    Test = "Model 0 vs Model 1",
    Chisq = lrtest_result$Chisq[2],
    DF = lrtest_result$Df[2],
    P_value = lrtest_result$`Pr(>Chisq)`[2]
  ))

  # Write to Excel
  addWorksheet(wb, "mdl0"); writeData(wb, "mdl0", mdl0_results)
  addWorksheet(wb, "mdl1"); writeData(wb, "mdl1", mdl1_results)
  addWorksheet(wb, "mdl2"); writeData(wb, "mdl2", mdl2_results)
  addWorksheet(wb, "Likelihood_Ratio_Test"); writeData(wb, "Likelihood_Ratio_Test", lrtest_result)

  saveWorkbook(wb, "cox_threshold_results.xlsx", overwrite = TRUE)
  cat("\nResults saved to: cox_threshold_results.xlsx\n")

  # ======================
  # PLOTTING SECTION
  # ======================
  # Modified plot to match logistic_threshold style
  p <- ggplot(pred_df, aes(x = .data[[var_name]])) +
    # HR curve (blue solid line)
    geom_line(aes(y = fit), color = "#2C7FBF", linewidth = 1.2) +
    # 95% CI ribbon (light blue with transparency)
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
         y = "Hazard Ratio (HR)",
         title = paste("Threshold Effect Analysis for", var_name),
         subtitle = paste("Optimal cutoff:", cut_off),
         caption = "Blue line: HR curve; Shaded area: 95% CI; Dashed line: Threshold") +
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
      breaks = scales::extended_breaks(n = 6)  # Show 6 ticks on y-axis
    )

  # Print the plot
  print(p)

  # ======================
  # RETURN RESULTS
  # ======================
  list(
    cutoff = cut_off,
    cox_model = cox_model,
    mdl0 = mdl0,
    mdl1 = mdl1,
    mdl2 = mdl2,
    data = data,
    pred_data = pred_df
  )
}
