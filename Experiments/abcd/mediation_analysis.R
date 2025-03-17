require(tidyverse)
require(dplyr)
require(lme4)
require(lmerTest)
require(broom.mixed)
require(boot)
require(parallel)
require(lavaan)


#' Robust Parameter Analysis
#'
#' @param data Data frame containing all variables
#' @param formula Model formula
#' @param param_name Name of parameter to analyze
#' @param n_boot Number of bootstrap iterations
#' @param n_perm Number of permutation iterations
#' @param ncores Number of CPU cores to use (default: all but one)
#' @return Data frame with parameter estimates and statistical tests
robust_parameter_analysis <- function(data, formula, param_name, n_boot = 1000, 
                                      n_perm = 1000, ncores = detectCores() - 1) {
  # Fit original model and get estimate
  original_model <- lmerTest::lmer(formula, data = data, REML = TRUE)
  
  # Get all coefficient names
  coef_names <- names(fixef(original_model))
  
  # Find coefficients that match the param_name
  matching_coefs <- coef_names[grepl(param_name, coef_names)]
  
  # If no matches found or if param_name is an exact match, use it directly
  if (length(matching_coefs) == 0 || param_name %in% coef_names) {
    matching_coefs <- param_name
  }
  
  # Process each matching coefficient
  results_list <- do.call(rbind, lapply(matching_coefs, function(current_param) {
    original_estimate <- fixef(original_model)[current_param]
    
    # Function for bootstrap resampling
    boot_fn <- function(data, indices) {
      boot_data <- data[indices, ]
      boot_model <- try(lmer(formula, data = boot_data, REML = TRUE), silent = TRUE)
      if (inherits(boot_model, "try-error")) return(NA)
      return(fixef(boot_model)[current_param])
    }
    
    # Get bootstrap distribution
    boot_results <- boot(data, boot_fn, R = n_boot, parallel = "multicore", ncpus = ncores)
    
    # Calculate percentile confidence interval
    boot_ci <- tryCatch({
      boot.ci(boot_results, type = "perc", index = 1)
    }, error = function(e) {
      return(list(percent = c(NA, NA, NA, NA, NA)))
    })
    
    # Check if CI excludes zero
    sig_ci <- (boot_ci$percent[4] > 0 & boot_ci$percent[5] > 0) | 
      (boot_ci$percent[4] < 0 & boot_ci$percent[5] < 0)
    
    # Find which variable to permute - use the base variable name derived from the coefficient
    model_terms <- attr(terms(original_model), "term.labels")
    
    # Find base variable using lapply and Filter
    base_vars <- Filter(function(term) grepl(paste0("^", term), current_param), model_terms)
    var_name <- if (length(base_vars) > 0) base_vars[1] else param_name
    
    # Run permutation test
    perm_estimates <- unlist(mclapply(1:n_perm, function(i) {
      tryCatch({
        perm_data <- data
        # Only permute if the variable exists in the data
        if (var_name %in% names(perm_data)) {
          perm_data[[var_name]] <- sample(perm_data[[var_name]])
        }
        
        perm_model <- try(lmer(formula, data = perm_data, REML = TRUE), silent = TRUE)
        
        as.numeric(fixef(perm_model)[current_param])
      }, error = function(e) {
        return(NA)
      })
    }, mc.cores = ncores))
    
    # Calculate two-sided p-value with +1 adjustment
    perm_estimates <- perm_estimates[!is.na(perm_estimates)]
    n_perm_valid <- length(perm_estimates)
    p_lower <- (sum(perm_estimates <= original_estimate) + 1) / (n_perm_valid + 1)
    p_upper <- (sum(perm_estimates >= original_estimate) + 1) / (n_perm_valid + 1)
    perm_p <- 2 * min(p_lower, p_upper)
    perm_p <- min(perm_p, 1)
    
    data.frame(
      parameter = current_param,
      estimate = original_estimate,
      lower.ci = boot_ci$percent[4],
      upper.ci = boot_ci$percent[5],
      significant.ci = sig_ci,
      p.value = perm_p,
      stringsAsFactors = FALSE
    )
  }))
  
  return(results_list)
}


parallel_med_analysis <- function(data, covariates, mediators, outcome, 
                                          predictor, random_effect = "family_id", 
                                          bootstrap = TRUE, n_boot = 1000) {
  # Step 1: Handle random effects
  # For random effects in lavaan, we can use clustered standard errors or two-level modeling
  # Here we'll use a cluster-robust approach
  
  null_formula <- as.formula(paste0(outcome, " ~ 1 + (1|", random_effect, ")"))
  null_model <- lmer(null_formula, data = data)
  icc <- performance::icc(null_model)$ICC_adjusted
  
  cluster_var <- random_effect
  
  # paths from predictor to mediators (a paths)
  a_paths <- paste0(mediators, " ~ a", seq_along(mediators), "*", predictor)
  
  # Add paths from mediators to outcome (b paths)
  b_paths <- paste0(outcome, " ~ b", seq_along(mediators), "*", mediators)
  
  # Add direct effect
  c_prime_path <- paste0(outcome, " ~ c_prime*", predictor)
  
  # Add covariate paths if needed
  covar_paths <- c()
  if (length(covariates) > 0) {
    # Covariates affect mediators
    for (med in mediators) {
      covar_paths <- c(covar_paths, paste0(med, " ~ ", paste0(covariates, collapse = " + ")))
    }
    # Covariates affect outcome
    covar_paths <- c(covar_paths, paste0(outcome, " ~ ", paste0(covariates, collapse = " + ")))
  }
  
  # Allow mediators to correlate
  med_corr <- c()
  if (length(mediators) > 1) {
    med_pairs <- combn(mediators, 2)
    for (i in 1:ncol(med_pairs)) {
      med_corr <- c(med_corr, paste0(med_pairs[1,i], " ~~ ", med_pairs[2,i]))
    }
  }
  
  # Define indirect effects
  indirect_effects <- paste0("indirect", seq_along(mediators), " := a", 
                             seq_along(mediators), "*b", seq_along(mediators))
  
  # Total indirect effect
  if (length(mediators) > 1) {
    total_indirect <- paste0("total_indirect := ", 
                             paste0("indirect", seq_along(mediators), collapse = " + "))
  } else {
    total_indirect <- "total_indirect := indirect1"
  }
  
  # Total effect
  total_effect <- "total := c_prime + total_indirect"
  
  # Combine all paths
  model <- c(a_paths, b_paths, c_prime_path, covar_paths, med_corr, 
             indirect_effects, total_indirect, total_effect)
  
  # Convert to lavaan syntax
  model_syntax <- paste(model, collapse = "\n")
  
  # Step 3: Fit the model using lavaan
  if (bootstrap) {
    fit <- sem(model_syntax, data = data, 
               estimator = "MLR", 
               se = "bootstrap", 
               bootstrap = n_boot,
               cluster = cluster_var)
  } else {
    fit <- sem(model_syntax, data = data, 
               estimator = "MLR", 
               cluster = cluster_var)
  }
  
  # Step 4: Format results similar to the original function
  # Extract parameter estimates
  params <- parameterEstimates(fit, standardized = TRUE)
  
  # Create result dataframe
  result_list <- list()
  
  # Total effect
  total_row <- params[params$label == "total", ]
  total_df <- data.frame(
    roi = outcome,
    label = "TOTAL_EFFECT",
    a_path = NA,
    b_path = NA,
    effect_size = total_row$est[1],
    lower_ci = total_row$ci.lower[1],
    upper_ci = total_row$ci.upper[1],
    p.value = total_row$pvalue[1],
    prop_total = 1,
    stringsAsFactors = FALSE
  )
  result_list[[1]] <- total_df
  
  # Total indirect effect
  total_ind_row <- params[params$label == "total_indirect", ]
  total_ind_df <- data.frame(
    roi = outcome,
    label = "TOTAL_INDIRECT",
    a_path = NA,
    b_path = NA,
    effect_size = total_ind_row$est[1],
    lower_ci = total_ind_row$ci.lower[1],
    upper_ci = total_ind_row$ci.upper[1],
    p.value = total_ind_row$pvalue[1],
    prop_total = total_ind_row$est[1] / total_row$est[1],
    stringsAsFactors = FALSE
  )
  result_list[[2]] <- total_ind_df
  
  # Individual mediator effects
  for (i in seq_along(mediators)) {
    med <- mediators[i]
    a_path_row <- params[params$label == paste0("a", i), ][1, ]
    b_path_row <- params[params$label == paste0("b", i), ][1, ]
    indirect_row <- params[params$label == paste0("indirect", i), ][1, ]
    
    med_df <- data.frame(
      roi = outcome,
      label = med,
      a_path = a_path_row$est,
      b_path = b_path_row$est,
      effect_size = indirect_row$est,
      lower_ci = indirect_row$ci.lower,
      upper_ci = indirect_row$ci.upper,
      p.value = indirect_row$pvalue,
      prop_total = indirect_row$est / total_row$est[1],
      stringsAsFactors = FALSE
    )
    result_list[[i + 2]] <- med_df
  }
  
  # Direct effect
  direct_row <- params[params$label == "c_prime", ][1, ]
  direct_df <- data.frame(
    roi = outcome,
    label = "DIRECT_EFFECT",
    a_path = NA,
    b_path = NA,
    effect_size = direct_row$est,
    lower_ci = direct_row$ci.lower,
    upper_ci = direct_row$ci.upper,
    p.value = direct_row$pvalue,
    prop_total = direct_row$est / total_row$est[1],
    stringsAsFactors = FALSE
  )
  result_list[[length(result_list) + 1]] <- direct_df
  
  # Combine all results
  final_df <- do.call(rbind, result_list)
  
  attr(final_df, "model") <- fit
  
  return(final_df)
}