
robust_parameter_analysis <- function(data, formula, param_name, n_boot = 1000, n_perm = 1000, ncores = detectCores() - 1) {
  # Fit original model and get estimate
  original_model <- lmerTest::lmer(formula, data = data, REML = TRUE)
  
  # Get all coefficient names
  coef_names <- names(fixef(original_model))
  
  # Find coefficients that match the param_name
  matching_coefs <- coef_names[grepl(param_name, coef_names)]
  
  # If no matches found or if param_name is an exact match, use it directly
  if(length(matching_coefs) == 0 || param_name %in% coef_names) {
    matching_coefs <- param_name
  }
  
  # Process each matching coefficient
  results_list <- do.call(rbind, lapply(matching_coefs, function(current_param) {
    original_estimate <- fixef(original_model)[current_param]
    
    # Function for bootstrap resampling
    boot_fn <- function(data, indices) {
      boot_data <- data[indices, ]
      boot_model <- try(lmer(formula, data = boot_data, REML = TRUE), silent = TRUE)
      if(inherits(boot_model, "try-error")) return(NA)
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
    var_name <- if(length(base_vars) > 0) base_vars[1] else param_name
    
    # Run permutation test
    perm_estimates <- unlist(mclapply(1:n_perm, function(i) {
      tryCatch({
        perm_data <- data
        # Only permute if the variable exists in the data
        if(var_name %in% names(perm_data)) {
          perm_data[[var_name]] <- sample(perm_data[[var_name]])
        }
        
        perm_model <- try(lmer(formula, data = perm_data, REML = TRUE), silent = TRUE)
        
        as.numeric(fixef(perm_model)[current_param])
      }, error=function(e) {return(NA)})
    }, mc.cores=ncores))
    
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