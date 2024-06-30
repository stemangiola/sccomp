#' @importFrom tidyr complete
#' @importFrom tidyr nesting
#' @importFrom tidyr replace_na
sccomp_glm_data_frame_raw = function(.data,
                                     formula_composition = ~ 1 ,
                                     formula_variability = ~ 1,
                                     .sample,
                                     .cell_group,
                                     .count = NULL,
                                     
                                     # Secondary arguments
                                     contrasts = NULL,
                                     prior_mean = list(intercept = c(0,1), coefficients = c(0,1)),                        
                                     prior_overdispersion_mean_association = list(intercept = c(5, 2), slope = c(0,  0.6), standard_deviation = c(20, 40)),
                                     percent_false_positive =  5,
                                     check_outliers = TRUE,
                                     variational_inference = NULL,
                                     inference_method = "variational",
                                     test_composition_above_logit_fold_change = 0.1, .sample_cell_group_pairs_to_exclude = NULL,
                                     verbose = FALSE,
                                     exclude_priors = FALSE,
                                     bimodal_mean_variability_association = FALSE,
                                     enable_loo = FALSE,
                                     use_data = TRUE,
                                     cores = 4,
                                     mcmc_seed = sample(1e5, 1),
                                     max_sampling_iterations = 20000,
                                     pass_fit = TRUE ) {
  
  # See https://community.rstudio.com/t/how-to-make-complete-nesting-work-with-quosures-and-tidyeval/16473
  # See https://github.com/tidyverse/tidyr/issues/506
  
  
  # Prepare column same enquo
  .sample = enquo(.sample)
  .cell_group = enquo(.cell_group)
  .sample_cell_group_pairs_to_exclude = enquo(.sample_cell_group_pairs_to_exclude)
  
  # Check if columns exist
  check_columns_exist(.data, c(
    quo_name(.sample),
    quo_name(.cell_group),
    parse_formula(formula_composition)
  ))
  
  # Check if any column is NA or null
  check_if_any_NA(.data, c(
    quo_name(.sample),
    quo_name(.cell_group),
    parse_formula(formula_composition)
  ))
  
  .grouping_for_random_intercept = parse_formula_random_intercept(formula_composition) |> pull(grouping) |> unique()
  
  # Make counts
  .data %>%

    class_list_to_counts(!!.sample, !!.cell_group) %>%
    
    # Add formula_composition information
    add_formula_columns(.data, !!.sample,formula_composition) |>
    
    # Attach possible exclusion of data points
    left_join(.data %>%
                as_tibble() |>
                select(
                  !!.sample, !!.cell_group, any_of(quo_name(.sample_cell_group_pairs_to_exclude))
                ) %>%
                distinct(),
              by = c(quo_name(.sample), quo_name(.cell_group))
    ) |>
    mutate(
      across(!!.sample_cell_group_pairs_to_exclude, ~replace_na(.x, 0))
    ) |>
    
    # Return
    sccomp_glm_data_frame_counts(
      formula_composition = formula_composition,
      formula_variability = formula_variability,
      
      .sample = !!.sample,
      .cell_group = !!.cell_group,
      .count = count,
      contrasts = contrasts,
      #.grouping_for_random_intercept = !! .grouping_for_random_intercept,
      prior_mean = prior_mean,
      prior_overdispersion_mean_association = prior_overdispersion_mean_association,
      percent_false_positive =  percent_false_positive,
      check_outliers = check_outliers,
      inference_method = inference_method,
      exclude_priors = exclude_priors,
      bimodal_mean_variability_association = bimodal_mean_variability_association,
      enable_loo = enable_loo,
      use_data = use_data,
      cores = cores,
      test_composition_above_logit_fold_change = test_composition_above_logit_fold_change, .sample_cell_group_pairs_to_exclude = !!.sample_cell_group_pairs_to_exclude,
      verbose = verbose,
      mcmc_seed = mcmc_seed,
      max_sampling_iterations = max_sampling_iterations,
      pass_fit = pass_fit
    )
}

sccomp_glm_data_frame_counts = function(.data,
                                        formula_composition = ~ 1 ,
                                        formula_variability = ~ 1,
                                        .sample,
                                        .cell_group,
                                        .count = NULL,
                                        
                                        # Secondary arguments
                                        contrasts = NULL,
                                        #.grouping_for_random_intercept = NULL,
                                        prior_mean = list(intercept = c(0,1), coefficients = c(0,1)),                        
                                        prior_overdispersion_mean_association = list(intercept = c(5, 2), slope = c(0,  0.6), standard_deviation = c(20, 40)),
                                        percent_false_positive = 5,
                                        check_outliers = TRUE,
                                        variational_inference = NULL,
                                        inference_method = "variational",
                                        test_composition_above_logit_fold_change = 0.1, .sample_cell_group_pairs_to_exclude = NULL,
                                        verbose = FALSE,
                                        exclude_priors = FALSE,
                                        bimodal_mean_variability_association = FALSE,
                                        enable_loo = FALSE,
                                        use_data = TRUE,
                                        cores = 4,
                                        mcmc_seed = sample(1e5, 1),
                                        max_sampling_iterations = 20000,
                                        pass_fit = TRUE) {
  
  # Prepare column same enquo
  .sample = enquo(.sample)
  .cell_group = enquo(.cell_group)
  .count = enquo(.count)
  .sample_cell_group_pairs_to_exclude = enquo(.sample_cell_group_pairs_to_exclude)
  #.grouping_for_random_intercept = enquo(.grouping_for_random_intercept)
  
  
  #Check column class
  check_if_columns_right_class(.data, !!.sample, !!.cell_group)
  
  # Check that count is integer
  if(.data %>% pull(!!.count) %>% is("integer")) message(sprintf("sccomp says: %s column is an integer. The sum-constrained beta binomial model will be used", quo_name(.count)))
  else if(.data %>% pull(!!.count) %>% is("integer") |> not() & .data %>% pull(!!.count) |> dplyr::between(0, 1) |> all()) message(sprintf("sccomp says: %s column is a proportion. The sum-constrained beta model will be used. When possible using counts is preferred as the binomial noise component is often dominating for rare groups (e.g. rare cell types).", quo_name(.count)))
  else stop(sprintf("sccomp: %s column must be an integer or a proportion", quo_name(.count)))
  
  # Check if columns exist
  check_columns_exist(.data, c(
    quo_name(.sample),
    quo_name(.cell_group),
    quo_name(.count),
    parse_formula(formula_composition)
  ))
  
  # Check that I have rectangular data frame
  if(
    .data |> count(!!.sample) |> distinct(n) |> nrow() > 1 || 
    .data |> count(!!.cell_group) |> distinct(n) |> nrow() > 1 
  ){
    warning(sprintf("sccomp says: the input data frame does not have the same number of `%s`, for all `%s`. We have made it so, adding 0s for the missing sample/feature pairs.", quo_name(.cell_group), quo_name(.sample)))
    .data = .data |> 
      
      # I need renaming trick because complete list(...) cannot accept quosures
      select(!!.sample, !!.cell_group, count := !!.count) |> 
      distinct() |> 
      complete( !!.sample, !!.cell_group, fill = list(count = 0) ) %>%
      rename(!!.count := count) |> 
      mutate(!!.count := as.integer(!!.count)) |> 
      
      # Add formula_composition information
      add_formula_columns(.data, !!.sample, formula_composition) 
  }
  
  
  # Check if test_composition_above_logit_fold_change is 0, as the Bayesian FDR does not allow it
  if(test_composition_above_logit_fold_change <= 0)
    stop("sccomp says: test_composition_above_logit_fold_change should be > 0 for the FDR to be calculated in the Bayesian context (doi: 10.1093/biostatistics/kxw041). Also, testing for > 0 differences avoids significant but meaningless (because of the small magnitude) estimates.")
  
  
  # Check if any column is NA or null
  check_if_any_NA(.data, c(
    quo_name(.sample),
    quo_name(.cell_group),
    quo_name(.count),
    parse_formula(formula_composition)
  ))
  
  # Return
  
  # Credible interval
  CI = 1 - (percent_false_positive/100)
  
  # Produce data list
  factor_names = parse_formula(formula_composition)
  
  # Random intercept
  random_intercept_elements = parse_formula_random_intercept(formula_composition)
  
  # Variational only if no random intercept
  if(inference_method=="variational" & random_intercept_elements |> filter(factor != "(Intercept)") |> nrow() > 0)
    stop("sccomp says: for random effect modelling plese use `inference_method` = \"hmc\", for the full Bayes HMC inference.")
  
  # If no random intercept fake it
  if(nrow(random_intercept_elements)>0){
    .grouping_for_random_intercept = random_intercept_elements |> pull(grouping) |> unique() |>   map(~ quo(!! sym(.x)))
    
  } else{
    .grouping_for_random_intercept = list(quo(!! sym("random_intercept")))
    .data = .data |> mutate(!!.grouping_for_random_intercept[[1]] := "1")
  }
  
  # If .sample_cell_group_pairs_to_exclude 
  if(quo_is_symbolic(.sample_cell_group_pairs_to_exclude)){
    
    # Error if not logical
    if(.data |> pull(!!.sample_cell_group_pairs_to_exclude) |> is("logical") |> not())
      stop(glue("sccomp says: {quo_name(.sample_cell_group_pairs_to_exclude)} must be logical"))
    
    # Error if not consistent to sample/cell group
    if(.data |> count(!!.sample, !!.cell_group, name = "n") |> filter(n>1) |> nrow() |> gt(0))
      stop(glue("sccomp says: {quo_name(.sample_cell_group_pairs_to_exclude)} must be unique with .sample/.cell_group pairs. You might have a .sample/.cell_group pair with both TRUE and FALSE {quo_name(.sample_cell_group_pairs_to_exclude)}."))
    
  } else {
    
    # If no .sample_cell_group_pairs_to_exclude fake it
    .sample_cell_group_pairs_to_exclude = quo(!! sym(".sample_cell_group_pairs_to_exclude"))
    .data = .data |> mutate(!!.sample_cell_group_pairs_to_exclude := FALSE)
  }
  
  
  # Original - old
  # prec_sd ~ normal(0,2);
  # prec_coeff ~ normal(0,5);
  
  message("sccomp says: estimation")
  
  data_for_model =
    .data %>%
    data_to_spread ( formula_composition, !!.sample, !!.cell_group, !!.count, .grouping_for_random_intercept) %>%
    data_spread_to_model_input(
      formula_composition, !!.sample, !!.cell_group, !!.count,
      truncation_ajustment = 1.1,
      approximate_posterior_inference = inference_method %in% c("variational", "pathfinder"),
      formula_variability = formula_variability,
      contrasts = contrasts,
      bimodal_mean_variability_association = bimodal_mean_variability_association,
      use_data = use_data,
      random_intercept_elements
    )
  
  # Print design matrix
  message(sprintf("sccomp says: the composition design matrix has columns: %s", data_for_model$X %>% colnames %>% paste(collapse=", ")))
  message(sprintf("sccomp says: the variability design matrix has columns: %s", data_for_model$Xa %>% colnames %>% paste(collapse=", ")))
  
  # Force outliers, Get the truncation index
  data_for_model$user_forced_truncation_not_idx = 
    .data |> 
    select(!!.sample, !!.cell_group, !!.sample_cell_group_pairs_to_exclude) |>
    left_join( data_for_model$y |> rownames() |> enframe(name="N", value=quo_name(.sample)), by = join_by(!!.sample) ) |>  
    left_join( data_for_model$y |> colnames() |> enframe(name="M", value=quo_name(.cell_group)), by = join_by(!!.cell_group) ) |> 
    select(!!.sample_cell_group_pairs_to_exclude, N, M) |>
    arrange(N, M) |>
    pull(!!.sample_cell_group_pairs_to_exclude) |>
    not() |>
    which()
  
  data_for_model$truncation_not_idx = data_for_model$user_forced_truncation_not_idx
  data_for_model$TNS = length(data_for_model$truncation_not_idx)
  
  # Prior
  data_for_model$prior_prec_intercept = prior_overdispersion_mean_association$intercept
  data_for_model$prior_prec_slope  = prior_overdispersion_mean_association$slope
  data_for_model$prior_prec_sd = prior_overdispersion_mean_association$standard_deviation
  data_for_model$prior_mean_intercept = prior_mean$intercept
  data_for_model$prior_mean_coefficients = prior_mean$coefficients
  data_for_model$exclude_priors = exclude_priors
  data_for_model$enable_loo = TRUE & enable_loo
  
  # # Check that design matrix is not too big
  # if(ncol(data_for_model$X)>20)
  #   message("sccomp says: the design matrix has more than 20 columns. Possibly some numerical factors are erroneously of type character/factor.")
  
  fit =
    data_for_model %>%
    
    # Run the first discovery phase with permissive false discovery rate
    fit_model(
      stanmodels$glm_multi_beta_binomial,
      cores= cores,
      quantile = CI,
      inference_method = inference_method,
      verbose = verbose,
      seed = mcmc_seed,
      max_sampling_iterations = max_sampling_iterations,
      pars = c("beta", "alpha", "prec_coeff","prec_sd",   "alpha_normalised", "beta_random_intercept", "log_lik")
    )
  

  
  estimate_tibble = 
    # Create a dummy tibble
    tibble() |>
    # Attach association mean concentration
    add_attr(fit, "fit") %>%
    add_attr(data_for_model, "model_input") |>
    add_attr(.data, "truncation_df2") |>
    add_attr(.sample, ".sample") |>
    add_attr(.cell_group, ".cell_group") |>
    add_attr(.count, ".count") |>
    add_attr(check_outliers, "check_outliers") |>
    add_attr(formula_composition, "formula_composition") |>
    add_attr(formula_variability, "formula_variability") |>
    add_attr(parse_formula(formula_composition), "factors" ) |> 
    
    # Add class to the tbl
    add_class("sccomp_tbl") |> 
    
    # Print estimates
    sccomp_test() |>
    
    # drop hypothesis testing as the estimation exists without probabilities.
    # For hypothesis testing use sccomp_test
    select(-contains("_FDR"), -contains("_pH0")) 
  
  
  if(inference_method %in% c("variational") && max(na.omit(estimate_tibble$c_R_k_hat)) > 4)
    warning("sccomp says: using variational inference, c_R_k_hat resulted too high for some parameters, indicating lack of convergence of the model. We reccomend using inference_method = \"hmc\" to use the state-of-the-art (although slower) HMC sampler.")
  
  estimate_tibble
  
}




#' @importFrom stats model.matrix
get_mean_precision = function(fit, data_for_model){
  
  # Define the variables as NULL to avoid CRAN NOTES
  M <- NULL
  `2.5%` <- NULL
  `97.5%` <- NULL
  
  fit %>%
    summary_to_tibble("alpha", "C", "M") %>%

    # Just grub intercept of alpha
    filter(C==1) %>%
    select( M, mean, `2.5%` , `97.5%`) %>%
    nest(concentration = -M)

  # WRONG ATTEMPT TO PLOT ALPHA FROM MULTIPLE factorS
  # fit %>%
  #   draws_to_tibble_x_y("alpha", "C", "M") %>%
  #   nest(data = -c(M)) %>%
  #
  #   # Add Design matrix
  #   left_join(
  #     as_tibble(data_for_model$X, rownames="M") %>%
  #       nest(X = -M) %>%
  #       mutate(M = as.integer(M)),
  #     by="M"
  #   ) %>%
  #   mutate(concentration = map2(
  #     data, X,
  #     ~ .x %>%
  #           select(.draw, C, .value) %>%
  #           spread(C, .value) %>%
  #           as_matrix(rownames=".draw") %*%
  #
  #         ( .y %>% as_matrix() %>% t() ) %>%
  #       quantile(c(0.025, 0.5, 0.975)) %>%
  #       enframe() %>%
  #       spread(name, value) %>%
  #       rename(mean = `50%`)
  #   )) %>%
  #   select(-data, -X)
}

get_mean_precision_association = function(fit){
  c(
    fit$summary("prec_coeff")[,2] |> set_names("prec_coeff") ,
    fit$summary("prec_sd")[,2] |> set_names("prec_sd")
  )
}
