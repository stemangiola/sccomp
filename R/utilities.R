
sccomp_stan_models_cache_dir = file.path(path.expand("~"), ".sccomp_models")

# Greater than
gt = function(a, b){	a > b }

# Smaller than
st = function(a, b){	a < b }

# Negation
not = function(is){	!is }

#' Add attribute to abject
#'
#' @keywords internal
#' @noRd
#'
#'
#' @param var A tibble
#' @param attribute An object
#' @param name A character name of the attribute
#'
#' @return A tibble with an additional attribute
add_attr = function(var, attribute, name) {
  attr(var, name) <- attribute
  var
}


#' Formula parser
#'
#' @param fm A formula
#'
#' @importFrom stringr str_subset
#' @importFrom magrittr extract2
#' @importFrom stats terms
#'
#' @return A character vector
#'
#' @keywords internal
#' @noRd
parse_formula <- function(fm) {
  stopifnot("The formula must be of the kind \"~ factors\" " = attr(terms(fm), "response") == 0)

    as.character(attr(terms(fm), "variables")) |>
    str_subset("\\|", negate = TRUE) %>%

      # Does not work the following
      # |>
      # extract2(-1)
      .[-1]
}


#' Formula parser
#'
#' @param fm A formula
#'
#' @importFrom stringr str_subset
#' @importFrom stringr str_split
#' @importFrom stringr str_remove_all
#' @importFrom rlang set_names
#' @importFrom purrr map_dfr
#' @importFrom stringr str_trim
#'
#' @importFrom magrittr extract2
#'
#' @return A character vector
#'
#' @keywords internal
#' @noRd
formula_to_random_effect_formulae <- function(fm) {

  # Define the variables as NULL to avoid CRAN NOTES
  formula <- NULL
  
  stopifnot("The formula must be of the kind \"~ factors\" " = attr(terms(fm), "response") == 0)

  random_intercept_elements =
    as.character(attr(terms(fm), "variables")) |>

    # Select random intercept part
    str_subset("\\|")

  if(length(random_intercept_elements) > 0){

    random_intercept_elements |>

      # Divide grouping from factors
      str_split("\\|") |>

      # Set name
      map_dfr(~ .x |> set_names(c("formula", "grouping"))) |>

      # Create formula
      mutate(formula = map(formula, ~ formula(glue("~ {.x}")))) |>
      mutate(grouping = grouping |> str_trim())

  }

  else
    tibble(`formula` = list(), grouping = character())

}

#' Formula parser
#'
#' @param fm A formula
#'
#' @importFrom stringr str_subset
#' @importFrom stringr str_split
#' @importFrom stringr str_remove_all
#' @importFrom rlang set_names
#' @importFrom purrr map_dfr
#'
#' @importFrom magrittr extract2
#'
#' @return A character vector
#'
#' @keywords internal
#' @noRd
parse_formula_random_intercept <- function(fm) {

  # Define the variables as NULL to avoid CRAN NOTES
  formula <- NULL
  
  stopifnot("The formula must be of the kind \"~ factors\" " = attr(terms(fm), "response") == 0)

  random_intercept_elements =
    as.character(attr(terms(fm), "variables")) |>

    # Select random intercept part
    str_subset("\\|")

  if(length(random_intercept_elements) > 0){

    formula_to_random_effect_formulae(fm) |>

      # Divide factors
      mutate(factor = map(
        formula,
        ~
          # Attach intercept
          .x |>
          terms() |>
          attr("intercept") |>
          str_replace("^1$", "(Intercept)") |>
          str_subset("0", negate = TRUE) |>

          # Attach variables
          c(
            .x |>
            terms() |>
            attr("variables") |>
            as.character() |>
            str_split("\\+") |>
            as.character() %>%
            .[-1]
          )
      )) |>
      unnest(factor)

  }

  else
    tibble(factor = character(), grouping = character())

}

#' Get matrix from tibble
#'
#' @import dplyr
#'
#' @keywords internal
#' @noRd
#'
#' @import dplyr
#' @importFrom purrr as_mapper
#'
#' @param .x A tibble
#' @param .p A boolean
#' @param .f1 A function
#' @param .f2 A function
#'
#' @return A tibble
ifelse_pipe = function(.x, .p, .f1, .f2 = NULL) {
  switch(.p %>% `!` %>% sum(1),
         as_mapper(.f1)(.x),
         if (.f2 %>% is.null %>% `!`)
           as_mapper(.f2)(.x)
         else
           .x)

}

#' @importFrom tidyr gather
#' @importFrom magrittr set_rownames
#'
#' @keywords internal
#' @noRd
#'
#' @param tbl A tibble
#' @param rownames A character string of the rownames
#'
#' @return A matrix
as_matrix <- function(tbl, rownames = NULL) {
  
  # Define the variables as NULL to avoid CRAN NOTES
  variable <- NULL
  
  tbl %>%

    ifelse_pipe(
      tbl %>%
        ifelse_pipe(!is.null(rownames),		~ .x %>% dplyr::select(-contains(rownames))) %>%
        summarise_all(class) %>%
        gather(variable, class) %>%
        pull(class) %>%
        unique() %>%
        `%in%`(c("numeric", "integer")) %>% not() %>% any(),
      ~ {
        warning("to_matrix says: there are NON-numerical columns, the matrix will NOT be numerical")
        .x
      }
    ) %>%
    as.data.frame() %>%

    # Deal with rownames column if present
    ifelse_pipe(!is.null(rownames),
                ~ .x  %>%
                  set_rownames(tbl %>% pull(!!rownames)) %>%
                  select(-!!rownames)) %>%

    # Convert to matrix
    as.matrix()
}

#' vb_iterative
#'
#' @description Runs iteratively variational bayes until it suceeds
#'
#'
#' @keywords internal
#' @noRd
#'
#' @param model A Stan model
#' @param output_samples An integer of how many samples from posteriors
#' @param iter An integer of how many max iterations
#' @param tol_rel_obj A real
#' @param additional_parameters_to_save A character vector
#' @param data A data frame
#' @param seed An integer
#' @param ... List of paramaters for vb function of Stan
#'
#' @return A Stan fit object
#'
vb_iterative = function(model,
                        output_samples,
                        iter,
                        tol_rel_obj,
                        additional_parameters_to_save = c(),
                        data,
                        output_dir = output_dir,
                        seed, 
                        init = "random",
                        inference_method,
                        cores = 1,
                        ...) {
  res = NULL
  i = 0
  while  (is.null(res) & i < 5) {
    res = tryCatch({

      if(inference_method=="pathfinder")
        my_res = model$pathfinder(
        	data = data,
        	tol_rel_obj = tol_rel_obj,
        	output_dir = output_dir,
        	seed = seed+i,
        	# init = init,
        	num_paths=50, 
          num_threads = cores,
        	single_path_draws = output_samples / 50 ,
        	max_lbfgs_iters=100, 
        	history_size = 100, 
        	...

        )
    
      else if(inference_method=="variational")
        my_res = model$variational(
          data = data,
          output_samples = output_samples,
          iter = iter,
          tol_rel_obj = tol_rel_obj,
          output_dir = output_dir,
          seed = seed+i,
          init = init,
          ...
        )

      boolFalse <- TRUE
      return(my_res)
    },
    error = function(e) {
      writeLines(sprintf("Further attempt with Variational Bayes: %s", e))
      return(NULL)
    },
    finally = {
    })
    i = i + 1
  }

  if(is.null(res)) stop(sprintf("sccomp says: variational Bayes did not converge after %s attempts. Please use variational_inference = FALSE for a HMC fitting.", i))
  
  return(res)
}




#' draws_to_tibble_x_y
#'
#' @importFrom tidyr pivot_longer
#' @importFrom rlang :=
#'
#' @param fit A fit object
#' @param par A character vector. The parameters to extract.
#' @param x A character. The first index.
#' @param y A character. The first index.
#'
#' @keywords internal
#' @noRd
draws_to_tibble_x_y = function(fit, par, x, y, number_of_draws = NULL) {

  # Define the variables as NULL to avoid CRAN NOTES
  dummy <- NULL
  .variable <- NULL
  .chain <- NULL
  .iteration <- NULL
  .draw <- NULL
  .value <- NULL
  
  par_names =
    fit$metadata()$stan_variables %>% grep(sprintf("%s", par), ., value = TRUE)

  fit$draws(variables = par, format = "draws_df") %>%
    mutate(.iteration = seq_len(n())) %>%

    pivot_longer(
      names_to = "parameter", # c( ".chain", ".variable", x, y),
      cols = contains(par),
      #names_sep = "\\.?|\\[|,|\\]|:",
      # names_ptypes = list(
      #   ".variable" = character()),
      values_to = ".value"
    ) %>%
    tidyr::extract(parameter, c(".chain", ".variable", x, y), "([1-9]+)?\\.?([a-zA-Z0-9_\\.]+)\\[([0-9]+),([0-9]+)") |> 

    # Warning message:
    # Expected 5 pieces. Additional pieces discarded
    suppressWarnings() %>%

    mutate(
      !!as.symbol(x) := as.integer(!!as.symbol(x)),
      !!as.symbol(y) := as.integer(!!as.symbol(y))
    ) %>%
    arrange(.variable, !!as.symbol(x), !!as.symbol(y), .chain) %>%
    group_by(.variable, !!as.symbol(x), !!as.symbol(y)) %>%
    mutate(.draw = seq_len(n())) %>%
    ungroup() %>%
    select(!!as.symbol(x), !!as.symbol(y), .chain, .iteration, .draw ,.variable ,     .value) %>%
    filter(.variable == par)

}


#' @importFrom tidyr separate
#' @importFrom purrr when
#' 
#' @param fit A fit object from a statistical model, from the 'rstan' package.
#' @param par A character vector specifying the parameters to extract from the fit object.
#' @param x A character string specifying the first index in the parameter names.
#' @param y A character string specifying the second index in the parameter names (optional).
#' @param probs A numerical vector specifying the quantiles to extract.
#' 
#' @keywords internal
#' @noRd
summary_to_tibble = function(fit, par, x, y = NULL, probs = c(0.025, 0.25, 0.50, 0.75, 0.975)) {
  
  # Extract parameter names from the fit object that match the 'par' argument
  par_names = names(fit) %>% grep(sprintf("%s", par), ., value = TRUE)

  # Avoid bug
  #if(fit@stan_args[[1]]$method %>% is.null) fit@stan_args[[1]]$method = "hmc"

  summary = 
    fit$summary(variables = par, "mean", ~quantile(.x, probs = probs,  na.rm=TRUE)) %>%
    rename(.variable = variable ) %>%

    when(
      is.null(y) ~ (.) %>% tidyr::separate(col = .variable,  into = c(".variable", x, y), sep="\\[|,|\\]", convert = TRUE, extra="drop"),
      ~ (.) %>% tidyr::separate(col = .variable,  into = c(".variable", x, y), sep="\\[|,|\\]", convert = TRUE, extra="drop")
    )
  
  # summaries are returned only for HMC
  if(!"n_eff" %in% colnames(summary)) summary = summary |> mutate(n_eff = NA)
  if(!"R_k_hat" %in% colnames(summary)) summary = summary |> mutate(R_k_hat = NA)
  
  summary

}

#' @importFrom rlang :=
label_deleterious_outliers = function(.my_data){

  # Define the variables as NULL to avoid CRAN NOTES
  .count <- NULL
  `95%` <- NULL
  `5%` <- NULL
  X <- NULL
  iteration <- NULL
  outlier_above <- NULL
  slope <- NULL
  is_group_right <- NULL
  outlier_below <- NULL
  
  .my_data %>%

    # join CI
    mutate(outlier_above = !!.count > `95%`) %>%
    mutate(outlier_below = !!.count < `5%`) %>%

    # Mark if on the right of the factor scale
    mutate(is_group_right = !!as.symbol(colnames(X)[2]) > mean( !!as.symbol(colnames(X)[2]) )) %>%

    # Check if outlier might be deleterious for the statistics
    mutate(
      !!as.symbol(sprintf("deleterious_outlier_%s", iteration)) :=
        (outlier_above & slope > 0 & is_group_right)  |
        (outlier_below & slope > 0 & !is_group_right) |
        (outlier_above & slope < 0 & !is_group_right) |
        (outlier_below & slope < 0 & is_group_right)
    ) %>%

    select(-outlier_above, -outlier_below, -is_group_right)

}

#' @importFrom readr write_file
fit_model = function(
  data_for_model, model, censoring_iteration = 1, cores = detectCores(), quantile = 0.95,
  warmup_samples = 300, approximate_posterior_inference = NULL, inference_method, verbose = FALSE,
  seed , pars = c("beta", "alpha", "prec_coeff","prec_sd"), output_samples = NULL, chains=NULL, max_sampling_iterations = 20000
)
{


  # # if analysis approximated
  # # If posterior analysis is approximated I just need enough
  # how_many_posterior_draws_practical = ifelse(approximate_posterior_analysis, 1000, how_many_posterior_draws)
  # additional_parameters_to_save = additional_parameters_to_save %>% c("lambda_log_param", "sigma_raw") %>% unique


  # Find number of draws
  draws_supporting_quantile = 50
  if(is.null(output_samples)){
    
    output_samples =
      (draws_supporting_quantile/((1-quantile)/2)) %>% # /2 because I have two tails
      max(4000) 
    
    if(output_samples > max_sampling_iterations) {
      # message("sccomp says: the number of draws used to defined quantiles of the posterior distribution is capped to 20K.") # This means that for very low probability threshold the quantile could become unreliable. We suggest to limit the probability threshold between 0.1 and 0.01")
      output_samples = max_sampling_iterations
    
  }}
    
  # Find optimal number of chains
  if(is.null(chains))
    chains =
      find_optimal_number_of_chains(
        how_many_posterior_draws = output_samples,
        warmup = warmup_samples,
        parallelisation_start_penalty = 100
      ) %>%
      min(cores)

  init_list=list(
    prec_coeff = c(5,0),
    prec_sd = 1,
    alpha = matrix(c(rep(5, data_for_model$M), rep(0, (data_for_model$A-1) *data_for_model$M)), nrow = data_for_model$A, byrow = TRUE),
    beta_raw_raw = matrix(0, data_for_model$C , data_for_model$M-1) ,
    mix_p = 0.1 
   )

  if(data_for_model$N_random_intercepts>0){
    init_list$random_intercept_raw = matrix(0, data_for_model$N_grouping  , data_for_model$M-1) |> as.data.frame()  
    init_list$random_intercept_sigma_mu = 0.5 |> as.array()
    init_list$random_intercept_sigma_sigma = 0.2 |> as.array()
    init_list$random_intercept_sigma_raw = matrix(0, data_for_model$M-1 , data_for_model$how_many_factors_in_random_design)
    init_list$sigma_correlation_factor = matrix(0, data_for_model$how_many_factors_in_random_design  , data_for_model$how_many_factors_in_random_design )
    init_list$zero_random_intercept = rep(0, size = 1) |> as.array()
    
  }
 
  
  init = map(1:chains, ~ init_list) %>%
    setNames(as.character(1:chains))

  output_directory = "sccomp_draws_files"
  dir.create(output_directory, showWarnings = FALSE)

  # Fit
  mod = load_model("glm_multi_beta_binomial")
  
  
  if(inference_method == "hmc"){
#mod$compile()
      mod$sample(
        data = data_for_model ,
        chains = chains,
        parallel_chains = chains,
        threads_per_chain = 1,
        iter_warmup = warmup_samples,
        iter_sampling = as.integer(output_samples /chains),
        #refresh = ifelse(verbose, 1000, 0),
        seed = seed,
        save_warmup = FALSE,
        init = init,
        output_dir = output_directory
      ) %>%
      suppressWarnings()

}
  else
    vb_iterative(
      mod,
      output_samples = output_samples ,
      iter = 10000,
      tol_rel_obj = 0.01,
      data = data_for_model, refresh = ifelse(verbose, 1000, 0),
      seed = seed,
      output_dir = output_directory,
      init = pf , #list(init_list),
      inference_method = inference_method, 
      cores = cores,
      psis_resample = FALSE
    ) %>%
      suppressWarnings()


}


#' @importFrom purrr map2_lgl
#' @importFrom tidyr pivot_wider
#' @importFrom rlang :=
#'
#' @keywords internal
#' @noRd
parse_fit = function(data_for_model, fit, censoring_iteration = 1, chains){

  # Define the variables as NULL to avoid CRAN NOTES
  M <- NULL
  
  
  fit %>%
    draws_to_tibble_x_y("beta", "C", "M") %>%
    left_join(tibble(C=seq_len(ncol(data_for_model$X)), C_name = colnames(data_for_model$X)), by = "C") %>%
    nest(!!as.symbol(sprintf("beta_posterior_%s", censoring_iteration)) := -M)

}

#' @importFrom purrr map2_lgl
#' @importFrom tidyr pivot_wider
#' @importFrom tidyr spread
#' @importFrom stats C
#' @importFrom rlang :=
#' @importFrom tibble enframe
#'
#' @keywords internal
#' @noRd
beta_to_CI = function(fitted, censoring_iteration = 1, false_positive_rate, factor_of_interest){

  # Define the variables as NULL to avoid CRAN NOTES
  M <- NULL
  C_name <- NULL
  .lower <- NULL
  .median <- NULL
  .upper <- NULL
  
  effect_column_name = sprintf("composition_effect_%s", factor_of_interest) %>% as.symbol()

  CI = fitted %>%
    unnest(!!as.symbol(sprintf("beta_posterior_%s", censoring_iteration))) %>%
    nest(data = -c(M, C, C_name)) %>%
    # Attach beta
    mutate(!!as.symbol(sprintf("beta_quantiles_%s", censoring_iteration)) := map(
      data,
      ~ quantile(
        .x$.value,
        probs = c(false_positive_rate/2,  0.5,  1-(false_positive_rate/2))
      ) %>%
        enframe() %>%
        mutate(name = c(".lower", ".median", ".upper")) %>%
        spread(name, value)
    )) %>%
    unnest(!!as.symbol(sprintf("beta_quantiles_%s", censoring_iteration))) %>%
    select(-data, -C) %>%
    pivot_wider(names_from = C_name, values_from=c(.lower , .median ,  .upper)) 
  
  # Create main effect if exists
  if(!is.na(factor_of_interest) )
    CI |>
    mutate(!!effect_column_name := !!as.symbol(sprintf(".median_%s", factor_of_interest))) %>%
    nest(composition_CI = -c(M, !!effect_column_name))
  
  else 
    CI |> nest(composition_CI = -c(M))

}

#' @importFrom purrr map2_lgl
#' @importFrom tidyr pivot_wider
#' @importFrom stats C
#' @importFrom rlang :=
#'
#' @keywords internal
#' @noRd
alpha_to_CI = function(fitted, censoring_iteration = 1, false_positive_rate, factor_of_interest){

  # Define the variables as NULL to avoid CRAN NOTES
  M <- NULL
  C_name <- NULL
  .lower <- NULL
  .median <- NULL
  .upper <- NULL
  
  effect_column_name = sprintf("variability_effect_%s", factor_of_interest) %>% as.symbol()

  fitted %>%
    unnest(!!as.symbol(sprintf("alpha_%s", censoring_iteration))) %>%
    nest(data = -c(M, C, C_name)) %>%
    # Attach beta
    mutate(!!as.symbol(sprintf("alpha_quantiles_%s", censoring_iteration)) := map(
      data,
      ~ quantile(
        .x$.value,
        probs = c(false_positive_rate/2,  0.5,  1-(false_positive_rate/2))
      ) %>%
        enframe() %>%
        mutate(name = c(".lower", ".median", ".upper")) %>%
        spread(name, value)
    )) %>%
    unnest(!!as.symbol(sprintf("alpha_quantiles_%s", censoring_iteration))) %>%
    select(-data, -C) %>%
    pivot_wider(names_from = C_name, values_from=c(.lower , .median ,  .upper)) %>%
    mutate(!!effect_column_name := !!as.symbol(sprintf(".median_%s", factor_of_interest))) %>%
    nest(variability_CI = -c(M, !!effect_column_name))



}


#' Get Random Intercept Design 2
#'
#' This function processes the formula composition elements in the data and creates design matrices
#' for random intercept models.
#'
#' @param .data_ A data frame containing the data.
#' @param .sample A quosure representing the sample variable.
#' @param formula_composition A data frame containing the formula composition elements.
#' 
#' @return A data frame with the processed design matrices for random intercept models.
#' 
#' @importFrom glue glue
#' @importFrom magrittr subtract
#' @importFrom purrr map2
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr pull
#' @importFrom dplyr filter
#' @importFrom dplyr left_join
#' @importFrom dplyr mutate_all
#' @importFrom dplyr mutate_if
#' @importFrom dplyr as_tibble
#' @importFrom tidyr pivot_longer
#' @importFrom rlang enquo
#' @importFrom rlang quo_name
#' @importFrom tidyselect all_of
#' @importFrom readr type_convert
#' @noRd
get_random_intercept_design2 = function(.data_, .sample, formula_composition ){

  # Define the variables as NULL to avoid CRAN NOTES
  formula <- NULL
  
  .sample = enquo(.sample)

 grouping_table =
   formula_composition |>
   formula_to_random_effect_formulae() |>

   mutate(design = map2(
     formula, grouping,
     ~ {

       mydesign = .data_ |> get_design_matrix(.x, !!.sample)

       mydesign_grouping = .data_ |> select(all_of(.y)) |> pull(1) |> rep(ncol(mydesign)) |> matrix(ncol = ncol(mydesign))
       mydesign_grouping[mydesign==0L] = NA
       colnames(mydesign_grouping) = colnames(mydesign)
       rownames(mydesign_grouping) = rownames(mydesign)

       mydesign_grouping |>
         as_tibble(rownames = quo_name(.sample)) |>
         pivot_longer(-!!.sample, names_to = "factor", values_to = "grouping") |>
         filter(!is.na(grouping)) |>

          mutate("mean_idx" = glue("{factor}___{grouping}") |> as.factor() |> as.integer() )|>
          with_groups(factor, ~ ..1 |> mutate(mean_idx = if_else(mean_idx == max(mean_idx), 0L, mean_idx))) |>
         mutate(minus_sum = if_else(mean_idx==0, factor |> as.factor() |> as.integer(), 0L)) |>

         # Make right rank
         mutate(mean_idx = mean_idx |> as.factor() |> as.integer() |> subtract(1)) |>

         # drop minus_sum if we just have one grouping per factor
         with_groups(factor, ~ {
           if(length(unique(..1$grouping)) == 1) ..1 |> mutate(., minus_sum = 0)
             else ..1
         }) |>

         # Add value
        left_join(

          mydesign |>
            as_tibble(rownames = quo_name(.sample)) |>
            mutate_all(as.character) |>
            readr::type_convert(guess_integer = TRUE ) |>
            suppressMessages() |>
            mutate_if(is.integer, ~1) |>
            pivot_longer(-!!.sample, names_to = "factor"),

          by = join_by(!!.sample, factor)
        ) |>

         # Create unique name
         mutate(group___label = glue("{factor}___{grouping}")) |>
         mutate(group___numeric = group___label |> as.factor() |> as.integer()) |>
         mutate(factor___numeric = `factor` |> as.factor() |> as.integer())



     }))

 }


#' Get Random Intercept Design
#'
#' This function processes random intercept elements in the data and creates design matrices
#' for random intercept models.
#'
#' @param .data_ A data frame containing the data.
#' @param .sample A quosure representing the sample variable.
#' @param random_intercept_elements A data frame containing the random intercept elements.
#' 
#' @return A data frame with the processed design matrices for random intercept models.
#' 
#' @importFrom glue glue
#' @importFrom magrittr subtract
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom dplyr pull
#' @importFrom dplyr if_else
#' @importFrom dplyr distinct
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr set_names
#' @importFrom purrr map_lgl
#' @importFrom purrr pmap
#' @importFrom purrr map_int
#' @importFrom tidyr with_groups
#' @importFrom rlang enquo
#' @importFrom rlang quo_name
#' @importFrom tidyselect all_of
#' @noRd
get_random_intercept_design = function(.data_, .sample, random_intercept_elements ){

  # Define the variables as NULL to avoid CRAN NOTES
  is_factor_continuous <- NULL
  design <- NULL
  max_mean_idx <- NULL
  max_minus_sum <- NULL
  max_factor_numeric <- NULL
  max_group_numeric <- NULL
  min_mean_idx <- NULL
  min_minus_sum <- NULL
  
  .sample = enquo(.sample)

  # If intercept is not defined create it
  if(nrow(random_intercept_elements) == 0 )
    return(
      random_intercept_elements |>
        mutate(
          design = list(),
          is_factor_continuous = logical()
        )
    )

  # Otherwise process
  random_intercept_elements |>
    mutate(is_factor_continuous = map_lgl(
      `factor`,
      ~ .x != "(Intercept)" && .data_ |> select(all_of(.x)) |> pull(1) |> is("numeric")
    )) |>
    mutate(design = pmap(
      list(grouping, `factor`, is_factor_continuous),
      ~ {

        # Make exception for random intercept
        if(..2 == "(Intercept)")
          .data_ = .data_ |> mutate(`(Intercept)` = 1)

        .data_ =
          .data_ |>
          select(!!.sample, ..1, ..2) |>
          set_names(c(quo_name(.sample), "group___", "factor___")) |>
          mutate(group___numeric = group___, factor___numeric = factor___) |>

          mutate(group___label := glue("{group___}___{.y}")) |>
          mutate(factor___ = ..2)


        # If factor is continuous
        if(..3)
          .data_ %>%

          # Mutate random intercept grouping to number
          mutate(group___numeric = factor(group___numeric) |> as.integer()) |>

          # If intercept is not defined create it
          mutate(., factor___numeric = 1L) |>

          # If categorical make sure the group is independent for factors
          mutate(mean_idx = glue("{group___numeric}") |> as.factor() |> as.integer()) |>
          mutate(mean_idx = if_else(mean_idx == max(mean_idx), 0L, mean_idx)) |>
          mutate(mean_idx = as.factor(mean_idx) |> as.integer() |> subtract(1L)) |>
          mutate(minus_sum = if_else(mean_idx==0, 1L, 0L))

        #|>
        #  distinct()

        # If factor is discrete
        else
          .data_ %>%

          # Mutate random intercept grouping to number
          mutate(group___numeric = factor(group___numeric) |> as.integer()) |>

          # If categorical make sure the group is independent for factors
          mutate(mean_idx = glue("{factor___numeric}{group___numeric}") |> as.factor() |> as.integer()) |>
          with_groups(factor___numeric, ~ ..1 |> mutate(mean_idx = if_else(mean_idx == max(mean_idx), 0L, mean_idx))) |>
          mutate(mean_idx = as.factor(mean_idx) |> as.integer() |> subtract(1L)) |>
          mutate(minus_sum = if_else(mean_idx==0, as.factor(factor___numeric) |> as.integer(), 0L)) |>

          # drop minus_sum if we just have one group___numeric per factor
          with_groups(factor___numeric, ~ {
            if(length(unique(..1$group___numeric)) == 1) ..1 |> mutate(., minus_sum = 0)
            else ..1
          }) |>
          mutate(factor___numeric = as.factor(factor___numeric) |> as.integer())

        #|>
        #  distinct()
      }
    )) |>

    # Make indexes unique across parameters
    mutate(
      max_mean_idx = map_int(design, ~ ..1 |> pull(mean_idx) |> max()),
      max_minus_sum = map_int(design, ~ ..1 |> pull(minus_sum) |> max()),
      max_factor_numeric = map_int(design, ~ ..1 |> pull(factor___numeric) |> max()),
      max_group_numeric = map_int(design, ~ ..1 |> pull(group___numeric) |> max())
    ) |>
    mutate(
      min_mean_idx = cumsum(max_mean_idx) - max_mean_idx ,
      min_minus_sum = cumsum(max_minus_sum) - max_minus_sum,
      max_factor_numeric = cumsum(max_factor_numeric) - max_factor_numeric,
      max_group_numeric = cumsum(max_group_numeric) - max_group_numeric
    ) |>
    mutate(design = pmap(
      list(design, min_mean_idx, min_minus_sum, max_factor_numeric, max_group_numeric),
      ~ ..1 |>
        mutate(
          mean_idx = if_else(mean_idx>0, mean_idx + ..2, mean_idx),
          minus_sum = if_else(minus_sum>0, minus_sum + ..3, minus_sum),
          factor___numeric = factor___numeric + ..4,
          group___numeric = group___numeric + ..5

        )
    ))

}

#' @importFrom glue glue
#' @noRd
get_design_matrix = function(.data_spread, formula, .sample){

  .sample = enquo(.sample)

  design_matrix =
  	.data_spread %>%

    select(!!.sample, parse_formula(formula)) |>
  	mutate(across(where(is.numeric),  scale)) |>
    model.matrix(formula, data=_)

  rownames(design_matrix) = .data_spread |> pull(!!.sample)

  design_matrix
}


#' Check Random Intercept Design
#'
#' This function checks the validity of the random intercept design in the data.
#'
#' @param .data A data frame containing the data.
#' @param factor_names A character vector of factor names.
#' @param random_intercept_elements A data frame containing the random intercept elements.
#' @param formula The formula used for the model.
#' @param X The design matrix.
#' 
#' @return A data frame with the checked random intercept elements.
#' 
#' @importFrom dplyr nest
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr pull
#' @importFrom dplyr filter
#' @importFrom dplyr distinct
#' @importFrom dplyr set_names
#' @importFrom tidyr unite
#' @importFrom purrr map2
#' @importFrom stringr str_subset
#' @importFrom readr type_convert
#' @noRd
check_random_intercept_design = function(.data, factor_names, random_intercept_elements, formula, X){

  # Define the variables as NULL to avoid CRAN NOTES
  factors <- NULL
  groupings <- NULL
  
  
  .data_ = .data

  # Loop across groupings
  random_intercept_elements |>
    nest(factors = `factor` ) |>
    mutate(checked = map2(
      grouping, factors,
      ~ {

        .y = unlist(.y)

        # Check that the group column is categorical
        stopifnot("sccomp says: the grouping column should be categorical (not numeric)" =
                    .data_ |>
                    select(all_of(.x)) |>
                    pull(1) |>
                    class() %in%
                    c("factor", "logical", "character")
        )


        # # Check sanity of the grouping if only random intercept
        # stopifnot(
        #   "sccomp says: the random intercept completely confounded with one or more discrete factors" =
        #     !(
        #       !.y |> equals("(Intercept)") &&
        #         .data_ |> select(any_of(.y)) |> suppressWarnings() |>  pull(1) |> class() %in% c("factor", "character") |> any() &&
        #         .data_ |>
        #         select(.x, any_of(.y)) |>
        #         select_if(\(x) is.character(x) | is.factor(x) | is.logical(x)) |>
        #         distinct() %>%
        #
        #         # TEMPORARY FIX
        #         set_names(c(colnames(.)[1], 'factor___temp')) |>
        #
        #         count(factor___temp) |>
        #         pull(n) |>
        #         equals(1) |>
        #         any()
        #     )
        # )

        # # Check if random intercept with random continuous slope. At the moment is not possible
        # # Because it would require I believe a multivariate prior
        # stopifnot(
        #   "sccomp says: continuous random slope is not supported yet" =
        #     !(
        #       .y |> str_subset("1", negate = TRUE) |> length() |> gt(0) &&
        #         .data_ |>
        #         select(
        #           .y |> str_subset("1", negate = TRUE)
        #         ) |>
        #         map_chr(class) %in%
        #         c("integer", "numeric")
        #     )
        # )

        # Check if random intercept with random continuous slope. At the moment is not possible
        # Because it would require I believe a multivariate prior
        stopifnot(
          "sccomp says: currently, discrete random slope is only supported in a intercept-free model. For example ~ 0 + treatment + (treatment | group)" =
            !(
              # If I have both random intercept and random discrete slope

                .y |> equals("(Intercept)") |> any() &&
                  length(.y) > 1 &&
                # If I have random slope and non-intercept-free model
                .data_ |> select(any_of(.y)) |> suppressWarnings() |>  pull(1) |> class() %in% c("factor", "character") |> any()

            )
        )


        # I HAVE TO REVESIT THIS
        #  stopifnot(
        #   "sccomp says: the groups in the formula (factor | group) should not be shared across factor groups" =
        #     !(
        #       # If I duplicated groups
        #       .y  |> identical("(Intercept)") |> not() &&
        #       .data_ |> select(.y |> setdiff("(Intercept)")) |> lapply(class) != "numeric" &&
        #         .data_ |>
        #         select(.x, .y |> setdiff("(Intercept)")) |>
        #
        #         # Drop the factor represented by the intercept if any
        #         mutate(`parameter` = .y |> setdiff("(Intercept)")) |>
        #         unite("factor_name", c(parameter, factor), sep = "", remove = FALSE) |>
        #         filter(factor_name %in% colnames(X)) |>
        #
        #         # Count
        #         distinct() %>%
        #         set_names(as.character(1:ncol(.))) |>
        #         count(`1`) |>
        #         filter(n>1) |>
        #         nrow() |>
        #         gt(1)
        #
        #     )
        # )

      }
    ))

  random_intercept_elements |>
    nest(groupings = grouping ) |>
    mutate(checked = map2(`factor`, groupings, ~{
      # Check the same group spans multiple factors
      stopifnot(
        "sccomp says: the groups in the formula (factor | group) should be present in only one factor, including the intercept" =
          !(
              # If I duplicated groups
            .y |> unlist() |> length() |> gt(1)

          )
      )


    }))




}

#' @importFrom purrr when
#' @importFrom stats model.matrix
#' @importFrom tidyr expand_grid
#' @importFrom stringr str_detect
#' @importFrom stringr str_remove_all
#' @importFrom purrr reduce
#' @importFrom stats as.formula
#'
#' @keywords internal
#' @noRd
#'
data_spread_to_model_input =
  function(
    .data_spread, formula, .sample, .cell_type, .count,
    truncation_ajustment = 1, approximate_posterior_inference ,
    formula_variability = ~ 1,
    contrasts = NULL,
    bimodal_mean_variability_association = FALSE,
    use_data = TRUE,
    random_intercept_elements){

    # Define the variables as NULL to avoid CRAN NOTES
    exposure <- NULL
    design <- NULL
    mat <- NULL
    factor___numeric <- NULL
    mean_idx <- NULL
    design_matrix <- NULL
    minus_sum <- NULL
    group___numeric <- NULL
    idx <- NULL
    group___label <- NULL
    parameter <- NULL
    group <- NULL
    design_matrix_col <- NULL
    
    # Prepare column same enquo
    .sample = enquo(.sample)
    .cell_type = enquo(.cell_type)
    .count = enquo(.count)
    .grouping_for_random_intercept =
      random_intercept_elements |>
      pull(grouping) |>
      unique() 
    
    if (length(.grouping_for_random_intercept)==0 ) .grouping_for_random_intercept = "random_intercept"


    X  =

    .data_spread |>
      get_design_matrix(
      # Drop random intercept
      formula |>
      as.character() |>
      str_remove_all("\\+ ?\\(.+\\|.+\\)") |>
      paste(collapse="") |>
      as.formula(),
       !!.sample
    )

    Xa  =
      .data_spread |>
      get_design_matrix(
      # Drop random intercept
      formula_variability |>
      as.character() |>
      str_remove_all("\\+ ?\\(.+\\|.+\\)") |>
      paste(collapse="") |>
      as.formula() ,
      !!.sample
    )

    XA = Xa %>%
      as_tibble() %>%
      distinct()

    A = ncol(XA);
    Ar = nrow(XA);





    factor_names = parse_formula(formula)
    factor_names_variability = parse_formula(formula_variability)
    cell_cluster_names = .data_spread %>% select(-!!.sample, -any_of(factor_names), -exposure, -!!.grouping_for_random_intercept) %>% colnames()

    # Random intercept
    if(nrow(random_intercept_elements)>0 ) {

      #check_random_intercept_design(.data_spread, any_of(factor_names), random_intercept_elements, formula, X)
      random_intercept_grouping = get_random_intercept_design2(.data_spread, !!.sample,  formula )

      # Actual parameters, excluding for the sum to one parameters
      N_random_intercepts = random_intercept_grouping |> mutate(n = map_int(design, ~ .x |> filter(mean_idx>0) |> distinct(mean_idx) |> nrow())) |> pull(n) |> sum()

      # Number of sum to one
      N_minus_sum = random_intercept_grouping |> mutate(n = map_int(design, ~ .x |> filter(minus_sum>0) |> distinct(minus_sum) |> nrow())) |> pull(n) |> sum()

      paring_cov_random_intercept =
        random_intercept_grouping |>
        mutate(mat = map(design, ~ .x |> distinct(factor___numeric, mean_idx) |> filter(mean_idx>0) )) |>
        select(mat) |>
        unnest(mat) |>
        arrange(factor___numeric, mean_idx) |>
        as_matrix()

      X_random_intercept =
        random_intercept_grouping |>
        mutate(design_matrix = map(
          design,
          ~ ..1 |>
            select(!!.sample, group___label, value) |>
            pivot_wider(names_from = group___label, values_from = value) |>
            mutate(across(everything(), ~ .x |> replace_na(0)))
        )) |>

        # Merge
        pull(design_matrix) |>
      	reduce(left_join, by = join_by(!!.sample)) |>
        as_matrix(rownames = quo_name(.sample))

    idx_group_random_intercepts =
      random_intercept_grouping |>
      mutate(design = map(design, ~ .x |> select(mean_idx, minus_sum, group___numeric, group___label))) |>
      select(design) |>
      unnest(design) |>

      mutate(minus_sum = -minus_sum) |>
      mutate(idx = mean_idx + minus_sum) |>
      distinct(group___numeric, idx, group___label) |>
      as_matrix(rownames = "group___label")


    N_grouping =
      random_intercept_grouping |>
      mutate(n = map_int(design, ~.x |> distinct(group___numeric) |> nrow())) |>
      pull(n) |> sum()

    
    # TEMPORARY
    group_factor_indexes_for_covariance = 
    	X_random_intercept |> 
    	colnames() |> 
    	enframe(value = "parameter", name = "order")  |> 
    	separate(parameter, c("factor", "group"), "___", remove = FALSE) |> 
    	complete(factor, group, fill = list(order=0)) |> 
    	select(-parameter) |> 
    	pivot_wider(names_from = group, values_from = order)  |> 
    	as_matrix(rownames = "factor")
    
    how_many_groups = ncol(group_factor_indexes_for_covariance )
    how_many_factors_in_random_design = nrow(group_factor_indexes_for_covariance )
    
    
    } else {
      X_random_intercept = matrix(rep(1, nrow(.data_spread)))[,0]
      N_random_intercepts = 0
      N_minus_sum = 0
      N_grouping =0
      paring_cov_random_intercept = matrix(c(1, 1), ncol = 2)[0,]
      idx_group_random_intercepts = matrix(c(1, 1), ncol = 2)[0,]
      how_many_groups = 0
      how_many_factors_in_random_design = 0
      group_factor_indexes_for_covariance = matrix()[0,0]
    }
    
    
    
    data_for_model =
      list(
        N = .data_spread %>% nrow(),
        M = .data_spread %>% select(-!!.sample, -any_of(factor_names), -exposure, -!!.grouping_for_random_intercept) %>% ncol(),
        exposure = .data_spread$exposure,
        y = .data_spread %>% select(-any_of(factor_names), -exposure, -!!.grouping_for_random_intercept) %>% as_matrix(rownames = quo_name(.sample)),
        X = X,
        XA = XA,
        Xa = Xa,
        C = ncol(X),
        A = A,
        Ar = Ar,
        truncation_ajustment = truncation_ajustment,
        is_vb = as.integer(approximate_posterior_inference),
        bimodal_mean_variability_association = bimodal_mean_variability_association,
        use_data = use_data,

        # Random intercept
        N_random_intercepts = N_random_intercepts,
        N_minus_sum = N_minus_sum,
        paring_cov_random_intercept = paring_cov_random_intercept,
        N_grouping = N_grouping,
        X_random_intercept = X_random_intercept,
        idx_group_random_intercepts = idx_group_random_intercepts,
        group_factor_indexes_for_covariance = group_factor_indexes_for_covariance,
        how_many_groups = how_many_groups,
        how_many_factors_in_random_design = how_many_factors_in_random_design,

        # For parallel chains
        grainsize = 1,
        
        ## LOO
        enable_loo = FALSE
      )

    # Add censoring
    data_for_model$is_truncated = 0
    data_for_model$truncation_up = matrix(rep(-1, data_for_model$M * data_for_model$N), ncol = data_for_model$M)
    data_for_model$truncation_down = matrix(rep(-1, data_for_model$M * data_for_model$N), ncol = data_for_model$M)
    data_for_model$truncation_not_idx = seq_len(data_for_model$M*data_for_model$N)
    data_for_model$TNS = length(data_for_model$truncation_not_idx)

    # Add parameter factor dictionary
    data_for_model$factor_parameter_dictionary = tibble()

    if(.data_spread  |> select(any_of(parse_formula(formula))) |> lapply(class) %in% c("factor", "character") |> any())
      data_for_model$factor_parameter_dictionary =
      data_for_model$factor_parameter_dictionary |> bind_rows(
        # For discrete
        .data_spread  |>
          select(any_of(parse_formula(formula)))  |>
          distinct()  |>

          # Drop numerical
          select_if(function(x) !is.numeric(x)) |>
          pivot_longer(everything(), names_to =  "factor", values_to = "parameter") %>%
          unite("design_matrix_col", c(`factor`, parameter), sep="", remove = FALSE)  |>
          select(-parameter) |>
          filter(design_matrix_col %in% colnames(data_for_model$X)) %>%
          distinct()

      )

 # For continuous
    if(.data_spread  |> select(all_of(parse_formula(formula))) |> lapply(class) |> equals("numeric") |> any())
      data_for_model$factor_parameter_dictionary =
      data_for_model$factor_parameter_dictionary |>
          bind_rows(
            tibble(
              design_matrix_col =  .data_spread  |>
                select(all_of(parse_formula(formula)))  |>
                distinct()  |>

                # Drop numerical
                select_if(function(x) is.numeric(x)) |>
                names()
            ) |>
              mutate(`factor` = design_matrix_col)
)

    # If constrasts is set it is a bit more complicated
    if(! is.null(contrasts))
      data_for_model$factor_parameter_dictionary =
        data_for_model$factor_parameter_dictionary |>
        distinct() |>
        expand_grid(parameter=contrasts) |>
        filter(str_detect(parameter, design_matrix_col )) |>
        select(-design_matrix_col) |>
        rename(design_matrix_col = parameter) |>
        distinct()

    data_for_model$intercept_in_design = X[,1] |> unique() |> identical(1)

    
    if (data_for_model$intercept_in_design | length(factor_names_variability) == 0) {
      data_for_model$A_intercept_columns = 1
    } else {
      data_for_model$A_intercept_columns = 
        .data_spread |> 
        select(any_of(factor_names[1])) |> 
        distinct() |> 
        nrow()
    }
    
    
    if (data_for_model$intercept_in_design ) {
      data_for_model$B_intercept_columns = 1
    } else {
      data_for_model$B_intercept_columns = 
        .data_spread |> 
        select(any_of(factor_names[1])) |> 
        distinct() |> 
        nrow()
    }
    
    # Return
    data_for_model
  }

data_to_spread = function(.data, formula, .sample, .cell_type, .count, .grouping_for_random_intercept){

  .sample = enquo(.sample)
  .cell_type = enquo(.cell_type)
  .count = enquo(.count)
  .grouping_for_random_intercept = .grouping_for_random_intercept |> map(~ .x |> quo_name() ) |> unlist()

  .data %>%
    nest(data = -!!.sample) %>%
    mutate(exposure = map_int(data, ~ .x %>% pull(!!.count) %>% sum() )) %>%
    unnest(data) %>%
    select(!!.sample, !!.cell_type, exposure, !!.count, parse_formula(formula), any_of(.grouping_for_random_intercept)) %>%
    spread(!!.cell_type, !!.count)

}

#' @importFrom purrr when
#' @importFrom stats model.matrix
#'
#' @keywords internal
#' @noRd
#'
data_simulation_to_model_input =
  function(.data, formula, .sample, .cell_type, .exposure, .coefficients, truncation_ajustment = 1, approximate_posterior_inference ){

    # Define the variables as NULL to avoid CRAN NOTES
    sd <- NULL
    . <- NULL
    
    
    # Prepare column same enquo
    .sample = enquo(.sample)
    .cell_type = enquo(.cell_type)
    .exposure = enquo(.exposure)
    .coefficients = enquo(.coefficients)

    factor_names = parse_formula(formula)

    sample_data =
      .data %>%
      select(!!.sample, any_of(factor_names)) %>%
      distinct() %>%
      arrange(!!.sample)
    X =
      sample_data %>%
      model.matrix(formula, data=.) %>%
      apply(2, function(x) {
        
        if(sd(x)==0 ) x
        else x |> scale(scale=FALSE)
        
      } ) %>%
      {
        .x = (.)
        rownames(.x) = sample_data %>% pull(!!.sample)
        .x
      }

    if(factor_names == "1") XA = X[,1, drop=FALSE]
    else XA = X[,c(1,2), drop=FALSE]
    
    XA = XA |> 
      as_tibble()  |> 
      distinct()

    cell_cluster_names =
      .data %>%
      distinct(!!.cell_type) %>%
      arrange(!!.cell_type) %>%
      pull(!!.cell_type)

    coefficients =
      .data %>%
      select(!!.cell_type, !!.coefficients) %>%
      unnest(!!.coefficients) %>%
      distinct() %>%
      arrange(!!.cell_type) %>%
      as_matrix(rownames = quo_name(.cell_type)) %>%
      t()

    list(
      N = .data %>% distinct(!!.sample) %>% nrow(),
      M = .data %>% distinct(!!.cell_type) %>% nrow(),
      exposure = .data %>% distinct(!!.sample, !!.exposure) %>% arrange(!!.sample) %>% pull(!!.exposure),
      X = X,
      XA = XA,
      C = ncol(X),
      A =  ncol(XA),
      beta = coefficients
    )

  }


#' Choose the number of chains baed on how many draws we need from the posterior distribution
#' Because there is a fix cost (warmup) to starting a new chain,
#' we need to use the minimum amount that we can parallelise
#' @param how_many_posterior_draws A real number of posterior draws needed
#' @param max_number_to_check A sane upper plateau
#'
#' @keywords internal
#' @noRd
#'
#' @return A Stan fit object
find_optimal_number_of_chains = function(how_many_posterior_draws = 100,
                                         max_number_to_check = 100, warmup = 200, parallelisation_start_penalty = 100) {



  # Define the variables as NULL to avoid CRAN NOTES
  chains <- NULL


  chains_df =
    tibble(chains = seq_len(max_number_to_check)) %>%
    mutate(tot = (how_many_posterior_draws / chains) + warmup + (parallelisation_start_penalty * chains))

  d1 <- diff(chains_df$tot) / diff(seq_len(nrow(chains_df))) # first derivative
  abs(d1) %>% order() %>% .[1] # Find derivative == 0


}


get.elbow.points.indices <- function(x, y, threshold) {
  # From https://stackoverflow.com/questions/41518870/finding-the-elbow-knee-in-a-curve
  d1 <- diff(y) / diff(x) # first derivative
  d2 <- diff(d1) / diff(x[-1]) # second derivative
  indices <- which(abs(d2) > threshold)
  return(indices)
}

#' @importFrom magrittr divide_by
#' @importFrom magrittr multiply_by
#' @importFrom stats C
#'
#' @keywords internal
#' @noRd
#'
get_probability_non_zero_OLD = function(.data, prefix = "", test_above_logit_fold_change = 0){

  # Define the variables as NULL to avoid CRAN NOTES
  .draw <- NULL
  M <- NULL
  C_name <- NULL
  bigger_zero <- NULL
  smaller_zero <- NULL
  
  probability_column_name = sprintf("%s_prob_H0", prefix) %>% as.symbol()

  total_draws = .data %>% pull(2) %>% .[[1]] %>% distinct(.draw) %>% nrow()

  .data %>%
    unnest(2 ) %>%
    filter(C ==2) %>%
    nest(data = -c(M, C_name)) %>%
    mutate(
      bigger_zero = map_int(data, ~ .x %>% filter(.value>test_above_logit_fold_change) %>% nrow),
      smaller_zero = map_int(data, ~ .x %>% filter(.value< -test_above_logit_fold_change) %>% nrow)
    ) %>%
    rowwise() %>%
    mutate(
      !!probability_column_name :=
        1 - (
          max(bigger_zero, smaller_zero) %>%
            #max(1) %>%
            divide_by(total_draws)
        )
    )  %>%
    ungroup() %>%
    select(M, !!probability_column_name)
  # %>%
  # mutate(false_discovery_rate = cummean(prob_non_zero))

}

#' @importFrom magrittr divide_by
#' @importFrom magrittr multiply_by
#' @importFrom stats C
#' @importFrom stats setNames
#'
#' @keywords internal
#' @noRd
#'
get_probability_non_zero_ = function(fit, parameter, prefix = "", test_above_logit_fold_change = 0){

  # Define the variables as NULL to avoid CRAN NOTES
  M <- NULL
  C_name <- NULL
  bigger_zero <- NULL
  smaller_zero <- NULL
  
  

  draws = fit$draws(
      variables = parameter,
      inc_warmup = FALSE,
      format = getOption("cmdstanr_draws_format", "draws_matrix")
    )

  total_draws = dim(draws)[1]

  bigger_zero =
    draws %>%
      apply(2, function(x) (x>test_above_logit_fold_change) %>% which %>% length)

  smaller_zero =
    draws %>%
        apply(2, function(x) (x< -test_above_logit_fold_change) %>% which %>% length)

  (1 - (pmax(bigger_zero, smaller_zero) / total_draws)) %>%
    enframe() %>%
    tidyr::extract(name, c("C", "M"), ".+\\[([0-9]+),([0-9]+)\\]") %>%
    mutate(across(c(C, M), ~ as.integer(.x))) %>%
    tidyr::spread(C, value)

}

get_probability_non_zero = function(draws, test_above_logit_fold_change = 0, probability_column_name){

  draws %>%
    with_groups(c(M, C_name), ~ .x |> summarise(
      bigger_zero = which(.value>test_above_logit_fold_change) |> length(),
      smaller_zero = which(.value< -test_above_logit_fold_change) |> length(),
      n=n()
    )) |>
    mutate(!!as.symbol(probability_column_name) :=  (1 - (pmax(bigger_zero, smaller_zero) / n)))

}

#' @keywords internal
#' @noRd
#'
parse_generated_quantities = function(rng, number_of_draws = 1){

  # Define the variables as NULL to avoid CRAN NOTES
  .draw <- NULL
  N <- NULL
  .value <- NULL
  generated_counts <- NULL
  M <- NULL
  generated_proportions <- NULL
  
  draws_to_tibble_x_y(rng, "counts", "N", "M", number_of_draws) %>%
    with_groups(c(.draw, N), ~ .x %>% mutate(generated_proportions = .value/max(1, sum(.value)))) %>%
    filter(.draw<= number_of_draws) %>%
    rename(generated_counts = .value, replicate = .draw) %>%

    mutate(generated_counts = as.integer(generated_counts)) %>%
    select(M, N, generated_proportions, generated_counts, replicate)

}

#' design_matrix_and_coefficients_to_simulation
#'
#' @description Create simulation from design matrix and coefficient matrix
#'
#' @importFrom dplyr left_join
#' @importFrom tidyr expand_grid
#'
#' @keywords internal
#' @noRd
#'
#' @param design_matrix A matrix
#' @param coefficient_matrix A matrix
#'
#' @return A data frame
#'
#'
#'
design_matrix_and_coefficients_to_simulation = function(

  design_matrix, coefficient_matrix, .estimate_object

){

  # Define the variables as NULL to avoid CRAN NOTES
  cell_type <- NULL
  beta_1 <- NULL
  beta_2 <- NULL
  
  design_df = as.data.frame(design_matrix)
  coefficient_df = as.data.frame(coefficient_matrix)

  rownames(design_df) = sprintf("sample_%s", seq_len(nrow(design_df)))
  colnames(design_df) = sprintf("factor_%s", seq_len(ncol(design_df)))

  rownames(coefficient_df) = sprintf("cell_type_%s", seq_len(nrow(coefficient_df)))
  colnames(coefficient_df) = sprintf("beta_%s", seq_len(ncol(coefficient_df)))

  input_data =
    expand_grid(
      sample = rownames(design_df),
      cell_type = rownames(coefficient_df)
    ) |>
    left_join(design_df |> as_tibble(rownames = "sample") , by = "sample") |>
    left_join(coefficient_df |>as_tibble(rownames = "cell_type"), by = "cell_type")

  simulate_data(.data = input_data,

                .estimate_object = .estimate_object,

                formula_composition = ~ factor_1 ,
                .sample = sample,
                .cell_group = cell_type,
                .coefficients = c(beta_1, beta_2),
                mcmc_seed = sample(1e5, 1)
  )


}

design_matrix_and_coefficients_to_dir_mult_simulation =function(design_matrix, coefficient_matrix, precision = 100, seed = sample(1:100000, size = 1)){

  # Define the variables as NULL to avoid CRAN NOTES
  cell_type <- NULL
  generated_counts <- NULL
  factor_1 <- NULL

  # design_df = as.data.frame(design_matrix)
  # coefficient_df = as.data.frame(coefficient_matrix)
  #
  # rownames(design_df) = sprintf("sample_%s", 1:nrow(design_df))
  # colnames(design_df) = sprintf("factor_%s", 1:ncol(design_df))
  #
  # rownames(coefficient_df) = sprintf("cell_type_%s", 1:nrow(coefficient_df))
  # colnames(coefficient_df) = sprintf("beta_%s", 1:ncol(coefficient_df))

  exposure = 500

  prop.means =
    design_matrix %*%
    t(coefficient_matrix) %>%
    boot::inv.logit()

  extraDistr::rdirmnom(length(design_matrix), exposure, prop.means * precision) %>%
    as_tibble(.name_repair = "unique", rownames = "sample") %>%
    mutate(factor_1= design_matrix) %>%
    gather(cell_type, generated_counts, -sample, -factor_1) %>%
    mutate(generated_counts = as.integer(generated_counts))


}

#' @importFrom rlang ensym
#' @noRd
class_list_to_counts = function(.data, .sample, .cell_group){

  .sample_for_tidyr = ensym(.sample)
  .cell_group_for_tidyr = ensym(.cell_group)

  .sample = enquo(.sample)
  .cell_group = enquo(.cell_group)

  .data %>%
    count(!!.sample,
          !!.cell_group,
          name = "count") %>%

    complete(
      !!.sample_for_tidyr,!!.cell_group_for_tidyr,
      fill = list(count = 0)
    ) %>%
    mutate(count = as.integer(count))
}

#' @importFrom dplyr cummean
#' @noRd
get_FDR = function(x){
  
  # Define the variables as NULL to avoid CRAN NOTES
  value <- NULL
  name <- NULL
  FDR <- NULL
  
  
  enframe(x) %>%
    arrange(value) %>%
    mutate(FDR = cummean(value)) %>%
    arrange(name) %>%
    pull(FDR)
}

#' Plot 1D Intervals for Cell-group Effects
#'
#' This function creates a series of 1D interval plots for cell-group effects, highlighting significant differences based on a given significance threshold.
#'
#' @param .data Data frame containing the main data.
#' @param significance_threshold Numeric value specifying the significance threshold for highlighting differences. Default is 0.025.
#' @param test_composition_above_logit_fold_change A positive integer. It is the effect threshold used for the hypothesis test. A value of 0.2 correspond to a change in cell proportion of 10% for a cell type with baseline proportion of 50%. That is, a cell type goes from 45% to 50%. When the baseline proportion is closer to 0 or 1 this effect thrshold has consistent value in the logit uncontrained scale.
#' @importFrom patchwork wrap_plots
#' @importFrom forcats fct_reorder
#' @importFrom tidyr drop_na
#' 
#' @export
#' 
#' @return A combined plot of 1D interval plots.
#' @examples
#' # Example usage:
#' # plot_1D_intervals(.data, "cell_group", 0.025, theme_minimal())
plot_1D_intervals = function(.data, significance_threshold = 0.05, test_composition_above_logit_fold_change = .data |> attr("test_composition_above_logit_fold_change")){
  
  # Define the variables as NULL to avoid CRAN NOTES
  parameter <- NULL
  estimate <- NULL
  value <- NULL
  
  .cell_group = attr(.data, ".cell_group")
  
  # Check if test have been done
  if(.data |> select(ends_with("FDR")) |> ncol() |> equals(0))
    stop("sccomp says: to produce plots, you need to run the function sccomp_test() on your estimates.")
  
  
  plot_list = 
    .data |>
    filter(parameter != "(Intercept)") |>
    
    # Reshape data
    select(-contains("n_eff"), -contains("R_k_hat")) |> 
    pivot_longer(c(contains("c_"), contains("v_")), names_sep = "_", names_to = c("which", "estimate")) |>
    pivot_wider(names_from = estimate, values_from = value) |>
    
    # Nest data by parameter and which
    nest(data = -c(parameter, which)) |>
    mutate(plot = pmap(
      list(data, which, parameter),
      ~  {

        # Check if there are any statistics to plot
        if(..1 |> filter(!effect |> is.na()) |> nrow() |> equals(0))
          return(NA)
        
        # Create ggplot for each nested data
        ggplot(..1, aes(x = effect, y = fct_reorder(!!.cell_group, effect))) +
          geom_vline(xintercept = test_composition_above_logit_fold_change, colour = "grey") +
          geom_vline(xintercept = -test_composition_above_logit_fold_change, colour = "grey") +
          geom_errorbar(aes(xmin = lower, xmax = upper, color = FDR < significance_threshold)) +
          geom_point() +
          scale_color_brewer(palette = "Set1") +
          xlab("Credible interval of the slope") +
          ylab("Cell group") +
          ggtitle(sprintf("%s %s", ..2, ..3)) +
          multipanel_theme +
          theme(legend.position = "bottom") 
      }
    )) %>%
    
    # Filter out NA plots
    filter(!plot |> is.na()) |> 
    pull(plot) 
  
  # Combine all individual plots into one plot
  plot_list |>
    wrap_plots(ncol = plot_list |> length() |> sqrt() |> ceiling())
}


#' Plot 2D Intervals for Mean-Variance Association
#'
#' This function creates a 2D interval plot for mean-variance association, highlighting significant differences based on a given significance threshold.
#'
#' @param .data Data frame containing the main data.
#' @param significance_threshold Numeric value specifying the significance threshold for highlighting differences. Default is 0.025.
#' @param test_composition_above_logit_fold_change A positive integer. It is the effect threshold used for the hypothesis test. A value of 0.2 correspond to a change in cell proportion of 10% for a cell type with baseline proportion of 50%. That is, a cell type goes from 45% to 50%. When the baseline proportion is closer to 0 or 1 this effect thrshold has consistent value in the logit uncontrained scale.
#' 
#' 
#' @importFrom dplyr filter arrange mutate if_else row_number
#' @importFrom ggplot2 ggplot geom_vline geom_hline geom_errorbar geom_point annotate geom_text_repel aes facet_wrap
#' @importFrom scales trans_new
#' @importFrom stringr str_replace
#' @importFrom stats quantile
#' @importFrom magrittr equals
#' 
#' @export
#' 
#' @return A ggplot object representing the 2D interval plot.
#' @examples
#' # Example usage:
#' # plot_2D_intervals(.data, "cell_group", theme_minimal(), 0.025)
plot_2D_intervals = function(.data, significance_threshold = 0.05, test_composition_above_logit_fold_change = .data |> attr("test_composition_above_logit_fold_change")){
  
  # Define the variables as NULL to avoid CRAN NOTES
  v_effect <- NULL
  parameter <- NULL
  c_effect <- NULL
  c_lower <- NULL
  c_upper <- NULL
  c_FDR <- NULL
  v_lower <- NULL
  v_upper <- NULL
  v_FDR <- NULL
  cell_type_label <- NULL
  multipanel_theme <- NULL
  
  
  .cell_group = attr(.data, ".cell_group")
  
  # Check if test have been done
  if(.data |> select(ends_with("FDR")) |> ncol() |> equals(0))
    stop("sccomp says: to produce plots, you need to run the function sccomp_test() on your estimates.")
  
  
  # Mean-variance association
  .data %>%
    
    # Filter where variance is inferred
    filter(!is.na(v_effect)) %>%
    
    # Add labels for significant cell groups
    with_groups(
      parameter,
      ~ .x %>%
        arrange(c_FDR) %>%
        mutate(cell_type_label = if_else(row_number() <= 3 & c_FDR < significance_threshold & parameter != "(Intercept)", !!.cell_group, ""))
    ) %>%
    with_groups(
      parameter,
      ~ .x %>%
        arrange(v_FDR) %>%
        mutate(cell_type_label = if_else((row_number() <= 3 & v_FDR < significance_threshold & parameter != "(Intercept)"), !!.cell_group, cell_type_label))
    ) %>%
    
    {
      .x = (.)
      
      # Plot
      ggplot(.x, aes(c_effect, v_effect)) +
        
        # Add vertical and horizontal lines
        geom_vline(xintercept = c(-test_composition_above_logit_fold_change, test_composition_above_logit_fold_change), colour = "grey", linetype = "dashed", linewidth = 0.3) +
        geom_hline(yintercept = c(-test_composition_above_logit_fold_change, test_composition_above_logit_fold_change), colour = "grey", linetype = "dashed", linewidth = 0.3) +
        
        # Add error bars
        geom_errorbar(aes(xmin = `c_lower`, xmax = `c_upper`, color = `c_FDR` < significance_threshold, alpha = `c_FDR` < significance_threshold), linewidth = 0.2) +
        geom_errorbar(aes(ymin = v_lower, ymax = v_upper, color = `v_FDR` < significance_threshold, alpha = `v_FDR` < significance_threshold), linewidth = 0.2) +
        
        # Add points
        geom_point(size = 0.2) +
        
        # Add annotations
        annotate("text", x = 0, y = 3.5, label = "Variable", size = 2) +
        annotate("text", x = 5, y = 0, label = "Abundant", size = 2, angle = 270) +
        
        # Add text labels for significant cell groups
        geom_text_repel(aes(c_effect, -v_effect, label = cell_type_label), size = 2.5, data = .x %>% filter(cell_type_label != "")) +
        
        # Set color and alpha scales
        scale_color_manual(values = c("#D3D3D3", "#E41A1C")) +
        scale_alpha_manual(values = c(0.4, 1)) +
        
        # Facet by parameter
        facet_wrap(~parameter, scales = "free") +
        
        # Apply custom theme
        multipanel_theme 
    }
}


#' Plot Boxplot of Cell-group Proportion
#'
#' This function creates a boxplot of cell-group proportions, optionally highlighting significant differences based on a given significance threshold.
#'
#' @param .data Data frame containing the main data.
#' @param data_proportion Data frame containing proportions of cell groups.
#' @param factor_of_interest A factor indicating the biological condition of interest.
#' @param .cell_group The cell group to be analysed.
#' @param .sample The sample identifier.
#' @param significance_threshold Numeric value specifying the significance threshold for highlighting differences. Default is 0.025.
#' @param my_theme A ggplot2 theme object to be applied to the plot.
#' @importFrom scales trans_new
#' @importFrom stringr str_replace
#' @importFrom stats quantile
#' 
#' 
#' @return A ggplot object representing the boxplot.
#' @examples
#' # Example usage:
#' # plot_boxplot(.data, data_proportion, "condition", "cell_group", "sample", 0.025, theme_minimal())
plot_boxplot = function(
    .data, data_proportion, factor_of_interest, .cell_group,
    .sample, significance_threshold = 0.05, my_theme
){
  
  # Define the variables as NULL to avoid CRAN NOTES
  stats_name <- NULL
  parameter <- NULL
  stats_value <- NULL
  count_data <- NULL
  generated_proportions <- NULL
  proportion <- NULL
  name <- NULL
  outlier <- NULL
  
  # Function to calculate boxplot statistics
  calc_boxplot_stat <- function(x) {
    coef <- 1.5
    n <- sum(!is.na(x))
    
    # Calculate quantiles
    stats <- quantile(x, probs = c(0.0, 0.25, 0.5, 0.75, 1.0))
    names(stats) <- c("ymin", "lower", "middle", "upper", "ymax")
    iqr <- diff(stats[c(2, 4)])
    
    # Set whiskers
    outliers <- x < (stats[2] - coef * iqr) | x > (stats[4] + coef * iqr)
    if (any(outliers)) {
      stats[c(1, 5)] <- range(c(stats[2:4], x[!outliers]), na.rm = TRUE)
    }
    return(stats)
  }
  
  # Function to remove leading zero from labels
  dropLeadingZero <- function(l){  stringr::str_replace(l, '0(?=.)', '') }
  
  # Define square root transformation and its inverse
  S_sqrt <- function(x){sign(x)*sqrt(abs(x))}
  IS_sqrt <- function(x){x^2*sign(x)}
  S_sqrt_trans <- function() scales::trans_new("S_sqrt",S_sqrt,IS_sqrt)
  
  .cell_group = enquo(.cell_group)
  .sample = enquo(.sample)
  
  # Prepare significance colors
  significance_colors =
    .data %>%
    pivot_longer(
      c(contains("c_"), contains("v_")),
      names_pattern = "([cv])_([a-zA-Z0-9]+)",
      names_to = c("which", "stats_name"),
      values_to = "stats_value"
    ) %>%
    filter(stats_name == "FDR") %>%
    filter(parameter != "(Intercept)") %>%
    filter(stats_value < significance_threshold) %>%
    filter(`factor` == factor_of_interest) 
  
  if(nrow(significance_colors) > 0){
    
    if(.data |> attr("contrasts") |> is.null())
      significance_colors =
        significance_colors %>%
        unite("name", c(which, parameter), remove = FALSE) %>%
        distinct() %>%
        
        # Get clean parameter
        mutate(!!as.symbol(factor_of_interest) := str_replace(parameter, sprintf("^%s", `factor`), "")) %>%
        with_groups(c(!!.cell_group, !!as.symbol(factor_of_interest)), ~ .x %>% summarise(name = paste(name, collapse = ", ")))
    else
      significance_colors =
        significance_colors |>
        mutate(count_data = map(count_data, ~ .x |> select(all_of(factor_of_interest)) |> distinct())) |>
        unnest(count_data) |>
        
        # Filter relevant parameters
        mutate( !!as.symbol(factor_of_interest) := as.character(!!as.symbol(factor_of_interest) ) ) |>
        filter(str_detect(parameter, !!as.symbol(factor_of_interest) )) |>
        
        # Rename
        select(!!.cell_group, !!as.symbol(factor_of_interest), name = parameter) |>
        
        # Merge contrasts
        with_groups(c(!!.cell_group, !!as.symbol(factor_of_interest)), ~ .x %>% summarise(name = paste(name, collapse = ", ")))
  }
  
  my_boxplot = ggplot()
  
  if("fit" %in% names(attributes(.data))){
    
    simulated_proportion =
      .data |>
      sccomp_replicate(number_of_draws = 100) |>
      left_join(data_proportion %>% distinct(!!as.symbol(factor_of_interest), !!.sample, !!.cell_group))
    
    my_boxplot = my_boxplot +
      
      # Add boxplot for simulated proportions
      stat_summary(
        aes(!!as.symbol(factor_of_interest), (generated_proportions)),
        fun.data = calc_boxplot_stat, geom="boxplot",
        outlier.shape = NA, outlier.color = NA,outlier.size = 0,
        fatten = 0.5, lwd=0.2,
        data =
          simulated_proportion %>%
          inner_join(data_proportion %>% distinct(!!as.symbol(factor_of_interest), !!.cell_group)),
        color="blue"
      )
  }
  
  if(nrow(significance_colors) == 0 |
     length(intersect(
       significance_colors |> pull(!!as.symbol(factor_of_interest)),
       data_proportion |> pull(!!as.symbol(factor_of_interest))
     )) == 0){
    
    my_boxplot=
      my_boxplot +
      
      # Add boxplot without significance colors
      geom_boxplot(
        aes(!!as.symbol(factor_of_interest), proportion,  group=!!as.symbol(factor_of_interest), fill = NULL),
        outlier.shape = NA, outlier.color = NA,outlier.size = 0,
        data =
          data_proportion |>
          mutate(!!as.symbol(factor_of_interest) := as.character(!!as.symbol(factor_of_interest))) %>%
          left_join(significance_colors, by = c(quo_name(.cell_group), factor_of_interest)),
        fatten = 0.5,
        lwd=0.5,
      )
  } else {
    my_boxplot=
      my_boxplot +
      
      # Add boxplot with significance colors
      geom_boxplot(
        aes(!!as.symbol(factor_of_interest), proportion,  group=!!as.symbol(factor_of_interest), fill = name),
        outlier.shape = NA, outlier.color = NA,outlier.size = 0,
        data =
          data_proportion |>
          mutate(!!as.symbol(factor_of_interest) := as.character(!!as.symbol(factor_of_interest))) %>%
          left_join(significance_colors, by = c(quo_name(.cell_group), factor_of_interest)),
        fatten = 0.5,
        lwd=0.5,
      )
  }
  
  my_boxplot +
    
    # Add jittered points for individual data
    geom_jitter(
      aes(!!as.symbol(factor_of_interest), proportion, shape=outlier, color=outlier,  group=!!as.symbol(factor_of_interest)),
      data = data_proportion,
      position=position_jitterdodge(jitter.height = 0, jitter.width = 0.2),
      size = 0.5
    ) +
    
    # Facet wrap by cell group
    facet_wrap(
      vars(!!.cell_group),
      scales = "free_y",
      nrow = 4
    ) +
    scale_color_manual(values = c("black", "#e11f28")) +
    scale_y_continuous(trans=S_sqrt_trans(), labels = dropLeadingZero) +
    scale_fill_discrete(na.value = "white") +
    xlab("Biological condition") +
    ylab("Cell-group proportion") +
    guides(color="none", alpha="none", size="none") +
    labs(fill="Significant difference") +
    ggtitle("Note: Be careful judging significance (or outliers) visually for lowly abundant cell groups. \nVisualising proportion hides the uncertainty characteristic of count data, that a count-based statistical model can estimate.") +
    my_theme +
    theme(axis.text.x =  element_text(angle=20, hjust = 1), title = element_text(size = 3))
}

#' Plot Scatterplot of Cell-group Proportion
#'
#' This function creates a scatterplot of cell-group proportions, optionally highlighting significant differences based on a given significance threshold.
#'
#' @param .data Data frame containing the main data.
#' @param data_proportion Data frame containing proportions of cell groups.
#' @param factor_of_interest A factor indicating the biological condition of interest.
#' @param .cell_group The cell group to be analysed.
#' @param .sample The sample identifier.
#' @param significance_threshold Numeric value specifying the significance threshold for highlighting differences. Default is 0.025.
#' @param my_theme A ggplot2 theme object to be applied to the plot.
#' @importFrom scales trans_new
#' @importFrom stringr str_replace
#' @importFrom stats quantile
#' @importFrom magrittr equals
#' 
#' 
#' @return A ggplot object representing the scatterplot.
#' @examples
#' # Example usage:
#' # plot_scatterplot(.data, data_proportion, "condition", "cell_group", "sample", 0.025, theme_minimal())
plot_scatterplot = function(
    .data, data_proportion, factor_of_interest, .cell_group,
    .sample, significance_threshold = 0.05, my_theme
){
  
  # Define the variables as NULL to avoid CRAN NOTES
  stats_name <- NULL
  parameter <- NULL
  stats_value <- NULL
  count_data <- NULL
  generated_proportions <- NULL
  proportion <- NULL
  name <- NULL
  outlier <- NULL
  
  # Function to remove leading zero from labels
  dropLeadingZero <- function(l){  stringr::str_replace(l, '0(?=.)', '') }
  
  # Define square root transformation and its inverse
  S_sqrt <- function(x){sign(x)*sqrt(abs(x))}
  IS_sqrt <- function(x){x^2*sign(x)}
  S_sqrt_trans <- function() scales::trans_new("S_sqrt",S_sqrt,IS_sqrt)
  
  .cell_group = enquo(.cell_group)
  .sample = enquo(.sample)
  
  # Prepare significance colors
  significance_colors =
    .data %>%
    pivot_longer(
      c(contains("c_"), contains("v_")),
      names_pattern = "([cv])_([a-zA-Z0-9]+)",
      names_to = c("which", "stats_name"),
      values_to = "stats_value"
    ) %>%
    filter(stats_name == "FDR") %>%
    filter(parameter != "(Intercept)") %>%
    filter(stats_value < significance_threshold) %>%
    filter(`factor` == factor_of_interest) 
  
  if(nrow(significance_colors) > 0){
    
    if(.data |> attr("contrasts") |> is.null())
      significance_colors =
        significance_colors %>%
        unite("name", c(which, parameter), remove = FALSE) %>%
        distinct() %>%
        
        # Get clean parameter
        mutate(!!as.symbol(factor_of_interest) := str_replace(parameter, sprintf("^%s", `factor`), "")) %>%
        with_groups(c(!!.cell_group, !!as.symbol(factor_of_interest)), ~ .x %>% summarise(name = paste(name, collapse = ", ")))
    else
      significance_colors =
        significance_colors |>
        mutate(count_data = map(count_data, ~ .x |> select(all_of(factor_of_interest)) |> distinct())) |>
        unnest(count_data) |>
        
        # Filter relevant parameters
        mutate( !!as.symbol(factor_of_interest) := as.character(!!as.symbol(factor_of_interest) ) ) |>
        filter(str_detect(parameter, !!as.symbol(factor_of_interest) )) |>
        
        # Rename
        select(!!.cell_group, !!as.symbol(factor_of_interest), name = parameter) |>
        
        # Merge contrasts
        with_groups(c(!!.cell_group, !!as.symbol(factor_of_interest)), ~ .x %>% summarise(name = paste(name, collapse = ", ")))
  }
  
  my_scatterplot = ggplot()
  
  if("fit" %in% names(attributes(.data))){
    
    simulated_proportion =
      .data |>
      sccomp_replicate(number_of_draws = 1000) |>
      left_join(data_proportion %>% distinct(!!as.symbol(factor_of_interest), !!.sample, !!.cell_group))
    
    my_scatterplot = 
      my_scatterplot +
      
      # Add smoothed line for simulated proportions
      geom_smooth(
        aes(!!as.symbol(factor_of_interest), (generated_proportions)),
        lwd=0.2,
        data =
          simulated_proportion %>%
          inner_join(data_proportion %>% distinct(!!as.symbol(factor_of_interest), !!.cell_group, !!.sample)) ,
        color="blue", fill="blue",
        span = 1
      )
  }
  
  if(
    nrow(significance_colors)==0 ||
    
    significance_colors |> 
    pull(!!as.symbol(factor_of_interest)) |> 
    intersect(
      data_proportion |> 
      pull(!!as.symbol(factor_of_interest))
    ) |> 
    length() |> 
    equals(0)
  ) {
    
    my_scatterplot=
      my_scatterplot +
      
      # Add smoothed line without significance colors
      geom_smooth(
        aes(!!as.symbol(factor_of_interest), proportion, fill = NULL),
        data =
          data_proportion ,
        lwd=0.5,
        color = "black",
        span = 1
      )
  } else {
    my_scatterplot=
      my_scatterplot +
      
      # Add smoothed line with significance colors
      geom_smooth(
        aes(!!as.symbol(factor_of_interest), proportion, fill = name),
        outlier.shape = NA, outlier.color = NA,outlier.size = 0,
        data = data_proportion ,
        fatten = 0.5,
        lwd=0.5,
        color = "black",
        span = 1
      )
  }
  
  my_scatterplot +
    
    # Add jittered points for individual data
    geom_point(
      aes(!!as.symbol(factor_of_interest), proportion, shape=outlier, color=outlier),
      data = data_proportion,
      position=position_jitterdodge(jitter.height = 0, jitter.width = 0.2),
      size = 0.5
    ) +
    
    # Facet wrap by cell group
    facet_wrap(
      vars(!!.cell_group),
      scales = "free_y",
      nrow = 4
    ) +
    scale_color_manual(values = c("black", "#e11f28")) +
    scale_y_continuous(trans=S_sqrt_trans(), labels = dropLeadingZero) +
    scale_fill_discrete(na.value = "white") +
    xlab("Biological condition") +
    ylab("Cell-group proportion") +
    guides(color="none", alpha="none", size="none") +
    labs(fill="Significant difference") +
    ggtitle("Note: Be careful judging significance (or outliers) visually for lowly abundant cell groups. \nVisualising proportion hides the uncertainty characteristic of count data, that a count-based statistical model can estimate.") +
    my_theme +
    theme(axis.text.x =  element_text(angle=20, hjust = 1), title = element_text(size = 3))
}

draws_to_statistics = function(draws, false_positive_rate, test_composition_above_logit_fold_change, .cell_group, prefix = ""){

  # Define the variables as NULL to avoid CRAN NOTES
  M <- NULL
  parameter <- NULL
  bigger_zero <- NULL
  smaller_zero <- NULL
  lower <- NULL
  effect <- NULL
  upper <- NULL
  pH0 <- NULL
  FDR <- NULL
  n_eff <- NULL
  R_k_hat <- NULL
  
  .cell_group = enquo(.cell_group)

  draws =
    draws |>
    with_groups(c(!!.cell_group, M, parameter), ~ .x |> summarise(
      lower = quantile(.value, false_positive_rate/2),
      effect = quantile(.value, 0.5),
      upper = quantile(.value, 1-(false_positive_rate/2)),
      bigger_zero = which(.value>test_composition_above_logit_fold_change) |> length(),
      smaller_zero = which(.value< -test_composition_above_logit_fold_change) |> length(),
      R_k_hat = unique(R_k_hat),
      n_eff = unique(n_eff),
      n=n()
    )) |>

    # Calculate probability non 0
    mutate(pH0 =  (1 - (pmax(bigger_zero, smaller_zero) / n))) |>
    with_groups(parameter, ~ mutate(.x, FDR = get_FDR(pH0))) |>

    select(!!.cell_group, M, parameter, lower, effect, upper, pH0, FDR, any_of(c("n_eff", "R_k_hat"))) |>
    suppressWarnings()

  # Setting up names separately because |> is not flexible enough
  draws |>
    setNames(c(colnames(draws)[1:3], sprintf("%s%s", prefix, colnames(draws)[4:ncol(draws)])))
}

enquos_from_list_of_symbols <- function(...) {
  enquos(...)
}

contrasts_to_enquos = function(contrasts){
  
  # Define the variables as NULL to avoid CRAN NOTES
  . <- NULL
  
  contrasts |> enquo() |> quo_names() |> syms() %>% do.call(enquos_from_list_of_symbols, .)
}

#' Mutate Data Frame Based on Expression List
#'
#' @description
#' `mutate_from_expr_list` takes a data frame and a list of formula expressions, 
#' and mutates the data frame based on these expressions. It allows for ignoring 
#' errors during the mutation process.
#'
#' @param x A data frame to be mutated.
#' @param formula_expr A named list of formula expressions used for mutation.
#' @param ignore_errors Logical flag indicating whether to ignore errors during mutation.
#'
#' @return A mutated data frame with added or modified columns based on `formula_expr`.
#'
#' @details
#' The function performs various checks and transformations on the formula expressions,
#' ensuring that the specified transformations are valid and can be applied to the data frame.
#' It supports advanced features like handling special characters in column names and intelligent
#' parsing of formulas.
#'
#' @importFrom purrr map2_dfc
#' @importFrom tibble add_column
#' @importFrom tidyselect last_col
#' @importFrom dplyr mutate
#' @importFrom stringr str_subset
#' 
#' @noRd
#' 
mutate_from_expr_list = function(x, formula_expr, ignore_errors = TRUE){

  if(formula_expr |> names() |> is.null())
    names(formula_expr) = formula_expr

  
  
  # Check if all elements of contrasts are in the parameter
  parameter_names = x |> colnames()

  # Creating a named vector where the names are the strings to be replaced
  # and the values are empty strings
  
  # Using str_replace_all to replace each instance of the strings in A with an empty string in B
  contrasts_elements <- 
    formula_expr |> 
    
    # Remove fractions
    str_remove_all_ignoring_if_inside_backquotes("[0-9]+/[0-9]+ ?\\*") |>  
    
    # Remove decimals
    str_remove_all_ignoring_if_inside_backquotes("[-+]?[0-9]+\\.[0-9]+ ?\\*") |> 
    
    str_split_ignoring_if_inside_backquotes("\\+|-|\\*") |> 
    unlist() |> 
    str_remove_all_ignoring_if_inside_backquotes("[\\(\\) ]") 
  
  
  # Check is backquoted are not used
  require_back_quotes = !contrasts_elements |>  str_remove_all("`") |> contains_only_valid_chars_for_column() 
  has_left_back_quotes = contrasts_elements |>  str_detect("^`") 
  has_right_back_quotes = contrasts_elements |>  str_detect("`$") 
  if_true_not_good = require_back_quotes & !(has_left_back_quotes & has_right_back_quotes)
  
  if(any(if_true_not_good))
    warning(sprintf("sccomp says: for columns which have special characters e.g. %s, you need to use surrounding backquotes ``.", paste(contrasts_elements[!if_true_not_good], sep=", ")))
  
  # Check if columns exist
  contrasts_not_in_the_model = 
    contrasts_elements |> 
    str_remove_all("`") |> 
    setdiff(parameter_names)
  
  contrasts_not_in_the_model = contrasts_not_in_the_model[contrasts_not_in_the_model!=""]
  
  if(length(contrasts_not_in_the_model) > 0 & !ignore_errors)
    warning(sprintf("sccomp says: These components of your contrasts are not present in the model as parameters: %s. Factors including special characters, e.g. \"(Intercept)\" require backquotes e.g. \"`(Intercept)`\" ", paste(contrasts_not_in_the_model, sep = ", ")))
  
  # Calculate
  if(ignore_errors) my_mutate = mutate_ignore_error
  else my_mutate = mutate
  
  map2_dfc(
    formula_expr,
    names(formula_expr),
    ~  x |>
      my_mutate(!!.y := eval(rlang::parse_expr(.x))) |>
        # mutate(!!column_name := eval(rlang::parse_expr(.x))) |>
        select(any_of(.y))
  ) |>

  	# I could drop this to just result contrasts
    add_column(x |> select(-any_of(names(formula_expr))), .before = 1)

}

mutate_ignore_error = function(x, ...){
  tryCatch(
    {  x |> mutate(...) },
    error=function(cond) {  x  }
  )
}

simulate_multinomial_logit_linear = function(model_input, sd = 0.51){

  mu = model_input$X %*% model_input$beta

  proportions =
    rnorm(length(mu), mu, sd) %>%
    matrix(nrow = nrow(model_input$X)) %>%
    boot::inv.logit()
  apply(1, function(x) x/sum(x)) %>%
    t()

  rownames(proportions) = rownames(model_input$X)
  colnames(proportions) = colnames(model_input$beta )
}

compress_zero_one = function(y){
  # https://stats.stackexchange.com/questions/48028/beta-regression-of-proportion-data-including-1-and-0

  n = length(y)
  (y * (n-1) + 0.5) / n
}

# this can be helpful if we want to draw PCA with uncertainty
get_abundance_contrast_draws = function(.data, contrasts){

  # Define the variables as NULL to avoid CRAN NOTES
  X <- NULL
  .value <- NULL
  N_random_intercepts <- NULL
  X_random_intercept <- NULL
  .variable <- NULL
  y <- NULL
  M <- NULL
  khat <- NULL
  parameter <- NULL
  n_eff <- NULL
  R_k_hat <- NULL
  
  .cell_group = .data |>  attr(".cell_group")

  # Beta
  beta_factor_of_interest = .data |> attr("model_input") %$% X |> colnames()
  beta =
    .data |>
    attr("fit") %>%
    draws_to_tibble_x_y("beta", "C", "M") |>
    pivot_wider(names_from = C, values_from = .value) %>%
    setNames(colnames(.)[1:5] |> c(beta_factor_of_interest))

  # Random intercept
  is_random_intercept =
    .data |>
    attr("model_input") %$%
    N_random_intercepts |>
    equals(0) |>
    not()

  if(is_random_intercept){
    beta_random_intercept_factor_of_interest = .data |> attr("model_input") %$% X_random_intercept |> colnames()
    beta_random_intercept =
      .data |>
      attr("fit") %>%
      draws_to_tibble_x_y("beta_random_intercept", "C", "M") |>
      pivot_wider(names_from = C, values_from = .value) %>%
      setNames(colnames(.)[1:5] |> c(beta_random_intercept_factor_of_interest))
  } else {
    beta_random_intercept_factor_of_interest = ""
  }


  # Abundance
  draws = select(beta, -.variable)
  
  # Random intercept
  if(is_random_intercept) 
    draws = draws |> 
    left_join(select(beta_random_intercept, -.variable),
              by = c("M", ".chain", ".iteration", ".draw")
    )
  
  # If I have constrasts calculate
  if(!is.null(contrasts))
    draws = 
      draws |> 
      mutate_from_expr_list(contrasts, ignore_errors = FALSE) |>
      select(- any_of(c(beta_factor_of_interest, beta_random_intercept_factor_of_interest) |> setdiff(contrasts)) ) 
  
  # Add cell name
  draws = draws |> 
    left_join(
      .data |>
        attr("model_input") %$%
        y %>%
        colnames() |>
        enframe(name = "M", value  = quo_name(.cell_group)),
      by = "M"
    ) %>%
    select(!!.cell_group, everything())


  # If no contrasts of interest just return an empty data frame
  if(ncol(draws)==5) return(draws |> distinct(M, !!.cell_group))

  # Get convergence
  convergence_df =
    .data |>
      attr("fit") |>
      summary_to_tibble("beta", "C", "M") |>

      # Add cell name
      left_join(
        .data |>
          attr("model_input") %$%
          y %>%
          colnames() |>
          enframe(name = "M", value  = quo_name(.cell_group)),
        by = "M"
      ) |>

      # factor names
      left_join(
        beta_factor_of_interest |>
          enframe(name = "C", value = "parameter"),
        by = "C"
      )

  if ("Rhat" %in% colnames(convergence_df)) {
    convergence_df <- rename(convergence_df, R_k_hat = Rhat)
  } else if ("khat" %in% colnames(convergence_df)) {
    convergence_df <- rename(convergence_df, R_k_hat = khat)
  }
  

    convergence_df =
      convergence_df |> 
      select(!!.cell_group, parameter, any_of(c("n_eff", "R_k_hat"))) |>
      suppressWarnings()

  draws |>
    pivot_longer(-c(1:5), names_to = "parameter", values_to = ".value") |>

    # Attach convergence if I have no contrasts
    left_join(convergence_df, by = c(quo_name(.cell_group), "parameter")) |>

    # Reorder because pivot long is bad
    mutate(parameter = parameter |> fct_relevel(colnames(draws)[-c(1:5)])) |>
    arrange(parameter)

}

#' @importFrom forcats fct_relevel
#' @noRd
get_variability_contrast_draws = function(.data, contrasts){

  # Define the variables as NULL to avoid CRAN NOTES
  XA <- NULL
  .value <- NULL
  y <- NULL
  M <- NULL
  khat <- NULL
  parameter <- NULL
  n_eff <- NULL
  R_k_hat <- NULL
  
  .cell_group = .data |>  attr(".cell_group")

  variability_factor_of_interest = .data |> attr("model_input") %$% XA |> colnames()

  draws =

  .data |>
    attr("fit") %>%
    draws_to_tibble_x_y("alpha_normalised", "C", "M") |>

    # We want variability, not concentration
    mutate(.value = -.value) |>

    pivot_wider(names_from = C, values_from = .value) %>%
    setNames(colnames(.)[1:5] |> c(variability_factor_of_interest)) |>

    select( -.variable) 
  
  # If I have constrasts calculate
  if (!is.null(contrasts)) 
    draws <- mutate_from_expr_list(draws, contrasts, ignore_errors = TRUE)
    
  draws =  draws |>

    # Add cell name
    left_join(
      .data |> attr("model_input") %$%
        y %>%
        colnames() |>
        enframe(name = "M", value  = quo_name(.cell_group)),
      by = "M"
    ) %>%
    select(!!.cell_group, everything())

  # If no contrasts of interest just return an empty data frame
  if(ncol(draws)==5) return(draws |> distinct(M, !!.cell_group))

  # Get convergence
  convergence_df =
    .data |>
    attr("fit") |>
    summary_to_tibble("alpha_normalised", "C", "M") |>

    # Add cell name
    left_join(
      .data |>
        attr("model_input") %$%
        y %>%
        colnames() |>
        enframe(name = "M", value  = quo_name(.cell_group)),
      by = "M"
    ) |>

    # factor names
    left_join(
      variability_factor_of_interest |>
        enframe(name = "C", value = "parameter"),
      by = "C"
    )


  if ("Rhat" %in% colnames(convergence_df)) {
    convergence_df <- rename(convergence_df, R_k_hat = Rhat)
  } else if ("khat" %in% colnames(convergence_df)) {
    convergence_df <- rename(convergence_df, R_k_hat = khat)
  }
  

    convergence_df =
    convergence_df |> 
      select(!!.cell_group, parameter, any_of(c("n_eff", "R_k_hat"))) |>
      suppressWarnings()


  draws |>
    pivot_longer(-c(1:5), names_to = "parameter", values_to = ".value") |>

    # Attach convergence if I have no contrasts
    left_join(convergence_df, by = c(quo_name(.cell_group), "parameter")) |>

    # Reorder because pivot long is bad
    mutate(parameter = parameter |> fct_relevel(colnames(draws)[-c(1:5)])) |>
    arrange(parameter)

}

#' @importFrom tibble deframe
#'
#' @noRd
replicate_data = function(.data,
          formula_composition = NULL,
          formula_variability = NULL,
          new_data = NULL,
          number_of_draws = 1,
          mcmc_seed = sample(1e5, 1)){


  # Select model based on noise model
  noise_model = attr(.data, "noise_model")

  model_input = attr(.data, "model_input")
  .sample = attr(.data, ".sample")
  .cell_group = attr(.data, ".cell_group")

  # Composition
  if(is.null(formula_composition)) formula_composition =  .data |> attr("formula_composition")

  # New data
  if(new_data |> is.null())
    new_data =
    .data |>
    select(count_data) |>
    unnest(count_data) |>
    distinct()

  # If seurat
  else if(new_data |> is("Seurat")) new_data = new_data[[]]

  # Just subset
  new_data = new_data |> .subset(!!.sample)


  # Check if the input new data is not suitable
  if(!parse_formula(formula_composition) %in% colnames(new_data) |> all())
    stop("sccomp says: your `new_data` might be malformed. It might have the covariate columns with multiple values for some element of the \"%s\" column. As a generic example, a sample identifier (\"Sample_123\") might be associated with multiple treatment values, or age values.")


  # Match factors with old data
  nrow_new_data = nrow(new_data)
  new_exposure = new_data |>
    nest(data = -!!.sample) |>
    mutate(exposure = map_dbl(
      data,
      ~{
        if ("count" %in% colnames(.x))  sum(.x$count)
        else 5000
      })) |>
    select(!!.sample, exposure) |>
    deframe() |>
    as.array()

  # Update data, merge with old data because
  # I need the same ordering of the design matrix
  new_data =

    # Old data
    .data |>
    select(count_data) |>
    unnest(count_data) |>
    select(-count) |>
    select(new_data |> as_tibble() |> colnames() |>  any_of()) |>
    distinct() |>

    # Change sample names to make unique
    mutate(dummy = "OLD") |>
    tidyr::unite(!!.sample, c(!!.sample, dummy), sep="___") |>

    # New data
    bind_rows(
      new_data |> as_tibble()
    )

  new_X =
    new_data |>
    get_design_matrix(

      # Drop random intercept
      formula_composition |>
        as.character() |>
        str_remove_all("\\+ ?\\(.+\\|.+\\)") |>
        paste(collapse="") |>
        as.formula(),
      !!.sample
    ) |>
    tail(nrow_new_data) %>%

    # Remove columns that are not in the original design matrix
    .[,colnames(.) %in% colnames(model_input$X), drop=FALSE]

  X_which =
    colnames(new_X) |>
    match(
      model_input$X %>%
        colnames()
    ) |>
    na.omit() |>
    as.array()

  # Variability
  if(is.null(formula_variability)) formula_variability =  .data |> attr("formula_variability")

  new_Xa =
    new_data |>
    get_design_matrix(

      # Drop random intercept
      formula_variability |>
        as.character() |>
        str_remove_all("\\+ ?\\(.+\\|.+\\)") |>
        paste(collapse="") |>
        as.formula(),
      !!.sample
    ) |>
    tail(nrow_new_data) %>%

    # Remove columns that are not in the original design matrix
    .[,colnames(.) %in% colnames(model_input$Xa), drop=FALSE]

  XA_which =
    colnames(new_Xa) |>
    match(
      model_input %$%
        Xa %>%
        colnames()
    ) |>
    na.omit() |>
    as.array()

  # If I want to replicate data with intercept and I don't have intercept in my fit
  create_intercept =
    model_input %$% intercept_in_design |> not() &
    "(Intercept)" %in% colnames(new_X)
  if(create_intercept) warning("sccomp says: your estimated model is intercept free, while your desired replicated data do have an intercept term. The intercept estimate will be calculated averaging your first factor in your formula ~ 0 + <factor>. If you don't know the meaning of this warning, this is likely undesired, and please reconsider your formula for replicate_data()")

  # Random intercept
  random_intercept_elements = parse_formula_random_intercept(formula_composition)
  if(random_intercept_elements |> nrow() |> equals(0)) {
    X_random_intercept_which = array()[0]
    new_X_random_intercept = matrix(rep(0, nrow_new_data))[,0, drop=FALSE]


  }
  else {

    random_intercept_grouping =
      new_data %>%

        get_random_intercept_design2(
        !!.sample,
        formula_composition
      )

    new_X_random_intercept =
      random_intercept_grouping |>
      mutate(design_matrix = map(
        design,
        ~ ..1 |>
          select(!!.sample, group___label, value) |>
          pivot_wider(names_from = group___label, values_from = value) |>
          mutate(across(everything(), ~ .x |> replace_na(0)))
      )) |>

      # Merge
      pull(design_matrix) |>
      bind_cols() |>
      as_matrix(rownames = quo_name(.sample))  |>

      tail(nrow_new_data)

    # Check if I have column in the new design that are not in the old one
    missing_columns = new_X_random_intercept |> colnames() |> setdiff(colnames(model_input$X_random_intercept))
    if(missing_columns |> length() > 0)
    	stop(glue("sccomp says: the columns in the design matrix {paste(missing_columns, collapse= ' ,')} are missing from the design matrix of the estimate-input object. Please make sure your new model is a sub-model of your estimated one."))

    # I HAVE TO KEEP GROUP NAME IN COLUMN NAME
    X_random_intercept_which =
      colnames(new_X_random_intercept) |>
      match(
        model_input %$%
          X_random_intercept %>%
          colnames()
      ) |>
      as.array()
  }

  # New X
  model_input$X_original = model_input$X
  model_input$X = new_X
  model_input$Xa = new_Xa
  model_input$N_original = model_input$N
  model_input$N = nrow_new_data
  model_input$exposure = new_exposure

  model_input$X_random_intercept = new_X_random_intercept
  model_input$N_grouping_new = ncol(new_X_random_intercept)


  
  number_of_draws_in_the_fit = attr(.data, "fit") |>  get_output_samples()
  
  # To avoid error in case of a NULL posterior sample
  number_of_draws = min(number_of_draws, number_of_draws_in_the_fit)
  
  # Load model
  mod_rng = load_model("glm_multi_beta_binomial_generate_data")
  
  
  # Generate quantities
  mod_rng$generate_quantities(
    attr(.data, "fit")$draws(format = "matrix")[
      sample(seq_len(number_of_draws_in_the_fit), size=number_of_draws),, drop=FALSE
    ],
    data = model_input |> c(list(

      # Add subset of coefficients
      length_X_which = length(X_which),
      length_XA_which = length(XA_which),
      X_which = X_which,
      XA_which = XA_which,

      # Random intercept
      X_random_intercept_which = X_random_intercept_which,
      length_X_random_intercept_which = length(X_random_intercept_which),

      # Should I create intercept for generate quantities
      create_intercept = create_intercept

    )),
    seed = mcmc_seed, 
    threads_per_chain = 1
  )



  
  
  
  
}

get_model_from_data = function(file_compiled_model, model_code){
  if(file.exists(file_compiled_model))
    readRDS(file_compiled_model)
  else {
    model_generate = stan_model(model_code = model_code)
    model_generate  %>% saveRDS(file_compiled_model)
    model_generate

  }
}

add_formula_columns = function(.data, .original_data, .sample,  formula_composition){

  .sample = enquo(.sample)

  formula_elements = parse_formula(formula_composition)

  # If no formula return the input
  if(length(formula_elements) == 0) return(.data)

  # Get random intercept
  .grouping_for_random_intercept = parse_formula_random_intercept(formula_composition) |> pull(grouping) |> unique()

  data_frame_formula =
    .original_data %>%
    as_tibble() |>
    select( !!.sample, formula_elements, any_of(.grouping_for_random_intercept) ) %>%
    distinct()

  .data |>
    left_join(data_frame_formula, by = quo_name(.sample) )

}

#' chatGPT - Remove Specified Regex Pattern from Each String in a Vector
#'
#' This function takes a vector of strings and a regular expression pattern.
#' It removes occurrences of the pattern from each string, except where the pattern
#' is found inside backticks. The function returns a vector of cleaned strings.
#'
#' @param text_vector A character vector with the strings to be processed.
#' @param regex A character string containing a regular expression pattern to be removed
#' from the text.
#'
#' @return A character vector with the regex pattern removed from each string.
#' Occurrences of the pattern inside backticks are not removed.
#'
#' @examples
#' texts <- c("A string with (some) parentheses and `a (parenthesis) inside` backticks",
#'            "Another string with (extra) parentheses")
#' cleaned_texts <- str_remove_all_ignoring_if_inside_backquotes(texts, "\\(")
#' print(cleaned_texts)
#' 
#' @noRd
str_remove_all_ignoring_if_inside_backquotes <- function(text_vector, regex) {
  # Nested function to handle regex removal for a single string
  remove_regex_chars <- function(text, regex) {
    inside_backticks <- FALSE
    result <- ""
    skip <- 0
    
    chars <- strsplit(text, "")[[1]]
    for (i in seq_along(chars)) {
      if (skip > 0) {
        skip <- skip - 1
        next
      }
      
      char <- chars[i]
      if (char == "`") {
        inside_backticks <- !inside_backticks
        result <- paste0(result, char)
      } else if (!inside_backticks) {
        # Check the remaining text against the regex
        remaining_text <- paste(chars[i:length(chars)], collapse = "")
        match <- regexpr(regex, remaining_text)
        
        if (attr(match, "match.length") > 0 && match[1] == 1) {
          # Skip the length of the matched text
          skip <- attr(match, "match.length") - 1
          next
        } else {
          result <- paste0(result, char)
        }
      } else {
        result <- paste0(result, char)
      }
    }
    
    return(result)
  }
  
  # Apply the function to each element in the vector
  sapply(text_vector, remove_regex_chars, regex)
}


#' chatGPT - Split Each String in a Vector by a Specified Regex Pattern
#'
#' This function takes a vector of strings and a regular expression pattern. It splits
#' each string based on the pattern, except where the pattern is found inside backticks.
#' The function returns a list, with each element being a vector of the split segments
#' of the corresponding input string.
#'
#' @param text_vector A character vector with the strings to be processed.
#' @param regex A character string containing a regular expression pattern used for splitting
#' the text.
#'
#' @return A list of character vectors. Each list element corresponds to an input string
#' from `text_vector`, split according to `regex`, excluding occurrences inside backticks.
#'
#' @examples
#' texts <- c("A string with, some, commas, and `a, comma, inside` backticks",
#'            "Another string, with, commas")
#' split_texts <- split_regex_chars_from_vector(texts, ",")
#' print(split_texts)
#' 
#' @noRd
str_split_ignoring_if_inside_backquotes <- function(text_vector, regex) {
  # Nested function to handle regex split for a single string
  split_regex_chars <- function(text, regex) {
    inside_backticks <- FALSE
    result <- c()
    current_segment <- ""
    
    chars <- strsplit(text, "")[[1]]
    for (i in seq_along(chars)) {
      char <- chars[i]
      if (char == "`") {
        inside_backticks <- !inside_backticks
        current_segment <- paste0(current_segment, char)
      } else if (!inside_backticks) {
        # Check the remaining text against the regex
        remaining_text <- paste(chars[i:length(chars)], collapse = "")
        match <- regexpr(regex, remaining_text)
        
        if (attr(match, "match.length") > 0 && match[1] == 1) {
          # Add current segment to result and start a new segment
          result <- c(result, current_segment)
          current_segment <- ""
          # Skip the length of the matched text
          skip <- attr(match, "match.length") - 1
          i <- i + skip
        } else {
          current_segment <- paste0(current_segment, char)
        }
      } else {
        current_segment <- paste0(current_segment, char)
      }
    }
    
    # Add the last segment to the result
    result <- c(result, current_segment)
    return(result)
  }
  
  # Apply the function to each element in the vector
  lapply(text_vector, split_regex_chars, regex)
}


#' chatGPT - Check for Valid Column Names in Tidyverse Context
#'
#' This function checks if each given column name in a vector contains only valid characters 
#' (letters, numbers, periods, and underscores) and does not start with a digit 
#' or an underscore, which are the conditions for a valid column name in `tidyverse`.
#'
#' @param column_names A character vector representing the column names to be checked.
#'
#' @return A logical vector: `TRUE` for each column name that contains only valid characters 
#' and does not start with a digit or an underscore; `FALSE` otherwise.
#'
#' @examples
#' contains_only_valid_chars_for_column(c("valid_column", "invalid column", "valid123", 
#' "123startWithNumber", "_startWithUnderscore"))
#'
#' @noRd
contains_only_valid_chars_for_column <- function(column_names) {
  # Function to check a single column name
  check_validity <- function(column_name) {
    # Regex pattern for valid characters (letters, numbers, periods, underscores)
    valid_char_pattern <- "[A-Za-z0-9._]"
    
    # Check if all characters in the string match the valid pattern
    all_chars_valid <- stringr::str_detect(column_name, paste0("^", valid_char_pattern, "+$"))
    
    # Check for leading digits or underscores
    starts_with_digit_or_underscore <- stringr::str_detect(column_name, "^[0-9_]")
    
    return(all_chars_valid && !starts_with_digit_or_underscore)
  }
  
  # Apply the check to each element of the vector
  sapply(column_names, check_validity)
}


#' chatGPT - Intelligently Remove Surrounding Brackets from Each String in a Vector
#'
#' This function processes each string in a vector and removes surrounding brackets if the content
#' within the brackets includes any of '+', '-', or '*', and if the brackets are not 
#' within backticks. This is particularly useful for handling formula-like strings.
#'
#' @param text A character vector with strings from which the brackets will be removed based on
#' specific conditions.
#'
#' @return A character vector with the specified brackets removed from each string.
#'
#' @examples
#' str_remove_brackets_from_formula_intelligently(c("This is a test (with + brackets)", "`a (test) inside` backticks", "(another test)"))
#'
#' @noRd
str_remove_brackets_from_formula_intelligently <- function(text) {
  # Function to remove brackets from a single string
  remove_brackets_single <- function(s) {
    inside_backticks <- FALSE
    bracket_depth <- 0
    valid_bracket_content <- FALSE
    result <- ""
    bracket_content <- ""
    
    chars <- strsplit(s, "")[[1]]
    
    for (i in seq_along(chars)) {
      char <- chars[i]
      
      if (char == "`") {
        inside_backticks <- !inside_backticks
      }
      
      if (!inside_backticks) {
        if (char == "(") {
          bracket_depth <- bracket_depth + 1
          if (bracket_depth > 1) {
            bracket_content <- paste0(bracket_content, char)
          }
          next
        } else if (char == ")") {
          bracket_depth <- bracket_depth - 1
          if (bracket_depth == 0) {
            if (grepl("[\\+\\-\\*]", bracket_content)) {
              result <- paste0(result, bracket_content)
            } else {
              result <- paste0(result, "(", bracket_content, ")")
            }
            bracket_content <- ""
            next
          }
        }
        
        if (bracket_depth >= 1) {
          bracket_content <- paste0(bracket_content, char)
        } else {
          result <- paste0(result, char)
        }
      } else {
        result <- paste0(result, char)
      }
    }
    
    return(result)
  }
  
  # Apply the function to each element in the vector
  sapply(text, remove_brackets_single)
}



# Negation
not = function(is){	!is }

#' Convert array of quosure (e.g. c(col_a, col_b)) into character vector
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom rlang quo_name
#' @importFrom rlang quo_squash
#'
#' @param v A array of quosures (e.g. c(col_a, col_b))
#'
#' @return A character vector
quo_names <- function(v) {

	v = quo_name(quo_squash(v))
	gsub('^c\\(|`|\\)$', '', v) |>
		strsplit(', ') |>
		unlist()
}

#' Add class to abject
#'
#' @keywords internal
#' @noRd
#'
#' @param var A tibble
#' @param name A character name of the attribute
#'
#' @return A tibble with an additional attribute
add_class = function(var, name) {

  if(!name %in% class(var)) class(var) <- append(class(var),name, after = 0)

  var
}


#' Get Output Samples from a Stan Fit Object
#'
#' This function retrieves the number of output samples from a Stan fit object, 
#' supporting different methods (MHC and Variational) based on the available data within the object.
#'
#' @param fit A `stanfit` object, which is the result of fitting a model via Stan.
#' @return The number of output samples used in the Stan model. 
#'         Returns from MHC if available, otherwise from Variational inference.
#' @examples
#' # Assuming 'fit' is a stanfit object obtained from running a Stan model
#' samples_count = get_output_samples(fit)
#'
get_output_samples = function(fit){
  
  # Check if the output_samples field is present in the metadata of the fit object
  # This is generally available when the model is fit using MHC (Markov chain Monte Carlo)
  if(!is.null(fit$metadata()$output_samples)) {
    # Return the output_samples from the metadata
    fit$metadata()$output_samples
  }
  
  # If the output_samples field is not present, check for iter_sampling
  # This occurs typically when the model is fit using Variational inference methods
  else if(!is.null(fit$metadata()$iter_sampling)) {
    # Return the iter_sampling from the metadata
    fit$metadata()$iter_sampling
  }
  else
    fit$metadata()$num_psis_draws
}





#' Load, Compile, and Cache a Stan Model
#'
#' This function attempts to load a precompiled Stan model using the `instantiate` package.
#' If the model is not found, it will locate the Stan model file within the `sccomp` package,
#' compile it using `cmdstanr`, and save the compiled model to the cache directory.
#'
#' @param name A character string representing the name of the Stan model.
#' @param cache_dir A character string representing the path to the cache directory.
#' 
#' @return A compiled Stan model object.
#' 
#' @importFrom instantiate stan_package_model
#' @importFrom cmdstanr cmdstan_model
#' @export
#' @examples
#' \dontrun{
#'   model <- load_model("glm_multi_beta_binomial_", "~/cache")
#' }
load_model <- function(name, cache_dir = sccomp_stan_models_cache_dir) {
  
  
  # tryCatch({
  #   # Attempt to load a precompiled Stan model using the instantiate package
  #   instantiate::stan_package_model(
  #     name = name,
  #     package = "sccomp"
  #   )
  # }, error = function(e) {
    # Try to load the model from cache
  cache_dir |> dir.create(showWarnings = FALSE, recursive = TRUE)
    cache_file <- file.path(cache_dir, paste0(name, ".rds"))
    if (file.exists(cache_file)) {
      message("Loading model from cache...")
      return(readRDS(cache_file))
    }
    
    # If loading the precompiled model fails, find the Stan model file within the package
    message("Precompiled model not found. Compiling the model...")
    stan_model_path <- system.file("stan", paste0(name, ".stan"), package = "sccomp")
    
    # Compile the Stan model using cmdstanr with threading support enabled
    mod <- cmdstan_model(
      stan_model_path, 
      cpp_options = list(stan_threads = TRUE),
      force_recompile = TRUE
    )
    
    # Save the compiled model object to cache
    saveRDS(mod, file = cache_file)
    message("Model compiled and saved to cache successfully.")
    
    return(mod)
  # })

}

#' Check and Install cmdstanr and CmdStan
#'
#' This function checks if the `cmdstanr` package and CmdStan are installed. 
#' If they are not installed, it installs them automatically in non-interactive sessions
#' or asks for permission to install them in interactive sessions.
#'
#' @importFrom instantiate stan_cmdstan_exists
#' @importFrom utils install.packages
#' @importFrom utils menu
#' @return NULL
#' 
#' @examples
#' \dontrun{
#'   check_and_install_cmdstanr()
#' }
check_and_install_cmdstanr <- function() {
  # Check if cmdstanr is installed
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    message("The 'cmdstanr' package is not installed.")
    if (interactive()) {
      install <- menu(c("yes", "no"), title = "Do you want to install 'cmdstanr'?")
      if (install == 1) {
        install.packages(pkgs = "cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
        library(cmdstanr)
      } else {
        stop("cmdstanr is required to proceed.")
      }
    } else {
      message("Installing 'cmdstanr' package...")
      install.packages(pkgs = "cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
      library(cmdstanr)
    }
  }
  
  # Check if CmdStan is installed
  if (!stan_cmdstan_exists()) {
    message("CmdStan is not installed.")
    if (interactive()) {
      install <- menu(c("yes", "no"), title = "Do you want to install CmdStan?")
      if (install == 1) {
        cmdstanr::install_cmdstan()
      } else {
        stop("CmdStan is required to proceed.")
      }
    } else {
      message("Installing CmdStan...")
      cmdstanr::install_cmdstan()
    }
  }
}


