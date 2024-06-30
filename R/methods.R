

#' Main Function for SCCOMP Estimate
#'
#' @description
#' The `sccomp_estimate` function performs linear modeling on a table of cell counts, 
#' which includes a cell-group identifier, sample identifier, integer count, and factors 
#' (continuous or discrete). The user can define a linear model with an input R formula, 
#' where the first factor is the factor of interest. Alternatively, `sccomp` accepts 
#' single-cell data containers (e.g., Seurat, SingleCellExperiment, cell metadata, or 
#' group-size) and derives the count data from cell metadata.
#'
#' @import dplyr
#' @importFrom magrittr %$%
#' @importFrom magrittr divide_by
#' @importFrom magrittr multiply_by
#' @importFrom magrittr equals
#' @importFrom rlang quo_is_null
#' @importFrom SingleCellExperiment colData
#' @importFrom parallel detectCores
#' @importFrom rlang inform
#' @importFrom lifecycle is_present
#' @importFrom lifecycle deprecate_warn
#'
#' @param .data A tibble including cell_group name column, sample name column, 
#'              read counts column (optional depending on the input class), and factor columns.
#' @param formula_composition A formula describing the model for differential abundance.
#' @param formula_variability A formula describing the model for differential variability.
#' @param .sample A column name as symbol for the sample identifier.
#' @param .cell_group A column name as symbol for the cell_group identifier.
#' @param .count A column name as symbol for the cell_group abundance (read count).
#' @param cores Number of cores to use for parallel calculations.
#' @param bimodal_mean_variability_association Boolean for modeling mean-variability as bimodal.
#' @param prior_mean List with prior knowledge about mean distribution, they are the intercept and coefficient.
#' @param prior_overdispersion_mean_association List with prior knowledge about mean/variability association.
#' @param percent_false_positive Real number between 0 and 100 for outlier identification.
#' @param variational_inference Boolean for using variational Bayes for posterior inference. It is faster and convenient. Setting this argument to FALSE runs the full Bayesian (Hamiltonian Monte Carlo) inference, slower but it is the gold standard.
#' @param enable_loo Boolean to enable model comparison using the LOO package.
#' @param exclude_priors Boolean to run a prior-free model.
#' @param use_data Boolean to run the model data-free.
#' @param max_sampling_iterations Integer to limit maximum iterations for large datasets.
#' @param pass_fit Boolean to include the Stan fit as attribute in the output.
#' @param .sample_cell_group_pairs_to_exclude Column name with boolean for sample/cell-group pairs exclusion.
#' @param verbose Boolean to print progression.
#' @param noise_model Character string for the noise model (e.g., 'multi_beta_binomial').
#' @param mcmc_seed Integer for MCMC reproducibility.
#' @param approximate_posterior_inference DEPRECATED please use the `variational_inference` argument.
#' 
#' 
#' @return A nested tibble `tbl`, with the following columns
#' \itemize{
#'   \item cell_group - column including the cell groups being tested
#'   \item parameter - The parameter being estimated, from the design matrix dscribed with the input formula_composition and formula_variability
#'   \item factor - The factor in the formula corresponding to the covariate, if exists (e.g. it does not exist in case og Intercept or contrasts, which usually are combination of parameters)
#'
#'   \item c_lower - lower (2.5%) quantile of the posterior distribution for a composition (c) parameter.
#'   \item c_effect - mean of the posterior distribution for a composition (c) parameter.
#'   \item c_upper - upper (97.5%) quantile of the posterior distribution fo a composition (c)  parameter.
#'   \item c_pH0 - Probability of the null hypothesis (no difference) for  a composition (c). This is not a p-value.
#'   \item c_FDR - False-discovery rate of the null hypothesis (no difference) for  a composition (c).
#'   \item c_n_eff - Effective sample size - the number of independent draws in the sample, the higher the better (mc-stan.org/docs/2_25/cmdstan-guide/stansummary.html).
#'   \item c_R_k_hat - R statistic, a measure of chain equilibrium, should be within 0.05 of 1.0 (mc-stan.org/docs/2_25/cmdstan-guide/stansummary.html).
#'
#'   \item v_lower - Lower (2.5%) quantile of the posterior distribution for a variability (v) parameter
#'   \item v_effect - Mean of the posterior distribution for a variability (v) parameter
#'   \item v_upper - Upper (97.5%) quantile of the posterior distribution for a variability (v) parameter
#'   \item v_pH0 - Probability of the null hypothesis (no difference) for a variability (v). This is not a p-value.
#'   \item v_FDR - False-discovery rate of the null hypothesis (no difference), for a variability (v).
#'   \item v_n_eff - Effective sample size for a variability (v) parameter - the number of independent draws in the sample, the higher the better (mc-stan.org/docs/2_25/cmdstan-guide/stansummary.html).
#'   \item v_R_k_hat - R statistic for a variability (v) parameter, a measure of chain equilibrium, should be within 0.05 of 1.0 (mc-stan.org/docs/2_25/cmdstan-guide/stansummary.html).
#'
#'   \item count_data Nested input count data.
#'
#' }
#'
#' @examples
#'
#' data("counts_obj")
#'
#' estimate =
#'   sccomp_estimate(
#'   counts_obj ,
#'    ~ type,
#'    ~1,
#'    sample,
#'    cell_group,
#'    count,
#'     cores = 1
#'   )
#'
#' @export
#'
#'
sccomp_estimate <- function(.data,
                       formula_composition = ~ 1 ,
                       formula_variability = ~ 1,
                       .sample,
                       .cell_group,
                       .count = NULL,

                       # Secondary arguments
                       cores = detectCores(),
                       bimodal_mean_variability_association = FALSE,
                       percent_false_positive = 5,
                       inference_method = "pathfinder",
                       prior_mean = list(intercept = c(0,1), coefficients = c(0,1)),
                       prior_overdispersion_mean_association = list(intercept = c(5, 2), slope = c(0,  0.6), standard_deviation = c(10, 20)),
                       .sample_cell_group_pairs_to_exclude = NULL,
                       verbose = TRUE,
                       enable_loo = FALSE,
                       noise_model = "multi_beta_binomial",
                       exclude_priors = FALSE,
                       use_data = TRUE,
                       mcmc_seed = sample(1e5, 1),
                       max_sampling_iterations = 20000,
                       pass_fit = TRUE,
                       
                       # DEPRECATED
                       approximate_posterior_inference = NULL,
                       variational_inference = NULL
                       
                       ) {
  if(inference_method == "variational") 
    rlang::inform(
      message = "sccomp says: From version 1.7.7 the model by default is fit with the variational inference method (inference_method = “variational”; much faster). For a full Bayesian inference (HMC method; the gold standard) use inference_method = \"hmc\".", 
      .frequency = "once", 
      .frequency_id = "variational_message"
    )
  
  rlang::inform(
      message = "sccomp says: From version 1.7.12 the logit fold change threshold for significance has be changed from 0.2 to 0.1.", 
      .frequency = "once", 
      .frequency_id = "new_logit_fold_change_threshold"
  )
  
  # Run the function
  check_and_install_cmdstanr()
  
  UseMethod("sccomp_estimate", .data)
}

#' @export
sccomp_estimate.Seurat = function(.data,
                                  formula_composition = ~ 1 ,
                                  formula_variability = ~ 1,
                                  .sample,
                                  .cell_group,
                                  .count = NULL,
                                  
                                  # Secondary arguments
                                  cores = detectCores(),
                                  bimodal_mean_variability_association = FALSE,
                                  percent_false_positive = 5,
                                  inference_method = "pathfinder",
                                  prior_mean = list(intercept = c(0,1), coefficients = c(0,1)),                        
                                  prior_overdispersion_mean_association = list(intercept = c(5, 2), slope = c(0,  0.6), standard_deviation = c(10, 20)),
                                  .sample_cell_group_pairs_to_exclude = NULL,
                                  verbose = TRUE,
                                  enable_loo = FALSE,
                                  noise_model = "multi_beta_binomial",
                                  exclude_priors = FALSE,
                                  use_data = TRUE,
                                  mcmc_seed = sample(1e5, 1),
                                  max_sampling_iterations = 20000,
                                  pass_fit = TRUE,
                                  
                                  # DEPRECATED
                                  approximate_posterior_inference = NULL,
                                  variational_inference = NULL){

  if(!is.null(.count)) stop("sccomp says: .count argument can be used only for data frame input")

  # DEPRECATION OF approximate_posterior_inference
  if (is_present(approximate_posterior_inference) & !is.null(approximate_posterior_inference)) {
    deprecate_warn("1.7.7", "sccomp::sccomp_estimate(approximate_posterior_inference = )", details = "The argument approximate_posterior_inference is now deprecated please use inference_method By default variational_inference value is inferred from approximate_posterior_inference.")
    
     inference_method = ifelse(approximate_posterior_inference == "all", "variational","hmc")
  }
  
  # DEPRECATION OF variational_inference
  if (is_present(variational_inference) & !is.null(variational_inference)) {
    deprecate_warn("1.7.11", "sccomp::sccomp_estimate(variational_inference = )", details = "The argument variational_inference is now deprecated please use variational_inference. By default inference_method value is inferred from variational_inference")
    
    inference_method = ifelse(variational_inference, "variational","hmc")
  }
  
  # DEPRECATION OF variational_inference
  if (is_present(variational_inference) & !is.null(variational_inference)) {
    deprecate_warn("1.7.11", "sccomp::sccomp_estimate(variational_inference = )", details = "The argument variational_inference is now deprecated please use variational_inference. By default inference_method value is inferred from variational_inference")
    
    inference_method = ifelse(variational_inference, "variational","hmc")
  }
  
  # Prepare column same enquo
  .sample = enquo(.sample)
  .cell_group = enquo(.cell_group)
  .sample_cell_group_pairs_to_exclude = enquo(.sample_cell_group_pairs_to_exclude)

  .data[[]] %>%
    sccomp_estimate(
      formula_composition = formula_composition,
      formula_variability = formula_variability,
      !!.sample,
      !!.cell_group,
      
      # Secondary arguments
      cores = cores,
      bimodal_mean_variability_association = bimodal_mean_variability_association,
      percent_false_positive = percent_false_positive,
      inference_method = inference_method,
      prior_mean = prior_mean, 
      prior_overdispersion_mean_association = prior_overdispersion_mean_association,
      .sample_cell_group_pairs_to_exclude = !!.sample_cell_group_pairs_to_exclude,
      verbose = verbose,
      enable_loo = enable_loo,
      noise_model = noise_model,
      exclude_priors = exclude_priors,
      use_data = use_data,
      mcmc_seed = mcmc_seed,
      max_sampling_iterations = max_sampling_iterations,
      pass_fit = pass_fit
    )


}

#' @export
sccomp_estimate.SingleCellExperiment = function(.data,
                                                formula_composition = ~ 1 ,
                                                formula_variability = ~ 1,
                                                .sample,
                                                .cell_group,
                                                .count = NULL,
                                                
                                                # Secondary arguments
                                                cores = detectCores(),
                                                bimodal_mean_variability_association = FALSE,
                                                percent_false_positive = 5,
                                                inference_method = "pathfinder",
                                                prior_mean = list(intercept = c(0,1), coefficients = c(0,1)),                        
                                                prior_overdispersion_mean_association = list(intercept = c(5, 2), slope = c(0,  0.6), standard_deviation = c(10, 20)),
                                                .sample_cell_group_pairs_to_exclude = NULL,
                                                verbose = TRUE,
                                                enable_loo = FALSE,
                                                noise_model = "multi_beta_binomial",
                                                exclude_priors = FALSE,
                                                use_data = TRUE,
                                                mcmc_seed = sample(1e5, 1),
                                                max_sampling_iterations = 20000,
                                                pass_fit = TRUE,
                                                
                                                # DEPRECATED
                                                approximate_posterior_inference = NULL,
                                                variational_inference = NULL) {


  if(!is.null(.count)) stop("sccomp says: .count argument can be used only for data frame input")

  
  # DEPRECATION OF approximate_posterior_inference
  if (is_present(approximate_posterior_inference) & !is.null(approximate_posterior_inference)) {
    deprecate_warn("1.7.7", "sccomp::sccomp_estimate(approximate_posterior_inference = )", details = "The argument approximate_posterior_inference is now deprecated please use variational_inference. By default variational_inference value is inferred from approximate_posterior_inference.")
     inference_method = ifelse(approximate_posterior_inference == "all", "variational","hmc")
  }

  # DEPRECATION OF variational_inference
  if (is_present(variational_inference) & !is.null(variational_inference)) {
    deprecate_warn("1.7.11", "sccomp::sccomp_estimate(variational_inference = )", details = "The argument variational_inference is now deprecated please use variational_inference. By default inference_method value is inferred from variational_inference")
    
    inference_method = ifelse(variational_inference, "variational","hmc")
  }
  
  # Prepare column same enquo
  .sample = enquo(.sample)
  .cell_group = enquo(.cell_group)
  .sample_cell_group_pairs_to_exclude = enquo(.sample_cell_group_pairs_to_exclude)


  .data %>%
    colData() %>%
    sccomp_estimate(
      formula_composition = formula_composition,
      formula_variability = formula_variability,
      !!.sample,
      !!.cell_group,
      
      # Secondary arguments
      cores = cores,
      bimodal_mean_variability_association = bimodal_mean_variability_association,
      percent_false_positive = percent_false_positive,
      inference_method = inference_method,
      prior_mean = prior_mean, 
      prior_overdispersion_mean_association = prior_overdispersion_mean_association,
      .sample_cell_group_pairs_to_exclude = !!.sample_cell_group_pairs_to_exclude,
      verbose = verbose,
      enable_loo = enable_loo,
      noise_model = noise_model,
      exclude_priors = exclude_priors,
      use_data = use_data,
      mcmc_seed = mcmc_seed,
      max_sampling_iterations = max_sampling_iterations,
      pass_fit = pass_fit
    )


}

#' @export
sccomp_estimate.DFrame = function(.data,
                                  formula_composition = ~ 1 ,
                                  formula_variability = ~ 1,
                                  .sample,
                                  .cell_group,
                                  .count = NULL,
                                  
                                  # Secondary arguments
                                  cores = detectCores(),
                                  bimodal_mean_variability_association = FALSE,
                                  percent_false_positive = 5,
                                  inference_method = "pathfinder",
                                  prior_mean = list(intercept = c(0,1), coefficients = c(0,1)),                        
                                  prior_overdispersion_mean_association = list(intercept = c(5, 2), slope = c(0,  0.6), standard_deviation = c(10, 20)),
                                  .sample_cell_group_pairs_to_exclude = NULL,
                                  verbose = TRUE,
                                  enable_loo = FALSE,
                                  noise_model = "multi_beta_binomial",
                                  exclude_priors = FALSE,
                                  use_data = TRUE,
                                  mcmc_seed = sample(1e5, 1),
                                  max_sampling_iterations = 20000,
                                  pass_fit = TRUE,
                                  
                                  # DEPRECATED
                                  approximate_posterior_inference = NULL,
                                  variational_inference = NULL) {


  if(!is.null(.count)) stop("sccomp says: .count argument can be used only for data frame input")

  
  # DEPRECATION OF approximate_posterior_inference
  if (is_present(approximate_posterior_inference) & !is.null(approximate_posterior_inference)) {
    deprecate_warn("1.7.7", "sccomp::sccomp_estimate(approximate_posterior_inference = )", details = "The argument approximate_posterior_inference is now deprecated please use variational_inference. By default variational_inference value is inferred from approximate_posterior_inference.")
     inference_method = ifelse(approximate_posterior_inference == "all", "variational","hmc")
  }
  
  # DEPRECATION OF variational_inference
  if (is_present(variational_inference) & !is.null(variational_inference)) {
    deprecate_warn("1.7.11", "sccomp::sccomp_estimate(variational_inference = )", details = "The argument variational_inference is now deprecated please use variational_inference. By default inference_method value is inferred from variational_inference")
    
    inference_method = ifelse(variational_inference, "variational","hmc")
  }
  
  # Prepare column same enquo
  .sample = enquo(.sample)
  .cell_group = enquo(.cell_group)
  .count = enquo(.count)
  .sample_cell_group_pairs_to_exclude = enquo(.sample_cell_group_pairs_to_exclude)

  .data %>%
    as.data.frame %>%
    sccomp_estimate(
      formula_composition = formula_composition,
      formula_variability = formula_variability,
      !!.sample,
      !!.cell_group,
      
      # Secondary arguments
      cores = cores,
      bimodal_mean_variability_association = bimodal_mean_variability_association,
      percent_false_positive = percent_false_positive,
      inference_method = inference_method,
      prior_mean = prior_mean, 
      prior_overdispersion_mean_association = prior_overdispersion_mean_association,
      .sample_cell_group_pairs_to_exclude = !!.sample_cell_group_pairs_to_exclude,
      verbose = verbose,
      enable_loo = enable_loo,
      noise_model = noise_model,
      exclude_priors = exclude_priors,
      use_data = use_data,
      mcmc_seed = mcmc_seed,
      max_sampling_iterations = max_sampling_iterations,
      pass_fit = pass_fit
    )
}

#' @importFrom purrr when
#' @export
sccomp_estimate.data.frame = function(.data,
                                      formula_composition = ~ 1 ,
                                      formula_variability = ~ 1,
                                      .sample,
                                      .cell_group,
                                      .count = NULL,
                                      
                                      # Secondary arguments
                                      cores = detectCores(),
                                      bimodal_mean_variability_association = FALSE,
                                      percent_false_positive = 5,
                                      inference_method = "pathfinder",
                                      prior_mean = list(intercept = c(0,1), coefficients = c(0,1)),                        
                                      prior_overdispersion_mean_association = list(intercept = c(5, 2), slope = c(0,  0.6), standard_deviation = c(10, 20)),
                                      .sample_cell_group_pairs_to_exclude = NULL,
                                      verbose = TRUE,
                                      enable_loo = FALSE,
                                      noise_model = "multi_beta_binomial",
                                      exclude_priors = FALSE,
                                      use_data = TRUE,
                                      mcmc_seed = sample(1e5, 1),
                                      max_sampling_iterations = 20000,
                                      pass_fit = TRUE,
                                      
                                      # DEPRECATED
                                      approximate_posterior_inference = NULL,
                                      variational_inference = NULL) {

  
  # DEPRECATION OF approximate_posterior_inference
  if (is_present(approximate_posterior_inference) & !is.null(approximate_posterior_inference)) {
    deprecate_warn("1.7.7", "sccomp::sccomp_estimate(approximate_posterior_inference = )", details = "The argument approximate_posterior_inference is now deprecated please use variational_inference. By default variational_inference value is inferred from approximate_posterior_inference.")
     inference_method = ifelse(approximate_posterior_inference == "all", "variational","hmc")
  }
  
  # DEPRECATION OF variational_inference
  if (is_present(variational_inference) & !is.null(variational_inference)) {
    deprecate_warn("1.7.11", "sccomp::sccomp_estimate(variational_inference = )", details = "The argument variational_inference is now deprecated please use variational_inference. By default inference_method value is inferred from variational_inference")
    
    inference_method = ifelse(variational_inference, "variational","hmc")
  }
  

  # Prepare column same enquo
  .sample = enquo(.sample)
  .cell_group = enquo(.cell_group)
  .count = enquo(.count)
  .sample_cell_group_pairs_to_exclude = enquo(.sample_cell_group_pairs_to_exclude)

  if( quo_is_null(.count))
    res = sccomp_glm_data_frame_raw(
      .data,
      formula_composition = formula_composition,
      formula_variability = formula_variability,
      !!.sample,
      !!.cell_group,
      
      # Secondary arguments
      cores = cores,
      bimodal_mean_variability_association = bimodal_mean_variability_association,
      percent_false_positive = percent_false_positive,
      inference_method = inference_method,
      prior_mean = prior_mean, 
      prior_overdispersion_mean_association = prior_overdispersion_mean_association,
      .sample_cell_group_pairs_to_exclude = !!.sample_cell_group_pairs_to_exclude,
      verbose = verbose,
      enable_loo = enable_loo,
      exclude_priors = exclude_priors,
      use_data = use_data,
      mcmc_seed = mcmc_seed,
      max_sampling_iterations = max_sampling_iterations,
      pass_fit = pass_fit
    )
  
  else 
    res = sccomp_glm_data_frame_counts(
      .data,
      formula_composition = formula_composition,
      formula_variability = formula_variability,
      !!.sample,
      !!.cell_group,
      !!.count,
      
      # Secondary arguments
      cores = cores,
      bimodal_mean_variability_association = bimodal_mean_variability_association,
      percent_false_positive = percent_false_positive,
      inference_method = inference_method,
      prior_mean = prior_mean, 
      prior_overdispersion_mean_association = prior_overdispersion_mean_association,
      .sample_cell_group_pairs_to_exclude = !!.sample_cell_group_pairs_to_exclude,
      verbose = verbose,
      enable_loo = enable_loo,
      exclude_priors = exclude_priors,
      use_data = use_data,
      mcmc_seed = mcmc_seed,
      max_sampling_iterations = max_sampling_iterations,
      pass_fit = pass_fit
    )
  
  res  |> 

    # Track input parameters
    add_attr(noise_model, "noise_model")  |> 
    add_attr(.sample, ".sample")  |> 
    add_attr(.cell_group, ".cell_group")
}


#' sccomp_remove_outliers main
#'
#' @description The function for linear modelling takes as input a table of cell counts with three columns containing a cell-group identifier, sample identifier, integer count and the factors (continuous or discrete). The user can define a linear model with an input R formula, where the first factor is the factor of interest. Alternatively, sccomp accepts single-cell data containers (Seurat, SingleCellExperiment44, cell metadata or group-size). In this case, sccomp derives the count data from cell metadata.
#'
#' @import dplyr
#' @importFrom magrittr %$%
#' @importFrom magrittr divide_by
#' @importFrom magrittr multiply_by
#' @importFrom magrittr equals
#' @importFrom rlang quo_is_null
#' @importFrom SingleCellExperiment colData
#' @importFrom parallel detectCores
#' @importFrom rlang quo_is_symbolic
#' @importFrom rlang inform
#' @importFrom tidyr unnest
#' @importFrom tidyr nest
#'
#' @param .estimate A tibble including a cell_group name column | sample name column | read counts column (optional depending on the input class) | factor columns.
#' @param percent_false_positive A real between 0 and 100 non included. This used to identify outliers with a specific false positive rate.
#' @param variational_inference Boolean for using variational Bayes for posterior inference. It is faster and convenient. Setting this argument to FALSE runs the full Bayesian (Hamiltonian Monte Carlo) inference, slower but it is the gold standard.
#' @param cores An integer. How many cored to be used with parallel calculations.
#' @param verbose A boolean. Prints progression.
#' @param mcmc_seed An integer. Used for Markov-chain Monte Carlo reproducibility. By default a random number is sampled from 1 to 999999. This itself can be controlled by set.seed()
#' @param max_sampling_iterations An integer. This limit the maximum number of iterations in case a large dataset is used, for limiting the computation time.
#' @param enable_loo A boolean. Enable model comparison by the R package LOO. This is helpful when you want to compare the fit between two models, for example, analogously to ANOVA, between a one factor model versus a interceot-only model.
#' @param approximate_posterior_inference DEPRECATED please use the `variational_inference` argument.
#'
#' @return A nested tibble `tbl`, with the following columns
#' \itemize{
#'   \item cell_group - column including the cell groups being tested
#'   \item parameter - The parameter being estimated, from the design matrix dscribed with the input formula_composition and formula_variability
#'   \item factor - The factor in the formula corresponding to the covariate, if exists (e.g. it does not exist in case og Intercept or contrasts, which usually are combination of parameters)
#'
#'   \item c_lower - lower (2.5%) quantile of the posterior distribution for a composition (c) parameter.
#'   \item c_effect - mean of the posterior distribution for a composition (c) parameter.
#'   \item c_upper - upper (97.5%) quantile of the posterior distribution fo a composition (c)  parameter.
#'   \item c_n_eff - Effective sample size - the number of independent draws in the sample, the higher the better (mc-stan.org/docs/2_25/cmdstan-guide/stansummary.html).
#'   \item c_R_k_hat - R statistic, a measure of chain equilibrium, should be within 0.05 of 1.0 (mc-stan.org/docs/2_25/cmdstan-guide/stansummary.html).
#'
#'   \item v_lower - Lower (2.5%) quantile of the posterior distribution for a variability (v) parameter
#'   \item v_effect - Mean of the posterior distribution for a variability (v) parameter
#'   \item v_upper - Upper (97.5%) quantile of the posterior distribution for a variability (v) parameter
#'   \item v_n_eff - Effective sample size for a variability (v) parameter - the number of independent draws in the sample, the higher the better (mc-stan.org/docs/2_25/cmdstan-guide/stansummary.html).
#'   \item v_R_k_hat - R statistic for a variability (v) parameter, a measure of chain equilibrium, should be within 0.05 of 1.0 (mc-stan.org/docs/2_25/cmdstan-guide/stansummary.html).
#'
#'   \item count_data Nested input count data.
#'
#' }
#'
#' @examples
#'
#' data("counts_obj")
#'
#' estimate =
#'   sccomp_estimate(
#'   counts_obj ,
#'    ~ type,
#'    ~1,
#'    sample,
#'    cell_group,
#'    count,
#'     cores = 1
#'   ) |>
#'   sccomp_remove_outliers(cores = 1)
#'
#' @export
#'
#'
sccomp_remove_outliers <- function(.estimate,
                                   percent_false_positive = 5,
                                   cores = detectCores(),
                                   inference_method = "pathfinder",
                                   verbose = TRUE,
                                   mcmc_seed = sample(1e5, 1),
                                   max_sampling_iterations = 20000,
                                   enable_loo = FALSE,
                                   
                                   # DEPRECATED
                                   approximate_posterior_inference = NULL,
                                   variational_inference = NULL
) {
  if(inference_method == "variational") 
    rlang::inform(
      message = "sccomp says: From version 1.7.7 the model by default is fit with the variational inference method (inference_method = “variational”; much faster). For a full Bayesian inference (HMC method; the gold standard) use inference_method = “hmc”.", 
      .frequency = "once", 
      .frequency_id = "variational_message"
    )
    
  UseMethod("sccomp_remove_outliers", .estimate)
}

#' @importFrom readr write_file
#' @export
sccomp_remove_outliers.sccomp_tbl = function(.estimate,
                                             percent_false_positive = 5,
                                             cores = detectCores(),
                                             inference_method = "pathfinder",
                                             verbose = TRUE,
                                             mcmc_seed = sample(1e5, 1),
                                             max_sampling_iterations = 20000,
                                             enable_loo = FALSE,
                                             
                                             # DEPRECATED
                                             approximate_posterior_inference = NULL,
                                             variational_inference = NULL
) {
  
  
  # DEPRECATION OF approximate_posterior_inference
  if (is_present(approximate_posterior_inference) & !is.null(approximate_posterior_inference)) {
    deprecate_warn("1.7.7", "sccomp::sccomp_estimate(approximate_posterior_inference = )", details = "The argument approximate_posterior_inference is now deprecated please use variational_inference. By default variational_inference value is inferred from approximate_posterior_inference.")
     inference_method = ifelse(approximate_posterior_inference == "all", "variational","hmc")
  }
  
  # DEPRECATION OF variational_inference
  if (is_present(variational_inference) & !is.null(variational_inference)) {
    deprecate_warn("1.7.11", "sccomp::sccomp_estimate(variational_inference = )", details = "The argument variational_inference is now deprecated please use variational_inference. By default inference_method value is inferred from variational_inference")
    
    inference_method = ifelse(variational_inference, "variational","hmc")
  }
  
  # Prepare column same enquo
  .sample = .estimate |>  attr(".sample")
  .cell_group = .estimate |>  attr(".cell_group")
  .count = .estimate |>  attr(".count")
  .sample_cell_group_pairs_to_exclude = .estimate |> attr(".sample_cell_group_pairs_to_exclude")
  
  # Formulae
  formula_composition = .estimate |> attr("formula_composition")
  formula_variability = .estimate |> attr("formula_variability")
  
  noise_model = .estimate |> attr("noise_model")
  
  # Get model input
  data_for_model = .estimate |> attr("model_input")
  
  # Credible interval
  CI = 1 - (percent_false_positive/100)
  
  # Count data
  .data = 
    .estimate |> 
    select(!!.cell_group, count_data) |> 
    unnest(count_data) |> 
    distinct() |> 
    
    # Drop previous outlier estimation for the new one
    select(-any_of(c(
      ".lower" ,  ".median" ,  ".upper" , "Rhat" ,  "truncation_down" , 
      "truncation_up"   , "outlier"    ,   "contains_outliers"
    )))
  
  # Random intercept
  random_intercept_elements = .estimate |> attr("formula_composition") |> parse_formula_random_intercept()
  
  # Load model
  mod_rng = load_model("glm_multi_beta_binomial_generate_data")

  rng = mod_rng$generate_quantities(
    attr(.estimate , "fit")$draws(format = "matrix"),
    
    # This is for the new data generation with selected factors to do adjustment
    data = 
      .estimate |>
      attr("model_input") |> 
      c(list(
        
        # Add subset of coefficients
        X_original = data_for_model$X,
        N_original = data_for_model$N,
        length_X_which = ncol(data_for_model$X),
        length_XA_which = ncol(data_for_model$XA),
        X_which = seq_len(ncol(data_for_model$X)) |> as.array(),
        XA_which = seq_len(ncol(data_for_model$Xa)) |> as.array(),
        
        # Random intercept
        N_grouping_new = ncol(data_for_model$X_random_intercept), # I could put this in the intial data
        length_X_random_intercept_which = ncol(data_for_model$X_random_intercept),
        X_random_intercept_which = seq_len(ncol(data_for_model$X_random_intercept)) |> as.array(),
        create_intercept = FALSE
      )),
    parallel_chains = ifelse(data_for_model$is_vb, 1, attr(.estimate , "fit")$num_chains()), 
    threads_per_chain = cores
    
  )
  
  # Free memory
  rm(.estimate)
  
  # Detect outliers
  truncation_df =
    .data %>%
    left_join(
      summary_to_tibble(rng, "counts", "N", "M", probs = c(0.05, 0.95)) %>%
        nest(data = -N) %>%
        mutate(!!.sample := rownames(data_for_model$y)) %>%
        unnest(data) %>%
        nest(data = -M) %>%
        mutate(!!.cell_group := colnames(data_for_model$y)) %>%
        unnest(data) ,
      
      by = c(quo_name(.sample), quo_name(.cell_group))
    ) %>%
    
    # Add truncation
    mutate(   truncation_down = `5%`,   truncation_up =  `95%`) %>%
    
    # Add outlier stats
    mutate( outlier = !(!!.count >= `5%` & !!.count <= `95%`) ) %>%
    nest(data = -M) %>%
    mutate(contains_outliers = map_lgl(data, ~ .x %>% filter(outlier) %>% nrow() %>% `>` (0))) %>%
    unnest(data) %>%
    
    mutate(
      truncation_down = case_when( outlier ~ -1, TRUE ~ truncation_down),
      truncation_up = case_when(outlier ~ -1, TRUE ~ truncation_up),
    )
  
  # Add censoring
  data_for_model$is_truncated = 1
  data_for_model$truncation_up = truncation_df %>% select(N, M, truncation_up) %>% spread(M, truncation_up) %>% as_matrix(rownames = "N") %>% apply(2, as.integer)
  data_for_model$truncation_down = truncation_df %>% select(N, M, truncation_down) %>% spread(M, truncation_down) %>% as_matrix(rownames = "N") %>% apply(2, as.integer)
  data_for_model$truncation_not_idx = 
    (data_for_model$truncation_down >= 0) %>% 
    t() %>% 
    as.vector()  %>% 
    which() |>
    intersect(data_for_model$user_forced_truncation_not_idx) |>
    sort()
  data_for_model$TNS = length(data_for_model$truncation_not_idx)
  
  message("sccomp says: outlier identification - step 1/2")
  
  my_quantile_step_2 = 1 - (0.1 / data_for_model$N)
  
  # This step gets the credible interval to control for within-category false positive rate
  # We want a category-wise false positive rate of 0.1, and we have to correct for how many samples we have in each category
  CI_step_2 = (1-my_quantile_step_2) / 2 * 2
  
  fit2 =
    data_for_model %>%
    fit_model(
      stanmodels$glm_multi_beta_binomial,
      cores = cores,
      quantile = my_quantile_step_2,
      inference_method = inference_method,
      verbose = verbose,
      seed = mcmc_seed,
      max_sampling_iterations = max_sampling_iterations,
      pars = c("beta", "alpha", "prec_coeff", "prec_sd",   "alpha_normalised", "beta_random_intercept")
    )

  
  #fit_model(stan_model("inst/stan/glm_multi_beta_binomial.stan"), chains= 4, output_samples = 500, approximate_posterior_inference = FALSE, verbose = TRUE)
  
  rng2 =  mod_rng$generate_quantities(
    fit2$draws(format = "matrix"),
    
    # This is for the new data generation with selected factors to do adjustment
    data = data_for_model |> c(list(
      
      # Add subset of coefficients
      X_original = data_for_model$X,
      N_original = data_for_model$N,
      length_X_which = ncol(data_for_model$X),
      length_XA_which = ncol(data_for_model$XA),
      X_which = seq_len(ncol(data_for_model$X)) |> as.array(),
      XA_which = seq_len(ncol(data_for_model$Xa)) |> as.array(),
      
      # Random intercept
      N_grouping_new = ncol(data_for_model$X_random_intercept), # I could put this in the intial data
      length_X_random_intercept_which = ncol(data_for_model$X_random_intercept),
      X_random_intercept_which = seq_len(ncol(data_for_model$X_random_intercept)) |> as.array(),
      create_intercept = FALSE
      
    )),
    parallel_chains = ifelse(data_for_model$is_vb, 1, fit2$num_chains()), 
    threads_per_chain = cores
    
  )
  
  rng2_summary = 
    summary_to_tibble(rng2, "counts", "N", "M", probs = c(CI_step_2, 0.5, 1-CI_step_2)) 
  
  column_quantile_names = rng2_summary |> colnames() |> str_subset("\\%") |> _[c(1,3)]
  
  rng2_summary = 
    rng2_summary %>%
    
    # !!! THIS COMMAND RELIES ON POSITION BECAUSE IT'S NOT TRIVIAL TO MATCH
    # !!! COLUMN NAMES BASED ON LIMITED PRECISION AND/OR PERIODICAL QUANTILES
    rename(
      .lower := !!as.symbol(column_quantile_names[1]) ,
      .median = `50%`,
      .upper := !!as.symbol(column_quantile_names[2])
    ) %>%
    nest(data = -N) %>%
    mutate(!!.sample := rownames(data_for_model$y)) %>%
    unnest(data) %>%
    nest(data = -M) %>%
    mutate(!!.cell_group := colnames(data_for_model$y)) %>%
    unnest(data) 
  
  # Detect outliers
  truncation_df2 =
    .data %>%
    left_join(rng2_summary,
      by = c(quo_name(.sample), quo_name(.cell_group))
    ) %>%
    
    # Add truncation
    mutate(   truncation_down = .lower,   truncation_up =  .upper) %>%
    
    # Add outlier stats
    mutate( outlier = !(!!.count >= .lower & !!.count <= .upper) ) %>%
    nest(data = -M) %>%
    mutate(contains_outliers = map_lgl(data, ~ .x %>% filter(outlier) %>% nrow() %>% `>` (0))) %>%
    unnest(data) %>%
    
    mutate(
      truncation_down = case_when( outlier ~ -1, TRUE ~ truncation_down),
      truncation_up = case_when(outlier ~ -1, TRUE ~ truncation_up)
    )
  
  data_for_model$truncation_up = truncation_df2 %>% select(N, M, truncation_up) %>% spread(M, truncation_up) %>% as_matrix(rownames = "N") %>% apply(2, as.integer)
  data_for_model$truncation_down = truncation_df2 %>% select(N, M, truncation_down) %>% spread(M, truncation_down) %>% as_matrix(rownames = "N") %>% apply(2, as.integer)
  data_for_model$truncation_not_idx = 
    (data_for_model$truncation_down >= 0) %>% 
    t() %>% 
    as.vector()  %>% 
    which() |>
    intersect(data_for_model$user_forced_truncation_not_idx) |>
    sort()
  data_for_model$TNS = length(data_for_model$truncation_not_idx)
  
  # LOO
  data_for_model$enable_loo = TRUE & enable_loo
  
  message("sccomp says: outlier-free model fitting - step 2/2")
  
  # Print design matrix
  message(sprintf("sccomp says: the composition design matrix has columns: %s", data_for_model$X %>% colnames %>% paste(collapse=", ")))
  message(sprintf("sccomp says: the variability design matrix has columns: %s", data_for_model$Xa %>% colnames %>% paste(collapse=", ")))
  
  fit3 =
    data_for_model %>%
    # Run the first discovery phase with permissive false discovery rate
    fit_model(
      stanmodels$glm_multi_beta_binomial,
      cores = cores,
      quantile = CI,
      inference_method = inference_method,
      verbose = verbose, 
      seed = mcmc_seed,
      max_sampling_iterations = max_sampling_iterations,
      pars = c("beta", "alpha", "prec_coeff","prec_sd",   "alpha_normalised", "beta_random_intercept", "log_lik")
    )
  
  # Create a dummy tibble
  tibble() |>
    # Attach association mean concentration
    add_attr(fit3, "fit") %>%
    add_attr(data_for_model, "model_input") |>
    add_attr(truncation_df2, "truncation_df2") |>
    add_attr(.sample, ".sample") |>
    add_attr(.cell_group, ".cell_group") |>
    add_attr(.count, ".count") |>
    
    add_attr(formula_composition, "formula_composition") |>
    add_attr(formula_variability, "formula_variability") |>
    add_attr(parse_formula(formula_composition), "factors" ) |> 
    
    # Add class to the tbl
    add_class("sccomp_tbl") |> 
    
    # Print estimates
    sccomp_test() |>
    
    # drop hypothesis testing as the estimation exists without probabilities.
    # For hypothesis testing use sccomp_test
    select(-contains("_FDR"), -contains("_pH0")) |> 
    
    add_attr(noise_model, "noise_model") 
  
  
}

#' sccomp_test
#'
#' @description This function test contrasts from a sccomp result.
#'
#'
#' @param .data A tibble. The result of sccomp_estimate.
#' @param contrasts A vector of character strings. For example if your formula is `~ 0 + treatment` and the factor treatment has values `yes` and `no`, your contrast could be "constrasts = c(treatmentyes - treatmentno)".
#' @param percent_false_positive A real between 0 and 100 non included. This used to identify outliers with a specific false positive rate.
#' @param test_composition_above_logit_fold_change A positive integer. It is the effect threshold used for the hypothesis test. A value of 0.2 correspond to a change in cell proportion of 10% for a cell type with baseline proportion of 50%. That is, a cell type goes from 45% to 50%. When the baseline proportion is closer to 0 or 1 this effect thrshold has consistent value in the logit uncontrained scale.
#' @param pass_fit A boolean. Whether to pass the Stan fit as attribute in the output. Because the Stan fit can be very large, setting this to FALSE can be used to lower the memory imprint to save the output.
#'
#' @return A nested tibble `tbl` with cell_group-wise statistics
#'
#' @export
#'
#' @examples
#'
#' data("counts_obj")
#'
#'   estimates =
#'   sccomp_estimate(
#'   counts_obj ,
#'    ~ 0 + type, ~1,  sample, cell_group, count,
#'     cores = 1
#'   ) |>
#'
#'   sccomp_test("typecancer - typebenign")
#'
sccomp_test <- function(.data,
                           contrasts = NULL,
                           percent_false_positive = 5,
                           test_composition_above_logit_fold_change = 0.1,
                           pass_fit = TRUE) {
  UseMethod("sccomp_test", .data)
}


#'
#' @importFrom dplyr any_of
#'
#' @export
#'
#'
sccomp_test.sccomp_tbl = function(.data,
                                     contrasts = NULL,
                                     percent_false_positive = 5,
                                     test_composition_above_logit_fold_change = 0.1,
                                     pass_fit = TRUE){


  .sample = .data |>  attr(".sample")
  .cell_group = .data |>  attr(".cell_group")
  .count = .data |>  attr(".count")
  model_input = .data |> attr("model_input")
  truncation_df2 =  .data |>  attr("truncation_df2")

  # Abundance
  abundance_CI =
    get_abundance_contrast_draws(.data, contrasts)

  # If my contrasts do not match my model. I have to do something more elegant.
  if ("parameter" %in% colnames(abundance_CI))
    abundance_CI =
      abundance_CI |>
      draws_to_statistics(
        percent_false_positive/100,
        test_composition_above_logit_fold_change,
        !!.cell_group,
        "c_"
      )

  variability_CI =
    get_variability_contrast_draws(.data, contrasts)
  
  # Variability
  if ("parameter" %in% colnames(variability_CI))
  variability_CI = 
    variability_CI |>
    draws_to_statistics(
      percent_false_positive / 100,
      test_composition_above_logit_fold_change,
      !!.cell_group,
      "v_"
    )

  # If I don't have factors (~1)
  if (!"factor" %in% colnames(model_input$factor_parameter_dictionary))
    factor_parameter_dictionary = tibble(`factor` = character(), design_matrix_col = character())
  else
    factor_parameter_dictionary =
    model_input$factor_parameter_dictionary |>
    select(`factor`, design_matrix_col)

  # Merge and parse
  result =
    abundance_CI |>

    # Add ALPHA
    left_join(variability_CI) |>
    suppressMessages() |>

    # Add easy to understand factor labels
    left_join(factor_parameter_dictionary,
              by = c("parameter" = "design_matrix_col")) %>%
    select(parameter, `factor`, everything()) %>%

    select(!!.cell_group, everything(),-M)


  result =
    result |>
    left_join(
      truncation_df2 |>
        select(-any_of(c(
          "M",
          "N",
          ".variable",
          "mean",
          "se_mean",
          "sd",
          "n_eff",
          "R_hat",
          "k_hat",
          "Rhat",
          ".lower", ".median", ".upper"
        ))) |>
        nest(count_data = -!!.cell_group),
      by = quo_name(.cell_group)
    )

  result =
    result |>

    # Add back attributes
    add_attr(
      .data |> attr("fit") |> get_mean_precision_association(),
      "mean_concentration_association"
    )

  if(pass_fit)
    result =
      result |>
      add_attr(.data |> attr("fit") , "fit")

  result |>
    
    add_attr(test_composition_above_logit_fold_change, "test_composition_above_logit_fold_change") |>
    
    # Attach association mean concentration
    add_attr(.data |> attr("model_input") , "model_input") |>
    add_attr(.data |> attr("truncation_df2"), "truncation_df2") |>
    add_attr(.data |> attr("noise_model") , "noise_model") |>

    add_attr(.sample, ".sample") |>
    add_attr(.cell_group, ".cell_group") |>
    add_attr(.count, ".count") |>

    add_attr(.data |> attr("formula_composition"), "formula_composition") |>
    add_attr(.data |> attr("formula_variability"), "formula_variability") |>
    
    # Add class to the tbl
    add_class("sccomp_tbl") 
}


#' sccomp_replicate
#'
#' @description This function replicates counts from a real-world dataset.
#'
#'
#' @param fit The result of sccomp_estimate.
#' @param formula_composition A formula. The formula describing the model for differential abundance, for example ~treatment. This formula can be a sub-formula of your estimated model; in this case all other factor will be factored out.
#' @param formula_variability A formula. The formula describing the model for differential variability, for example ~treatment. In most cases, if differentially variability is of interest, the formula should only include the factor of interest as a large anount of data is needed to define variability depending to each factors. This formula can be a sub-formula of your estimated model; in this case all other factor will be factored out.
#' @param number_of_draws An integer. How may copies of the data you want to draw from the model joint posterior distribution.
#' @param mcmc_seed An integer. Used for Markov-chain Monte Carlo reproducibility. By default a random number is sampled from 1 to 999999. This itself can be controlled by set.seed()
#'
#' @return A nested tibble `tbl` with cell_group-wise statistics
#'
#' @export
#'
#' @examples
#'
#' data("counts_obj")
#'
#' if(.Platform$OS.type == "unix")
#'   sccomp_estimate(
#'   counts_obj ,
#'    ~ type, ~1,  sample, cell_group, count,
#'     cores = 1
#'   ) |>
#'
#'   sccomp_replicate()
#'
sccomp_replicate <- function(fit,
                             formula_composition = NULL,
                             formula_variability = NULL,
                             number_of_draws = 1,
                             mcmc_seed = sample(1e5, 1)) {
  UseMethod("sccomp_replicate", fit)
}

#' @export
#'
sccomp_replicate.sccomp_tbl = function(fit,
                                       formula_composition = NULL,
                                       formula_variability = NULL,
                                       number_of_draws = 1,
                                       mcmc_seed = sample(1e5, 1)){

  .sample = attr(fit, ".sample")
  .cell_group = attr(fit, ".cell_group")

  rng =
    replicate_data(
      fit,
      formula_composition = formula_composition,
      formula_variability = formula_variability,
      number_of_draws = number_of_draws,
      mcmc_seed = mcmc_seed
    )

  model_input = attr(fit, "model_input")

  # mean generated
  rng |>

    # Parse
    parse_generated_quantities(number_of_draws = number_of_draws) %>%

    # Get sample name
    nest(data = -N) %>%
    arrange(N) %>%
    mutate(!!.sample := rownames(model_input$y)) %>%
    unnest(data) %>%

    # get cell type name
    nest(data = -M) %>%
    mutate(!!.cell_group := colnames(model_input$y)) %>%
    unnest(data) %>%

    select(-N, -M) |>
    select(!!.cell_group, !!.sample, everything())

}

#' sccomp_predict
#'
#' @description This function replicates counts from a real-world dataset.
#'
#'
#' @param fit The result of sccomp_estimate.
#' @param formula_composition A formula. The formula describing the model for differential abundance, for example ~treatment. This formula can be a sub-formula of your estimated model; in this case all other factor will be factored out.
#' @param new_data A sample-wise data frame including the column that represent the factors in your formula. If you want to predict proportions for 10 samples, there should be 10 rows. T
#' @param number_of_draws An integer. How may copies of the data you want to draw from the model joint posterior distribution.
#' @param mcmc_seed An integer. Used for Markov-chain Monte Carlo reproducibility. By default a random number is sampled from 1 to 999999. This itself can be controlled by set.seed()
#'
#' @return A nested tibble `tbl` with cell_group-wise statistics
#'
#' @export
#'
#' @examples
#'
#' data("counts_obj")
#'
#' if(.Platform$OS.type == "unix")
#'   sccomp_estimate(
#'   counts_obj ,
#'    ~ type, ~1,  sample, cell_group, count,
#'     cores = 1
#'   ) |>
#'
#'   sccomp_predict()
#'
sccomp_predict <- function(fit,
                           formula_composition = NULL,
                           new_data = NULL,
                           number_of_draws = 500,
                           mcmc_seed = sample(1e5, 1)) {
  UseMethod("sccomp_predict", fit)
}

#' @export
#'
sccomp_predict.sccomp_tbl = function(fit,
                                     formula_composition = NULL,
                                     new_data = NULL,
                                     number_of_draws = 500,
                                     mcmc_seed = sample(1e5, 1)){



  model_input = attr(fit, "model_input")
  .sample = attr(fit, ".sample")
  .cell_group = attr(fit, ".cell_group")
  
  rng =
    replicate_data(
      fit,
      formula_composition = formula_composition,
      formula_variability = ~ 1,
      new_data = new_data,
      number_of_draws = number_of_draws,
      mcmc_seed = mcmc_seed
    )

  # New data
  if(new_data |> is.null())
    sample_names =
      fit |>
      select(count_data) |>
      unnest(count_data) |>
      distinct(!!.sample) |> 
      pull(!!.sample)
  
  # If seurat
  else if(new_data |> is("Seurat")) 
    sample_names = 
      new_data[[]] |> 
      distinct(!!.sample) |> 
      pull(!!.sample)
  
  # Just subset
  else 
    sample_names = 
      new_data |> 
      distinct(!!.sample) |> 
      pull(!!.sample)
  
  # mean generated
  rng |>
    summary_to_tibble("mu", "M", "N") |>
    select(M, N, proportion_mean = mean, proportion_lower = `2.5%`, proportion_upper = `97.5%`) |>

    # Get sample name
    nest(data = -N) %>%
    arrange(N) %>%
    mutate(!!.sample := sample_names) %>%
    unnest(data) %>%

    # get cell type name
    nest(data = -M) %>%
    mutate(!!.cell_group := colnames(model_input$y)) %>%
    unnest(data) %>%

    select(-N, -M) |>
    select(!!.cell_group, !!.sample, everything())

}


#' sccomp_remove_unwanted_variation
#'
#' @description This function uses the model to remove unwanted variation from a dataset using the estimated of the model. For example if you fit your data with this formula `~ factor_1 + factor_2` and use this formula to remove unwanted variation `~ factor_1`, the `factor_2` will be factored out.
#'
#'
#' @param .data A tibble. The result of sccomp_estimate.
#' @param formula_composition A formula. The formula describing the model for differential abundance, for example ~treatment. This formula can be a sub-formula of your estimated model; in this case all other factor will be factored out.
#' @param formula_variability A formula. The formula describing the model for differential variability, for example ~treatment. In most cases, if differentially variability is of interest, the formula should only include the factor of interest as a large anount of data is needed to define variability depending to each factors. This formula can be a sub-formula of your estimated model; in this case all other factor will be factored out.
#'
#' @return A nested tibble `tbl` with cell_group-wise statistics
#'
#' @export
#'
#' @examples
#'
#' data("counts_obj")
#'
#'   estimates = sccomp_estimate(
#'   counts_obj ,
#'    ~ type, ~1,  sample, cell_group, count,
#'     cores = 1
#'   )
#'
#'   sccomp_remove_unwanted_variation(estimates)
#'
sccomp_remove_unwanted_variation <- function(.data,
                                      formula_composition = ~1,
                                      formula_variability = NULL) {
  UseMethod("sccomp_remove_unwanted_variation", .data)
}

#' @importFrom readr write_file
#' @export
#'
sccomp_remove_unwanted_variation.sccomp_tbl = function(.data,
                                                formula_composition = ~1,
                                                formula_variability = NULL){


  model_input = attr(.data, "model_input")
  .sample = attr(.data, ".sample")
  .cell_group = attr(.data, ".cell_group")
  .grouping_for_random_intercept = attr(.data, ".grouping_for_random_intercept")
  .count = attr(.data, ".count")

  fit = attr(.data, "fit")

  # Load model
  mod = load_model("glm_multi_beta_binomial_generate_data")
  


  message("sccomp says: calculating residuals")

  # Residuals
  residuals =
    .data |>
    sccomp_predict(
      number_of_draws = attr(.data, "fit") |>  get_output_samples() |> min(500)
    ) |>
    distinct(!!.sample, !!.cell_group, proportion_mean) |>
    mutate( proportion_mean =
              proportion_mean |>
             #compress_zero_one() |>
             boot::logit()
    )|>
  # Join counts
  left_join(
    .data |>
      attr("model_input") %$%
      y |>
      as_tibble(rownames = quo_name(.sample)) |>
      pivot_longer(-!!.sample, names_to = quo_name(.cell_group), values_to = quo_name(.count)) |>
      with_groups(!!.sample,  ~ .x |> mutate(observed_proportion := !!.count / sum(!!.count ))) |>

      with_groups(!!.sample,  ~ .x |>  mutate(exposure := sum(!!.count))  ) |>


      mutate(observed_proportion =
               observed_proportion |>
               compress_zero_one() |>
               boot::logit()
      ),
    by = c(quo_name(.sample), quo_name(.cell_group))
  ) |>
  mutate(logit_residuals = observed_proportion - proportion_mean) |>
  select(!!.sample, !!.cell_group, logit_residuals, exposure)

  message("sccomp says: regressing out unwanted factors")


  # Generate quantities
  .data |>
    sccomp_predict(
      formula_composition = formula_composition,
      number_of_draws = attr(.data, "fit") |>  get_output_samples() |> min(500)
    ) |>
    distinct(!!.sample, !!.cell_group, proportion_mean) |>
    mutate(proportion_mean =
             proportion_mean |>
             # compress_zero_one() |>
             boot::logit()
    ) |>
    left_join(residuals,  by = c(quo_name(.sample), quo_name(.cell_group))) |>
    mutate(adjusted_proportion = proportion_mean + logit_residuals) |>
    mutate(adjusted_proportion = adjusted_proportion |> boot::inv.logit()) |>
    with_groups(!!.sample,  ~ .x |> mutate(adjusted_proportion := adjusted_proportion / sum(adjusted_proportion ))) |>

    # Recostituite counts
    mutate(adjusted_counts = adjusted_proportion * exposure) |>

    select(!!.sample, !!.cell_group, adjusted_proportion, adjusted_counts, logit_residuals)



}

#' simulate_data
#'
#' @description This function simulates counts from a linear model.
#'
#' @import dplyr
#' @importFrom magrittr %$%
#' @importFrom magrittr divide_by
#' @importFrom magrittr multiply_by
#' @importFrom magrittr equals
#' @importFrom rlang quo_is_null
#' @importFrom SingleCellExperiment colData
#' @importFrom parallel detectCores
#'
#' @param .data A tibble including a cell_group name column | sample name column | read counts column | factor columns | Pvalue column | a significance column
#' @param .estimate_object The result of sccomp_estimate execution. This is used for sampling from real-data properties.
#' @param formula_composition A formula. The sample formula used to perform the differential cell_group abundance analysis
#' @param formula_variability A formula. The formula describing the model for differential variability, for example ~treatment. In most cases, if differentially variability is of interest, the formula should only include the factor of interest as a large anount of data is needed to define variability depending to each factors.
#' @param .sample A column name as symbol. The sample identifier
#' @param .cell_group A column name as symbol. The cell_group identifier
#' @param .coefficients The column names for coefficients, for example, c(b_0, b_1)
#' @param variability_multiplier A real scalar. This can be used for artificially increasing the variability of the simulation for benchmarking purposes.
#' @param number_of_draws An integer. How may copies of the data you want to draw from the model joint posterior distribution.
#' @param mcmc_seed An integer. Used for Markov-chain Monte Carlo reproducibility. By default a random number is sampled from 1 to 999999. This itself can be controlled by set.seed()
#'
#' @return A nested tibble `tbl` with cell_group-wise statistics
#'
#' @export
#'
#' @examples
#'
#' data("counts_obj")
#' library(dplyr)
#'
#' estimate =
#'  sccomp_estimate(
#'  counts_obj ,
#'   ~ type, ~1,  sample, cell_group, count,
#'    cores = 1
#'  )
#'
#' # Set coefficients for cell_groups. In this case all coefficients are 0 for simplicity.
#' counts_obj = counts_obj |> mutate(b_0 = 0, b_1 = 0)
#' # Simulate data
#' simulate_data(counts_obj, estimate, ~type, ~1, sample, cell_group, c(b_0, b_1))
#'
simulate_data <- function(.data,
                          .estimate_object,
                          formula_composition,
                          formula_variability = NULL,
                       .sample = NULL,
                       .cell_group = NULL,
                       .coefficients = NULL,
                       variability_multiplier = 5,
                       number_of_draws = 1,
                       mcmc_seed = sample(1e5, 1)) {
  UseMethod("simulate_data", .data)
}

#' @export
#' @importFrom magrittr divide_by
#' @importFrom magrittr multiply_by
#' @importFrom SingleCellExperiment counts
#' @importFrom purrr map_dbl
#' @importFrom purrr pmap
#' @importFrom readr read_file
#'
simulate_data.tbl = function(.data,
                                    .estimate_object,
                                    formula_composition,
                                    formula_variability = NULL,
                                    .sample = NULL,
                                    .cell_group = NULL,
                                    .coefficients = NULL,
                                    variability_multiplier = 5,
                                    number_of_draws = 1,
                                    mcmc_seed = sample(1e5, 1)){


  .sample = enquo(.sample)
  .cell_group = enquo(.cell_group)
  .coefficients = enquo(.coefficients)

  #Check column class
  check_if_columns_right_class(.data, !!.sample, !!.cell_group)

  model_data = attr(.estimate_object, "model_input")

  # Select model based on noise model
  if(attr(.estimate_object, "noise_model") == "multi_beta_binomial") my_model = stanmodels$glm_multi_beta_binomial_simulate_data
  else if(attr(.estimate_object, "noise_model") == "dirichlet_multinomial") my_model = get_model_from_data("model_glm_dirichlet_multinomial_generate_quantities.rds", glm_dirichlet_multinomial_generate_quantities)
  else if(attr(.estimate_object, "noise_model") == "logit_normal_multinomial") my_model = get_model_from_data("glm_multinomial_logit_linear_simulate_data.stan", read_file("~/PostDoc/sccomp/dev/stan_models/glm_multinomial_logit_linear_simulate_data.stan"))


  model_input =
    .data %>%
    nest(data___ = -!!.sample) %>%
    mutate(.exposure = sample(model_data$exposure, size = n(), replace = TRUE )) %>%
    unnest(data___) %>%
    data_simulation_to_model_input(
      formula_composition,
      #formula_variability,
      !!.sample, !!.cell_group, .exposure, !!.coefficients
    )

    # [1]  5.6260004 -0.6940178
    # prec_sd  = 0.816423129

  fit =
    rstan::gqs(
    my_model,
    draws =  as.matrix(attr(.estimate_object, "fit") ),
    data = model_input %>% c(list(variability_multiplier = variability_multiplier)),
    seed = mcmc_seed
  )

  parsed_fit =
    fit %>%
    parse_generated_quantities(number_of_draws = number_of_draws) %>%

    # Get sample name
    nest(data = -N) %>%
    arrange(N) %>%
    mutate(!!.sample := rownames(model_input$X)) %>%
    unnest(data) %>%

    # get cell type name
    nest(data = -M) %>%
    mutate(!!.cell_group := colnames(model_input$beta)) %>%
    unnest(data) %>%

    select(-N, -M)

  .data %>%
    left_join(
      parsed_fit,
      by = c(quo_name(.sample), quo_name(.cell_group))
    )


}

#' sccomp_boxplot
#'
#' @description This function plots a boxplot of the results of the model.
#'
#' @importFrom ggrepel geom_text_repel
#' @importFrom tidyr pivot_longer
#' @import ggplot2
#' @importFrom tidyr unite
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr with_groups
#'
#' @param .data A tibble including a cell_group name column | sample name column | read counts column | factor columns | Pvalue column | a significance column
#' @param factor A character string for a factor of interest included in the model
#' @param significance_threshold A real. FDR threshold for labelling significant cell-groups.
#' @param test_composition_above_logit_fold_change A positive integer. It is the effect threshold used for the hypothesis test. A value of 0.2 correspond to a change in cell proportion of 10% for a cell type with baseline proportion of 50%. That is, a cell type goes from 45% to 50%. When the baseline proportion is closer to 0 or 1 this effect thrshold has consistent value in the logit uncontrained scale.
#'
#' @return A `ggplot`
#'
#' @export
#'
#' @examples
#'
#' data("counts_obj")
#'
#' estimate =
#'   sccomp_estimate(
#'   counts_obj ,
#'    ~ type, ~1, sample, cell_group, count,
#'     cores = 1
#'   ) |>
#'   sccomp_test()
#'
#' # estimate |> sccomp_boxplot()
sccomp_boxplot = function(.data, factor, significance_threshold = 0.05, test_composition_above_logit_fold_change = .data |> attr("test_composition_above_logit_fold_change")){


  .cell_group = attr(.data, ".cell_group")
  .count = attr(.data, ".count")
  .sample = attr(.data, ".sample")

  # Check if test have been done
  if(.data |> select(ends_with("FDR")) |> ncol() |> equals(0))
    stop("sccomp says: to produce plots, you need to run the function sccomp_test() on your estimates.")
  
  
  data_proportion =
    .data %>%

    # Otherwise does not work
    select(-`factor`) %>%

    pivot_wider(names_from = parameter, values_from = c(contains("c_"), contains("v_"))) %>%
    unnest(count_data) %>%
    with_groups(!!.sample, ~ mutate(.x, proportion = (!!.count)/sum(!!.count)) ) 
  
  # If I don't have outliers add them
  if(!"outlier" %in% colnames(data_proportion)) data_proportion = data_proportion |> mutate(outlier = FALSE) 

  plot_boxplot(
        .data,
        data_proportion,
        factor,
        !!.cell_group,
        !!.sample,
        significance_threshold = significance_threshold,
        multipanel_theme
      ) +
    ggtitle(sprintf("Grouped by %s (for multi-factor models, associations could be hardly observable with unidimensional data stratification)", factor))


}

#' plot
#'
#' @description This function plots a summary of the results of the model.
#'
#' @importFrom ggrepel geom_text_repel
#' @importFrom tidyr pivot_longer
#' @import ggplot2
#' @importFrom tidyr unite
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr with_groups
#' @importFrom magrittr equals
#'
#' @param x A tibble including a cell_group name column | sample name column | read counts column | factor columns | Pvalue column | a significance column
#' @param significance_threshold Numeric value specifying the significance threshold for highlighting differences. Default is 0.025.
#' @param test_composition_above_logit_fold_change A positive integer. It is the effect threshold used for the hypothesis test. A value of 0.2 correspond to a change in cell proportion of 10% for a cell type with baseline proportion of 50%. That is, a cell type goes from 45% to 50%. When the baseline proportion is closer to 0 or 1 this effect thrshold has consistent value in the logit uncontrained scale.
#' @param ... For internal use 
#'
#' @return A `ggplot`
#'
#' @export
#'
#' @examples
#'
#' data("counts_obj")
#'
#' estimate =
#'   sccomp_estimate(
#'   counts_obj ,
#'    ~ type, ~1, sample, cell_group, count,
#'     cores = 1
#'   )
#'
#' # estimate |> plot()
#'
plot.sccomp_tbl <- function(x,  significance_threshold = 0.05, test_composition_above_logit_fold_change = .data |> attr("test_composition_above_logit_fold_change"), ...) {

  # Define the variables as NULL to avoid CRAN NOTES
  parameter <- NULL
  count_data <- NULL
  v_effect <- NULL
  
  .cell_group = attr(x, ".cell_group")
  .count = attr(x, ".count")
  .sample = attr(x, ".sample")

  plots = list()

  # Check if test have been done
  if(x |> select(ends_with("FDR")) |> ncol() |> equals(0))
    stop("sccomp says: to produce plots, you need to run the function sccomp_test() on your estimates.")

  data_proportion =
    x %>%
  
    # Otherwise does not work
    select(-`factor`) %>%
  
    pivot_wider(names_from = parameter, values_from = c(contains("c_"), contains("v_"))) %>%
    unnest(count_data) %>%
    with_groups(!!.sample, ~ mutate(.x, proportion = (!!.count)/sum(!!.count)) )
  
  
  # If I don't have outliers add them
  if(!"outlier" %in% colnames(data_proportion)) data_proportion = data_proportion |> mutate(outlier = FALSE) 

# Select the factors to plot  
my_factors =
  x %>%
  filter(!is.na(`factor`)) |>
  distinct(`factor`) %>%
  pull(`factor`) 

if(length(my_factors)==0)
  message("sccomp says: the contrasts you have tested do not represent factors. Therefore, plot of the the posterior predictive check will be omitted.")
else {
  # Boxplot
  plots$boxplot =
    
    my_factors |>
    map(
      ~ {
        
        # If variable is continuous
        if(data_proportion |> select(all_of(.x)) |> pull(1) |> is("numeric"))
          my_plot = 
            plot_scatterplot(
              .data = x,
              data_proportion = data_proportion,
              factor_of_interest = .x,
              .cell_group = !!.cell_group,
              .sample =  !!.sample,
              my_theme = multipanel_theme,
              significance_threshold = significance_threshold
            )
        
        # If discrete
        else 
          my_plot = 
            plot_boxplot(
              .data = x,
              data_proportion = data_proportion,
              factor_of_interest = .x,
              .cell_group = !!.cell_group,
              .sample =  !!.sample,
              my_theme = multipanel_theme,
              significance_threshold = significance_threshold
            ) 
        
        # Return
        my_plot +
          ggtitle(sprintf("Grouped by %s (for multi-factor models, associations could be hardly observable with unidimensional data stratification)", .x))
      } )
  
}

# 1D intervals
plots$credible_intervals_1D = plot_1D_intervals(.data = x, significance_threshold = significance_threshold)

# 2D intervals
if("v_effect" %in% colnames(x) && (x |> filter(!is.na(v_effect)) |> nrow()) > 0)  plots$credible_intervals_2D = plot_2D_intervals(.data = x, significance_threshold = significance_threshold)

plots

}

#' Clear Stan Model Cache
#'
#' This function attempts to delete the Stan model cache directory and its contents.
#' If the cache directory does not exist, it prints a message indicating this.
#'
#' @param cache_dir A character string representing the path of the cache directory to delete. Defaults to `sccomp_stan_models_cache_dir`.
#' 
#' @return NULL
#' 
#' @examples
#' \dontrun{
#'   clear_stan_model_cache("path/to/cache_dir")
#' }
#' @export
clear_stan_model_cache <- function(cache_dir = sccomp_stan_models_cache_dir) {
  
  # Check if the directory exists
  if (dir.exists(cache_dir)) {
    # Attempt to delete the directory and its contents
    unlink(cache_dir, recursive = TRUE)
    message("Cache deleted: ", cache_dir)
  } else {
    message("Cache does not exist: ", cache_dir)
  }
}
