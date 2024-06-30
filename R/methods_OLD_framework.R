#' DEPRECATED - sccomp_glm main
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
#'
#' @param .data A tibble including a cell_group name column | sample name column | read counts column (optional depending on the input class) | factor columns.
#' @param formula_composition A formula. The formula describing the model for differential abundance, for example ~treatment.
#' @param formula_variability A formula. The formula describing the model for differential variability, for example ~treatment. In most cases, if differentially variability is of interest, the formula should only include the factor of interest as a large anount of data is needed to define variability depending to each factors.
#' @param .sample A column name as symbol. The sample identifier
#' @param .cell_group A column name as symbol. The cell_group identifier
#' @param .count A column name as symbol. The cell_group abundance (read count). Used only for data frame count output. The variable in this column should be of class integer.
#'
#' @param contrasts A vector of character strings. For example if your formula is `~ 0 + treatment` and the factor treatment has values `yes` and `no`, your contrast could be constrasts = c("treatmentyes - treatmentno").
#' @param prior_mean_variable_association A list of the form list(intercept = c(5, 2), slope = c(0,  0.6), standard_deviation = c(20, 40)). Where for intercept and slope parameters, we specify mean and standard deviation, while for standard deviation, we specify shape and rate. This is used to incorporate prior knowledge about the mean/variability association of cell-type proportions.
#' @param check_outliers A boolean. Whether to check for outliers before the fit.
#' @param bimodal_mean_variability_association A boolean. Whether to model the mean-variability as bimodal, as often needed in the case of single-cell RNA sequencing data, and not usually for CyTOF and microbiome data. The plot summary_plot()$credible_intervals_2D can be used to assess whether the bimodality should be modelled.
#' @param enable_loo A boolean. Enable model comparison by the R package LOO. This is helpful when you want to compare the fit between two models, for example, analogously to ANOVA, between a one factor model versus a interceot-only model.
#'
#' @param percent_false_positive A real between 0 and 100 non included. This used to identify outliers with a specific false positive rate.
#' @param exclude_priors A boolean. Whether to run a prior-free model, for benchmarking purposes.
#' @param use_data A booelan. Whether to sun the model data free. This can be used for prior predictive check.
#' @param max_sampling_iterations An integer. This limit the maximum number of iterations in case a large dataset is used, for limiting the computation time.
#' @param pass_fit A boolean. Whether to pass the Stan fit as attribute in the output. Because the Stan fit can be very large, setting this to FALSE can be used to lower the memory imprint to save the output.
#' @param approximate_posterior_inference A boolean. Whether the inference of the joint posterior distribution should be approximated with variational Bayes. It confers execution time advantage.
#' @param test_composition_above_logit_fold_change A positive integer. It is the effect threshold used for the hypothesis test. A value of 0.2 correspond to a change in cell proportion of 10% for a cell type with baseline proportion of 50%. That is, a cell type goes from 45% to 50%. When the baseline proportion is closer to 0 or 1 this effect thrshold has consistent value in the logit uncontrained scale.
#' @param .sample_cell_group_pairs_to_exclude A column name that includes a boolean variable for the sample/cell-group pairs to be ignored in the fit. This argument is for pro-users.
#' @param verbose A boolean. Prints progression.
#' @param noise_model A character string. The two noise models available are multi_beta_binomial (default) and dirichlet_multinomial.
#' @param cores An integer. How many cored to be used with parallel calculations.
#' @param mcmc_seed An integer. Used for Markov-chain Monte Carlo reproducibility. By default a random number is sampled from 1 to 999999. This itself can be controlled by set.seed()
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
#'   sccomp_glm(
#'   counts_obj ,
#'    ~ type,
#'    ~1,
#'    sample,
#'    cell_group,
#'    count,
#'     check_outliers = FALSE,
#'     cores = 1
#'   )
#'
#' @export
#'
#'
sccomp_glm <- function(.data,
                       formula_composition = ~ 1 ,
                       formula_variability = ~ 1,
                       .sample,
                       .cell_group,
                       .count = NULL,
                       
                       # Secondary arguments
                       contrasts = NULL,
                       prior_mean_variable_association = list(intercept = c(5, 2), slope = c(0,  0.6), standard_deviation = c(20, 40)),
                       check_outliers = TRUE,
                       bimodal_mean_variability_association = FALSE,
                       enable_loo = FALSE,
                       
                       # Tertiary arguments
                       cores = detectCores(),
                       percent_false_positive = 5,
                       approximate_posterior_inference = "none",
                       test_composition_above_logit_fold_change = 0.1, .sample_cell_group_pairs_to_exclude = NULL,
                       verbose = FALSE,
                       noise_model = "multi_beta_binomial",
                       exclude_priors = FALSE,
                       use_data = TRUE,
                       mcmc_seed = sample(1e5, 1),
                       max_sampling_iterations = 20000,
                       pass_fit = TRUE) {
  UseMethod("sccomp_glm", .data)
}

#' @export
sccomp_glm.Seurat = function(.data,
                             formula_composition = ~ 1 ,
                             formula_variability = ~ 1,
                             .sample,
                             .cell_group,
                             .count = NULL,
                             
                             # Secondary arguments
                             contrasts = NULL,
                             prior_mean_variable_association = list(intercept = c(5, 2), slope = c(0,  0.6), standard_deviation = c(20, 40)),
                             check_outliers = TRUE,
                             bimodal_mean_variability_association = FALSE,
                             enable_loo = FALSE,
                             
                             # Tertiary arguments
                             cores = detectCores(),
                             percent_false_positive = 5,
                             approximate_posterior_inference = "none",
                             test_composition_above_logit_fold_change = 0.1, .sample_cell_group_pairs_to_exclude = NULL,
                             verbose = FALSE,
                             noise_model = "multi_beta_binomial",
                             exclude_priors = FALSE,
                             use_data = TRUE,
                             mcmc_seed = sample(1e5, 1),
                             max_sampling_iterations = 20000,
                             pass_fit = TRUE) {
  
  if(!is.null(.count)) stop("sccomp says: .count argument can be used only for data frame input")
  
  # Prepare column same enquo
  .sample = enquo(.sample)
  .cell_group = enquo(.cell_group)
  .sample_cell_group_pairs_to_exclude = enquo(.sample_cell_group_pairs_to_exclude)
  
  .data[[]] %>%
    sccomp_glm(
      formula_composition = formula_composition,
      formula_variability = formula_variability,
      
      !!.sample,!!.cell_group,
      contrasts = contrasts,
      prior_mean_variable_association = prior_mean_variable_association,
      percent_false_positive = percent_false_positive ,
      check_outliers = check_outliers,
      approximate_posterior_inference = approximate_posterior_inference,
      test_composition_above_logit_fold_change = test_composition_above_logit_fold_change, .sample_cell_group_pairs_to_exclude = !!.sample_cell_group_pairs_to_exclude,
      verbose = verbose,
      noise_model = noise_model,
      exclude_priors = exclude_priors,
      bimodal_mean_variability_association = bimodal_mean_variability_association,
      enable_loo = enable_loo,
      use_data = use_data,
      cores = cores,
      mcmc_seed = mcmc_seed,
      max_sampling_iterations = max_sampling_iterations,
      pass_fit = pass_fit
    )
  
  
}

#' @export
sccomp_glm.SingleCellExperiment = function(.data,
                                           formula_composition = ~ 1 ,
                                           formula_variability = ~ 1,
                                           .sample,
                                           .cell_group,
                                           .count = NULL,
                                           
                                           # Secondary arguments
                                           contrasts = NULL,
                                           prior_mean_variable_association = list(intercept = c(5, 2), slope = c(0,  0.6), standard_deviation = c(20, 40)),
                                           check_outliers = TRUE,
                                           bimodal_mean_variability_association = FALSE,
                                           enable_loo = FALSE,
                                           
                                           # Tertiary arguments
                                           cores = detectCores(),
                                           percent_false_positive = 5,
                                           approximate_posterior_inference = "none",
                                           test_composition_above_logit_fold_change = 0.1, .sample_cell_group_pairs_to_exclude = NULL,
                                           verbose = FALSE,
                                           noise_model = "multi_beta_binomial",
                                           exclude_priors = FALSE,
                                           use_data = TRUE,
                                           mcmc_seed = sample(1e5, 1),
                                           max_sampling_iterations = 20000,
                                           pass_fit = TRUE) {
  
  if(!is.null(.count)) stop("sccomp says: .count argument can be used only for data frame input")
  
  
  # Prepare column same enquo
  .sample = enquo(.sample)
  .cell_group = enquo(.cell_group)
  .sample_cell_group_pairs_to_exclude = enquo(.sample_cell_group_pairs_to_exclude)
  
  
  .data %>%
    colData() %>%
    sccomp_glm(
      formula_composition = formula_composition,
      formula_variability = formula_variability,
      
      !!.sample,!!.cell_group,
      check_outliers = check_outliers,
      contrasts = contrasts,
      prior_mean_variable_association = prior_mean_variable_association,
      percent_false_positive = percent_false_positive ,
      approximate_posterior_inference = approximate_posterior_inference,
      test_composition_above_logit_fold_change = test_composition_above_logit_fold_change, .sample_cell_group_pairs_to_exclude = !!.sample_cell_group_pairs_to_exclude,
      verbose = verbose,
      noise_model = noise_model,
      exclude_priors = exclude_priors,
      bimodal_mean_variability_association = bimodal_mean_variability_association,
      enable_loo = enable_loo,
      use_data = use_data,
      cores = cores,
      mcmc_seed = mcmc_seed,
      max_sampling_iterations = max_sampling_iterations,
      pass_fit = pass_fit
    )
  
  
}

#' @export
sccomp_glm.DFrame = function(.data,
                             formula_composition = ~ 1 ,
                             formula_variability = ~ 1,
                             .sample,
                             .cell_group,
                             .count = NULL,
                             
                             # Secondary arguments
                             contrasts = NULL,
                             prior_mean_variable_association = list(intercept = c(5, 2), slope = c(0,  0.6), standard_deviation = c(20, 40)),
                             check_outliers = TRUE,
                             bimodal_mean_variability_association = FALSE,
                             enable_loo = FALSE,
                             
                             # Tertiary arguments
                             cores = detectCores(),
                             percent_false_positive = 5,
                             approximate_posterior_inference = "none",
                             test_composition_above_logit_fold_change = 0.1, .sample_cell_group_pairs_to_exclude = NULL,
                             verbose = FALSE,
                             noise_model = "multi_beta_binomial",
                             exclude_priors = FALSE,
                             use_data = TRUE,
                             mcmc_seed = sample(1e5, 1),
                             max_sampling_iterations = 20000,
                             pass_fit = TRUE) {
  
  if(!is.null(.count)) stop("sccomp says: .count argument can be used only for data frame input")
  
  
  # Prepare column same enquo
  .sample = enquo(.sample)
  .cell_group = enquo(.cell_group)
  .count = enquo(.count)
  .sample_cell_group_pairs_to_exclude = enquo(.sample_cell_group_pairs_to_exclude)
  
  .data %>%
    as.data.frame %>%
    sccomp_glm(
      formula_composition = formula_composition,
      formula_variability = formula_variability,
      !!.sample,!!.cell_group,
      contrasts = contrasts,
      prior_mean_variable_association = prior_mean_variable_association,
      percent_false_positive = percent_false_positive ,
      check_outliers = check_outliers,
      approximate_posterior_inference = approximate_posterior_inference,
      test_composition_above_logit_fold_change = test_composition_above_logit_fold_change, .sample_cell_group_pairs_to_exclude = !!.sample_cell_group_pairs_to_exclude,
      verbose = verbose,
      noise_model = noise_model,
      exclude_priors = exclude_priors,
      bimodal_mean_variability_association = bimodal_mean_variability_association,
      enable_loo = enable_loo,
      use_data = use_data,
      cores = cores,
      mcmc_seed = mcmc_seed,
      max_sampling_iterations = max_sampling_iterations,
      pass_fit = pass_fit
    )
}

#' @importFrom purrr when
#' @export
sccomp_glm.data.frame = function(.data,
                                 formula_composition = ~ 1 ,
                                 formula_variability = ~ 1,
                                 .sample,
                                 .cell_group,
                                 .count = NULL,
                                 
                                 # Secondary arguments
                                 contrasts = NULL,
                                 prior_mean_variable_association = list(intercept = c(5, 2), slope = c(0,  0.6), standard_deviation = c(20, 40)),
                                 check_outliers = TRUE,
                                 bimodal_mean_variability_association = FALSE,
                                 enable_loo = FALSE,
                                 
                                 # Tertiary arguments
                                 cores = detectCores(),
                                 percent_false_positive = 5,
                                 approximate_posterior_inference = "none",
                                 test_composition_above_logit_fold_change = 0.1, .sample_cell_group_pairs_to_exclude = NULL,
                                 verbose = FALSE,
                                 noise_model = "multi_beta_binomial",
                                 exclude_priors = FALSE,
                                 use_data = TRUE,
                                 mcmc_seed = sample(1e5, 1),
                                 max_sampling_iterations = 20000,
                                 pass_fit = TRUE) {
  
  # Prepare column same enquo
  .sample = enquo(.sample)
  .cell_group = enquo(.cell_group)
  .count = enquo(.count)
  .sample_cell_group_pairs_to_exclude = enquo(.sample_cell_group_pairs_to_exclude)
  
  # DEPRECATE
  deprecate_warn(
    when="1.7.1",
    what="sccomp_glm()",
    details="sccomp says: sccomp_glm() is soft-deprecated. Please use the new modular framework instead, which includes sccomp_estimate(), sccomp_test(), sccomp_remove_outliers(), among other functions."
  )

  # DEPRECATION OF approximate_posterior_inference
  if (is_present(approximate_posterior_inference) & !is.null(approximate_posterior_inference)) {
    deprecate_warn("1.7.7", "sccomp::sccomp_estimate(approximate_posterior_inference = )", details = "The argument approximate_posterior_inference is now deprecated please use inference_method By default variational_inference value is inferred from approximate_posterior_inference.")
    
    inference_method = ifelse(approximate_posterior_inference == "all", "variational","hmc")
  }
  
  if(quo_is_null(.count) )
  result =    sccomp_glm_data_frame_raw(
        .data,
        formula_composition = formula_composition,
        formula_variability = formula_variability,
        
        !!.sample,
        !!.cell_group,
        contrasts = contrasts,
        prior_overdispersion_mean_association = prior_mean_variable_association,
        percent_false_positive = percent_false_positive ,
        check_outliers = check_outliers,
        inference_method = inference_method,
        test_composition_above_logit_fold_change = test_composition_above_logit_fold_change, 
        .sample_cell_group_pairs_to_exclude = !!.sample_cell_group_pairs_to_exclude,
        verbose = verbose,
        # my_glm_model = my_glm_model,
        exclude_priors = exclude_priors,
        bimodal_mean_variability_association = bimodal_mean_variability_association,
        enable_loo = enable_loo,
        use_data = use_data,
        cores = cores,
        mcmc_seed = mcmc_seed,
        max_sampling_iterations = max_sampling_iterations,
        pass_fit = pass_fit
      )
  
  # If the dataframe does includes counts
  else
    result =  sccomp_glm_data_frame_counts(
        .data,
        formula_composition = formula_composition,
        formula_variability = formula_variability,
        
        !!.sample,
        !!.cell_group,
        !!.count,
        contrasts = contrasts,
        prior_overdispersion_mean_association = prior_mean_variable_association,
        percent_false_positive = percent_false_positive ,
        check_outliers = check_outliers,
        variational_inference = approximate_posterior_inference == "all",
        test_composition_above_logit_fold_change = test_composition_above_logit_fold_change, .sample_cell_group_pairs_to_exclude = !!.sample_cell_group_pairs_to_exclude,
        verbose = verbose,
        # my_glm_model = my_glm_model,
        exclude_priors = exclude_priors,
        bimodal_mean_variability_association = bimodal_mean_variability_association,
        enable_loo = enable_loo,
        use_data = use_data,
        cores = cores,
        mcmc_seed = mcmc_seed,
        max_sampling_iterations = max_sampling_iterations,
        pass_fit = pass_fit
      )
    
  result = result |> 
    
    # Track input parameters
    add_attr(noise_model, "noise_model") |> 
    add_attr(.sample, ".sample") |> 
    add_attr(.cell_group, ".cell_group") 
  
  
  # Remove outliers
  if(check_outliers)
    result = result |> 
    sccomp_remove_outliers(
      percent_false_positive = percent_false_positive,
      cores = cores, 
      approximate_posterior_inference = approximate_posterior_inference,
      verbose = verbose,
      mcmc_seed = mcmc_seed,
      max_sampling_iterations = max_sampling_iterations,
      enable_loo = enable_loo
    )
  
  
  result |> 
    sccomp_test(
      contrasts = contrasts, 
      percent_false_positive = percent_false_positive, 
      test_composition_above_logit_fold_change = test_composition_above_logit_fold_change, 
      pass_fit = pass_fit
    )
}



#' test_contrasts
#'
#' @description This function test ocntrasts from a sccomp result.
#'
#'
#' @param .data A tibble. The result of sccomp_glm.
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
#'   sccomp_glm(
#'   counts_obj ,
#'    ~ 0 + type, ~1,  sample, cell_group, count,
#'     check_outliers = FALSE,
#'     cores = 1
#'   ) |>
#'
#'   test_contrasts("typecancer - typebenign")
#'
test_contrasts <- function(.data,
                           contrasts = NULL,
                           percent_false_positive = 5,
                           test_composition_above_logit_fold_change = 0.1,
                           pass_fit = TRUE) {
  
  # DEPRECATE
  deprecate_warn(
    when="1.7.1",
    what="test_contrasts()",
    details="sccomp says: test_contrasts() is soft-deprecated. Please use sccomp_test()."
  )
  
  UseMethod("sccomp_test", .data)
}

#' plot_summary
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
#' @param .data A tibble including a cell_group name column | sample name column | read counts column | factor columns | Pvalue column | a significance column
#' @param significance_threshold A real. FDR threshold for labelling significant cell-groups.
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
#'     approximate_posterior_inference = "all",
#'     cores = 1
#'   )
#'
#' # estimate |> plot_summary()
#'
plot_summary <- function(.data, significance_threshold = 0.025) {
  
  # DEPRECATE
  deprecate_warn(
    when="1.7.1",
    what="plot_summary()",
    details="sccomp says: plot_summary() is soft-deprecated. Please use sccomp_test()."
  )
  
  UseMethod("plot", .data)
 
  
}
