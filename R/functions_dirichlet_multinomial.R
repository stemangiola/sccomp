#' dirichlet_multinomial_glm main
#'
#' @description This function runs the data modelling and statistical test for the hypothesis that a cell_type includes outlier biological replicate.
#'
#' @importFrom tibble as_tibble
#' @import dplyr
#' @importFrom tidyr spread
#' @import tidybayes
#' @importFrom magrittr %$%
#' @importFrom magrittr divide_by
#' @importFrom magrittr multiply_by
#' @importFrom purrr map2
#' @importFrom purrr map_int
#' @importFrom tidybulk scale_abundance
#' @importFrom benchmarkme get_ram
#' @importFrom magrittr multiply_by
#' @importFrom magrittr equals
#' @importFrom purrr map
#' @importFrom tibble rowid_to_column
#' @importFrom furrr future_map
#' @importFrom purrr map_lgl
#'
#' @param .data A tibble including a cell_type name column | sample name column | read counts column | covariate columns | Pvaue column | a significance column
#' @param formula A formula. The sample formula used to perform the differential cell_type abundance analysis
#' @param .sample A column name as symbol. The sample identifier
#' @param .cell_type A column name as symbol. The cell_type identifier
#' @param .count A column name as symbol. The cell_type abundance (read count)
#' @param approximate_posterior_inference A boolean. Whether the inference of the joint posterior distribution should be approximated with variational Bayes. It confers execution time advantage.
#' @param cores An integer. How many cored to be used with parallel calculations.
#' @param seed An integer. Used for development and testing purposes
#'
#' @return A nested tibble `tbl` with cell_type-wise information: `sample wise data` | plot | `ppc samples failed` | `exposure deleterious outliers`
#'
#' @export
#'
dirichlet_multinomial_glm = function(.data,
                                     formula = ~ 1,
                                     .sample,
                                     .cell_type,
                                     .count,
                                     check_outliers = FALSE,
                                     approximate_posterior_inference = T,
                                     cores = detect_cores(), # For development purpose,
                                     seed = sample(1:99999, size = 1)
) {
  # Prepare column same enquo
  .sample = enquo(.sample)
  .cell_type = enquo(.cell_type)
  .count = enquo(.count)

  data_for_model =
    data_to_spread (.data, formula, !!.sample, !!.cell_type, !!.count) %>%
    data_spread_to_model_input(formula, !!.sample, !!.cell_type, !!.count)

  # Produce data list
  if(!check_outliers){

    data_for_model %>%
      # Run the first discovery phase with permissive false discovery rate
      fit_model(stanmodels$glm_dirichlet_multinomial, censoring_iteration, chains= 4) %>%
      parse_fit(data_for_model, ., censoring_iteration = 1, chains) %>%
      beta_to_CI(censoring_iteration = 1 ) %>%

      # Join filtered
      mutate(
        significant =
          !!as.symbol(sprintf(".lower_%s", colnames(data_for_model$X)[2])) *
          !!as.symbol(sprintf(".upper_%s", colnames(data_for_model$X)[2])) > 0
      ) %>%

      # Clesn
      select(-M) %>%
      mutate(!!.cell_type := data_for_model$y %>% colnames()) %>%
      select(!!.cell_type, everything())

  }

  else{

    .data_1 =
      .data %>%
      fit_model_and_parse_out_no_missing_data(data_for_model, formula, !!.sample, !!.cell_type, !!.count, iteration = 1, chains = 4)

    .data_2 =
      .data_1 %>%
      select(-contains("posterior")) %>%
      fit_model_and_parse_out_missing_data(formula, !!.sample, !!.cell_type, !!.count, iteration = 2)

    .data_2 %>%

      # Join filtered
      mutate(
        significant =
          !!as.symbol(sprintf(".lower_%s", colnames(data_for_model$X)[2])) *
          !!as.symbol(sprintf(".upper_%s", colnames(data_for_model$X)[2])) > 0
      ) %>%

      # #Join unfiltered
      # mutate(significant_pre_filtering = map_lgl(
      #   beta_quantiles_1,
      #   ~ .x$`2.5%` * .x$`97.5%` > 0
      # )) %>%

      # Define outlier
      rename(outlier = outlier_2 ) %>%

      # Clean
      select(-N, -M, -contains("posterior"))

  }


}

#' @importFrom purrr map2_lgl
fit_model_and_parse_out_no_missing_data = function(.data, data_for_model, formula, .sample, .cell_type, .count, iteration = 1, chains){


  # Prepare column same enquo
  .sample = enquo(.sample)
  .cell_type = enquo(.cell_type)
  .count = enquo(.count)

  # Run the first discovery phase with permissive false discovery rate
  fit_and_generated  = fit_and_generate_quantities(data_for_model, stanmodels$glm_dirichlet_multinomial, iteration, chains= 4, output_samples = 5000)

  # Integrate
  .data %>%
    left_join(tibble(.sample = rownames(data_for_model$y), N = 1:nrow(data_for_model$y)), by = ".sample" %>% setNames(quo_name(.sample))) %>%
    left_join(tibble(.cell_type = colnames(data_for_model$y), M = 1:ncol(data_for_model$y)), by = ".cell_type" %>% setNames(quo_name(.cell_type))) %>%

    # Add covariate from design
    left_join(fit_and_generated, by = c("M", "N")) %>%

    # Add theoretical data quantiles
    mutate(!!as.symbol(sprintf("generated_data_quantiles_%s", iteration)) := map(
      !!as.symbol(sprintf("generated_data_posterior_%s", iteration)),
      ~ quantile(
        .x$generated_quantity,
        probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)
      ) %>%
        enframe() %>%
        spread(name, value)
    )) %>%

    # # Attach beta
    # mutate(!!as.symbol(sprintf("beta_quantiles_%s", iteration)) := map(
    #   !!as.symbol(sprintf("beta_posterior_%s", iteration)),
    #   ~ quantile(
    #     .x$.value,
    #     probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)
    #   ) %>%
    #     enframe() %>%
    #     spread(name, value)
    # ))   %>%

  #     # Join slope
  # left_join(  beta_posterior %>% select(M, !!as.symbol(sprintf("slope_%s", iteration)) := `50%`)  ) %>%

  mutate(!!as.symbol(sprintf("outlier_%s", iteration)) := map2_lgl(
    !!.count, !!as.symbol(sprintf("generated_data_quantiles_%s", iteration)),
    ~ .x < .y$`5%` | .x > .y$`95%`)
  ) %>%

    # Add precision as attribute
    add_attr( attr(fit_and_generated, "precision"), "precision" )



}

fit_model_and_parse_out_missing_data = function(.data, formula, .sample, .cell_type, .count, iteration){


  # Prepare column same enquo
  .sample = enquo(.sample)
  .cell_type = enquo(.cell_type)
  .count = enquo(.count)

  # Produce data list
  data_for_model =
    .data %>%
    data_to_spread (formula, !!.sample, !!.cell_type, !!.count) %>%
    data_spread_to_model_input(formula, !!.sample, !!.cell_type, !!.count)

  # .count = enquo(.count)
  #
  # .data_wide =
  #   .my_data %>%
  #   select(N, M, !!.count, parse_formula(formula)) %>%
  #   distinct() %>%
  #   spread(M, !!.count)

  # .data_wide_no_covariates = .data_wide %>% select(-parse_formula(formula))

  to_exclude =
    .data %>%
    filter(!!as.symbol(sprintf("outlier_%s", iteration - 1))  ) %>%
    distinct(N, M)

  to_include =
    .data %>%
    filter(!!as.symbol(sprintf("outlier_%s", iteration - 1)) %>% `!`  ) %>%
    distinct(N, M)

  # To mix with imputed data
  .data_parsed_inliers =
    .data %>%
    anti_join(to_exclude,by = c("N", "M")) %>%
    select(.value = count, N, M)

  # Dirichlet with missing data
  fit_imputation =
    data_for_model %>%
    do_inference_imputation(
      formula,
      approximate_posterior_inference,
      approximate_posterior_analysis,
      C,
      X,
      cores,
      additional_parameters_to_save,
      pass_fit = pass_fit,
      to_include = to_include,
      tol_rel_obj = tol_rel_obj,
      #truncation_compensation = 0.7352941, # Taken by approximation study
      seed = seed,
      precision = .data %>% attr("precision")
    )

  beta_posterior_corrected =
    fit_imputation %>%
    draws_to_tibble_x_y("counts", "N", "M") %>%
    rename(.draw_imputation = .draw) %>%
    nest(data = -c(.chain ,.iteration, .draw_imputation ,.variable)) %>%
    sample_n(100) %>%
    mutate(fit = future_map(
      data,
      ~ {
        data_for_model$y =
          .x %>%
          anti_join(to_include,by = c("N", "M")) %>%
          bind_rows(.data_parsed_inliers) %>%
          spread(M, .value) %>%
          as_matrix(rownames = "N")

        # Run model
        fit_and_generate_quantities(data_for_model, stanmodels$glm_dirichlet_multinomial, iteration, chains=1, output_samples = 200)
      }
    )) %>%

    # Add precision
    mutate(precision = map(
      fit,
      ~ attr(.x, "precision")
    ))

  beta_posterior_corrected %>%

    select(fit) %>%
    unnest(fit) %>%
    nanny::nest_subset(data = -c(N, M)) %>%

    # Merge posterior data
    mutate(!!as.symbol(sprintf("generated_data_posterior_%s", iteration)) := map(
      data,
      ~ .x %>%
        select( !!as.symbol(sprintf("generated_data_posterior_%s", iteration))) %>%
        unnest( !!as.symbol(sprintf("generated_data_posterior_%s", iteration)))
    )) %>%

    # Add theoretical data quantiles
    mutate(!!as.symbol(sprintf("generated_data_quantiles_%s", iteration)) := map(
      !!as.symbol(sprintf("generated_data_posterior_%s", iteration)),
      ~ quantile(
        .x$generated_quantity,
        probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)
      ) %>%
        enframe() %>%
        spread(name, value)
    )) %>%

    # Merge posterior data
    mutate(!!as.symbol(sprintf("beta_posterior_%s", iteration)) := map(
      data,
      ~ .x %>%
        select( !!as.symbol(sprintf("beta_posterior_%s", iteration))) %>%
        unnest( !!as.symbol(sprintf("beta_posterior_%s", iteration)))
    )) %>%

    left_join(
      (.) %>% select(N, M, beta_posterior_2) %>% beta_to_CI(iteration)
    ) %>%

    select(-data) %>%

    right_join( .data) %>%

    mutate(!!as.symbol(sprintf("outlier_%s", iteration)) := map2_lgl(
      !!.count, !!as.symbol(sprintf("generated_data_quantiles_%s", iteration)),
      ~ .x < .y$`5%` | .x > .y$`95%`)
    ) %>%

    # Add precision as attribute
    add_attr(
      beta_posterior_corrected %>%
        select(precision) %>%
        unnest(precision) %>%
        as.matrix(),
      "precision" )

}