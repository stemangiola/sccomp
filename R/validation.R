check_columns_exist = function(.data, columns){

  if((!columns %in% (.data %>% colnames)) %>% any)
    stop(
      sprintf(
        "The columns %s are not present in your tibble",
        paste(columns[(!columns %in% (.data %>% colnames))], collapse=" ")
      )
    )
}

check_if_count_integer = function(.data, .count){
  .count = enquo(.count)

  # Check if the counts column is an integer
  if (!.data %>% pull(!!.count) %>% is("integer"))
    stop(
      sprintf(
        "The column %s must be of class integer. You can do as mutate(`%s` = as.integer(`%s`))",
        quo_name(.count),
        quo_name(.count),
        quo_name(.count)
      )
    )
}

#' Check if NA
#'
#' @importFrom tidyr drop_na
#' @importFrom dplyr enquo
#'
#' @param .data A tibble including a gene name column | sample name column | read counts column | factors column
#' @param columns Columns to check
#'
#' @keywords internal
#' @noRd
check_if_any_NA = function(.data, columns){


  if(
    .data %>%
    drop_na(all_of(columns)) %>%
    nrow %>% st(      .data %>% nrow    )
  )
    stop(sprintf("There are NA values in you tibble for any of the column %s", paste(columns, collapse=", ")))
}

check_if_within_posterior = function(.data, my_df, .do_check, .count){

  # Define the variables as NULL to avoid CRAN NOTES
  .lower <- NULL
  .upper <- NULL
  ppc <- NULL
  
  writeLines(sprintf("executing %s", "check_if_within_posterior"))

  .data %>%
    left_join(my_df, by = c("S", "G")) %>%
    filter((!!.do_check)) %>% # Filter only DE genes
    rowwise() %>%
    mutate(`ppc` = !!.count %>% between(`.lower`, `.upper`)) %>%
    mutate(`is higher than mean` = (!`ppc`) &
             (!!.count > mean)) %>%
    ungroup
}

check_if_columns_right_class = function(.data, .sample, .cell_group){

  .sample = enquo(.sample)
  .cell_group = enquo(.cell_group)

  # Check that sample and cell_type is chr of factor
  if(
    !(
      .data %>% pull(!!.sample) %>% is("character") |
      .data %>% pull(!!.sample) %>% is("factor")
    ) |
    !(
      .data %>% pull(!!.cell_group) %>% is("character") |
      .data %>% pull(!!.cell_group) %>% is("factor")
    )
  ) stop(sprintf("sccomp says: the columns %s and %s must be of class character or factor", quo_name(.sample), quo_name(.cell_group)))


}
