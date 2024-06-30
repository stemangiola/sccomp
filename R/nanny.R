#' @keywords internal
#' @noRd
#'
.subset <- 		function(.data,	 .column)	{
  # Make col names
  .column <- enquo(.column)

  # Check if column present
  if(quo_names(.column) %in% colnames(.data) %>% all %>% `!`)
    stop("nanny says: some of the .column specified do not exist in the input data frame.")

  .data %>%

    # Selecting the right columns
    select(	!!.column,	all_of(get_specific_annotation_columns(.data, !!.column)	)) %>%
    distinct()

}

#' @importFrom dplyr distinct_at
#' @importFrom rlang enquo
#' @importFrom magrittr equals
#' @importFrom dplyr select
#'
#' @keywords internal
#' @noRd
#'
#'
get_specific_annotation_columns <- function(.data, .col){


  # Comply with CRAN NOTES
  . <- NULL

  # Make col names
  .col <- enquo(.col)

  # x-annotation df
  n_x <- .data %>% distinct_at(vars(!!.col)) %>% nrow

  # element wise columns
  .data %>%
    select(-!!.col) %>%
    colnames %>%
    map(
      ~
        .x %>%
        ifelse_pipe(
          .data %>%
            distinct_at(vars(!!.col, .x)) %>%
            nrow %>%
            equals(n_x),
          ~ .x,
          ~ NULL
        )
    ) %>%

    # Drop NULL
    {	(.)[lengths((.)) != 0]	} %>%
    unlist

}

#' Convert array of quosure (e.g. c(col_a, col_b)) into character vector
#'
#' @keywords internal
#'
#' @importFrom rlang quo_name
#' @importFrom rlang quo_squash
#'
#' @keywords internal
#' @noRd
#'
#' @param v A array of quosures (e.g. c(col_a, col_b))
#'
#' @return A character vector
quo_names <- function(v) {

  v <- quo_name(quo_squash(v))
  gsub('^c\\(|`|\\)$', '', v) %>%
    strsplit(', ') %>%
    unlist
}

# From tidyr
strip_names <- function(df, base, names_sep) {
  base <- paste0(base, names_sep)
  names <- names(df)

  has_prefix <- regexpr(base, names, fixed = TRUE) == 1L
  names[has_prefix] <- substr(names[has_prefix], nchar(base) + 1, nchar(names[has_prefix]))

  set_names(df, names)
}

# Set internal
#' @importFrom rlang enquo
#' @importFrom rlang set_names
#' @importFrom rlang :=
#' @importFrom purrr map
#' @importFrom purrr imap
#' @importFrom tidyselect eval_select
#'
#' @keywords internal
#' @noRd
#'
.nest_subset <- 		function(.data, ..., .exclude = NULL, .names_sep = NULL)	{

  # Make col names - from tidyr
  cols <- enquos(...)
  .exclude <- enquo(.exclude)

  # Name of the new data column
  col_name_data  <- names(cols)

  # Column names
  cols <- map(cols, ~ names(eval_select(.x, .data)))
  cols <- map(cols, set_names)
  if (!is.null(.names_sep)) cols <- imap(cols, strip_names, .names_sep)
  asis <- setdiff(names(.data), unlist(cols))

  # Check if column present
  if(asis %in% colnames(.data) %>% all %>% `!`)
    stop("nanny says: some of the .column specified do not exist in the input data frame.")

  # Get my subset columns
  asis_subset <- asis %>%
    c(get_specific_annotation_columns(.data, asis)) %>%

    # Exclude custom columns
    setdiff(quo_names(.exclude))

  # Apply nest on those
  tidyr::nest(.data, !!col_name_data := -c(asis_subset))

}
