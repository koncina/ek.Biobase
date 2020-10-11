NULL

#' @export
format_publication_date <- function(x) {
  assertr::verify(x, str_detect(status, "^Public on")) %>%
    rename(publication_date = status) %>%
    mutate_at("publication_date", str_remove, "^Public on ") %>%
    mutate_at(vars(ends_with("date")), parse_date, format = "%b %d %Y")
}

#' Tidy the GEO characteristics columns
#'
#' @details
#' Characterstics columns might contain gathered and united key value pairs in form of "key: value".
#' This function will separate the key value pairs and spread the variables into separate columns.
#'
#' @importFrom tidyselect vars_select
#' @importFrom rlang enquos
#' @importFrom assertr verify
#' @importFrom tidyr separate_rows separate unite gather
#'
#' @param data A tibble.
#' @param ... columns containing the key value pairs.
#' @export
tidy_key_values <- function(data, ...) {

  data <- rowid_to_column(data, ".rowid")
  columns <- c(".rowid", vars_select(colnames(data), !!! enquos(...)))
  x <- select(data, columns)
  x <- pivot_longer(x, names_to = "characteristic", values_to = "key_value", c(-.rowid))
  x <- separate_rows(x, key_value, sep = ";")
  x <- separate(x, "key_value", c("key", "value"), ":[ ]?", extra = "merge", fill = "right")
  x <- mutate_at(x, "key", str_trim, "both")
  x <- unite(x, "key", characteristic, key, sep = "_")
  x <- pivot_wider(x, names_from = "key", values_from = "value")
  x <- verify(x, n_distinct(.rowid) == nrow(x))

  data <- left_join(data, x, by = ".rowid")
  select(data, -one_of(columns))
}

#' Extract Metadata from ExpressionSet with tidy characteristics columns
#'
#' @importFrom dplyr mutate_if
#' @importFrom janitor clean_names
#' @importFrom tibble as_tibble
#' @importFrom tidyr starts_with
#'
#' @param x A GEO ExpressionSet
#' @export
get_eset_metadata <- function(x) {
  x <- as_tibble(Biobase::pData(x))
  x <- tidy_key_values(x, starts_with("characteristics_"))
  x <- janitor::clean_names(x)
  x <- mutate_if(x, is.factor, as.character)
  mutate_if(x, is.character, parse_character, na = c(NA, "na", "NA", ""))
}
