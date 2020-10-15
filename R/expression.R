NULL

globalVariables(c("probeset_id", "gene_symbol", ".id"))

#' @importFrom tibble as_tibble
#' @importFrom dplyr select left_join
join_metadata <- function(x, eset, cols = everything()) {
  metadata <- pData(eset)
  metadata <- as_tibble(metadata, rownames = ".id")
  metadata <- select(metadata, .id, {{cols}})
  left_join(x, metadata, by = c("sample_id" = ".id"))
}

#' Extract genes of interest out of an ExpressionSet object
#'
#' A helper function to get the expression data of genes of interest out of the downloaded GEO data.
#'
#' @importFrom janitor make_clean_names
#' @importFrom dplyr select everything mutate
#' @importFrom stringr str_split
#' @importFrom tibble as_tibble
#' @importFrom tidyr unnest
#' @import Biobase
#'
#' @param eset An ExpressionSet objects.
#' @param goi A character vector containing the gene symbols of the genes of interest
#' @param metadata metadata columns to join (tidyselect). If not set will return only the expression values.
#'
#' @export
get_eset_goi <- function(eset, goi, metadata) {
  x <- as_tibble(fData(eset),
                 rownames = "probeset_id",
                 .name_repair = make_clean_names)
  x <- select(x, "probeset_id", "gene_symbol")
  x <- mutate(x, across(gene_symbol, str_split, " /// "))
  x <- unnest(x, gene_symbol)
  x <- filter(x, gene_symbol %in% goi)
  x <- left_join(x, as_tibble(exprs(eset), rownames = "probeset_id"), by = "probeset_id")
  x <- pivot_longer(x, names_to = "sample_id", values_to = "expression", -c(probeset_id, gene_symbol))

  if (!missing(metadata)) {
    x <- join_metadata(x, eset, {{metadata}})
  }

  x
}
