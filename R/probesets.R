#' @export
select_probeset <- function(data) {

  if (nrow(data) == 0) return(tibble(gene_symbol = character(0),
                                     probeset_id = character(0),
                                     selected_probeset = logical(0)))

  group_by(data, gene_symbol, probeset_id) %>%
    summarise(var = var(expression), .groups = "drop_last") %>%
    mutate(is_unspecific = str_detect(probeset_id, "[xs]_at$"),
           only_unspecific = all(is_unspecific),
           is_candidate = !xor(is_unspecific, only_unspecific),
           selected_probeset = var == max(var[is_candidate])) %>%
    select(-is_unspecific, -is_candidate, -only_unspecific)
}

#' @export
filter_probeset <- function(data, selected_probesets) {

  if (missing(selected_probesets)) {
    selected_probesets <- select_probeset(data)
  }

  semi_join(data,
            filter(selected_probesets, selected_probeset),
            by = c("probeset_id", "gene_symbol"))
}
