#' This function computes the local entropy index H for all localities.
#'
#' @param self list with population intensity calculated (result from popIntensity)
#' @return array with scores for n groups (number of groups)

localIndexH <- function(self){

  local_ids <- self$locality$id
  pop <- self$pop %>%
    dplyr::filter(.data$id %in% local_ids) %>%
    dplyr::select(-id)
  pop_sum <- self$pop_sum %>%
    dplyr::filter(.data$id %in% local_ids) %>%
    dplyr::select(-id)

  local_entropy <- localEntropy(self)
  global_entropy <- globalEntropy(self)

  et <- global_entropy * sum(pop_sum, na.rm=T)
  eei <- global_entropy - local_entropy
  h_local <- pop_sum * eei / et

  return(h_local)
}
