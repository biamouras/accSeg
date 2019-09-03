#' This function computes the local entropy score for a unit area Ei (diversity).
#'
#' @param self list with population intensity calculated (result from popIntensity)
#' @return array with local indices of entropy

localEntropy <- function(self){

  local_ids <- self$locality$id
  pop <- self$pop %>%
    dplyr::filter(.data$id %in% local_ids) %>%
    dplyr::select(-id)
  pop_sum <- self$pop_sum %>%
    dplyr::filter(.data$id %in% local_ids) %>%
    dplyr::select(-id)

  if (nrow(self$locality) == 0){
    proportion <- pop / as.numeric(pop_sum)
  } else{
    polygon_sum <- apply(self$locality[,-1], 1, sum, na.rm=T)
    proportion <- self$locality[,-1] / polygon_sum
  }
  entropy <- proportion * log(1 / proportion)
  entropy$id <- local_ids

  entropy <- entropy %>%
    tidyr::gather('group', 'ent', 1:self$n_group) %>%
    dplyr::mutate(ent = ifelse(is.na(.data$ent) | is.infinite(.data$ent), 0, .data$ent)) %>%
    tidyr::spread(.data$group, .data$ent) %>%
    dplyr::select(-id)

  entropy <- apply(entropy, 1, sum, na.rm=T)
  entropy <- matrix(entropy, ncol=1)

  return(entropy)
}
