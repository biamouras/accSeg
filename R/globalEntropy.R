#' This function computes the global entropy score E (diversity).
#'
#' @param self list with population intensity calculated (result from popIntensity)
#' @return diversity score

globalEntropy <- function(self){

  local_ids <- self$locality$id
  pop <- self$pop %>%
    dplyr::filter(.data$id %in% local_ids) %>%
    dplyr::select(-id)
  pop_sum <- self$pop_sum %>%
    dplyr::filter(.data$id %in% local_ids) %>%
    dplyr::select(-id)

  pop_total <- sum(pop, na.rm=T)
  prop <- apply(pop, 2, sum, na.rm=T)
  n <- length(prop)
  group_score <- array(dim=n)

  group_score <- sapply(1:n, function(i){
    group <- prop[i]
    group_idx <- group / pop_total * log(1 / (group / pop_total))
    group_idx
  })

  entropy <- sum(group_score, na.rm=T)

  return(entropy)
}
