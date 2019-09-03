#' Compute local dissimilarity for all groups.
#'
#' @param self list with population intensity calculated (result from popIntensity)
#' @return array with dissimilarity for all groups, size of number of localities

localDissimilarity <- function(self){

  local_ids <- self$locality$id
  pop <- self$pop %>%
    dplyr::filter(.data$id %in% local_ids) %>%
    dplyr::select(-id)
  pop_sum <- self$pop_sum %>%
    dplyr::filter(.data$id %in% local_ids) %>%
    dplyr::select(-id)

  if(nrow(self$locality) == 0){
    lj <- as.numeric(pop_sum)
    tjm <- pop / lj
  } else{
    locality <- self$locality[,-1]
    lj <- apply(locality, 1, sum)
    tjm <- locality / lj
  }

  pop_total <- sum(pop, na.rm=T)
  tm <- apply(pop, 2, sum) / pop_total

  index_i <- sum(tm * (1 - tm), na.rm = T)

  local_diss <- apply(abs(sweep(tjm, 2, tm, '-')) * t(pop_sum) / (2 * pop_total * index_i), 1, sum, na.rm=T)
  local_diss <- matrix(local_diss, ncol=1)
  local_diss[is.nan(local_diss)] <- 0

  return (local_diss)
}
