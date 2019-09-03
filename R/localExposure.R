#' This function computes the local exposure index of group m to group n.
#' in situations where m=n, then the result is the isolation index.
#'
#' @param self list with population intensity calculated (result from popIntensity)
#' @return data frame with local exposure indexes, size n*n (n=number of groups) x number of localities

localExposure <- function(self){


  m <- self$n_group
  j <- self$n_location
  local_ids <- self$locality$id
  pop <- self$pop %>%
    dplyr::filter(.data$id %in% local_ids) %>%
    dplyr::select(-id)
  pop_sum <- self$pop_sum %>%
    dplyr::filter(.data$id %in% local_ids) %>%
    dplyr::select(-id)

  local_expo <- sweep(pop, 2, apply(pop, 2, sum, na.rm=T), '/')
  local_expo <- dplyr::bind_cols(id = local_ids, local_expo) %>%
    tidyr::gather('group_exp', 'pop', (1:m)+1)

  locality <- self$locality[,-1]
  if (nrow(self$locality) == 0){
    locality_rate <- pop / pop_sum
  } else{
    locality_rate <- locality / apply(locality,1,sum, na.rm=T)
  }
  locality_rate <- dplyr::bind_cols(id = local_ids, locality_rate)

  exposure_rs <- local_expo %>%
    dplyr::left_join(locality_rate, by='id') %>%
    tidyr::gather('group_rt', 'rate', (4:(3+m))) %>%
    dplyr::mutate(exposure = .data$pop*.data$rate,
                  exposure = ifelse(is.nan(.data$exposure), 0, .data$exposure),
                  m_n = paste0(.data$group_exp, '_', .data$group_rt)) %>%
    dplyr::select(.data$id, .data$exposure, .data$m_n) %>%
    tidyr::spread(.data$m_n, .data$exposure)

  exposure_rs <- exposure_rs[,-1]

  return(exposure_rs)
}
