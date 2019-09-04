#' Processes the potential accessibility with Gaussian function
#' 
#' @inheritParams cumop
#' @return data frame with origin, acc and threshold

potGaussian <- function(cutoff, oportunities, matrix){
  matrix %>% 
    dplyr::mutate(weight = exp(-0.5 * (travel_time/cutoff)^2)) %>% 
    dplyr::left_join(oportunities, by=c('destination' = 'id')) %>% 
    dplyr::group_by(origin) %>% 
    dplyr::summarise(acc = sum(job*weight, na.rm=T)) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(threshold = paste0(cutoff,'min'))
}