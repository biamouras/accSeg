#' Processes the job cumulative accessibility
#' 
#' @param cutoff numeric threshold of cumulative accessibility
#' @param oportunities data frame with at least id and job fields
#' @param matrix data frame with origin, destination and travel_time structure
#' @return data frame with origin, acc and threshold

cumop <- function(cutoff, oportunities, matrix){
  matrix %>% 
    dplyr::filter(.data$travel_time<=cutoff & .data$travel_time >=0) %>% 
    dplyr::left_join(oportunities, by=c('destination' = 'id')) %>% 
    dplyr::group_by(.data$origin) %>% 
    dplyr::summarise(acc = sum(.data$job, na.rm=T)) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(threshold = paste0(cutoff,'min'))
}