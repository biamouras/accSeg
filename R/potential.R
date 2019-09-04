#' Processes the potential (or gravity) job accessibility.
#' 
#' @param betas numeric 'willing' decay factor
#' @param oportunities data frame with at least id and job fields
#' @param matrix data frame with origin, destination and travel_time structure
#' @return data frame with origin and acc
#' @description Uses the exponential function to weight the oportunities.

potential <- function(betas, oportunities, matrix){
  matrix %>% 
    dplyr::left_join(oportunities, by = c('destination' = 'id')) %>% 
    dplyr::mutate(decay = exp(.data$travel_time*betas),
                  job = .data$job*.data$decay) %>% 
    dplyr::group_by(.data$origin) %>% 
    dplyr::summarise(acc = sum(.data$job, na.rm=T)) %>% 
    dplyr::ungroup()
}