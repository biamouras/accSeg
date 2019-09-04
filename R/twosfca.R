#' Processes the jobs 2SFCA accessibility
#' 
#' @inheritParams cumop
#' @param prestep logical parameter. Defaults to TRUE to process the prestep 
#' that divides the population by the number of counts in the catchment area of the first step
#' @return data frame with origin, acc and threshold

twosfca <- function(cutoff, oportunities, matrix, prestep = TRUE){
  
  matrix <- matrix %>% 
    dplyr::filter(.data$travel_time<=cutoff & .data$travel_time >=0)
  
  if(prestep){
    # pre-step 
    # counts how many times the pop is counted at destination
    message('pre-step')
    freq <- matrix %>% 
      dplyr::group_by(.data$destination) %>% 
      dplyr::count()
    # divides population by frequency
    pop <- oportunities %>% 
      dplyr::left_join(freq, by=c('id'='destination')) %>% 
      dplyr::mutate(pop = .data$pop/n) %>% 
      dplyr::select(.data$id, .data$pop) 
    
    rm(freq)
    gc()
  }

  job <- oportunities %>% 
    dplyr::select(.data$id, .data$job) 
  
  message('1st step')
  # 1st step - rate of jobs (origin) by population (destination)
  pop_rate <- matrix %>% 
    dplyr::left_join(pop, by=c('destination'='id')) %>% 
    dplyr::group_by(.data$origin) %>% 
    dplyr::mutate(rate = 1/sum(.data$pop, na.rm=T)) %>% 
    dplyr::ungroup() %>% 
    dplyr::distinct(.data$origin, .data$rate)
  
  first <- matrix %>% 
    dplyr::left_join(job, by=c('origin'='id')) %>% 
    dplyr::left_join(pop_rate, by=c('origin')) %>% 
    dplyr::group_by(.data$origin) %>% 
    dplyr::mutate(rate = .data$job*.data$rate) %>% 
    dplyr::ungroup() %>% 
    dplyr::rename(id = .data$origin) %>% 
    dplyr::distinct(.data$id, .data$rate) 
  
  rm(pop_rate)
  gc()
  
  message('2nd step')
  # 2nd step - sum at the origin the rates at the destination
  matrix %>% 
    dplyr::left_join(first, by = c('destination' = 'id')) %>% 
    dplyr::group_by(.data$origin) %>% 
    dplyr::summarise(acc = sum(.data$rate, na.rm=T)) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(threshold = paste0(cutoff,'min'))
}