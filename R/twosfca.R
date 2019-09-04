#' Processes the jobs 2SFCA accessibility
#' 
#' @inheritParams cumop
#' @param prestep logical parameter. Defaults to TRUE to process the prestep 
#' that divides the population by the number of counts in the catchment area of the first step
#' @return data frame with origin, acc and threshold

twosfca <- function(cutoff, oportunities, matrix, prestep = TRUE){
  
  matrix <- matrix %>% 
    dplyr::filter(travel_time<=cutoff & travel_time >=0)
  
  if(prestep){
    # pre-step 
    # counts how many times the pop is counted at destination
    message('pre-step')
    freq <- matrix %>% 
      dplyr::group_by(destination) %>% 
      dplyr::count()
    # divides population by frequency
    pop <- oportunities %>% 
      dplyr::left_join(freq, by=c('id'='destination')) %>% 
      dplyr::mutate(pop = pop/n) %>% 
      dplyr::select(id, pop) 
    
    rm(freq)
    gc()
  }

  job <- oportunities %>% 
    dplyr::select(id, job) 
  
  message('1st step')
  # 1st step - rate of jobs (origin) by population (destination)
  pop_rate <- matrix %>% 
    dplyr::left_join(pop, by=c('destination'='id')) %>% 
    dplyr::group_by(origin) %>% 
    dplyr::mutate(rate = 1/sum(pop, na.rm=T)) %>% 
    dplyr::ungroup() %>% 
    dplyr::distinct(origin, rate)
  
  first <- matrix %>% 
    dplyr::left_join(job, by=c('origin'='id')) %>% 
    dplyr::left_join(pop_rate, by=c('origin')) %>% 
    dplyr::group_by(origin) %>% 
    dplyr::mutate(rate = job*rate) %>% 
    dplyr::ungroup() %>% 
    dplyr::rename(id = origin) %>% 
    dplyr::distinct(id, rate) 
  
  rm(pop_rate)
  gc()
  
  message('2nd step')
  # 2nd step - sum at the origin the rates at the destination
  matrix %>% 
    dplyr::left_join(first, by = c('destination' = 'id')) %>% 
    dplyr::group_by(origin) %>% 
    dplyr::summarise(acc = sum(rate, na.rm=T)) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(threshold = paste0(cutoff,'min'))
}