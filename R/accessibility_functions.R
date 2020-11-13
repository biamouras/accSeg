#' Processes the cumulative accessibility.
#' 
#' @param cutoff numeric threshold of cumulative accessibility. It must be in the same unit of the matrix.
#' @param oportunities data frame with id and oportunity structure.
#' @param matrix data frame with origin, destination and travel_time/distance structure.
#' @return data frame with origin, acc and threshold.
#' @description The job cumulative accessibility adds up the number of oportunities inside a given travel time or distance (threshold). 
#' @details As discussed by Handy and Niemeier (1997), the cumulative accessibility is easier to understand than the gravity one, 
#'    giving "some sense of the range of choice avaiable" (Handy & Niemeier, 1997, p. 1177).
#' @author Beatriz Moura dos Santos
#' @references Handy & Niemeier (1997). Measuring accessibility: An exploration of issues and alternatives. \emph{Environment and Planning A: Economy and Space}, 29(7), 1175-1194.

cumop <- function(cutoff, oportunities, matrix){
  # renaming fields
  names(oportunities) <- c('id', 'oportunity')
  names(matrix) <- c('origin', 'destination', 'travel_time')
  # processing cumulative accessibility
  matrix %>% 
    dplyr::filter(.data$travel_time<=cutoff & .data$travel_time >=0) %>% 
    dplyr::left_join(oportunities, by=c('destination' = 'id')) %>% 
    dplyr::group_by(.data$origin) %>% 
    dplyr::summarise(acc = sum(.data$oportunity, na.rm=T)) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(threshold = cutoff)
}

#' Processes the potential (or gravity) accessibility.
#' 
#' @param oportunities data frame with id and oportunity structure.
#' @param matrix data frame with origin, destination and travel_time/distance structure.
#' @param type string to choose between \emph{exponential} and \emph{gaussian} function to decay the oportunities. Defaults to exponential.
#' @param betas numeric 'willing' decay factor. It must be defined when type "exponential" is chosen. Defaults to NULL.
#' @param cutoff numeric threshold. It must be in the same unit of the matrix and it must be defined when type "gaussian" is chosen. Defaults to NULL.
#' @return data frame with origin and acc. If \emph{gaussian} type is chosen, adds the threshold column.
#' @description Uses the negative exponential (\code{exponential type}) or the gaussian (\code{gaussian type}) functions to weight the oportunities.
#' @details The negative exponential function is one of the most used impedance functions in recent studies of gravity-based accessibility (Handy & Niemeier, 1997).
#'     In order to avoid the rapid decline of the negative exponential function, Ingram (1971) proposed the use of the normal (gaussian) function, as it declines 
#'     smoother than the former.
#' @author Beatriz Moura dos Santos
#' @references Hansen (1959). How accessibility shapes land use. \emph{Journal of the American Institute of planners}, 25(2), 73-76.
#' @references Ingram (1971). The concept of accessibility: a search for an operational form. \emph{Regional studies}, 5(2), 101-107.
#' @references Handy & Niemeier (1997). Measuring accessibility: An exploration of issues and alternatives. \emph{Environment and Planning A: Economy and Space}, 29(7), 1175-1194.

potential <- function(oportunities, matrix, type = 'exponential', betas=NULL, cutoff=NULL){
  # renaming fields
  names(oportunities) <- c('id', 'oportunity')
  names(matrix) <- c('origin', 'destination', 'travel_time')
  # processing potential accessibility
  if(type=='exponential' & !is.null(betas)){
    matrix %>% 
      dplyr::left_join(oportunities, by = c('destination' = 'id')) %>% 
      dplyr::mutate(decay = exp(.data$travel_time*betas),
                    oportunity = .data$oportunity*.data$decay) %>% 
      dplyr::group_by(.data$origin) %>% 
      dplyr::summarise(acc = sum(.data$oportunity, na.rm=T)) %>% 
      dplyr::ungroup()
  } else if(type=='gaussian' & !is.null(cutoff)){
    matrix %>% 
      dplyr::mutate(weight = exp(-0.5 * (.data$travel_time/cutoff)^2)) %>% 
      dplyr::left_join(oportunities, by=c('destination' = 'id')) %>% 
      dplyr::group_by(.data$origin) %>% 
      dplyr::summarise(acc = sum(.data$oportunity*.data$weight, na.rm=T)) %>% 
      dplyr::ungroup() %>% 
      dplyr::mutate(threshold = cutoff)
  } else {
    stop("The required parameter to the chosen type is not correct.")
  }
}

#' Processes the competition accessibility.
#' 
#' @param oportunities data frame with id, oportunity and population structure.
#' @param matrix data frame with origin, destination and travel_time/distance structure.
#' @param type string to choose between \emph{2sfca} and \emph{shen} competition accessibilities. Defaults to 2sfca.
#' @param betas numeric 'willing' decay factor. It must be defined when type "exponential" is chosen. Defaults to NULL.
#' @param cutoff numeric threshold. It must be in the same unit of the matrix and it must be defined when type "gaussian" is chosen. Defaults to NULL.
#' @param prestep logical parameter TRUE to process the prestep in 2SFCA that divides the population by the number 
#'     of counts in the catchment area of the first step, FALSE otherwise. Defaults to FALSE.
#' @return data frame with origin and acc. If \emph{2sfca} type is chosen, adds the threshold column.
#' @description Process the 2 Step Floating Catchment Area or the Shen (1998) to process an accessibility with competition.
#' @author Beatriz Moura dos Santos
#' @references Shen (1998). Location characteristics of inner-city neighborhoods and employment accessibility of low-wage workers. 
#'     \emph{Environment and planning B: Planning and Design}, 25(3), 345-365.
#' @references Neutens (2015). Accessibility, equity and health care: review and research directions for transport geographers. 
#'     \emph{Journal of Transport Geography}, 43, 14-27.
#' @references RESOLUTION

twosfca <- function(oportunities, matrix, type = '2sfca', betas = NULL, cutoff = NULL, prestep = FALSE){
  # renaming fields
  names(oportunities) <- c('id', 'oportunity', 'pop')
  names(matrix) <- c('origin', 'destination', 'travel_time')
  
  if(type == '2sfca'){
    matrix <- matrix %>% 
      dplyr::filter(.data$travel_time<=cutoff & .data$travel_time >=0)
    
    # processing 2SFCA accessibility
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
      dplyr::select(.data$id, .data$oportunity) 
    
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
      dplyr::mutate(rate = .data$oportunity*.data$rate) %>% 
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
  }else{
    # 1st step - demand potential
    d_pot <- matrix %>% 
      dplyr::left_join(oportunities, by = c('destination' = 'id')) %>% 
      dplyr::mutate(decay = exp(.data$travel_time*betas),
                    pop = .data$pop*.data$decay) %>% 
      dplyr::group_by(.data$origin) %>% 
      dplyr::summarise(d_pot = sum(.data$pop, na.rm=T)) %>% 
      dplyr::ungroup()
    
    # 2nd step - accessibility
    matrix %>% 
      dplyr::left_join(d_pot, by = c('destination' = 'origin')) %>% 
      dplyr::mutate(decay = exp(.data$travel_time*betas),
                    oportunity = .data$oportunity*.data$decay/.data$d_pot) %>% 
      dplyr::group_by(.data$origin) %>% 
      dplyr::summarise(acc = sum(.data$d_pot, na.rm=T)) %>% 
      dplyr::ungroup()
  }
}