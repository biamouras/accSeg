#' Processes the ids of data frame 1 and data frame 2
#' 
#' @param df1 A data frame with ids to be processed.
#' @param df2 A data frame with ids to be processed.
#' @return A list with the two data frames processed with common ids and with the same order.
#' #' @author Beatriz Moura dos Santos

idProcess <- function(df1, df2){
  ids <- intersect(df1$id, df2$id)
  
  df1 <- df1 %>%
    dplyr::filter(.data$id %in% ids) %>%
    dplyr::arrange(match(.data$id, ids)) 
  
  df2 <- df2 %>%
    dplyr::filter(.data$id %in% ids) %>%
    dplyr::arrange(match(.data$id, ids)) 
  
  list(df1, df2)
}

#' Computes the weights for neighborhood.
#'
#' @param matrix A data frame with travel time/distance matrix (origin, destination, travel_time).
#' @param bandwidth A numeric objetc bandwidth selected to perform neighborhood (same unit of distance).
#' @param weightmethod A string object to choose the method to be used: gaussian , bi-square and moving window. 
#'    Defaults to gaussian.
#' @return A data frame with weight value for internal use (process the population intensity).
#' @author Beatriz Moura dos Santos
#' @references Smith, Goodchild & Longley (2018). 
#'    \emph{Geospatial Analysis. A Comprehensive Guide to PrinciplesTechniques and Software Tools}. 
#'    \url{https://www.spatialanalysisonline.com/HTML/index.html}

getWeight <- function(matrix, bandwidth, weightmethod = 'gaussian'){
  matrix$rate <- matrix$travel_time/bandwidth
  if(weightmethod=='gaussian'){
    matrix$weight <- exp(-0.5 * matrix$rate^2)
  } else if(weightmethod=='bi-squared'){
    matrix$weight <- (1 - matrix$rate^2)^2
    matrix[matrix$rate > 1, 'weight'] <- 0
  } else if(weightmethod=='moving window'){
    matrix$weight <- 1
    matrix[matrix$rate > 1, 'weight'] <- 0
  } else{
    stop('Invalid weight method selected!')
  }
  
  matrix <- matrix[,c('origin', 'destination','weight')]
  return(matrix)
}

#' Processes the local population intensity for each locality
#'
#' @param matrix A data frame with the weight processed in \code{link{getWeight}}. 
#'    It must be structured as origin, destination and weight.
#' @param pop A data frame structured as id, group and number of population.
#' @return A data frame with local intensity for each group.
#' @import dplyr
#' @author Beatriz Moura dos Santos
#' @references Feitosa, Camara, Monteiro, Koschitzi & Silva (2007). 
#'    Global and local spatial indices of urban segregation. 
#'    \emph{International Journal of Geographical Information Science}, 21(3), 299-323.
#' @references Sousa (2017). Segregation Metrics \url{https://github.com/sandrofsousa/Resolution/blob/master/Pysegreg/segregationMetrics.py}

localIntensity <- function(matrix, pop){
  matrix %>%
    dplyr::right_join(pop, by=c('destination'='id')) %>%
    dplyr::select(-.data$destination) %>%
    dplyr::mutate(int = .data$weight * .data$population) %>%
    dplyr::group_by(.data$group) %>%
    dplyr::summarise(local_int = sum(.data$int, na.rm=T)) %>%
    dplyr::ungroup()
}

#' Computes the local population intensity for all groups based on a time/distance matrix.
#'
#' @param matrix A data frame with travel time/distance matrix (origin, destination, travel_time).
#' @param pop A data frame structured as id and columns with number of population for each group.
#' @param bandwidth A numeric bandwidth for neighborhood (same unit of matrix). Defaults to 60 (minutes).
#' @param weightmethod A string object to choose the method to be used: gaussian , bi-square and moving window. 
#'    Defaults to gaussian.
#' @return A data frame with population intensity for all groups.
#' @author Beatriz Moura dos Santos
#' @references Smith, Goodchild & Longley (2018). 
#'    \emph{Geospatial Analysis. A Comprehensive Guide to PrinciplesTechniques and Software Tools}. 
#'    \url{https://www.spatialanalysisonline.com/HTML/index.html}

popIntensity <- function(matrix, pop,
                         bandwidth=60,
                         weightmethod='gaussian'){
  
  names(matrix) <- c('origin', 'destination', 'travel_time')
  
  # getting attributes values
  ids <- intersect(pop$id, unique(matrix$origin))
  n_local <- length(ids)
  pop$id <- as.character(pop$id)
  
  message('Filtering origin IDs')
  # selecting the same population and matrix localities
  matrix <- matrix %>% 
    dplyr::mutate(origin = as.character(.data$origin),
                  destination = as.character(.data$destination)) %>% 
    dplyr::filter(.data$origin %in% ids)
  
  message('Processing weights')
  matrix <- getWeight(matrix, bandwidth, weightmethod)
  
  pop <- pop %>% 
    dplyr::filter(.data$id %in% ids) %>% 
    tidyr::gather('group', 'population', -id)
  
  message('Processing local intensity')
  locality_temp <- matrix %>%
    dplyr::group_by(.data$origin) %>%
    dplyr::group_modify(~ localIntensity(.,pop)) %>%
    dplyr::ungroup()
  
  locality_temp[is.na(locality_temp$local_int), 'local_int'] <- 0
  locality_temp[is.nan(locality_temp$local_int), 'local_int'] <- 0
  locality_temp[locality_temp$local_int < 0, 'local_int'] <- 0
  
  locality_temp <- locality_temp %>%
    tidyr::spread(.data$group, .data$local_int)
  names(locality_temp) <- c('id', names(locality_temp)[-1])
  
  return(locality_temp)
}
