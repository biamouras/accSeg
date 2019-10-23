#' Processes the ids of data frame 1 and data frame 2
#' 
#' @param df1 A data frame with ids to be processed.
#' @param df2 A data frame with ids to be processed.
#' @return A list with the two data frames processed with common ids and with the same order.
#' #' @author Beatriz Moura dos Santos

idProcess <- function(df1, df2){
  
  df1 <- data.table::data.table(df1)
  df2 <- data.table::data.table(df2)
  
  nm_df1 <- names(df1)
  names(df1) <- c('id', nm_df1[-1])
  nm_df2 <- names(df2)
  names(df2) <- c('id', nm_df2[-1])
  
  ids <- intersect(df1$id, df2$id)
  
  df1 <- data.table::setorder(df1, id)[id %in% ids]
  df2 <- data.table::setorder(df2, id)[id %in% ids]
  
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
  df <- merge(matrix, pop, by.x = 'destination', by.y = 'id')
  df$int <- df$weight * df$population
  
  data.table::data.table(df)[,.(local_int = sum(int, na.rm = T)),by=group]
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
  
  message('Filtering origin and population IDs')
  # selecting the same population and matrix localities
  matrix <- data.table::data.table(matrix)
  matrix <- matrix[origin %in% ids & destination %in% ids]
  
  pop <- data.table::data.table(pop)
  pop <- pop[id %in% ids]
  pop_names <- names(pop)
  pop <- data.table::melt(pop,
                          id.vars = pop_names[1],
                          measure.vars = pop_names[-1],
                          variable.name = 'group',
                          value.name = 'population')
  
  message('Processing weights')
  matrix <- getWeight(matrix, bandwidth, weightmethod)
  
  message('Processing local intensity')
  locality_temp <- matrix[,localIntensity(.SD,pop), by = origin]
  locality_temp <- locality_temp[,.(id=origin, group, local_int)]
  
  
  
  return(locality_temp)
}

#' Processes the population values.
#' 
#' @param pop A data table structured as id and columns with number of population for each group.
#' @return A list with the following items
#' \describe{
#'     \item{Njm}{A data table with id, group and population values by locality.}
#'     \item{Nj}{A data table with id and population values by locality.}
#'     \item{Nm}{A data table with group and population totals.}
#'     \item{N}{A data table with population total.}
#' }
#' @author Beatriz Moura dos Santos
#' 

popSummary <- function(pop){
  pop_names <- names(pop)
  Njm <- data.table::melt(pop,
                          id.vars = pop_names[1],
                          measure.vars = pop_names[-1],
                          variable.names = 'group',
                          value.names = 'Njm') 
  
  Nj <- Njm[,Nj:=sum(Njm, na.rm=T), by = id] 
  
  Nm <- Njm[,Nm:=sum(Njm, na.rm=T), by = group] 
  
  N <- Njm[,sum(Njm, na.rm=T)]
}