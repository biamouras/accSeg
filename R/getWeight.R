#' Compute the weights for neighborhood.
#'
#' @param matrix data frame (origin, destination, travel_time structure) to be considered for weighting
#' @param bandwidth numeric bandwidth selected to perform neighborhood (same unit of distance)
#' @param weightmethod numeric to choose the method to be used: gaussian , bi-square and moving window. Defaults to gaussian
#' @return data frame with weight value for internal use (population intensity)

getWeight <- function(matrix, bandwidth, weightmethod = 'gaussian'){

  if(weightmethod=='gaussian'){
    matrix <- matrix %>%
      dplyr::mutate(rate = .data$travel_time/bandwidth,
                    weight = exp(-0.5 * .data$rate^2)) %>%
      dplyr::select(-.data$travel_time, -.data$rate)
  } else if(weightmethod=='bi-squared'){
    matrix <- matrix %>%
      dplyr::mutate(rate = .data$travel_time/bandwidth,
                    weight = ifelse(.data$rate <= 1, (1 - .data$rate^2)^2, 0)) %>%
      dplyr::select(-.data$travel_time, -.data$rate)
  } else if(weightmethod=='moving window'){
    matrix <- matrix %>%
      dplyr::mutate(rate = .data$travel_time/bandwidth,
                    weight = ifelse(.data$rate <= 1, 1, 0))%>%
      dplyr::select(-.data$travel_time, -.data$rate)
  } else{
    stop('Invalid weight method selected!')
  }

  return(matrix)
}
