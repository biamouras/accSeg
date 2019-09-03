#' This function calculate the local population intensity for all groups based on a time matrix.
#'
#' @param self list from readAttributesFile result
#' @param matrix data frame with origin, destination and travel_time structure
#' @param bandwidth numeric bandwidth for neighborhood (same unit of matrix). Defaults to 60 (minutes)
#' @param weightmethod numeric to choose the method to be used: gaussian , bi-square and moving window.
#'                     Defaults to gaussian
#' @param matrixpath string path for input time matrix. The matrix format MUST be as 'gathered' data file: origin, destination, travel time/distance
#' @param filepath string path with file to be read
#' @return list with:
#'    attribute matrix: raw data,
#'    location: x_coord, y_coord,
#'    pop: population of each group,
#'    n_group: number of groups,
#'    n_location: number of tracts,
#'    pop_sum: total population of each tract,
#'    tract_id: tract id,
#'    locality: population intensity for all groups.

popIntensity_matrix <- function(self=NULL,
                                matrix=NULL,
                                bandwidth=60,
                                weightmethod='gaussian',
                                filepath=NULL,
                                matrixpath=NULL){

  # analysing inputs and reading attributes
  if(is.null(self) & !is.null(filepath)){
    self <- readAttributesFile(filepath)
  } else if(is.null(filepath) & is.null(self)){
    stop('No attributes or file path inserted!')
  }
  if(is.null(matrix) & !is.null(matrixpath)){
    message('Reading matrix file')
    matrix <- utils::read.csv(matrixpath, stringsAsFactors=FALSE)
  } else if(is.null(matrixpath) & is.null(matrix)){
    stop('No matrix path inserted!')
  }
  names(matrix) <- c('origin', 'destination', 'travel_time')

  # getting attributes values
  ids <- self$tract_id
  n_subgroup <- self$n_group
  n_local <- self$n_location
  pop <- self$pop %>%
    tidyr::gather(key='group', value='population', (1:n_subgroup)+1)
  pop$id <- as.character(pop$id)

  message('Filtering origin IDs')
  # selecting the same population and matrix localities
  matrix <- dplyr::filter(matrix, .data$origin %in% ids)

  message('Processing weights')
  matrix <- getWeight(matrix, bandwidth, weightmethod)

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

  self$locality <- locality_temp

  return(self)
}
