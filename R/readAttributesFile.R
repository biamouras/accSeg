#' This function reads the csv file and populate the class's attributes.
#' Data has to be exactly in the following format or results will be wrong:
#' area id,  x_coord, y_coord, population of G1, pop G2, pop G3, ..., pop Gn (n groups)
#'
#' @param filepath string path with file to be read
#' @return list with:
#'    attribute matrix: raw data,
#'    location: x_coord, y_coord,
#'    pop: population of each group,
#'    n_group: number of groups,
#'    n_location: number of tracts,
#'    pop_sum: total population of each tract,
#'    tract_id: tract id

readAttributesFile <-  function(filepath){

  raw_data <- utils::read.csv(filepath)
  ncol <- ncol(raw_data)
  if (ncol < 5){
    stop('The data has less than 2 groups or it lacks parameters.\n
         The expected format is (delimeter = ,): id, x_coord, y_coord, population of group 1 (G1), pop of G2, ... pop of Gn.')
  }

  attributeMatrix <- raw_data[,-1]
  n <-  as.numeric(ncol(attributeMatrix))
  pop <- trunc(attributeMatrix[, 3:n])
  id <- as.character(raw_data$id)

  self <- list(
    attributeMatrix = raw_data,
    location = raw_data[, 1:3],
    pop = cbind(id, pop),
    pop_sum = data.frame(id = id, pop_sum = apply(pop, 1, sum, na.rm=T)),
    tract_id = id,
    n_group = n-2,
    n_location = as.numeric(nrow(attributeMatrix))
  )

  return(self)
}
