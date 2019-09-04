#' Processes the local Moran's I bivariate
#'
#' @inheritParams globalMoran
#'
#' @return the local values of Moran's I

localMoran <- function(shp, x=NULL, y = NULL, W=NULL){
  
  # if shp is NULL then it must have x and W
  if(is.null(shp)){
    # if x is NULL then stop
    if(is.null(x)) stop('shp and x are NULL. Please add any of these parameters')
    # if W is NULL then stop
    if(is.null(W)) stop('shp and W are NULL. Please add any of these parameters')
    # if y is NULL the process is not bivariate
    if(is.null(y)) y <- x
  } else{
    # if class of shp is sf
    if(class(shp)[1] == "sf") {
      shp_sf <- shp
      shp <- sf::as_Spatial(shp)
    } else{
      shp_sf <- sf::st_as_sf(shp)
    }
    # process the W
    W <- nbMatrix(shp)
    # x
    nm <- names(shp_sf)
    x <- shp_sf[[nm[2]]]
    if(ncol(shp) == 2){
      y <- x
    } else if(ncol(shp) == 3){
      y <- shp_sf[[nm[3]]]
    } else {stop('shp is bigger or shorter than expected')}
  }
  
  # scaling variables
  xp <- as.numeric(scale.default(x))
  yp <- as.numeric(scale.default(y))
  
  # size of the sample
  n <- length(x)
  
  # standardizing the matrix, which line sum equals to 1
  Wn <- W/rowSums(W)
  Wn[which(is.na(Wn))] <- 0
  
  # processing local Moran
  local  <- (xp*Wn%*%yp)
  
  return(local)
}