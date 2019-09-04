#' Evaluates the p-value of the Moran's I
#' 
#' @inheritParams globalMoranSim
#' 
#' @return local_sim (local Moran's I simulations)
#'         
localMoranSim <- function(shp, x=NULL, y = NULL, W=NULL, nsims = 999){
  
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
  
  # treating the main and auxiliary variables
  n   <- nrow(y)
  local_sims  <- matrix(NA, nrow = n, ncol=nsims)
  y_s <- replicate(nsims, sample(y, size = n)) 
  y_s <- scale.default(y_s)
  
  # processing the simulations
  local_sims  <- apply(y_s, 2, function(s) s*Wn%*%s)
  
  return(local_sims)
}