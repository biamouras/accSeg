#' Processes the neighbourhood matrix
#' 
#' @param shp shapefile or simple feature of locations
#' @return numeric binary neighbourhood matrix

nbMatrix <- function(shp){
  if(class(shp)[1] == "sf") {
    shp <- sf::as_Spatial(shp)
  } 
  nb <- spdep::poly2nb(shp)
  W <- spdep::nb2mat(nb, style = "B", zero.policy = T)
  return(W)
}

