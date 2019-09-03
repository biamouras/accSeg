#' Processes the local simulations and local Moran to produce 
#' LISA map. IT STILL CAN'T PROCESS THE BIVARIATE
#' 
#' @param x numeric vector with values of variable 1
#' @param y numeric vector with values of variable 2. Defaults to NULL (not bivariate)
#' @param W numeric matrix of neighbourhood matrix (must be binary matrix)
#' @param shp shapefile or simple feature of locations 
#' @param local_sims numeric vector with local simulations of Moran (result from `simulaMoran`)
#' @param local_moran numeric vector with local values of Moran's I (result from `moranI`)
#' 
#' @output sf with 

LISAmaps <- function(x, y=NULL, W, shp, local_sims, local_moran) {
  
  shp_sf <- st_as_sf(shp)
  nrow <- nrow(local_sims)
  if(is.null(y)) y <- x
  
  # Identificando e classificando os valores significantes locais
  # para fazer o mapa de significancia
  localPseudo <- array(dim= nrow)
  pvalor <- array(dim= nrow)
  apply(1:nrow, 1, function(i) {
    p <- ecdf(local_sims[i,])
    localPseudo[i] <- p(local_moran[i])
    
    if (localPseudo[i] <= 0.0001 | localPseudo[i] >= 0.9999) {
      pvalor[i] <- "p = 0.0001"
    } else if (localPseudo[i] <= 0.001 | localPseudo[i] >= 0.999) {
      pvalor[i] <- "p = 0.001"
    } else if (localPseudo[i] <= 0.01 | localPseudo[i] >= 0.99) {
      pvalor[i] <- "p = 0.01"
    } else if (localPseudo[i] <= 0.05 | localPseudo[i] >= 0.95) {
      pvalor[i] <- "p = 0.05"
    } else {
      pvalor[i] <- "Not significant"
    }
  }
  
  shp_sf$local <- localPseudo
  shp_sf$pvalor <- pvalor
  
  xp <- scale.default(x)
  yp <- scale.default(y)
  
  # Identificando e classificando os padroes dos LISA
  # para o mapa de agrupamento
  patterns <- as.character(interaction(xp > 0, W%*%yp > 0, sep = " - ") ) 
  patterns <- patterns %>% 
    str_replace_all("TRUE","High") %>% 
    str_replace_all("FALSE","Low")
  patterns[shp_sf$pvalor=="Not significant"] <- "Not significant"
  shp_sf$patterns <- patterns
  
  shp_sf
}