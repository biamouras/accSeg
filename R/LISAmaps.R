#' Processes the local simulations and local Moran to produce 
#' LISA map. IT STILL CAN'T PROCESS THE BIVARIATE
#' 
#' @param local_sims numeric vector with local simulations of Moran (result from \code{\link{localMoranSim}}
#' @param local_moran numeric vector with local values of Moran's I (result from \code{\link{localMoran}}
#' @inheritParams localMoran
#' 
#' @return sf with localPseudo (p-value of local), pvalue (p-value classified), patterns (cluster LISA map input)

LISAmaps <- function(shp=NULL, x=NULL, y=NULL, W=NULL, local_sims=NULL, local_moran=NULL){
  
  # if shp is NULL then it must have x and W
  if(is.null(shp)){
    # if x is NULL then stop
    if(is.null(x)) stop('shp and x are NULL. Please add any of these parameters')
    # if W is NULL then stop
    if(is.null(W)) stop('shp and W are NULL. Please add any of these parameters')
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
    if(ncol(shp) == 3){
      y <- shp_sf[[nm[3]]]
    } #else {stop('shp is bigger or shorter than expected')}
  }
  
  message('Processing local simulations')
  # if local_sims is NULL calls localMoranSim
  if(is.null(local_sims)){
    local_sims <- localMoranSim(x=x, y=y, W=W)
  }
  message('Processing local Moran')
  # if local_moran is NULL calls localMoran
  if(is.null(local_moran)){
    local_moran <- localMoran(x=x, y=y, W=W)
  }
  nrow <- nrow(local_sims)
  
  # classifying the local p-values
  df_pvalue <- purrr::map_df(1:nrow, function(i) {
    p <- stats::ecdf(local_sims[i,])
    localPseudo <- p(local_moran[i])
    
    pvalue <- as.character(factor(cut(
      localPseudo,
      breaks = c(0, 0.0001, 0.001, 0.01, 0.05, 
                 0.95, 0.99, 0.999, 0.9999, 1),
      labels = c("p = 0.0001", "p = 0.001", "p = 0.01", "p = 0.05", "Not significant",
                 "p = 0.05","p = 0.01", "p = 0.001","p = 0.0001"),
      include.lowest = T, 
      right = T)))
    
    data.frame(local = localPseudo,
               pvalue = pvalue)
  })
  
  shp_sf <- dplyr::bind_cols(shp_sf, df_pvalue)
  lagged <- laggedMoran(shp, x, y, W)
  xp <- scale.default(x_aux)
  
  # classifying the LISA patterns (HH / HL / LH / LL)
  patterns <- as.character(interaction(xp > 0, lagged > 0, sep = " - ") ) 
  patterns <- patterns %>% 
    stringr::str_replace_all("TRUE","High") %>% 
    stringr::str_replace_all("FALSE","Low")
  patterns[shp_sf$pvalue=="Not significant"] <- "Not significant"
  
  shp_sf$patterns <- patterns
  
  return(shp_sf)
}
