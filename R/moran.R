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

#' Processes the global Moran's I bivariate
#'
#' @param shp shapefile or simple feature to process the neighbourhood. Expected structure: id, x, y(optional)
#' @param x numeric vector with values of variable 1. Defaults to NULL.
#' @param y numeric vector with values of variable 2. Defaults to NULL (not bivariate).
#' @param W numeric matrix of neighbourhood matrix from \code{\link{nbMatrix}}. Defaults to NULL.
#' @return the global value of Moran's I

globalMoran <- function(shp, x=NULL, y=NULL, W=NULL){
  
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
  
  # processing global Moran
  global <- (n/sum(Wn, na.rm=T)) * (xp%*%Wn%*%yp)/(xp%*%xp)
  
  return(global)
}


#' Evaluates the p-value of the Moran's I
#' 
#' @param nsims numeric value of number of simulations. Defaults to 999.
#' @inheritParams globalMoran
#' @return global_sims (global Moran's I simulations) 

globalMoranSim <- function(shp, x=NULL, y = NULL, W=NULL, nsims = 999){
  
  # if shp is NULL then it must have a x and a W
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
  n <- length(y)
  global_sims <- NULL
  y_s <- replicate(nsims, sample(y, size = n)) 
  y_s <- scale.default(y_s)
  
  xp <- scale.default(x)
  
  # size of the sample
  n <- length(x)
  
  # standardizing the matrix, which line sum equals to 1
  Wn <- W/rowSums(W)
  Wn[which(is.na(Wn))] <- 0
  
  # processing the simulations
  global_sims  <- apply(y_s, 2, function(s) (n/sum(Wn, na.rm=T)) * (s%*%Wn%*%s)/(s%*%s))
  
  return(global_sims)
}

#' Processes the lagged value to plot Moran scatter plot
#'
#' @inheritParams globalMoran
#' @return lagged values to plot the Moran scatter plot

laggedMoran <- function(x, y = NULL, W=NULL, shp=NULL){
  
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
  yp <- as.numeric(scale.default(y))
  
  # standardizing the matrix, which line sum equals to 1
  Wl <- spdep::mat2listw(W)
  
  # neighbours' mean
  lagged <- scale.default(spdep::lag.listw(Wl, yp))
  
  return(lagged)
}

#' Processes the local simulations and local Moran to produce 
#' LISA map. IT STILL CAN'T PROCESS THE BIVARIATE
#' 
#' @param local_sims numeric vector with local simulations of Moran (result from \code{\link{localMoranSim}}
#' @param local_moran numeric vector with local values of Moran's I (result from \code{\link{localMoran}}
#' @inheritParams localMoran
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
  
  #REVER PARA RODAR COM LINHAS IGUAIS A ZERO
  # classifying the local p-values
  df_pvalue <- purrr::map_df(1:nrow, function(i) {
    message(i)
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

#' Processes the local Moran's I bivariate
#'
#' @inheritParams globalMoran
#'
#' @return the local values of Moran's I

localMoran <- function(x=NULL, y = NULL, W=NULL){
  
  # if x is NULL then stop
  if(is.null(x)) stop('shp and x are NULL. Please add any of these parameters')
  # if W is NULL then stop
  if(is.null(W)) stop('shp and W are NULL. Please add any of these parameters')
  # if y is NULL the process is not bivariate
  if(is.null(y)) y <- x

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

#' Evaluates the p-value of the Moran's I
#' 
#' @inheritParams globalMoranSim
#' 
#' @return local_sim (local Moran's I simulations)
#'         
localMoranSim <- function(x=NULL, y = NULL, W=NULL, nsims = 999){
  
    # if x is NULL then stop
    if(is.null(x)) stop('x is NULL. Please add this parameter')
    # if W is NULL then stop
    if(is.null(W)) stop('W is NULL. Please add this parameter')
    # if y is NULL the process is not bivariate
    if(is.null(y)) y <- x
  
  # treating the main and auxiliary variables
  n <- length(y)
  
  # standardizing the matrix, which line sum equals to 1
  Wn <- W/rowSums(W)
  Wn[which(is.na(Wn))] <- 0
  
  # processing
  local_sims  <- matrix(NA, nrow = n, ncol=nsims)
  local_sims <- replicate(nsims, {
    y_s <- sample(y, size = n)
    y_s <- scale.default(y_s)
    apply(y_s, 2, function(s) s*Wn%*%s)
  })
  local_sims <- as.data.frame(local_sims)
  return(local_sims)
}

#' Evaluates the p-value of the Moran's I
#' 
#' @param nsims numeric value of number of simulations. Defaults to 999.
#' @inheritParams globalMoran
#' 
#' @return list with global_sim (global Moran's I simulations) and 
#'         local_sim (local Moran's I simulations)and lagged (neighbours' mean) 
#'         
simulaMoran <- function(x, y = NULL, W=NULL, shp=NULL, nsims = 999){
  
  # if y is NULL the process is not bivariate
  if(is.null(y)) y = x
  
  # if W is NULL calls the nbMatrix
  if(is.null(W)){
    if(is.null(shp)){
      stop('W and shp are NULL. Please add any of these parameters.')
    } else{
      W <- nbMatrix(shp)
    }
  }
  
  # treating the main and auxiliary variables
  n   <- nrow(W)
  IDs <- 1:n
  xp <- as.vector(scale.default(x))
  Wn <- W/rowSums(W)
  global_sims <- NULL
  local_sims  <- matrix(NA, nrow = n, ncol=nsims)
  ID_sample <- sample(IDs, size = n*nsims, replace = T)
  
  y_s <- y[ID_sample]
  y_s <- matrix(y_s, nrow = n, ncol = nsims)
  y_s <- apply(y_s, 1, scale.default)
  
  # if it is NULL 
  if(is.null(y)) {
    # processing the simulations
    global_sims  <- apply(y_s, 1, 
                          FUN = function(x) 
                            (n/sum(Wn))*(x%*%Wn%*%x)/(x%*%x)
    )
    local_sims  <- apply(y_s, 1,function(x) (x*Wn%*%x))
  } else {
    # processing the simulations
    global_sims  <- apply(y_s, 1, 
                          FUN = function(x) 
                            (n/sum(Wn))*(xp%*%Wn%*%x)/(xp%*%xp)
    )
    local_sims  <- apply(y_s, 1,function(x) (xp*Wn%*%x))
  }
  
  
  # results
  list(global_sims = global_sims,
       local_sims  = local_sims)
}