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

#' Processes the variables to calculate Moran and local Moran
#' 
#' @inheritParams globalMoran
#' @return A list with variables to process Moran

aux <- function(shp, x=NULL, y=NULL, W=NULL){
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
      df <- sf::st_drop_geometry(shp)
      shp <- sf::as_Spatial(shp)
    } else{
      df <- sf::st_drop_geometry(sf::st_as_sf(shp))
    }
    # process the W
    W <- nbMatrix(shp)
    # x
    nm <- names(df)
    x <- df[[nm[2]]]
    if(ncol(shp) == 2){
      y <- x
    } else if(ncol(shp) == 3){
      y <- df[[nm[3]]]
    } else {stop('shp is bigger or shorter than expected')}
  }
  
  # scaling variables
  xp <- as.numeric(scale.default(x))
  yp <- as.numeric(scale.default(y))
  # standardizing the matrix, which line sum equals to 1
  Wn <- W/rowSums(W)
  Wn[which(is.na(Wn))] <- 0
  
  # treating the main and auxiliary variables
  n <- length(y)
  
  aux <- list(n = n,
              xp = xp,
              yp = yp,
              Wn = Wn)
  return(aux)
}

#' Processes the global Moran's I bivariate
#'
#' @param shp shapefile or simple feature to process the neighbourhood. Expected structure: id, x, y(optional)
#' @param x numeric vector with values of variable 1. Defaults to NULL.
#' @param y numeric vector with values of variable 2. Defaults to NULL (not bivariate).
#' @param W numeric matrix of neighbourhood matrix from \code{\link{nbMatrix}}. Defaults to NULL.
#' @return the global value of Moran's I

globalMoran <- function(shp, x=NULL, y=NULL, W=NULL){
  
  aux <- aux(shp, x, y, W)
  
  # processing global Moran
  global <- (aux$n/sum(aux$Wn, na.rm=T)) * (aux$xp%*%aux$Wn%*%aux$yp)/(aux$xp%*%aux$xp)
  
  return(global)
}


#' Evaluates the p-value of the Moran's I
#' 
#' @param nsims numeric value of number of simulations. Defaults to 999.
#' @inheritParams globalMoran
#' @return global_sims (global Moran's I simulations) 

globalMoranSim <- function(shp, x=NULL, y = NULL, W = NULL, nsims = 999){
  
  aux <- aux(shp, x, y, W)
  global_sims <- array(dim=nsims)
  
  y_s <- replicate(nsims, sample(aux$yp, size = aux$n))
  
  # processing the simulations
  global_sims  <- apply(y_s, 2, function(s) (aux$n/sum(aux$Wn, na.rm=T)) * (aux$xp%*%aux$Wn%*%s)/(aux$xp%*%aux$xp))
  
  return(global_sims)
}

#' Processes the lagged value to plot Moran scatter plot
#'
#' @inheritParams globalMoran
#' @return lagged values to plot the Moran scatter plot

laggedMoran <- function(shp, x=NULL, y = NULL, W=NULL){
  
  aux <- aux(shp, x, y, W)
  
  # process the W
  W <- nbMatrix(shp)
  # standardizing the matrix, which line sum equals to 1
  Wl <- spdep::mat2listw(W)
  
  # neighbours' mean
  lagged <- data.frame(x = aux$xp,
                       y =scale.default(spdep::lag.listw(Wl, aux$yp)))
  
  return(lagged)
}

#' Processes the local Moran's I bivariate
#'
#' @inheritParams globalMoran
#'
#' @return the local values of Moran's I

localMoran <- function(shp, x=NULL, y = NULL, W=NULL){
  
  aux <- aux(shp, x, y, W)
  
  # processing local Moran
  local  <- (aux$xp*aux$Wn%*%aux$yp)
  
  return(local)
}

#' Evaluates the p-value of the Moran's I
#' 
#' @inheritParams globalMoranSim
#' 
#' @return local_sim (local Moran's I simulations)
#'         
localMoranSim <- function(shp, x=NULL, y = NULL, W=NULL, nsims = 999){
  
  aux <- aux(shp, x, y, W)
  # processing
  local_sims  <- matrix(NA, nrow = aux$n, ncol=nsims)
  local_sims <- replicate(nsims, {
    y_s <- sample(aux$yp, size = aux$n)
    aux$xp*aux$Wn%*%y_s
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
simulaMoran <- function(shp, x=NULL, y = NULL, W=NULL,  nsims = 999){
  
  global_sims <- globalMoranSim(shp, x, y, W)
  local_sims  <- localMoranSim(shp, x, y, W)
  # results
  list(global_sims = global_sims,
       local_sims  = local_sims)
}

#' Processes the local simulations and local Moran to produce 
#' LISA map. IT STILL CAN'T PROCESS THE BIVARIATE
#' 
#' @param local_sims numeric vector with local simulations of Moran (result from \code{\link{localMoranSim}}
#' @param local_moran numeric vector with local values of Moran's I (result from \code{\link{localMoran}}
#' @inheritParams localMoran
#' @return sf with localPseudo (p-value of local), pvalue (p-value classified), patterns (cluster LISA map input)

LISAmaps <- function(shp, x=NULL, y=NULL, W=NULL, nsims = 999){
  
  message('Processing local simulations')
  local_sims <- localMoranSim(shp, x, y, W)
  message('Processing local Moran')
  local_moran <- localMoran(shp, x, y, W)
  nrow <- nrow(local_sims)
  
  #REVER PARA RODAR COM LINHAS IGUAIS A ZERO
  # classifying the local p-values
  df_pvalue <- purrr::map_df(1:nrow, function(i) {
    message(i)
    p <- quantile(local_sims[i,], probs = c(0, 0.001, 0.01, 0.05, 
                                            0.95, 0.99, 0.999, 1))
    p[1] <- -Inf
    p[8] <- Inf
    
    pvalue <- as.character(factor(.bincode(
      local_moran[i],
      breaks = p,
      include.lowest = T, 
      right = F),
      levels = 1:7,
      labels = c("p = 0.001", "p = 0.01", "p = 0.05", "Not significant",
                 "p = 0.05","p = 0.01", "p = 0.001")))
    
    data.frame(pvalue = pvalue)
  })
  
  # if class of shp is sf
  if(class(shp)[1] == "sf") {
    shp_sf <- shp
  } else{
    shp_sf <- sf::st_as_sf(shp)
  }
  shp_sf <- dplyr::bind_cols(shp_sf, df_pvalue)
  
  lagged <- laggedMoran(shp, x, y, W)
  # classifying the LISA patterns (HH / HL / LH / LL)
  patterns <- as.character(interaction(lagged$x > 0, lagged$y > 0, sep = " - ") ) 
  patterns <- patterns %>% 
    stringr::str_replace_all("TRUE","High") %>% 
    stringr::str_replace_all("FALSE","Low")
  patterns[shp_sf$pvalue=="Not significant"] <- "Not significant"
  
  shp_sf$patterns <- patterns
  
  return(shp_sf)
}

