#' Evaluates the p-value of the Moran's I
#' 
#' @param x numeric vector with values of variable 1
#' @param y numeric vector with values of variable 2. Defaults to NULL (not bivariate)
#' @param W numeric matrix of neighbourhood matrix (must be binary matrix)
#' @param nsims numeric value of number of simulations. Defaults to 999.
#' 
#' @output list with global_sim (global Moran's I simulations) and 
#'         local_sim (local Moran's I simulations)and lagged (neighbours' mean) 
#'         
simulaMoran <- function(x, y = NULL, W, nsims = 999){
  
  if(is.null(y)) y <- x
  
  # treating the main and auxiliary variables
  n   <- nrow(W)
  IDs <- 1:n
  xp <- as.vector(scale.default(x))
  W[which(is.na(W))] <- 0
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