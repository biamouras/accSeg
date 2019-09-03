#' Processes the Moran's I bivariate
#'
#' @param x numeric vector with values of variable 1
#' @param y numeric vector with values of variable 2. Defaults to NULL (not bivariate)
#' @param W numeric matrix of neighbourhood matrix (must be binary matrix)
#'
#' @output list with global Moran's  I, local Moran's I and lagged (neighbours' mean) that allows to plot the
#'         dispersion and cluster maps

moranI <- function(x, y = NULL, W){

  # if y is NULL the process is not bivariate
  if(is.null(y)) y = x

  # scaling variables
  xp <- scale.default(x)
  yp <- scale.default(y)

  # size of the sample
  n <- length(x)

  # adding 0s to the matrix
  W[which(is.na(W))] <- 0

  # standardizing the matrix, which line sum equals to 1
  Wn <- W/rowSums(W)
  Wl <- mat2listw(W)

  # processing global Moran
  global <- (n/sum(Wn)) * (xp%*%Wn%*%yp)/(xp%*%xp)
  # processing local Moran
  local  <- (xp*Wn%*%yp)
  # neighbours' mean
  lagged <- scale.default(lag.listw(Wl, y))

  # como saira a resposta
  list(global = global, local  = as.numeric(local), lagged = lagged)
}
