#' This function call local dissimilarity and compute the sum from individual values.
#'
#' @param self list with population intensity calculated (result from popIntensity)
#' @return display global value of dissimilarity

globalDissimilarity <- function(self){

  local_diss <- localDissimilarity(self)
  global_diss <- sum(local_diss, na.rm=T)

  return(global_diss)
}
