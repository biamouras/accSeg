#' Function to compute global index H
#'
#' @param self list with population intensity calculated (result from popIntensity)
#' @return: values with global index for each group.

globalIndexH <- function(self){

  h_local <- localIndexH(self)
  h_global <- sum(h_local, na.rm=T)

  return(h_global)
}
