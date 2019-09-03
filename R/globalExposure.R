#' This function call local exposure function and sum the results for the global index.
#'
#' @param self list with population intensity calculated (result from popIntensity)
#' @return data frame of global exposure result, size n x n (n=number of groups)

globalExposure <- function(self){

  m <- self$n_group
  local_exp  <- localExposure(self)
  global_exp  <- apply(local_exp, 2, sum, na.rm=T)
  global_exp <- matrix(global_exp, nrow = m, ncol=m, byrow=T)

  global_exp <- as.data.frame(global_exp,
                            row.names=paste0('G',1:self$n_group))
  names(global_exp) <- paste0('G',1:self$n_group)

  return(global_exp)
}
