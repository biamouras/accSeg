#' Summary of local segregation indexes
#'
#' @param self list with population intensity calculated (result from popIntensity)
#' @return data frame with calculated indexes for each tract with the following columns
#'   id: tract id
#'   diss: dissimilarity index
#'   ent: entropy
#'   iH: index H
#'   exp_m_n: exposure inde for group m with group n

localSegreg <- function(self){

  results <- data.frame(
    id = self$locality$id,
    diss = localDissimilarity(self),
    ent = localEntropy(self),
    iH = localIndexH(self)
  )

  results <- cbind(results, localExposure(self))

  return(results)
}
