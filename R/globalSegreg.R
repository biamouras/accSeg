#' Summary of global segregation indexes
#'
#' @param self list with population intensity calculated (result from \code{popIntensity_matrix})
#' @param digits numeric number of digits of display. Defaults to 5
#' @return list with calculated global indexes with the following items
#'   diss: numeric dissimilarity index
#'   ent: numeric entropy
#'   iH: numeric index H
#'   exp: data frame exposure indexes for group m with group n (size n x n)

globalSegreg <- function(self, digits=5){

  diss <- globalDissimilarity(self)
  exp <- globalExposure(self)
  ent <- globalEntropy(self)
  iH <- globalIndexH(self)

  message(paste0('Global Segregation Indexes\n\n' ,
                 'Dissimilarity Index\n', round(diss, digits),'\n\n',
                 'Exposure Indexes\n', paste0(utils::capture.output(round(exp, digits)), collapse='\n'),'\n\n',
                 'Entropy Index\n', round(ent, digits), '\n\n',
                 'Index H\n', round(iH, digits)))

  results <- list(
    diss = diss,
    ent = ent,
    iH = iH,
    exp = exp
  )

  return(results)
}
