#' Computes the global dissimilarity
#'
#' @inheritParams localDissimilarity
#' @return Display global value of dissimilarity.
#' @description This function calls local dissimilarity function and computes the sum from individual values.
#' @author Beatriz Moura dos Santos
#' @references Feitosa, Camara, Monteiro, Koschitzi & Silva (2007). 
#'    Global and local spatial indices of urban segregation. 
#'    \emph{International Journal of Geographical Information Science}, 21(3), 299-323.
#' @references Sousa (2017). Segregation Metrics. 
#'    \url{https://github.com/sandrofsousa/Resolution/blob/master/Pysegreg/segregationMetrics.py}. 

globalDissimilarity <- function(pop, pop_intensity){
  
  local_diss <- localDissimilarity(pop, pop_intensity)
  global_diss <- sum(local_diss, na.rm=T)
  
  return(global_diss)
}

#' Computes the global entropy score E (diversity).
#'
#' @inheritParams localDissimilarity
#' @return Display the diversity score.
#' @author Beatriz Moura dos Santos
#' @references Iceland (2004. 
#'    The multigroup entropy index (also known as Theil’s H or the information theory index). 
#'    \emph{US Census Bureau}, Retrieved July, 31, 2006.
#' @references Sousa (2017). Segregation Metrics. 
#'    \url{https://github.com/sandrofsousa/Resolution/blob/master/Pysegreg/segregationMetrics.py}. 

globalEntropy <- function(pop, pop_intensity){
  
  local_ids <- pop_intensity$id
  
  pop <- pop %>%
    dplyr::filter(.data$id %in% local_ids) %>%
    dplyr::select(-id)
  
  pop_sum <- apply(pop, 1, sum, na.rm=T)
  
  pop_total <- sum(pop, na.rm=T)
  
  prop <- apply(pop, 2, sum, na.rm=T)
  n <- length(prop)
  group_score <- array(dim=n)
  
  group_score <- sapply(1:n, function(i){
    group <- prop[i]
    group_idx <- group / pop_total * log(1 / (group / pop_total))
    group_idx
  })
  
  entropy <- sum(group_score, na.rm=T)
  
  return(entropy)
}

#' Computes the global exposure measures
#'
#' @inheritParams localDissimilarity
#' @return A data frame of global exposure result, size n x n (n=number of groups).
#' @description This function call local exposure function and sum the results for the global index.
#' @author Beatriz Moura dos Santos
#' @references Feitosa, Camara, Monteiro, Koschitzi & Silva (2007). 
#'    Global and local spatial indices of urban segregation. 
#'    \emph{International Journal of Geographical Information Science}, 21(3), 299-323.
#' @references Sousa (2017). Segregation Metrics. 
#'    \url{https://github.com/sandrofsousa/Resolution/blob/master/Pysegreg/segregationMetrics.py}. 

globalExposure <- function(pop, pop_intensity){
  
  m <- ncol(pop)-1
  local_exp  <- localExposure(pop, pop_intensity)
  global_exp  <- apply(local_exp, 2, sum, na.rm=T)
  global_exp <- matrix(global_exp, nrow = m, ncol=m, byrow=T)
  
  global_exp <- as.data.frame(global_exp,
                              row.names=paste0('G',1:m))
  names(global_exp) <- paste0('G',1:m)
  
  return(global_exp)
}

#' Computes global index H
#'
#' @inheritParams localDissimilarity
#' @return: values with global index for each group.
#' @author Beatriz Moura dos Santos
#' @references Iceland (2004. 
#'    The multigroup entropy index (also known as Theil’s H or the information theory index). 
#'    \emph{US Census Bureau}, Retrieved July, 31, 2006.
#' @references Sousa (2017). Segregation Metrics. 
#'    \url{https://github.com/sandrofsousa/Resolution/blob/master/Pysegreg/segregationMetrics.py}. 

globalIndexH <- function(pop, pop_intensity){
  
  h_local <- localIndexH(pop, pop_intensity)
  h_global <- sum(h_local, na.rm=T)
  
  return(h_global)
}

#' Summary of global segregation indexes
#'
#' @inheritParams localDissimilarity
#' @param digits A numeric object of number of digits to display. Defaults to 5.
#' @return list with calculated global indexes with the following items:
#'   \describe{
#'     \item{diss} {Numeric dissimilarity index.}
#'     \item{ent} {Numeric entropy.}
#'     \item{iH} {Numeric index H.}
#'     \item{exp} {Data frame exposure indexes for group m with group n (size n x n).}
#'   }

globalSegreg <- function(pop, pop_intensity, digits=5){
  
  diss <- globalDissimilarity(pop, pop_intensity)
  exp <- globalExposure(pop, pop_intensity)
  ent <- globalEntropy(pop, pop_intensity)
  iH <- globalIndexH(pop, pop_intensity)
  
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
