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
  global_diss <- sum(local_diss[,-1], na.rm=T)
  
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
  
  # processing the ids
  ls <- idProcess(pop, pop_intensity)
  pop <- ls[[1]]
  pop_intensity <- ls[[2]]
  rm(ls)
  gc()
  
  # processing population
  Njm <- pop %>% 
    tidyr::gather('group', 'Njm', -.data$id)
  
  Nj <- Njm %>% 
    dplyr::group_by(.data$id) %>% 
    dplyr::summarise(Nj = sum(.data$Njm, na.rm=T))
  
  Nm <- Njm %>% 
    dplyr::group_by(.data$group) %>% 
    dplyr::summarise(Nm = sum(.data$Njm, na.rm=T))
  
  N <- sum(Njm$Njm, na.rm=T)
  
  n <- nrow(Nm)
  
  score <-  Nm %>% 
    dplyr::mutate(score = .data$Nm/N * log(.data$Nm/N)) 
  
  entropy <- sum(score$score, na.rm = T)
  
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
  local_exp  <- accSeg:::localExposure(pop, pop_intensity)
  global_exp  <- local_exp %>% 
    tidyr::gather('m_n', 'expo', -.data$id) %>% 
    dplyr::group_by(.data$m_n) %>% 
    dplyr::summarise(global_exp = sum(.data$expo, na.rm=T)) %>% 
    tidyr::separate(.data$m_n, into = c('m','n'),sep='_')
  
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
  h_global <- sum(h_local[,-1], na.rm=T)
  
  return(h_global)
}

#' Summary of global segregation indexes
#'
#' @inheritParams localDissimilarity
#' @param digits A numeric object of number of digits to display. Defaults to 5.
#' @return list with calculated global indexes with the following items:
#'   \describe{
#'     \item{diss}{Numeric dissimilarity index.}
#'     \item{ent}{Numeric entropy.}
#'     \item{iH}{Numeric index H.}
#'     \item{exp}{Data frame exposure indexes for group m with group n (size n x n).}
#'   }

globalSegreg <- function(pop, pop_intensity, digits=5){
  
  diss <- globalDissimilarity(pop, pop_intensity)
  exp <- globalExposure(pop, pop_intensity)
  ent <- globalEntropy(pop, pop_intensity)
  iH <- globalIndexH(pop, pop_intensity)
  
  exp_out <- as.matrix(xtabs(global_exp ~ m + n, data = exp))
  
  message(paste0('Global Segregation Indexes\n\n' ,
                 'Dissimilarity Index\n', round(diss, digits),'\n\n',
                 'Exposure Indexes\n', paste0(utils::capture.output(round(exp_out, digits)), collapse='\n'),'\n\n',
                 'Entropy Index\n', round(ent, digits), '\n\n',
                 'Index H\n', round(iH, digits)))
  
  exp <- dplyr::rename(exp, value = global_exp) %>% 
    dplyr::mutate(measure = paste0(.data$m, '_', .data$n)) %>% 
    dplyr::select(-m,-n
                  )
  results <- data.frame(diss = diss,
                        ent = ent,
                        iH = iH) %>% 
    tidyr::gather('measure', 'value', 1:3) %>% 
    dplyr::bind_rows(exp)
  
  return(results)
}
