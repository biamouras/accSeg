#' Computes local dissimilarity for all groups.
#'
#' @param pop A data frame with id, group, and number of population.
#' @param pop_intensity A data frame with id and population intensity for all groups. Result from \code{\link{popIntensity}}.
#' @return  A data frame with id and local indexes of dissimilarity.
#' @author Beatriz Moura dos Santos
#' @references Feitosa, Camara, Monteiro, Koschitzi & Silva (2007). 
#'    Global and local spatial indices of urban segregation. 
#'    \emph{International Journal of Geographical Information Science}, 21(3), 299-323.
#' @references Sousa (2017). Segregation Metrics. 
#'    \url{https://github.com/sandrofsousa/Resolution/blob/master/Pysegreg/segregationMetrics.py}.

localDissimilarity <- function(pop, pop_intensity){
  
  # processing the ids
  ls <- idProcess(pop, pop_intensity)
  pop <- ls[[1]]
  pop_intensity <- ls[[2]]
  rm(ls)
  gc()
  
  # processing population
  Njm <- pop %>% 
    tidyr::gather('group', 'Njm', -.data$id)
  
  Nm <- Njm %>% 
    dplyr::group_by(.data$group) %>% 
    dplyr::summarise(Nm = sum(.data$Njm, na.rm=T))
  
  Nj <- Njm %>% 
    dplyr::group_by(.data$id) %>% 
    dplyr::summarise(Nj = sum(.data$Njm, na.rm=T))
  
  # processing the local proportion of group m in locality j
  if(nrow(pop_intensity) == 0){
    ljm <- Njm %>% 
      dplyr::rename(ljm = Njm)
    lj <- Nj %>% 
      dplyr::rename(lj = Nj)
  } else{
    ljm <- pop_intensity %>% 
      tidyr::gather('group', 'ljm', -.data$id)
    
    lj <- ljm %>% 
      dplyr::group_by(.data$id) %>% 
      dplyr::summarise(lj = sum(.data$ljm, na.rm=T))
  }
  
  tjm <- ljm %>% 
    dplyr::left_join(lj, by = 'id') %>% 
    dplyr::mutate(tjm = .data$ljm/.data$lj) 
  
  # processing the proportion of group m in the city
  pop_total <- sum(Njm$Njm, na.rm=T)
  Nm <-  Nm %>% 
    dplyr::mutate(tm = .data$Nm/pop_total)
  
  # processing index i
  index_i <- sum(Nm$tm * (1 - Nm$tm), na.rm = T)
  
  # processing local dissimilarity
  local_diss <- Nj %>% 
    dplyr::left_join(tjm, by = 'id') %>% 
    dplyr::left_join(Nm, by = 'group') %>% 
    dplyr::mutate(diss_group = (.data$Nj/(2*pop_total*index_i))*abs(.data$tjm-.data$tm)) %>% 
    dplyr::group_by(.data$id) %>% 
    dplyr::summarise(diss = sum(.data$diss_group, na.rm=T))
  
  return (local_diss)
}

#' Computes the local entropy score for a unit area Ei (diversity).
#'
#' @inheritParams localDissimilarity
#' @return A data frame with id and local indexes of entropy.
#' @author Beatriz Moura dos Santos
#' @references Iceland (2004. 
#'    The multigroup entropy index (also known as Theil’s H or the information theory index). 
#'    \emph{US Census Bureau}, Retrieved July, 31, 2006.
#' @references Sousa (2017). Segregation Metrics. 
#'    \url{https://github.com/sandrofsousa/Resolution/blob/master/Pysegreg/segregationMetrics.py}.

localEntropy <- function(pop, pop_intensity){
  
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
  
  if(nrow(pop_intensity) == 0){
    ljm <- Njm %>% 
      dplyr::rename(ljm = Njm)
    lj <- Nj %>% 
      dplyr::rename(lj = Nj)
  } else{
    ljm <- pop_intensity %>% 
      tidyr::gather('group', 'ljm', -.data$id)
    
    lj <- ljm %>% 
      dplyr::group_by(.data$id) %>% 
      dplyr::summarise(lj = sum(.data$ljm, na.rm=T))
  }
  
  tjm <- ljm %>% 
    dplyr::left_join(lj, by = 'id') %>% 
    dplyr::mutate(tjm = .data$ljm/.data$lj) 
  
  entropy <- tjm %>% 
    dplyr::mutate(ent_group = .data$tjm * log(1/.data$tjm)) %>% 
    dplyr::group_by(.data$id) %>% 
    dplyr::summarise(ent = sum(.data$ent_group))
  
  return(entropy)
}

#' Computes the local exposure index of group m to group n.
#' in situations where m=n, then the result is the isolation index.
#'
#' @inheritParams localDissimilarity
#' @return data frame with local exposure indexes, size number of localities x n^2 (n=number of groups).
#' @author Beatriz Moura dos Santos
#' @references Feitosa, Camara, Monteiro, Koschitzi & Silva (2007). 
#'    Global and local spatial indices of urban segregation. 
#'    \emph{International Journal of Geographical Information Science}, 21(3), 299-323.
#' @references Sousa (2017). Segregation Metrics. 
#'    \url{https://github.com/sandrofsousa/Resolution/blob/master/Pysegreg/segregationMetrics.py}.

localExposure <- function(pop, pop_intensity){
  
  # processing the ids
  ls <- idProcess(pop, pop_intensity)
  pop <- ls[[1]]
  pop_intensity <- ls[[2]]
  rm(ls)
  gc()
  
  # processing local exposure
  Njm <- pop %>% 
    tidyr::gather('group', 'Njm', -.data$id)
  
  Nm <- Njm %>% 
    dplyr::group_by(.data$group) %>% 
    dplyr::summarise(Nm = sum(.data$Njm, na.rm=T))
  
  local_expo <- Njm %>% 
    dplyr::left_join(Nm, by = 'group') %>% 
    dplyr::mutate(l_expo = .data$Njm/.data$Nm)
  
  # processing local rate
  if (nrow(pop_intensity) == 0){
    Nj <- Njm %>% 
      dplyr::group_by(.data$id) %>% 
      dplyr::summarise(Nj = sum(.data$Njm, na.rm=T))
    
    tjm <- Njm %>% 
      dplyr::left_join(Nj, by = 'id') %>% 
      dplyr::mutate(tjm = .data$Njm/.data$Nj)
  } else{
    
    ljn <- pop_intensity %>% 
      tidyr::gather('group', 'ljn', -.data$id)
    
    lj <- ljn %>% 
      dplyr::group_by(.data$id) %>% 
      dplyr::summarise(lj = sum(.data$ljn, na.rm=T))
    
    tjm <- ljn %>% 
      dplyr::left_join(lj, by = 'id') %>% 
      dplyr::mutate(tjm = .data$ljn/.data$lj) %>% 
      dplyr::rename(gp = .data$group)
  }
  
  # processing exposure
  exposure <- local_expo %>% 
    dplyr::left_join(tjm, by = 'id') %>% 
    dplyr::mutate(expo = .data$l_expo*.data$tjm,
                  m_n = paste0(.data$group,'_',.data$gp)) %>% 
    dplyr::select(.data$id, .data$expo, .data$m_n) %>% 
    tidyr::spread(.data$m_n, .data$expo)
  
  return(exposure)
}

#' Computes the local entropy index H for all localities.
#'
#' @inheritParams localDissimilarity
#' @return array with scores for n groups (number of groups).
#' @author Beatriz Moura dos Santos
#' @references Iceland (2004. 
#'    The multigroup entropy index (also known as Theil’s H or the information theory index). 
#'    \emph{US Census Bureau}, Retrieved July, 31, 2006.
#' @references Sousa (2017). Segregation Metrics. 
#'    \url{https://github.com/sandrofsousa/Resolution/blob/master/Pysegreg/segregationMetrics.py}.

localIndexH <- function(pop, pop_intensity){
  
  local_entropy <- localEntropy(pop, pop_intensity)
  global_entropy <- globalEntropy(pop, pop_intensity)
  
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
  
  N <- sum(Njm$Njm, na.rm = T)
  
  # processing local Theil's
  h_local <- Nj %>% 
    dplyr::left_join(local_entropy, by = 'id') %>% 
    dplyr::mutate(iH = .data$Nj*(global_entropy - .data$ent)/(global_entropy*N)) %>% 
    dplyr::select(.data$id, .data$iH)
  
  return(h_local)
}

#' Summary of local segregation indexes
#'
#' @inheritParams localDissimilarity
#' @return data frame with calculated indexes for each tract with the following columns
#' \describe{
#'   \item{id}{tract id}
#'   \item{diss}{dissimilarity index}
#'   \item{ent}{entropy}
#'   \item{iH}{index H}
#'   \item{exp_m_n}{exposure inde for group m with group n}
#' }  
#' @author Beatriz Moura dos Santos]
#' @references Feitosa, Camara, Monteiro, Koschitzi & Silva (2007). 
#'    Global and local spatial indices of urban segregation. 
#'    \emph{International Journal of Geographical Information Science}, 21(3), 299-323.
#' @references Iceland (2004. 
#'    The multigroup entropy index (also known as Theil’s H or the information theory index). 
#'    \emph{US Census Bureau}, Retrieved July, 31, 2006.
#' @references Sousa (2017). Segregation Metrics. 
#'    \url{https://github.com/sandrofsousa/Resolution/blob/master/Pysegreg/segregationMetrics.py}.

localSegreg <- function(pop, pop_intensity){
  
  results <- localDissimilarity(pop, pop_intensity) %>% 
    dplyr::left_join(localEntropy(pop, pop_intensity), by = 'id') %>% 
    dplyr::left_join(localIndexH(pop, pop_intensity), by = 'id') %>% 
    dplyr::left_join(localExposure(pop, pop_intensity), by = 'id')
  
  return(results)
}
