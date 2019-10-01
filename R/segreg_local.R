#' Computes local dissimilarity for all groups.
#'
#' @param pop A data frame with id, group, and number of population.
#' @param pop_intensity data frame with id and population intensity for all groups. Result from \code{\link{popIntensity}}.
#' @return array with dissimilarity for all groups, size of number of localities.
#' @author Beatriz Moura dos Santos
#' @references Feitosa, Camara, Monteiro, Koschitzi & Silva (2007). 
#'    Global and local spatial indices of urban segregation. 
#'    \emph{International Journal of Geographical Information Science}, 21(3), 299-323.
#' @references Sousa (2017). Segregation Metrics. 
#'    \url{https://github.com/sandrofsousa/Resolution/blob/master/Pysegreg/segregationMetrics.py}.

localDissimilarity <- function(pop, pop_intensity){
  
  local_ids <- pop_intensity$id
  
  pop <- pop %>%
    dplyr::filter(.data$id %in% local_ids) %>%
    dplyr::select(-id)
  
  pop_sum <- apply(pop, 1, sum, na.rm=T)
  
  if(nrow(pop_intensity) == 0){
    lj <- as.numeric(pop_sum)
    tjm <- pop / lj
  } else{
    locality <- pop_intensity[,-1]
    lj <- apply(locality, 1, sum)
    tjm <- locality / lj
  }
  
  pop_total <- sum(pop, na.rm=T)
  tm <- apply(pop, 2, sum) / pop_total
  
  index_i <- sum(tm * (1 - tm), na.rm = T)
  
  local_diss <- apply(abs(sweep(tjm, 2, tm, '-')) * t(pop_sum) / (2 * pop_total * index_i), 1, sum, na.rm=T)
  local_diss <- matrix(local_diss, ncol=1)
  local_diss[is.nan(local_diss)] <- 0
  
  return (local_diss)
}

#' Computes the local entropy score for a unit area Ei (diversity).
#'
#' @inheritParams localDissimilarity
#' @return An array with local indices of entropy.
#' @author Beatriz Moura dos Santos
#' @references Iceland (2004. 
#'    The multigroup entropy index (also known as Theil’s H or the information theory index). 
#'    \emph{US Census Bureau}, Retrieved July, 31, 2006.
#' @references Sousa (2017). Segregation Metrics. 
#'    \url{https://github.com/sandrofsousa/Resolution/blob/master/Pysegreg/segregationMetrics.py}.

localEntropy <- function(pop, pop_intensity){
  
  local_ids <- pop_intensity$id
  
  pop <- pop %>%
    dplyr::filter(.data$id %in% local_ids) %>%
    dplyr::select(-id)
  
  pop_sum <- apply(pop, 1, sum, na.rm=T)
  
  if (nrow(pop_intensity) == 0){
    proportion <- pop / as.numeric(pop_sum)
  } else{
    polygon_sum <- apply(pop_intensity[,-1], 1, sum, na.rm=T)
    proportion <- pop_intensity[,-1] / polygon_sum
  }
  entropy <- proportion * log(1 / proportion)
  entropy$id <- local_ids
  
  entropy <- entropy %>%
    tidyr::gather('group', 'ent', 1:ncol(pop)) %>%
    dplyr::mutate(ent = ifelse(is.na(.data$ent) | is.infinite(.data$ent), 0, .data$ent)) %>%
    tidyr::spread(.data$group, .data$ent) %>%
    dplyr::select(-id)
  
  entropy <- apply(entropy, 1, sum, na.rm=T)
  entropy <- matrix(entropy, ncol=1)
  
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
  
  local_ids <- pop_intensity$id
  
  pop <- pop %>%
    dplyr::filter(.data$id %in% local_ids) %>%
    dplyr::select(-id)
  
  pop_sum <- apply(pop, 1, sum, na.rm=T)
  
  m <- ncol(pop)
  j <- nrow(pop)
  
  local_expo <- sweep(pop, 2, apply(pop, 2, sum, na.rm=T), '/')
  local_expo <- dplyr::bind_cols(id = local_ids, local_expo) %>%
    tidyr::gather('group_exp', 'pop', (1:m)+1)
  
  locality <- pop_intensity[,-1]
  if (nrow(pop_intensity) == 0){
    locality_rate <- pop / pop_sum
  } else{
    locality_rate <- locality / apply(locality,1,sum, na.rm=T)
  }
  locality_rate <- dplyr::bind_cols(id = local_ids, locality_rate)
  
  exposure_rs <- local_expo %>%
    dplyr::left_join(locality_rate, by='id') %>%
    tidyr::gather('group_rt', 'rate', (4:(3+m))) %>%
    dplyr::mutate(exposure = .data$pop*.data$rate,
                  exposure = ifelse(is.nan(.data$exposure), 0, .data$exposure),
                  m_n = paste0(.data$group_exp, '_', .data$group_rt)) %>%
    dplyr::select(.data$id, .data$exposure, .data$m_n) %>%
    tidyr::spread(.data$m_n, .data$exposure)
  
  exposure_rs <- exposure_rs[,-1]
  
  return(exposure_rs)
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
  
  local_ids <- pop_intensity$id
  
  pop <- pop %>%
    dplyr::filter(.data$id %in% local_ids) %>%
    dplyr::select(-id)
  
  pop_sum <- apply(pop, 1, sum, na.rm=T)
  
  et <- global_entropy * sum(pop_sum, na.rm=T)
  eei <- global_entropy - local_entropy
  h_local <- pop_sum * eei / et
  
  return(h_local)
}

#' Summary of local segregation indexes
#'
#' @inheritParams localDissimilarity
#' @return data frame with calculated indexes for each tract with the following columns
#' \describe{
#'   \item{id} {tract id}
#'   \item{diss} {dissimilarity index}
#'   \item{ent} {entropy}
#'   \item{iH} {index H}
#'   \item{exp_m_n} {exposure inde for group m with group n}
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
  
  results <- data.frame(
    id = pop_intensity$id,
    diss = localDissimilarity(pop, pop_intensity),
    ent = localEntropy(pop, pop_intensity),
    iH = localIndexH(pop, pop_intensity)
  )
  
  results <- cbind(results, localExposure(pop, pop_intensity))
  
  return(results)
}
