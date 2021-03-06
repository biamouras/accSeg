% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/segreg_global.R
\name{globalDissimilarity}
\alias{globalDissimilarity}
\title{Computes the global dissimilarity}
\usage{
globalDissimilarity(pop, pop_intensity)
}
\arguments{
\item{pop}{A data frame with id, group, and number of population.}

\item{pop_intensity}{A data frame with id and population intensity for all groups. Result from \code{\link{popIntensity}}.}
}
\value{
Display global value of dissimilarity.
}
\description{
This function calls local dissimilarity function and computes the sum from individual values.
}
\references{
Feitosa, Camara, Monteiro, Koschitzi & Silva (2007). 
   Global and local spatial indices of urban segregation. 
   \emph{International Journal of Geographical Information Science}, 21(3), 299-323.

Sousa (2017). Segregation Metrics. 
   \url{https://github.com/sandrofsousa/Resolution/blob/master/Pysegreg/segregationMetrics.py}.
}
\author{
Beatriz Moura dos Santos
}
