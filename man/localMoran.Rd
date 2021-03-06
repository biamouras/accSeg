% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/moran.R
\name{localMoran}
\alias{localMoran}
\title{Processes the local Moran's I bivariate}
\usage{
localMoran(x = NULL, y = NULL, W = NULL)
}
\arguments{
\item{x}{numeric vector with values of variable 1. Defaults to NULL.}

\item{y}{numeric vector with values of variable 2. Defaults to NULL (not bivariate).}

\item{W}{numeric matrix of neighbourhood matrix from \code{\link{nbMatrix}}. Defaults to NULL.}
}
\value{
the local values of Moran's I
}
\description{
Processes the local Moran's I bivariate
}
