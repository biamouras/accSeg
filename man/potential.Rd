% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/accessibility_functions.R
\name{potential}
\alias{potential}
\title{Processes the potential (or gravity) accessibility.}
\usage{
potential(
  oportunities,
  matrix,
  type = "exponential",
  betas = NULL,
  cutoff = NULL
)
}
\arguments{
\item{oportunities}{data frame with id and oportunity structure.}

\item{matrix}{data frame with origin, destination and travel_time/distance structure.}

\item{type}{string to choose between \emph{exponential} and \emph{gaussian} function to decay the oportunities. Defaults to exponential.}

\item{betas}{numeric 'willing' decay factor. It must be defined when type "exponential" is chosen. Defaults to NULL.}

\item{cutoff}{numeric threshold. It must be in the same unit of the matrix and it must be defined when type "gaussian" is chosen. Defaults to NULL.}
}
\value{
data frame with origin and acc. If \emph{gaussian} type is chosen, adds the threshold column.
}
\description{
Uses the negative exponential (\code{exponential type}) or the gaussian (\code{gaussian type}) functions to weight the oportunities.
}
\details{
The negative exponential function is one of the most used impedance functions in recent studies of gravity-based accessibility (Handy & Niemeier, 1997).
    In order to avoid the rapid decline of the negative exponential function, Ingram (1971) proposed the use of the normal (gaussian) function, as it declines 
    smoother than the former.
}
\references{
Hansen (1959). How accessibility shapes land use. \emph{Journal of the American Institute of planners}, 25(2), 73-76.

Ingram (1971). The concept of accessibility: a search for an operational form. \emph{Regional studies}, 5(2), 101-107.

Handy & Niemeier (1997). Measuring accessibility: An exploration of issues and alternatives. \emph{Environment and Planning A: Economy and Space}, 29(7), 1175-1194.
}
\author{
Beatriz Moura dos Santos
}
