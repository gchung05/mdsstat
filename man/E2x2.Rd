% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/additional.R
\name{E2x2}
\alias{E2x2}
\title{Calculate expected for 2x2 table of observed
Returns a vector of expected in the same order as input}
\usage{
E2x2(x)
}
\arguments{
\item{x}{Vector of length 4 indicating the values in the 4 cells of the 2x2
table in the order A, B, C, D (across, then down)}
}
\value{
Expected values in a vector of length 4, same order as \code{x}
}
\description{
Calculate expected for 2x2 table of observed
Returns a vector of expected in the same order as input
}
\keyword{internal}
