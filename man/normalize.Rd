\name{normalize}
\alias{normalize}
\title{ Normalize Weights}
\description{
 Normalize weights to sum to one within each stratum.}
\usage{
normalize(wts, strata)
}
\arguments{
  \item{wts}{ A vector of non-negative weights.  At least one of the
weights in each stratum must be strictly positive.}
  \item{strata}{ A vector of the same length as \code{wts} defining
the strata.  It should be a vector of integers or a factor with
integer levels. }
}
\value{
  A vector of the same length as \code{wts} such that the sum of
the output values within each stratum is 1.}
\author{ Angelo J. Canty}
\keyword{ manip } 