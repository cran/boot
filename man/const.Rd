\name{const}
\alias{const}
\title{Test for Constant}
\description{
  Function checks if a numeric vector or matrix is constant to
  within a specified tolerance.}
\usage{
const(w, eps=1e-08)
}
\arguments{
  \item{w}{ The object to be tested }
  \item{eps}{ The tolerance required for comparisons. }
}
\details{
  Subtracts the mean from each element and then checks that all
  the remainders are within \code{eps} of 0.
}
\value{
  Logical value
}
\author{ Angelo J. Canty }
\note{
  This function was intended for internal use with the functions
  in the boot package.  It may not be sufficient for other uses.
}
\keyword{misc}

