\name{boot-arrays}
\alias{antithetic.array}
\alias{anti.arr}
\alias{balanced.array}
\alias{extra.array}
\alias{importance.array}
\alias{importance.array.bal}
\alias{index.array}
\alias{ordinary.array}
\alias{permutation.array}

\title{ Bootstrap Arrays}
\description{
Construct bootstrap index arrays
}
\usage{
antithetic.array(n, R, L, strata)
anti.arr(n, R, L, inds=1:n)
balanced.array(n, R, strata)
extra.array(n, R, m, strata=rep(1,n))
importance.array(n, R, weights, strata)
importance.array.bal(n, R, weights, strata)
index.array(n, R, sim, strata=rep(1,n), m=0, L=NULL, weights=NULL)}
ordinary.array(n, R, strata)
permutation.array(n, R, strata)
}
\details{
These are internal functions used to construct bootstrap arrays
which are called by \code{boot}.
They should not be called by the user.}
\seealso{\code{\link{boot}}}
\author{Angelo J. Canty}
\keyword{internal}
