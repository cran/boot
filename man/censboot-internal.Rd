\name{censboot-internal}
\alias{cens.case}
\alias{cens.resamp}
\alias{cens.weird}
\title{ Censored Bootstrap Functions}
\description{Internal fucntions for use with \code{censboot}.}
\usage{
cens.case(n, strata, R)
cens.resamp(data, R, F.surv, G.surv, strata, index = c(1, 2),
            cox = NULL, sim = "model")
cens.weird(data, surv, strata)
}
\details{
These are internal functions used by \code{censboot} to construct
the bootstrap resamples.  They should not be called directly by
the user.
}
\author{Angelo J. Canty}
\seealso{ \code{\link{censboot}}}
\keyword{internal}
