\name{empinf-internal}
\alias{empinf.reg}
\alias{inf.jack}
\alias{positive.jack}
\alias{usual.jack}
\title{ Internal empinf Functions }
\description{Internal empinf functions.
}
\usage{
empinf.reg(boot.out, t = boot.out$t[, 1])
inf.jack(data, stat, index = 1, strata = rep(1, n), eps = 0.001, ...)
positive.jack(data, stat, stype = "w", index = 1,
              strata = rep(1, n), ...)
usual.jack(data, stat, stype = "w", index = 1,
           strata = rep(1, n), ...)
}
\details{
These functions are used by \code{empinf} to find the empirical
influence values.  They should not be called directly by the user.
}
\seealso{\code{\link{empinf}}}
\author{Angelo J. Canty}
\keyword{internal}
