\name{boot-internal}
\alias{isMatrix}
\alias{norm.inter}
\alias{simplex1}
\title{Internal boot functions}
\description{
  Internal \code{boot} functions.
}
\usage{
isMatrix(x)
norm.inter(t, alpha)  
simplex1(a, A, b, init, basic, val = 0, stage = 2, n1 = N, 
         eps = 1e-10, n.iter = n1) 
}
\details{These are internal functions used in the package 
\code{boot}.  They should not be called directly by the user.}
\author{Angelo J. Canty}
\keyword{internal}

