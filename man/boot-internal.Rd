\name{boot-internal}
\alias{isMatrix}
\alias{const}
\alias{iden}
\alias{zero}
\alias{simplex1}
\alias{scramble}
\alias{normalize}
\alias{norm.inter}
\alias{make.ends}
\alias{index.array}
\alias{antithetic.array}
\alias{anti.arr}
\alias{balanced.array}
\alias{extra.array}
\alias{ordinary.array}
\alias{permutation.array}
\alias{importance.array}
\alias{importance.array.bal}
\title{Internal boot functions}
\description{
  Internal boot functions.
}
\usage{
const(w, eps = 1e-08) 
isMatrix(x)
iden(n)
zero(n, m)
simplex1(a, A, b, init, basic, val = 0, stage = 2, n1 = N, eps = 1e-10, 
    n.iter = n1)
scramble(ts, norm = TRUE)
normalize(wts, strata) 
norm.inter(t, alpha) 
make.ends(a, n)

index.array(n, R, sim, strata=rep(1,n), m=0, L=NULL, weights=NULL)}
antithetic.array(n, R, L, strata)
anti.arr(n, R, L, inds=1:n)
balanced.array(n, R, strata)
extra.array(n, R, m, strata=rep(1,n))
ordinary.array(n, R, strata) 
permutation.array(n, R, strata)
importance.array(n, R, weights, strata)
importance.array.bal(n, R, weights, strata)
}
\details{
  These are not to be called by the user.
}
\keyword{internal}
