\name{rock}
\alias{rock}
\title{
Measurements on Petroleum Rock Samples
}
\description{
The \code{rock} data frame has 48 rows and 4 columns.


Twelve core samples from petroleum reservoirs were sampled by 4
cross-sections.  Each core sample was measured for permeability, and
each cross-section has total area of pores, total perimeter of pores,
and shape.
}
\usage{
data(rock)
}
\format{
This data frame contains the following columns:


\describe{
\item{\code{area}}{
The area of the pores space (in pixels out of 256  x  256).
}
\item{\code{peri}}{
The total perimeter of the pores (in pixels).
}
\item{\code{shape}}{
The shape of the pores space defined as \code{peri/sqrt(area)}.
}
\item{\code{perm}}{
The permeability in milli-Darcies.
}}}
\source{
The data were collected and kindly made available by R. Katz of the Department
of Statistics at the University of Oxford.


Katz, R. (1995) \emph{Spatial Analysis of Pore Images}. D.Phil. thesis,
Department of Statistics, University of Oxford.
}
\references{
Davison, A.C. and Hinkley, D.V. (1997)
\emph{Bootstrap Methods and Their Application}.
Cambridge University Press.


Venables, W.N. and Ripley, B.D. (1994)
\emph{Modern Applied Statistics with S-Plus}.
Springer-Verlag.
}
\keyword{datasets}
% Converted by Sd2Rd version 1.15.
