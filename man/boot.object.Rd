\name{boot.object}
\alias{boot.object}
\title{
Bootstrap Output Objects
}
\description{
Class of objects that result from running a bootstrap.
}
\section{Generation}{
These classes of objects are returned from calls to the functions \code{boot},
\code{censboot}, \code{tilt.boot} or \code{tsboot}.
}
\section{Methods}{
The class \code{"boot"} has methods for \code{plot} and \code{print}.
}
\section{Structure}{
Objects of class \code{"boot"} are implemented as a list.  The exact structure of
this list depends on the generating function and can be found in the
relevant helpfiles for those functions. 
}
\seealso{
\code{\link{boot}}, \code{\link{censboot}}, \code{\link{plot.boot}}, \code{\link{print.boot}}, \code{\link{tilt.boot}}, \code{\link{tsboot}}
}
\keyword{htest}
\keyword{nonparametric}
\keyword{methods}
% Converted by Sd2Rd version 1.15.
