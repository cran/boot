\name{bootci.object}
\alias{bootci.object}
\title{
Bootstrap Confidence Interval Objects
}
\description{
Class of objects that result from calculating confidence intervals from
an object of class \code{"boot"}.
}
\section{GENERATION}{
This class of objects is returned from calls to the function \code{boot.ci},
}
\section{METHODS}{
The class \code{"bootci"} has methods for \code{print}.
}
\section{STRUCTURE}{
Objects of class \code{"bootci"} are implemented as a list with the following
components.


\item{R}{
The number of bootstrap replicates on which the intervals were based.
}
\item{t0}{
The observed value of the statistic on the same scale as the intervals.
}
\item{call}{
The call to \code{boot.ci} which generated the object


It will also contain one or more of the following components depending
on the value of \code{type} used in the call to \code{bootci}.


}
\item{normal}{
A matrix of intervals calculated using the normal approximation.  It will
have 3 columns, the first being the level and the other two being the upper
and lower endpoints of the intervals.
}
\item{basic}{
The intervals calculated using the basic bootstrap method.
}
\item{student}{
The intervals calculated using the studentized bootstrap method.
}
\item{percent}{
The intervals calculated using the bootstrap percentile method.
}
\item{bca}{
The intervals calculated using the adjusted bootstrap percentile (BCa) method.


These latter four components will be matrices with 5 columns,  the first column
containing the level, the next two containing the indices of the order 
statistics used in the calculations and the final two the calculated endpoints
themselves
}}
\seealso{
\code{\link{boot.ci}}, \code{\link{print.bootci}}
}
\keyword{htest}
\keyword{nonparametric}
\keyword{methods}
% Converted by Sd2Rd version 0.3-1.
