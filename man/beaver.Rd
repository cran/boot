\name{beaver}
\alias{beaver}
\title{
Beaver Body Temperature Data
}
\section{SUMMARY}{
The \code{beaver} data frame has 100 rows and 4 columns.  It is a multivariate
time series of class \code{"rts"} and also inherits from class \code{"data.frame"}.


This data set is part of a long study into body temperature regulation in
beavers.  Four adult female beavers were live-trapped and had a 
temperature-sensitive radiotransmitter surgically implanted.  Readings were
taken every 10 minutes.  The location of the beaver was also recorded and
her activity level was dichotomized by whether she was in the retreat or 
outside of it since high-intensity activities only occur outside of the 
retreat.


The data in this data frame are those readings for one of the beavers on a day
in autumn.  
}
\section{DATA DESCRIPTION}{
This data frame contains the following columns:


}
\arguments{
\item{day}{
The day number.  The data includes only data from day 307 and early 308.
}
\item{time}{
The time of day formatted as hour-minute.
}
\item{temp}{
The body temperature in degrees Celsius.
}
\item{activ}{
The dichotomized activity indicator. \code{1} indicates that the beaver is outside
of the retreat and therefore engaged in high-intensity activity.
}}
\section{SOURCE}{
The data were obtained from


Reynolds, P.S. (1994) Time-series analyses of beaver body temperatures.
In \emph{Case Studies in Biometry}. N. Lange, L. Ryan, L. Billard, 
D. Brillinger, L. Conquest and J. Greenhouse (editors), 211-228. John Wiley.
}
\section{OTHER REFERENCES}{
Davison, A.C. and Hinkley, D.V. (1997) 
\emph{Bootstrap Methods and Their Application}. Cambridge University Press.
}
\keyword{datasets}
% Converted by Sd2Rd version 0.3-1.
