\name{channing}
\alias{channing}
\title{
Channing House Data
}
\section{SUMMARY}{
The \code{channing} data frame has 462 rows and 5 columns.


Channing House is a retirement centre in Palo Alto, California.  These data 
were
collected between the opening of the house in 1964 until July 1, 1975.  In that
time 97 men and 365 women passed through the centre.  For each of these, their
age on entry and also on leaving or death was recorded.  A large number of the
observations were censored mainly due to the resident being alive on July 1,
1975 when the data was collected.  Over the time of the study 130 women and
46 men died at Channing House.  Differences between the survival of the sexes,
taking age into account, was one of the primary concerns of this study.
}
\section{DATA DESCRIPTION}{
This data frame contains the following columns:


}
\arguments{
\item{sex}{
A factor for the sex of each resident (\code{"Male"} or \code{"Female"}).
}
\item{entry}{
The residents age (in months) on entry to the centre
}
\item{exit}{
The age (in months) of the resident on death, leaving the centre or July 1, 
1975 whichever event occurred first.
}
\item{time}{
The length of time (in months) that the resident spent at Channing House.
(\code{time=exit-entry})
}
\item{cens}{
The indicator of right censoring.  1 indicates that the resident died at
Channing House, 0 indicates that they left the house prior to July 1, 1975 or
that they were still alive and living in the centre at that date.
}}
\section{SOURCE}{
The data were obtained from


Hyde, J. (1980) Testing survival with incomplete observations. 
\emph{Biostatistics Casebook}.  
R.G. Miller, B. Efron, B.W. Brown and L.E. Moses (editors), 
31-46. John Wiley.
}
\section{OTHER REFERENCES}{
Davison, A.C. and Hinkley, D.V. (1997) 
\emph{Bootstrap Methods and Their Application}. Cambridge University Press.
}
\keyword{sysdata}
% Converted by Sd2Rd version 0.3-1.
