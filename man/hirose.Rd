\name{hirose}
\alias{hirose}
\title{
Failure Time of PET Film
}
\section{SUMMARY}{
The \code{hirose} data frame has 44 rows and 3 columns.


PET film is used in electrical insulation.  In this accelerated life test
the failure times for 44 samples in gas insulated transformers.  4 different
voltage levels were used.
}
\section{DATA DESCRIPTION}{
This data frame contains the following columns:


}
\arguments{
\item{volt}{
The voltage (in kV).
}
\item{time}{
The failure or censoring time in hours.
}
\item{cens}{
The censoring indicator; \code{1} means right-censored data.
}}
\section{SOURCE}{
The data were obtained from


Hirose, H. (1993) Estimation of threshold stress in accelerated life-testing.
\emph{IEEE Transactions on Reliability}, \bold{42}, 650-657.
}
\section{OTHER REFERENCES}{
Davison, A.C. and Hinkley, D.V. (1997) 
\emph{Bootstrap Methods and Their Application}. Cambridge University Press.
}
\keyword{sysdata}
% Converted by Sd2Rd version 0.3-1.
