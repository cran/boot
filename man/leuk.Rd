\name{leuk}
\alias{leuk}
\title{
Survival Times and White Blood Counts for Leukaemia Patients
}
\section{SUMMARY}{
The \code{leuk} data frame has 33 rows and 3 columns.


Survival times are given for 33 patients who died from acute myelogenous 
leukaemia.  Also measured was the patient's white blood cell count at the
time of diagnosis.  The patients were also factored into 2 groups according to
the presence or absence of a morphologic characteristic of white blood cells.  
Patients termed AG positive were identified by the presence of Auer rods 
and/or  significant granulation  of the leukaemic cells in the bone marrow at 
the time of diagnosis.
}
\section{DATA DESCRIPTION}{
This data frame contains the following columns:


}
\arguments{
\item{wbc}{
White blood cell count.
}
\item{ag}{
Factor describing the patients group according to the test described above.
Levels are \code{present} and \code{absent}.
}
\item{time}{
Survival time in weeks after diagnosis.
}}
\section{SOURCE}{
The data were obtained from


Feigl, P. and Zelen, M. (1965) Estimation of exponential survival 
probabilities with concomitant information. \emph{Biometrics}, \bold{21}, 
826-838.
}
\section{OTHER REFERENCES}{
Davison, A.C. and Hinkley, D.V. (1997) 
\emph{Bootstrap Methods and Their Application}. Cambridge University Press.


Venables, W.N. and Ripley, B.D. (1994) 
\emph{Modern Applied Statistics with S-Plus}. Springer-Verlag.
}
\keyword{sysdata}
% Converted by Sd2Rd version 0.3-1.
