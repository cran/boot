\name{melanoma}
\alias{melanoma}
\title{
Survival from Malignant Melanoma
}
\section{SUMMARY}{
The \code{melanoma} data frame has 205 rows and 7 columns.


The data consist of measurements made on patients with malignant melanoma.
Each patient had their tumour removed by surgery at the Department of Plastic
Surgery, University Hospital of Odense, Denmark during the period 1962 to 1977.
The surgery consisted of complete removal of the tumour together with about
2.5cm of the surrounding skin.  Among the measurements taken were the thickness
of the tumour and whether it was ulcerated or not.  These are thought to be
important prognostic variables in that patients with a thick and/or ulcerated
tumour have an increased chance of death from melanoma.
Patients were followed until the end of 1977.  
}
\section{DATA DESCRIPTION}{
This data frame contains the following columns:


}
\arguments{
\item{time}{
Survival time in days since the operation, possibly censored.
}
\item{status}{
The patients status at the end of the study. 1 indicates that they had died
from melanoma, 2 indicates that they were still alive and 3 indicates that they
had died from causes unrelated to their melanoma.
}
\item{sex}{
The patients sex; 1=male, 2=female.
}
\item{age}{
Age in years at the time of the operation.
}
\item{year}{
Year of operation.
}
\item{thickness}{
Tumour thickness in mm.
}
\item{ulcer}{
Indicator of ulceration; 1=present, 0=absent.
}}
\section{NOTE}{
This dataset is not related to the dataset in the \code{trellis} library
with the same name.
}
\section{SOURCE}{
The data were obtained from


Andersen, P.K., Borgan, O., Gill, R.D. and Keiding, N. (1993)
\emph{Statistical Models Based on Counting Processes}. Springer-Verlag.
}
\section{OTHER REFERENCES}{
Davison, A.C. and Hinkley, D.V. (1997) \emph{Bootstrap Methods and Their Application}. Cambridge University Press.


Venables, W.N. and Ripley, B.D. (1994) 
\emph{Modern Applied Statistics with S-Plus}. Springer-Verlag.
}
\keyword{sysdata}
% Converted by Sd2Rd version 0.3-1.
