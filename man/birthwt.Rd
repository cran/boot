\name{birthwt}
\alias{birthwt}
\title{
Risk Factors for Low Birth Weight
}
\description{
The \code{birthwt} data frame has 189 rows and 10 columns.


Data was collected at the Baystate Medical Center, Springfield, Massachusetts
during 1986 to attempt to identify which factors contributed to an increased
risk of low birth weight infants.  Information was recorded for 189 women of
whom 59 had low birth weight infants.
}
\usage{
data(birthwt)
}
\format{
This data frame contains the following columns:
\describe{
\item{\code{low}}{
Indicator of low birth weight; \code{1} indicates low birth weight (< 2.5kg), \code{0}
indicates normal birth weight.
}
\item{\code{age}}{
Mother's age in years at the time of birth.
}
\item{\code{lwt}}{
The mother's weight (in pounds) at the time of the last menstrual period.
}
\item{\code{race}}{
A factor for the mother's race. Levels are \code{1} (white), \code{2} (black) and 
\code{3} (other).
}
\item{\code{smoke}}{
Indicator of mother's smoking status; \code{1}=smoker, \code{0}=non-smoker.
}
\item{\code{ptl}}{
Number of previous premature labours.
}
\item{\code{ht}}{
Indicator of hypertension history; \code{1}=yes, \code{0}=no.
}
\item{\code{ui}}{
Indicator of the presence of uterine irritability; \code{1}=yes, \code{0}=no.
}
\item{\code{ftv}}{
Number of physician visits during the first trimester
}
\item{\code{bwt}}{
Infant's birth weight in grams
}}}
\source{
The data were obtained from


Hosmer, D.W. and Lemeshow, S. (1989) \emph{Applied Logistic Regression}.
John Wiley.
}
\references{
Davison, A.C. and Hinkley, D.V. (1997) 
\emph{Bootstrap Methods and Their Application}. Cambridge University Press.


Venables, W.N. and Ripley, B.D. (1994) 
\emph{Modern Applied Statistics with S-Plus}. Springer-Verlag.
}
\keyword{datasets}
% Converted by Sd2Rd version 1.15.
