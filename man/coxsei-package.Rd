\name{coxsei-package}
\alias{coxsei-package}
\docType{package}
\title{
Fit a Cox-type self-exciting intensity model (CoxSEI) to right-censored
counting process data
}
\description{
Fit the CoxSEI model using the partial likelihood method.
}
\details{
To use the package, the data needs to be prepared into a data frame
containing a column named \code{Y} for observed event times in ascending
order of each individual process, a column named \code{delta} indicating
if the event is 'death' (1) or 'censoring' (0), a column named
\code{id} indicating the process id of each event time, and one or more
columns giving the value of any covariate variable at the observed event
times of each process. Then call the \code{coxseiest} function or the
identical but much faster \code{coxseiest2} function to estimate the
parametric part of the model and then the \code{coxseiInt} function to
estimate the cumulative baseline intensity function.
}
\author{
Feng Chen <feng.chen@unsw.edu.au>

Maintainer: Feng Chen <feng.chen@unsw.edu.au>
}
\references{
  Feng Chen and Kani Chen. (2014). Modeling Event Clustering Using the m-Memory
  Cox-Type
  Self-Exciting Intensity Model. International Journal of Statistics and
  Probability. 3(3): 126-137. doi:10.5539/ijsp.v3n3p126  URL:
  http://dx.doi.org/10.5539/ijsp.v3n3p126

  Feng Chen and Kani Chen. (2014). Case-cohort analysis of clusters of
  recurrent events. 20(1): 1-15. doi: 10.1007/s10985-013-9275-3
}
\keyword{ package }
\keyword{ survival}

