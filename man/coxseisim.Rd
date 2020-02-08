\name{coxseisim}
\alias{coxseisim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  %%  ~~function to do ... ~~
  A function to simulate a CoxSEI process conditional on specified
  covariate values 
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  simulate the sample path of the CoxSEI model with given covariate
  process values, and excitation function and order of autodependence in
  the excitation term.
}
\usage{
coxseisim(parreg, parg, lmd0 = function(tt) (1 + 0.5 * cos(2 * pi *
tt)),
          g = function(x, parg) {
                 ifelse(x <= 0, 0, parg[1] * parg[2] * exp(-parg[2] * x))
              },
          censor = 1, m = 2, trace=TRUE,
          Z = function(x) matrix(0, length(x), length(parreg))
         )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{parreg}{
    the regression parameter
%%     ~~Describe \code{parreg} here~~
}
  \item{parg}{
    %%     ~~Describe \code{parg} here~~
    parameters of the excitation function
}
  \item{lmd0}{
    %%     ~~Describe \code{lmd0} here~~
    the baseline intensity function
}
  \item{g}{
    %%     ~~Describe \code{g} here~~
    the excitation function
}
  \item{censor}{
    %%     ~~Describe \code{censor} here~~
    the censoring time
  }
  \item{m}{
    %%     ~~Describe \code{m} here~~
    order of autoregression in the excitation component of the intensity
    process 
  }
  \item{trace}{
    whether to trace the data generation process; defaults to \code{TRUE}
  }
  \item{Z}{
    %%     ~~Describe \code{Z} here~~
    a function to calculate the covariate values at a specified event time
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A data frame with provided covariate values and the censoring time, and
the generated event times. 
}
\author{
  %%  ~~who you are~~
  Feng Chen <feng.chen@unsw.edu.au>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
    n.smp <- 100;
    z <- matrix(NA,n.smp,3)
    for(i in 1:n.smp)
    z[i,] <- round(c(runif(1,0.5,1.5),runif(1,1.5,2.5),rbinom(1,1,0.5)),2)
    dat <- coxseisim(1:3*0.2,c(0.07,10),censor=rlnorm(1,0,0.1),m=2,
    Z=function(x)matrix(z[1,],length(x),3,byrow=TRUE))
    dat$id <- 1;
    for(i in 2:n.smp){
      dattmp <- coxseisim(1:3*0.2,c(0.07,10),censor=rlnorm(1,0,0.1),m=2,
      Z=function(x)matrix(z[i,],length(x),3,byrow=TRUE))
      dattmp$id <- i;
      dat <- rbind(dat,dattmp)
    }

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~datagen }
\keyword{ ~survival }% __ONLY ONE__ keyword per line
