\name{RND}
\alias{RND}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  %%  ~~function to do ... ~~
  Random number generator
}
\description{
  %% ~~ A concise (1-5 lines) description of what the function does. ~~
  RND takes a vectorized positive R function defined on positive reals
  and returns a vector of \code{n} values of the random variable
  (survival time) with the specifed function as its hazard/intensity
  rate function. 
}
\usage{
RND(n, int, tol = .Machine$double.eps^0.5, epsabs = 1e-10, epsrel =
1e-10, limit = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
    %%     ~~Describe \code{n} here~~
    number of observations.
  }
  \item{int}{
    %%     ~~Describe \code{int} here~~
    hazard rate function of the survival variable, or the intensity
    function of the one-event point process counting the number (0 or
    1) of deaths by following a sample of the surviving suject.
  }
  \item{tol}{
    %%     ~~Describe \code{tol} here~~
    tolerance of the numerical error in calculating the inverse of the
    cumulative distribution function of the survival variable. Defaults
    to the square root of the machine epsilon.
  }
  \item{epsabs}{
    maximum absolute error to be tolerated by the integrator.
  }
  \item{epsrel}{
    maximum relative error to be tolerated by the integrator.
  }
  \item{limit}{
    maximum number of iterations permitted by the integrator.
  }
%%   \item{\dots}{
%%     %%     ~~Describe \code{\dots} here~~
%%     other variables to be passed to the underlying functions such as the
%%     quantile functions etc. 
%% }
}

\value{
  %%  ~Describe the value returned
  %%  If it is a LIST, use
  %%  \item{comp1 }{Description of 'comp1'}
  %%  \item{comp2 }{Description of 'comp2'}
  %% ...
  a vector of \code{n} observations of the survival variable with the
  supplied intensity/hazard function.
}

\author{
  %%  ~~who you are~~
  Feng Chen <feng.chen@unsw.edu.au>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
set.seed(1)
dat <- RND(100,int=function(x)x^2)
ks.test(dat,pweibull,shape=3) # p-value = 0.6058
qqplot(dat,rweibull(100,shape=3))  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~datagen }
\keyword{ ~distribution }% __ONLY ONE__ keyword per line
\keyword{ ~survival }
