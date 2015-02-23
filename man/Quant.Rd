\name{Quant}
\alias{Quant}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  %%  ~~function to do ... ~~
  Quantile function
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  calculates the value of the quantile function (inverse of the
  distribution function) of the survival variable with given
  intensity/hazard function.
}
\usage{
Quant(p, int, tolerance = .Machine$double.eps, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
    %%     ~~Describe \code{p} here~~
    the (probability) values to calculate the quantiles at
}
  \item{int}{
    %%     ~~Describe \code{int} here~~
    the intensity/hazard function. Has to be vectorized.
}
  \item{tolerance}{
    %%     ~~Describe \code{tolerance} here~~
    tolerated numerical error in inverting the distribution function. 
}
  \item{\dots}{
    %%     ~~Describe \code{\dots} here~~
    arguments to be passed to \code{CumInt} (eventually to \code{integrate})
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
a numerical value or vector giving the values of the quantile function
at \code{x}
}

\author{
  %%  ~~who you are~~
  Feng Chen <feng.chen@unsw.edu.au>
}

\examples{
curve(Quant(x,int=function(x)3*x^2),from=1e-3,to=1 - 1e-3)
curve(qweibull(x,shape=3),col=3,lty=3,add=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~distribution }
\keyword{ ~survival }% __ONLY ONE__ keyword per line
