\name{Surv}
\alias{Surv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  %%  ~~function to do ... ~~
  Survival function
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  Evaluate the survival function corresponding to the given
  intensity/hazard function.
}
\usage{
Surv(x, int, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    %%     ~~Describe \code{x} here~~
    value to calculate the value of the survival function for
}
  \item{int}{
    %%     ~~Describe \code{int} here~~
    the intensity/hazard function
}
  \item{\dots}{
    %%     ~~Describe \code{\dots} here~~
    further arguments to be passed to \code{CumInt}
}
}

\value{
  a numerical value or vector giving the value(s) of the survival
  function at \code{x}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\author{
  %%  ~~who you are~~
  Feng Chen <feng.chen@unsw.edu.au>
}

\examples{
curve(Surv(x, int=function(x)3*x^2), from=0, to=5)
curve(pweibull(x,shape=3,lower=FALSE), add=TRUE, col=2, lty=3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~distribution }
\keyword{ ~survival }% __ONLY ONE__ keyword per line
