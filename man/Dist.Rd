\name{Dist}
\alias{Dist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  %%  ~~function to do ... ~~
  Distribution function
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  Calculate the value at \code{x} of the distribution function associated with
  the intensity/hazard function probived through \code{int}. 
}
\usage{
Dist(x, int, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    %%     ~~Describe \code{x} here~~
    the value to evaluate the distribution function at.
}
  \item{int}{
    %%     ~~Describe \code{int} here~~
    vectorized function specifying the intensity/hazard function
}
  \item{\dots}{
    %%     ~~Describe \code{\dots} here~~
    arguments to be passed to the \code{integrate} function
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A number between 0 and 1 inclusive, that gives the value of the
distribution function at the specified \code{x} value.
}

\author{
  %%  ~~who you are~~
  Feng Chen <feng.chen@unsw.edu.au>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
curve(Dist(x,int=function(x)3*x^2),0,5)
curve(pweibull(x,shape=3),0,5,add=TRUE,col=3,lty=3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~distribution }
\keyword{ ~survival }% __ONLY ONE__ keyword per line
