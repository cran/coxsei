\name{Dens}
\alias{Dens}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  %% ~~function to do ... ~~ A function to calculate the value of
  Density function
}

\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  Evaluate the density function corresponding to the specified
  intensity/hazard function \code{int}.  
}

\usage{
Dens(x, int, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    %%     ~~Describe \code{x} here~~
    the value at which to evaluate the density function
}
  \item{int}{
    %%     ~~Describe \code{int} here~~
    the intensity/hazard function. Has to be vectorized.
}
  \item{\dots}{
    other arguments to be passed to the underlying integrator 
%%     ~~Describe \code{\dots} here~~
}
}

\value{
  A numerical value or vector giving the value(s) of the density function
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
\note{
  %%  ~~further notes~~
  Relies on R's \code{integrate} function
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
set.seed(1); dat <- RND(1000,int=function(x)3*x^2)
hist(dat,freq=FALSE); curve(Dens(x,int=function(x)3*x^2),add=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~distribution }
\keyword{ ~survival }% __ONLY ONE__ keyword per line
