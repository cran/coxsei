\name{CumInt}
\alias{CumInt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  %%  ~~function to do ... ~~
  Cumulative intensity function
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  Calculate the cumulative/integrated hazard/intensity function
}
\usage{
CumInt(x, int, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    %%     ~~Describe \code{x} here~~
    the value at which to calculate the cumulative function value
}
  \item{int}{
    %%     ~~Describe \code{int} here~~
    the intensity/hazard rate function. Has to be vectorized.
}
  \item{\dots}{
    %%     ~~Describe \code{\dots} here~~
    the arguments to be passed in to control the behavior of the
    underlying \code{integrate} function. 
}
}
\details{
  %%  ~~ If necessary, more details than the description above ~~
  Relies on the numerical integration routine of R.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
 The value(s) of the cumulative hazard function at the specified
 \code{x} value(s).
}
\author{
  %%  ~~who you are~~
  Feng Chen <feng.chen@unsw.edu.au>
}
\note{
  %%  ~~further notes~~
  Not intended to be called by the user directly.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\section{Warning}{
  The validity of the user supplied intensity function is not checked.   
}

\examples{

curve(CumInt(x,int=function(y)1*( y>=0 & y<2)+3*(y>=2 & y<3)+1*(y>=3)),
      0,5,xlab="t",ylab="H(t) of a piece-wise constant hazard fun h(t)")   
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~survival }
% __ONLY ONE__ keyword per line
