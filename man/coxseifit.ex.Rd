\name{coxseifit.ex}
\alias{coxseifit.ex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  %%  ~~function to do ... ~~
  CoxSEI model with exponential function
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  Fit a CoxSEI model with exponential function to right censored
  counting process data
}
\usage{
coxseifit.ex(dat, par.init, m = 2, mit = 1000, tr = TRUE,
             method = "L-BFGS-B",lower=c(rep(-Inf,ncol(dat)-3),-Inf,0),
             upper=rep(Inf,ncol(dat)-3 + 2),...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
    %%     ~~Describe \code{dat} here~~
    The data
}
  \item{par.init}{
    %%     ~~Describe \code{par.init} here~~
    initial value of the regression coefficients and coefficients in the
    excitation function
}
  \item{m}{
    %%     ~~Describe \code{m} here~~
    the lag parameter (the m-dependence parameter)
}
  \item{mit}{
    %%     ~~Describe \code{mit} here~~
    maximum number of iterations allowed in the optimizer
}
  \item{tr}{
    %%     ~~Describe \code{tr} here~~
    whether to trace the optimization or not
}
  \item{method}{
    %%     ~~Describe \code{method} here~~
    the method of optimization used by the \code{optim} routine
  }
  \item{lower}{
    vector of lower boundary values of the parameter space
  }
  \item{upper}{
    vector of upper boundary of the parameter space
  }
  \item{\dots}{
    other arguments to be passed to the optimization routine
  }
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list of some components with kind of self-evident meanings by their name
}
\author{
  %%  ~~who you are~~
  Feng Chen <feng.chen@unsw.edu.au>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  %% ~~objects to See Also as \code{\link{help}}, ~~~
  \code{\link{coxseiest}}, \code{\link{coxseiInt}}
}
\examples{
data("dat")
csfit <- coxseifit.ex(dat,c(1:3*0.2,0.7,10))
coef(csfit)
plot(csfit$cintfn,do.points=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
\keyword{survival}% __ONLY ONE__ keyword per line
