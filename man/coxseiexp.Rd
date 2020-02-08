\name{coxseiexp}
\alias{coxseiexp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  %%  ~~function to do ... ~~
  CoxSEI model with exponential function
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
   fit CoxSEI model using an exponential excitation function
}
\usage{
coxseiexp(Y, delta, id, Z, par.init, m = 2, mit = 1000, tr = TRUE,
          method = "L-BFGS-B",lower=c(rep(-Inf,ncol(Z)),-Inf,0),
          upper=rep(Inf,ncol(Z) + 2),...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
    %%     ~~Describe \code{Y} here~~
    the observed times (including censoring times)
}
  \item{delta}{
    %%     ~~Describe \code{delta} here~~
    indicator of event: 1=event, 0=censoring
}
  \item{id}{
    %%     ~~Describe \code{id} here~~
    the id of the individual/group the event/censoring corresponds to
}
  \item{Z}{
    %%     ~~Describe \code{Z} here~~
    covariate matrix
}
  \item{par.init}{
    %%     ~~Describe \code{par.init} here~~
    initial parameter value to start the iteration
}
  \item{m}{
    %%     ~~Describe \code{m} here~~
    the lag parameter as in M-dependence
}
  \item{mit}{
    %%     ~~Describe \code{mit} here~~
    maximum number of iteration allowed in maximizing the loag partial
    likelihood 
}
  \item{tr}{
    %%     ~~Describe \code{tr} here~~
    should the optimization process be 'tr'aced
}
  \item{method}{
    %%     ~~Describe \code{method} here~~
    method of optimization; defaults to "L-BFGS-B"
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
an object of class "coxsei", basically a list with components
\item{coefficients }{a named vector of coefficients }
\item{vcov}{a symmetric matrix which is supposed to be positive definite
when m>0, or with the (np-2)x(np-2) major submatrix positive definite
when m=0}
}
\author{
  %%  ~~who you are~~
  Feng Chen <feng.chen@unsw.edu.au>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~regression}
\keyword{ ~survival }% __ONLY ONE__ keyword per line
