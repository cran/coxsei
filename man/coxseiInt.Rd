\name{coxseiInt}
\alias{coxseiInt}
% \alias{coxseiInt2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  %%  ~~function to do ... ~~
  Calculate the estimator of the cumulative baseline intensity function
  in the CoxSEI model.
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  It takes the paramter of the parametric part (or its theorized value)
  and calculate the values of the estimator at the jump times; it also
  gives the values of the estimator for the variance of the intensity
  estimator. 
}
\usage{
coxseiInt(dat, parest, hessian=NULL, vcovmat=solve(hessian), m = 2,
          gfun = function(x, pa) {
            ifelse(x <= 0, 0, pa[1] * pa[2] * exp(-pa[2] * x))
          },
          gfungrd = function(x, pa){
            if(length(x)==0)return(matrix(0,2,0));
            rbind(pa[2]*exp(-pa[2]*x),
                  pa[1]*exp(-pa[2]*x)*(1-pa[2]*x)
                 )
          })
}
% coxseiInt2(dat, parest, hessian=NULL, vcovmat=solve(hessian), m = 2,
%            gfun = function(x, pa) {
%              ifelse(x <= 0, 0, pa[1] * pa[2] * exp(-pa[2] * x))
%            },
%            gfungrd = function(x, pa){
%              if(length(x)==0)return(matrix(0,2,0));
%              rbind(pa[2]*exp(-pa[2]*x),
%                    pa[1]*exp(-pa[2]*x)*(1-pa[2]*x)
%                   )
% 	   })
% }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
    %%     ~~Describe \code{dat} here~~
    a data frame containing the right-censored counting process data 
  }
  \item{parest}{
    %%     ~~Describe \code{param} here~~
    the estimate of parameter of the parametric part of the CoxSEI model
  }
  \item{hessian}{
    the hessian matrix returned by the optimization procedure in the
    estimation of the parametric part based on partial likelihood
  }
  \item{vcovmat}{
    the variance-covariance matrix of the estimator of the the
    parametric components; defaulted to the inverse of the hessian
    matrix
  }    
  \item{m}{
    %%     ~~Describe \code{m} here~~
    autoregressive order in the excitation part of the intensity
  }
  \item{gfun}{
    %%     ~~Describe \code{gfun} here~~
    the excitation function; defaults to the exponential decay function
  }
  \item{gfungrd}{
    derivative/gradient function of the excitation function
  }
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
a list giving the jump times and values at these of the estimator of the
cumulative baseline intensity function.
\item{x}{ the ordered death/event times}
\item{y}{ the value of the estimator of the intensity function at the
  observed death/event times}
\item{varest}{ the value of the estimator of the variance of the
  estimator of the intensity function, at the jump times
}
The step function can be obtained using \code{stepfun}, and plotted by setting
\code{type="s"} in the \code{plot} function.
}
\author{
  %%  ~~who you are~~
  Feng Chen <feng.chen@unsw.edu.au>  
}
\note{
  %%  ~~further notes~~
  Currently doesn't compute the standard error or variance estimator of
  the baseline cumulative intensity estimator. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data("dat")
est <- coxseiest3(dat,c(0.2,0.4,0.6,log(0.06),log(5)))
pe <- est$par; pe[4:5] <- exp(pe[4:5]);
ve <- diag(pe) \%*\% solve(est$hessian, diag(pe));
cintest <- coxseiInt(dat,pe,vcovmat=ve)
plot(cintest,type="s")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~survival }
\keyword{ ~regression }% __ONLY ONE__ keyword per line
