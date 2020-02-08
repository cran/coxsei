\name{coxsei}
\alias{coxsei}
\alias{coxsei.default}
\alias{print.coxsei}
\alias{plot.coxsei}
\alias{summary.coxsei}
\alias{print.summary.coxsei}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  CoxSEI model
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  Fit a CoxSEI model to counting process data
}
\usage{
coxsei(x,...)
\method{coxsei}{default}(x,y,delta,id,par.init,m=2,mit=1000,tr=TRUE,
       method="L-BFGS-B",lower=c(rep(-Inf,ncol(x)),-Inf,0),
       upper=rep(Inf,ncol(x) + 2),...)
\method{print}{coxsei}(x,...)
\method{plot}{coxsei}(x,...)
\method{summary}{coxsei}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a covariate matrix, or an object of class \code{coxsei}}
  \item{y}{a vector of observed times}
  \item{delta}{
    a vector of event indicators: 1=event, 0=censoring
    %%     ~~Describe \code{delta} here~~
    
  }
  \item{id}{
    %%     ~~Describe \code{id} here~~
    the individual/group id to which the event/censoring time
    correspond 
  }
  \item{par.init}{
    %%     ~~Describe \code{par.init} here~~
    initial parameter guess to start the iteration
  }
  \item{m}{
    %%     ~~Describe \code{m} here~~
    lag parameter as in m-dependence
  }
  \item{mit}{
    %%     ~~Describe \code{mit} here~~
    max number of iteration
  }
  \item{tr}{
    %%     ~~Describe \code{tr} here~~
    whether to trace the optimization or not
  }
  \item{method}{
    %%     ~~Describe \code{method} here~~
    method used in optimization
  }
  \item{lower}{
    the lower bound of the parameter space if the L-BFGS-B
    method of optimization is used.
  }
  \item{upper}{
    the upper bound of the paramter space if the L-BFGS-B methodof
    optimaization is used.
  }
  \item{...}{further arguments to plot.stepfun}
  \item{object}{an object of the class coxsei}
}

\value{
  %%  ~Describe the value returned
  %%  If it is a LIST, use
  %%  \item{comp1 }{Description of 'comp1'}
  %%  \item{comp2 }{Description of 'comp2'}
  %% ...
  an object of class \code{coxsei}, basically a list of the following
  components
  \item{coefficients}{a numeric vector of coefficients}
  \item{vcov}{the variance-covariance matrix}
  \item{zval}{the vector of z-value of the Wald test statistic}
  \item{pval}{the vector of p-values}
  \item{details.par}{a list returned by the\code{optim} routine}
  \item{cintfn}{a step function as the estimated cumulative baseline
    intensity function}
  \item{cintvar}{a step function as the variance of the cumulative
    baseline intensity function estimator}
  \item{details.cint}{a list containing more details about the \code{cint}}
}
\references{
    Feng Chen and Kani Chen. (2014). Modeling Event Clustering Using the m-Memory
  Cox-Type
  Self-Exciting Intensity Model. International Journal of Statistics and
  Probability. 3(3): 126-137. doi:10.5539/ijsp.v3n3p126  URL:
  http://dx.doi.org/10.5539/ijsp.v3n3p126

  Feng Chen and Kani Chen. (2014). Case-cohort analysis of clusters of
  recurrent events. 20(1): 1-15. doi: 10.1007/s10985-013-9275-3
%% ~put references to the literature/web site here ~
}
\author{
  %%  ~~who you are~~
  Feng Chen <feng.chen@unsw.edu.au>
}

\seealso{
  %% ~~objects to See Also as \code{\link{help}}, ~~~
  \code{\link{coxseifit.ex}}
}
\examples{
data(dat,package="coxsei")
acoxsei <- coxsei(dat[,3:5],dat[,1],dat[,2],dat[,6],
                  c(0.2*1:3,log(0.07),log(10)))
summary(acoxsei)
plot(acoxsei,do.points=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression }
\keyword{survival }% __ONLY ONE__ keyword per line
