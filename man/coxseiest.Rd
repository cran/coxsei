\name{coxseiest}
\alias{coxseiest}
\alias{coxseiest2}
\alias{coxseiest3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  %%  ~~function to do ... ~~
  Function to estimate the parametric part of the Cox (proportional
  intensity) self-exciting point process (CoxSEI) model
}
\description{
  Estimate the parametric part of the CoxSEI model using (conditionally)
  right-censored counting process data. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
coxseiest(dat, par.init, m = 2, mit = 1000, tr = TRUE,
          method = "L-BFGS-B", lower=c(rep(-Inf,ncol(dat)-3),-Inf,0),
          upper=rep(Inf,ncol(dat)-3 + 2),
          gfun = function(x, pa) {
           ifelse(x <= 0, rep(0, length(x)), pa[1] * exp(-pa[2] * x))
          })
coxseiest2(dat, par.init, m = 2, mit = 1000, tr = TRUE,
           method = "L-BFGS-B", lower=c(rep(-Inf,ncol(dat)-3),-Inf,0),
           upper=rep(Inf,ncol(dat)-3 + 2),
           gfun = function(x, pa) {
             ifelse(x <= 0, rep(0, length(x)), pa[1] * exp(-pa[2] * x))
           })
coxseiest3(dat, par.init, m = 2, mit = 1000, tr = TRUE,
           method = "L-BFGS-B", lower=c(rep(-Inf,ncol(dat)-3),-Inf,0),
           upper=rep(Inf,ncol(dat)-3 + 2))
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
    %%     ~~Describe \code{dat} here~~
    a data frame with columns \code{Y} containing the censored
    event times of each individual process arranged in ascending order
    with the last time always being the the censoring time, \code{delta}
    containing the event time indicator with value indicator an event
    time and 0 a censoring time, \code{id} specifying the id (process
    number) of each event time recorded, and the others giving the value
    of the associated covariate process at the corresponding event
    times. 
}
  \item{par.init}{
    %%     ~~Describe \code{par.init} here~~
    init guess of the value of the parameters to start the optimization
    iteration with.
}
  \item{m}{
    %%     ~~Describe \code{m} here~~
    order of "autoregression" of the excitation term.
}
  \item{mit}{
    %%     ~~Describe \code{mit} here~~
    maximum number of iteration in the optimization routine
}
  \item{tr}{
    %%     ~~Describe \code{tr} here~~
    if set to \code{TRUE}, print some summary information while the
    optimization routine is running.
}
  \item{method}{
    %%     ~~Describe \code{method} here~~
    method of optimization
  }
  \item{lower}{
    vector of lower boundary values of the parameter space
  }
  \item{upper}{
    vector of upper boundary of the parameter space
  }  
  \item{gfun}{
    %%     ~~Describe \code{gfun} here~~
    the excitation function. Defaults to the exponential decay function 
    \deqn{g(t;\gamma)=\gamma_1 \gamma_2 e^{-\gamma_2 t}}{
      g(t;\gamma)=\gamma_1 \gamma_2 \exp (-\gamma_2 t)
    }
}
}
\details{
  %%  ~~ If necessary, more details than the description above ~~
  \code{coxseiest} uses only R code; \code{coxseiest2} uses external C
  code, and is expected to be 3~4 times fasters than the former;
  \code{coxseiest3} assumes the excitation function is the exponential
  function as defaulted by the former two, and hardwares it in the C
  side of the code, and therefore is much faster than the former two
  when the exponential excitation function is desired. 
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list as that returned by the call to the optimizer routine. For
instance, 
\item{par}{gives the estimate of the parameters}
\item{hessian}{gives the inverse of the estimate of the
  variance-covariance matrix} 
}
\references{
  %% ~put references to the literature/web site here ~
  Feng Chen and Kani Chen. (2014). Modeling Event Clustering Using the m-Memory
  Cox-Type
  Self-Exciting Intensity Model. International Journal of Statistics and
  Probability. 3(3): 126-137. doi:10.5539/ijsp.v3n3p126  URL:
  http://dx.doi.org/10.5539/ijsp.v3n3p126

  Feng Chen and Kani Chen. (2014). Case-cohort analysis of clusters of
  recurrent events. 20(1): 1-15. doi: 10.1007/s10985-013-9275-3
}
\author{
  %%  ~~who you are~~
  Feng Chen <feng.chen@unsw.edu.au>
}
\note{
  %%  ~~further notes~~
  the excitation function has to contain exactly two parameters; a
  feature that does not seem desiable and might change later.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  %% ~~objects to See Also as \code{\link{help}}, ~~~
 See \code{\link{optim}} for the components of the returned value
}
\examples{
data("dat")
## this takes over 15 minutes
##est0 <- coxseiest(dat,par.init=c(0.2,0.4,0.6,0.6,5))
## this one takes about 4 minutes
##est1 <- coxseiest2(dat,par.init=c(0.2,0.4,0.6,0.6,5))
## this one takes about 10 seconds
est2 <- coxseiest3(dat,par.init=c(0.2,0.4,0.6,0.6,5))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~models }
\keyword{ ~regression }% __ONLY ONE__ keyword per line
\keyword{ ~survival}
