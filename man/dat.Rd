\name{dat}
\alias{dat}
\docType{data}
\title{
  %%   ~~ data name/kind ... ~~
  A simulated data set from a CoxSEI model
}
\description{
  %%  ~~ A concise (1-5 lines) description of the dataset. ~~
  Simulated from a CoxSEI model with an exponential excitation function
  and an AR order 2 for the self-excitation effects. Generated using the
  following code:
  \code{
    set.seed(1);
    n.smp <- 50;
    z <- matrix(NA,n.smp,3);
    for(i in 1:n.smp)
    z[i,] <- round(c(runif(1,0.5,1.5),runif(1,1.5,2.5),rbinom(1,1,0.5)),2);
    dat <- coxseisim(1:3*0.2,c(0.07,10),censor=rlnorm(1,0,0.1),m=2,
    Z=function(x)matrix(z[1,],length(x),3,byrow=T));
    dat$id <- 1;
    for(i in 2:n.smp){
      dattmp <- coxseisim(1:3*0.2,c(0.07,10),censor=rlnorm(1,0,0.1),m=2,
      Z=function(x)matrix(z[i,],length(x),3,byrow=T))
      dattmp$id <- i;
      dat <- rbind(dat,dattmp)
    }
  }
}
\usage{data(dat)}
\format{
  A data frame with 307 observations on the following 6 variables.
  \describe{
    \item{\code{Y}}{a numeric vector}
    \item{\code{delta}}{a numeric vector}
    \item{\code{Z.1}}{a numeric vector}
    \item{\code{Z.2}}{a numeric vector}
    \item{\code{Z.3}}{a numeric vector}
    \item{\code{id}}{a numeric vector}
  }
}
\examples{
data(dat)
## maybe str(dat) ; plot(dat) ...
}
\keyword{datasets}
