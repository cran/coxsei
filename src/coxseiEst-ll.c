#include<R.h>
#include<Rinternals.h>
/* Arguments to pass in: */
/* Y: the sequence of event/censoring times (locally in ascending order) */
/* Z: the vectorized matrix of covariates */
/* Zs: the vectorized array of covariates */
/* cens: the censoring variable values */
/* gs: the group sizes */
/* ng: the number of groups */
/* gofst: the offset needed to locate the groups */
/* par: the parameter value at which to evel likelihood */
/* npex: no. of parameters of the exc. fn. */
/* exf: the exc fn. */
/* rho: the environment to evaluate the ex. fn. */

SEXP ll(SEXP Y, SEXP Z, SEXP Zs, SEXP cens, SEXP gs, 
	SEXP gofst, SEXP par, SEXP exf, SEXP npex,
	SEXP m, SEXP rho) 
{
  SEXP res, R_fcall, parex, argex;
  int i,j,l,k,posij, np=length(par), n=length(Y);
  int *gsp=INTEGER(gs), ng=length(gs), npexv=INTEGER(npex)[0];
  int *gofstp=INTEGER(gofst), N;
  /* declare temporary value holder and the convenience pointers  */
  double tmp, tmp1, *resp, *Yp=REAL(Y), *Zp=REAL(Z), *Zsp=REAL(Zs),
    *parp=REAL(par), *censp=REAL(cens);

  if(np<=npexv)error("length of np not bigger than length of npex");
  if(n!=length(Z)/(np-npexv))error("length of Y not equal to nrow of Z");
		    
  /* allocate memory for the result to be returned andy protect it */
  PROTECT(res=allocVector(REALSXP,1)); 
  /* create the R function call to be evaluated */ 
  PROTECT(R_fcall=lang3(exf, R_NilValue, R_NilValue)); 
  /* prepare the second argument (or parameter/coefficient) to the
   excitation function to be passed over later*/
  PROTECT(parex=allocVector(REALSXP,npexv));
  for(i=0;i<npexv;i++){
    REAL(parex)[npexv-i-1]=parp[np-i-1];/* exp(parp[np-i-1]); */
  }
  /* set the second argument of excitation function to parex*/
  SETCADDR(R_fcall,parex);

  /* prepare the (first) argument of the excitation function*/
  PROTECT(argex=allocVector(REALSXP,1));/* value to be set later */

  resp=REAL(res);

  resp[0]=0.0;
  for(i=0; i< ng; i++){
    if(gsp[i]==1)continue;
    
    for(j=0; j< gsp[i]-1; j++){
      posij= gofstp[i]+j;
      for(l=0;l<np-npexv;l++)
	resp[0] += Zp[posij + n*l]*parp[l];
      for(l=1;l <= INTEGER(m)[0] && l <= j+1 - 1; l++){
	REAL(argex)[0]=Yp[posij]-Yp[posij-l];
	SETCADR(R_fcall,argex);
	resp[0] += REAL(eval(R_fcall,rho))[0];
      }
      tmp = 0.0;
      for(l=0; l< ng; l++){
	if(Yp[posij] > censp[l]) continue;
	/* now calculate the exp{...} terms*/
	tmp1 = 0.0;
	for(k=0; k<np-npexv;k++)
	  tmp1 += Zsp[posij+k*n + l*n*(np-npexv)] * parp[k];
	N = gsp[l]-1;/*calculate the number of events of N_l(.) that are
		       before time Y_{ij}*/
	while(N > 0 && Yp[gofstp[l] +N -1]>=Yp[posij]) N--;
	for(k=1; k<=INTEGER(m)[0] && k<=N; k++){
	  REAL(argex)[0] = Yp[posij]-Yp[gofstp[l]+N-1 - (k-1)]; 
	  SETCADR(R_fcall,argex);
	  tmp1 += REAL(eval(R_fcall,rho))[0];
	}
	tmp += exp(tmp1);
      }
      resp[0] -= log(tmp);
    }
  }
  resp[0] = -resp[0];
  UNPROTECT(4);
  return(res);
}
  
