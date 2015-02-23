#include <R.h>
#include<Rinternals.h>
#include <R_ext/Applic.h>
#include<Rmath.h>
/* typedef void integr_fn(double *x, int n, void *ex); */
/* /\* void Rdqags(integr_fn f, void *ex, double *a, double *b, *\/ */
/* /\*             double *epsabs, double *epsrel, *\/ */
/* /\*             double *result, double *abserr, int *neval, int *ier, *\/ */
/* /\*             int *limit, int *lenw, int *last, *\/ */
/* /\*             int *iwork, double *work); *\/ */
typedef struct int_struct
{
    SEXP f;    /* function */
    SEXP env;  /* where to evaluate the calls */
} int_struct, *IntStruct;
/* This is *the* ``integr_fn f'' used when called from R : */
static void Rhazfn(double *x, int n, void *ex)
{
    SEXP args, resultsxp, tmp;
    int i;
    IntStruct IS = (IntStruct) ex;

    PROTECT(args = allocVector(REALSXP, n));
    for(i = 0; i < n; i++) REAL(args)[i] = x[i];

    PROTECT(tmp = lang2(IS->f , args));
    PROTECT(resultsxp = eval(tmp, IS->env));

    if(length(resultsxp) != n)
	error("evaluation of function gave a result of wrong length");
    if(TYPEOF(resultsxp) == INTSXP) {
	resultsxp = coerceVector(resultsxp, REALSXP);
    } else if(TYPEOF(resultsxp) != REALSXP)
	error("evaluation of function gave a result of wrong type");
    for(i = 0; i < n; i++) {
	x[i] = REAL(resultsxp)[i];
	if(!R_FINITE(x[i]))
	    error("non-finite function value");
    }
    UNPROTECT(3);
    return;
}


SEXP rndhaz(SEXP args){
  int_struct is;
  SEXP ans;
  double lo, up, mid, LO, epsabs, epsrel, HAZ, abserr, *work, abstol;
  int neval, ier, limit, lenw, last, *iwork;

  int i,k,n_rnd;

  args = CDR(args);
  is.f = CAR(args); args = CDR(args);
  n_rnd = asInteger(CAR(args)); args = CDR(args);
  is.env = CAR(args); args = CDR(args);
  abstol = asReal(CAR(args)); args = CDR(args);
  epsabs = asReal(CAR(args)); args = CDR(args);
  epsrel = asReal(CAR(args)); args = CDR(args);
  limit = asInteger(CAR(args)); args = CDR(args);
  lenw = 4 * limit;
  iwork = (int *) R_alloc((size_t) limit, sizeof(int));
  work = (double *) R_alloc((size_t) lenw, sizeof(double));

  PROTECT(ans = allocVector(REALSXP, n_rnd));
  GetRNGstate();
  for(i=0;i<n_rnd;i++){
    REAL(ans)[i]= -log(runif(0.0,1.0));
  }
  PutRNGstate();
  /* for(i=0;i<n_rnd;i++)printf("%f\t",REAL(ans)[i]); */
  /* printf("\n"); */
  LO=0.0;
  for(i=0; i<n_rnd; i++){
    lo=0.0; k=0; up=pow(2,k); 
    Rdqags(Rhazfn, (void*)&is,
	   &LO, &up, &epsabs, &epsrel, &HAZ,
	   &abserr, &neval, &ier, &limit, &lenw, &last, iwork, work);
    /* make sure that integrate(intfn,0,up)>re[i] */
    while(HAZ<REAL(ans)[i]){
      k+=1;
      up+=pow(2.0,k);      
      Rdqags(Rhazfn, (void*)&is,
	     &LO, &up, &epsabs, &epsrel, &HAZ,
	     &abserr, &neval, &ier, &limit, &lenw, &last, iwork, work);
    }
    /* printf("up=%f;HAZ(up)=%f\n",up,HAZ);     */
    /* now do the bisecting search: */
    while(fabs(up-lo)>abstol){
      mid=(lo+up)/2.0;
      Rdqags(Rhazfn, (void*)&is,
	     &LO, &mid, &epsabs, &epsrel, &HAZ,
	     &abserr, &neval, &ier, &limit, &lenw, &last, iwork, work);
      if(HAZ>=REAL(ans)[i]) up=mid; else lo=mid;
      /* printf("up=%f,lo=%f,mid=%f\n",up,lo,mid); */
    }
    REAL(ans)[i]=(up+lo)/2.0;
  }
  UNPROTECT(1);
  return(ans);
}
