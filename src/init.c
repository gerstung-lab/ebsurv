#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/*
  The following symbols/expressions for .NAME have been omitted

    _ebsurv_convolute_semiMarkov
    _ebsurv_convolute_Markov

  Most likely possible values need to be added below.
*/

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void agmssurv(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"agmssurv", (DL_FUNC) &agmssurv, 20},
    {NULL, NULL, 0}
};

/* .Call calls */
extern SEXP _ebsurv_convolute_Markov(SEXP, SEXP, SEXP,SEXP);
extern SEXP _ebsurv_convolute_semiMarkov(SEXP, SEXP,SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_ebsurv_convolute_Markov",            (DL_FUNC) &_ebsurv_convolute_Markov,    4},
  {"_ebsurv_convolute_semiMarkov",    (DL_FUNC) &_ebsurv_convolute_semiMarkov,    3},
  {NULL, NULL, 0}
};

void R_init_ebsurv(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
