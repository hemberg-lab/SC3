// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// consmx
arma::mat consmx(const arma::mat dat, int K);
RcppExport SEXP _SC3_consmx(SEXP datSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type dat(datSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(consmx(dat, K));
    return rcpp_result_gen;
END_RCPP
}
// norm_laplacian
arma::mat norm_laplacian(arma::mat A);
RcppExport SEXP _SC3_norm_laplacian(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(norm_laplacian(A));
    return rcpp_result_gen;
END_RCPP
}
// tmult
arma::mat tmult(arma::mat x);
RcppExport SEXP _SC3_tmult(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(tmult(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SC3_consmx", (DL_FUNC) &_SC3_consmx, 2},
    {"_SC3_norm_laplacian", (DL_FUNC) &_SC3_norm_laplacian, 1},
    {"_SC3_tmult", (DL_FUNC) &_SC3_tmult, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_SC3(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
