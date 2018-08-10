// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// pava_main
List pava_main(Eigen::VectorXd& y_, int nblocks, std::vector<int> sidx, std::vector<int> size, bool warmstart);
RcppExport SEXP _simPU_pava_main(SEXP y_SEXP, SEXP nblocksSEXP, SEXP sidxSEXP, SEXP sizeSEXP, SEXP warmstartSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type y_(y_SEXP);
    Rcpp::traits::input_parameter< int >::type nblocks(nblocksSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type sidx(sidxSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< bool >::type warmstart(warmstartSEXP);
    rcpp_result_gen = Rcpp::wrap(pava_main(y_, nblocks, sidx, size, warmstart));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_simPU_pava_main", (DL_FUNC) &_simPU_pava_main, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_simPU(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}