// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// run_tumopp
std::string run_tumopp(Rcpp::CharacterVector args);
RcppExport SEXP tumorr_run_tumopp(SEXP argsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type args(argsSEXP);
    __result = Rcpp::wrap(run_tumopp(args));
    return __result;
END_RCPP
}
