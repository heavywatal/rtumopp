#include <Rcpp.h>

// [[Rcpp::export]]
double nth_element(const Rcpp::NumericVector& x, int n) {
    Rcpp::NumericVector copy(Rcpp::clone(x));
    auto it = copy.begin() + n;
    std::nth_element(copy.begin(), it, copy.end());
    return *it;
}
