#ifndef TRAPZRCPP_H
#define TRAPZRCPP_H

#include <Rcpp.h>

using namespace Rcpp;

template <class iter> bool is_sorted (iter begin, iter end);

double trapzRcpp(const Rcpp::NumericVector X, const Rcpp::NumericVector Y);


#endif