#include <Rcpp.h>
#include <map>          
#include <string>       
#include <vector>       
#include <algorithm>  

using namespace Rcpp;

using Rcpp::List;
using Rcpp::as;


// Obtain order of vectors (same as `order()` in R base)
// [[Rcpp::export]]
IntegerVector order_(NumericVector x) {
  // if (is_true(any(duplicated(x)))) {
  //   Rf_warning("There are duplicates in 'x'; order not guaranteed to match that of R's base::order");
  // }
  NumericVector sorted = clone(x).sort();
  return match(sorted, x);
}


// [[Rcpp::export]]
double weighted_quantile_cpp(NumericVector dat,
                             NumericVector weights,
                             double p) {
  IntegerVector ord = order_(dat);
  weights = weights[ord - 1];   // start index to 0
  dat = dat[ord - 1];
  double q = sum(weights)*p;
  double m = 0;
  double w = 0;
  for (int i = 0; i < weights.length(); i++) {
    w = w + weights[i];
    if (w >= q) {
      m = dat[i];
      break;
    }
  }
  
  return m;
}


// [[Rcpp::export]]
double weighted_quantile_top_cpp(NumericVector dat,
                                 NumericVector weights,
                                 double p) {
  IntegerVector ord = order_(dat);
  weights = weights[ord - 1];   // start index to 0
  dat = dat[ord - 1];
  double q = sum(weights)*p;
  double m = 0;
  double w = 0;
  for (int i = weights.length()-1; i >= 0; i--) {
    w = w + weights[i];
    if (w >= q) {
      m = dat[i];
      break;
    } 
  }
  
  return m;
}


// [[Rcpp::export]]
NumericVector wrm_fit(NumericVector xdat,
                      NumericVector ydat,
                      double x0,
                      NumericVector weight_vec) {
  // double x0 = xdat[xdat.length()];
  weight_vec = weight_vec / sum(weight_vec);
  xdat = xdat[weight_vec != 0];
  ydat = ydat[weight_vec != 0];
  weight_vec = weight_vec[weight_vec != 0];
  
  NumericVector sm(xdat.length());
  NumericVector x_sub;
  NumericVector y_sub;
  NumericVector arg_1;
  NumericVector arg_2;
  IntegerVector ind = seq(0, xdat.length()-1);
  IntegerVector ind_sub;
  for (int i = 0; i < xdat.length(); i++) {
    ind_sub = ind[ind != i];
    y_sub = ydat[ind_sub];
    x_sub = xdat[ind_sub];
    arg_1 = (ydat[i] - y_sub) / (xdat[i] - x_sub);
    arg_2 = weight_vec[ind_sub];
    sm[i] = (weighted_quantile_top_cpp(arg_1, arg_2, 0.5) + weighted_quantile_cpp(arg_1, arg_2, 0.5)) / 2;
  }
  
  // double x0 = xdat[xdat.length()];
  
  double bm = (weighted_quantile_top_cpp(sm, weight_vec, 0.5) + weighted_quantile_cpp(sm, weight_vec, 0.5)) / 2;
  double am = (weighted_quantile_top_cpp(ydat-bm*(xdat-x0), weight_vec, 0.5) +
               weighted_quantile_cpp(ydat-bm*(xdat-x0), weight_vec, 0.5)) / 2;
  
  NumericVector res(2);
  res[0] = am;
  res[1] = bm;
  
  return res;
}


// [[Rcpp::export]]
NumericVector get_kernel(NumericVector x_i,
                         double x,
                         double h,
                         String method) {
  NumericVector kernel_weights;
  
  if (method == "epanechnikov") {
    kernel_weights = 1/h * ifelse(abs(x_i - x) <= h, 3./4.*(1 - pow((x_i - x)/h, 2)), 0);
  } else if (method == "gauss") {
    kernel_weights = 1/h * 1/sqrt(2.*M_PI)*exp(-1./2.*pow((x_i - x)/h, 2));
  }
  
  return kernel_weights;
}


// [[Rcpp::export]]
List wrm_smooth_cpp(NumericVector x,
                    NumericVector y,
                    double h,
                    NumericVector xgrid,
                    String kernel) {
  int N = xgrid.length();
  NumericVector mu(N);
  NumericVector beta(N);
  
  NumericVector weights;
  NumericVector fitted;
  for (int i = 0; i < N; i++) {
    weights = get_kernel(x, xgrid[i], h, kernel);
    fitted = wrm_fit(x, y, xgrid[i], weights);
    mu[i] = fitted[0];
    beta[i] = fitted[1];
  }
  
  return List::create(_["mu"] = mu,
                      _["beta"] = beta,
                      _["bw"] = h,
                      _["xgrid"] = xgrid,
                      _["kernel"] = kernel);
}

// 
// /*** R
// set.seed(123)
// dat <- rnorm(100, 10)
// gr <- seq(0, 1, length.out = 100)
// weights <- dat / sum(dat)
// p <- 0.5
// 
// wrm_smooth_cpp(x = gr, y = dat, h = p, gr, "epanechnikov")
// wrm.smooth(gr, dat, p, gr, 2)$level
// */

