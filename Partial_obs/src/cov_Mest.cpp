#include <Rcpp.h>
#include <map>
#include <string>
#include <vector>
#include <algorithm>

using namespace Rcpp;

// Huber's location M-estimator
// [[Rcpp::export]]
double huber_cpp(NumericVector y,
                 double k = 1.5,
                 double tol = 1e-6) {
  y = y[!is_na(y)];
  int n = y.size();
  double mu = median(y);   // median
  
  // MAD
  NumericVector s_ = abs(y - mu);
  double s = 1.4826 * median(s_);
  
  // Obtain Huber's M-estimator
  NumericVector yy;
  double mu1;
  while(1) {   // infinite loop
    yy = pmin(pmax(mu - k * s, y), mu + k * s);
    mu1 = sum(yy) / n;
    
    if (s == 0) {
      Rcout << "Estimated MAD is 0." << "\n";
      return 0;
    }
    
    if (abs(mu - mu1) < tol*s) {
      break;
    }
    
    mu = mu1;
  }
  
  return mu;
}


// Obtain column-wise M-estimator
// [[Rcpp::export]]
NumericVector mean_Mest_cpp(NumericMatrix X) {
  int p = X.ncol();
  NumericVector mu(p);
  for (int i = 0; i < p; i++) {
    mu[i] = huber_cpp(X.column(i));
  }
  
  return mu;
}


// Obtain marginal M-estimator for covariance
// [[Rcpp::export]]
NumericMatrix cov_Mest_cpp(NumericMatrix X) {
  int n = X.nrow();
  int p = X.ncol();
  NumericMatrix rob_var(p, p);
  
  LogicalVector NA_ind(n);
  IntegerVector ind_full = seq(0, n-1);
  IntegerVector ind_st;
  IntegerVector ind;
  NumericVector X_i;
  NumericVector mu;
  NumericMatrix X_sub;
  NumericVector A;
  
  if (sum(is_na(X)) == 0) {
    // Complete curves
    mu = mean_Mest_cpp(X);
    for (int i = 0; i < X.nrow(); i++) {
      X(i, _) = X(i, _) - mu;   // X-mu
    }
    A = NumericVector(n);
    for (int s = 0; s < p; s++) {
      for (int t = 0; t < p; t++) {
        if (s <= t) {
          for (int i = 0; i < X.nrow(); i++) {
            A[i] = X(i, s) * X(i, t);
          }
          
          rob_var(s, t) = huber_cpp(A);
        } else {
          rob_var(s, t) = rob_var(t, s);
        }
      }
    }
  } else {
    // Partially observed curves having missing
    for (int s = 0; s < p; s++) {
      for (int t = 0; t < p; t++) {
        if (s <= t) {
          // Find the index of curve which have NA at time s or t
          ind_st = IntegerVector::create(s, t);
          for (int i = 0; i < n; i++) {
            X_i = X.row(i);
            X_i = X_i[ind_st];
            NA_ind[i] = sum(is_na(X_i)) > 0;
          }
          
          if (sum(NA_ind) > 0) {
            ind = ind_full[!NA_ind];
          } else {
            ind = ind_full;
          }
          
          X_sub = NumericMatrix(ind.size(), p);
          for (int i = 0; i < ind.size(); i++) {
            X_sub(i, _) = X(ind[i], _);
          }
          
          mu = mean_Mest_cpp(X_sub);
          
          A = NumericVector(ind.size());
          for (int i = 0; i < X_sub.nrow(); i++) {
            X_sub(i, _) = X_sub(i, _) - mu;
            A[i] = X_sub(i, s) * X_sub(i, t);
          }
          
          rob_var(s, t) = huber_cpp(A);
        } else {
          rob_var(s, t) = rob_var(t, s);
        }
      }
    }
  }
  
  // Rcout << rob_var << "\n";
    
  return rob_var;
}


// 
// /*** R
// # x <- matrix(rnorm(10000), nrow = 10)
// system.time({
//   cov1 <- cov_Mest_cpp(x)
// })
// system.time({
//   cov2 <- var.rob.missfd(x, make.pos.semidef=F)
// })
// all.equal(cov1, cov2)
// 
// # y <- matrix(rnorm(10000), nrow = 10)
// # system.time({
// #   mu1 <- mean_Mest_cpp(x)
// # })
// # 
// # system.time({
// #   mu2 <- as.numeric(mean.rob.missfd(x))
// # })
// # 
// # all.equal(mu1, mu2)
// # # # mu <- apply(x, 2, function(t){ MASS::huber(t)$mu })
// */
