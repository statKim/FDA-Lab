#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include <map>          
#include <string>       
#include <vector>       
#include <algorithm>  

typedef Eigen::Map<Eigen::MatrixXd> MapMatd;
typedef Eigen::Map<Eigen::VectorXd> MapVecd;

using Rcpp::List;
using Rcpp::as;

using Eigen::Map;                       // 'maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision

// Get median using a R function
// [[Rcpp::export]]
double median_r(Eigen::VectorXd x) {
  // Obtain environment containing function
  Rcpp::Environment base("package:stats"); 
  
  // Make function callable from C++
  Rcpp::Function median_r = base["median"];  
  
  // Call the function and receive its list output
  Rcpp::NumericVector res = median_r(Rcpp::_["x"] = x,
                                     Rcpp::_["na.rm"]  = true);
  
  // Return test object in list structure
  return res[0];
  
}

// Iteratively re-weighted least sqaures
// [[Rcpp::export]]
inline Rcpp::List IRLScpp(const Eigen::VectorXd Y, 
                          const Eigen::MatrixXd X,
                          const int maxit = 30,
                          Rcpp::Nullable<Rcpp::NumericVector> weight_ = R_NilValue,
                          const double tol = 0.0001,
                          const double k = 1.345) {
  int n = Y.size();   // number of observations (unlist(Lt))
  
  // initial value for beta (LSE estimator)
  Eigen::MatrixXd beta(maxit, X.cols());   // container of beta for iterations
  beta.setZero();   // initialize to 0
  Eigen::LDLT<Eigen::MatrixXd> ldlt_XTX(X.transpose() * X);
  beta.row(0) = ldlt_XTX.solve(X.transpose() * Y);   // LSE estimator
  
  Eigen::VectorXd Y_hat(n);   // Y_hat
  Eigen::VectorXd beta_hat(X.cols());   // beta_hat
  beta_hat = beta.row(0);
  Y_hat = X * beta_hat;
  Eigen::VectorXd resid = Y - Y_hat;
  Eigen::VectorXd resid_med = resid.array() - median_r(resid);
  double s = median_r(resid_med.cwiseAbs()) * 1.4826;  // re-scaled MAD by MAD*1.4826
  
  // if weight is NULL, set 1
  Eigen::VectorXd weight(n);
  if (weight_.isNotNull()) {
    weight = Rcpp::as<Eigen::VectorXd>(weight_);
  } else {
    weight.setOnes();
  }

  // variables for iterations
  Eigen::VectorXd tmp(n);
  Eigen::VectorXd psi(n);
  Eigen::VectorXd w(n);
  Eigen::MatrixXd XTWX(X.cols(), X.cols());   // X'WX
  Eigen::LDLT<Eigen::MatrixXd> ldlt_XTWX(XTWX);   // LDLT obj for X'WX
  int iter = maxit;
  
  // iterate until beta converged
  for (int i = 0; i < maxit; i++) {
    beta_hat = beta.row(i);
    Y_hat = X * beta_hat;
    tmp = (Y - Y_hat).array() / s;

    // psi function of Huber loss
    for (int j = 0; j < n; j++) {
      if (abs(tmp(j)) >= k) {
        if (tmp(j) < 0) {
          psi(j) = -1 * k;
        } else {
          psi(j) = k;
        }
      } else {
        psi(j) = tmp(j);
      }
    }

    // weight matrix for WLS
    w = (weight.array() * psi.array()).array() / (Y - Y_hat).array();

    // WLS
    XTWX = X.transpose() * w.asDiagonal() * X;
    ldlt_XTWX.compute(XTWX);
    beta.row(i+1) = ldlt_XTWX.solve(X.transpose() * w.asDiagonal() * Y);
    beta_hat = beta.row(i+1);

    // if beta converges before maxit, break.
    // conv = (beta.row(i+1) - beta.row(i)).cwiseAbs().sum();
    if ((beta.row(i+1) - beta.row(i)).cwiseAbs().sum() < tol) {
      iter = i+1;
      break;
    }
  }
  
  return Rcpp::List::create(Rcpp::_("beta") = beta_hat, 
                            Rcpp::_["iter"] = iter,
                            Rcpp::_["scale.est"] = s);
}



// 
// /*** R
// library(MASS) 
// IRLScpp(Y = stackloss$stack.loss,
//         X = as.matrix(stackloss[, 1:3]))
// IRLS(Y = stackloss$stack.loss, 
//      X = stackloss[, 1:3],
//      method = "huber")
// */
