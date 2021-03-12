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
using Eigen::SelfAdjointEigenSolver;    // one of the eigenvalue solvers

// [[Rcpp::export]]
VectorXd getEigenValues(Map<MatrixXd> M) {
  SelfAdjointEigenSolver<MatrixXd> es(M);
  return es.eigenvalues();
}

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
                          Rcpp::Nullable<Rcpp::NumericVector> weight_ = R_NilValue,
                          const int maxit = 30,
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


// [[Rcpp::export]]
inline Eigen::VectorXd locpolysmooth(Eigen::VectorXd Lt,
                                     Eigen::VectorXd Ly,
                                     Eigen::VectorXd newt,
                                     std::string kernel = "gauss",
                                     const double bw = 0.1,
                                     const double k = 1.345,
                                     const int deg = 1) {
  int n_newt = newt.size();   // number of grid which is predicted
  int n = Lt.size();   // number of Lt
  Eigen::VectorXd w(n);
  w.setOnes();
  Eigen::VectorXd weig = w.array() / n;   // 1/length(Lt)
  
  Eigen::VectorXd tmp(n);
  Eigen::VectorXd kern(n);
  Rcpp::NumericVector W(n);
  Eigen::MatrixXd X(n, deg+1);
  X.setOnes();
  Eigen::VectorXd Y(n);
  Rcpp::List fit;
  Eigen::VectorXd beta_hat(deg+1);
  Eigen::VectorXd mu_hat(n_newt);
  
  for (int t = 0; t < n_newt; t++) {
    tmp = (Lt.array() - newt(t)) / bw;
    
    if (kernel == "epanechnikov") {
      kern = (3./4.) * (1 - tmp.array().pow(2));   // Epanechnikov kernel
    } else if (kernel == "gauss") {
      kern = 1./sqrt(2.*M_PI) * (-1./2. * tmp.array().pow(2)).exp();   // gaussian kernel
    }
    
    // weight vector for w_i * K_h(x)
    W = (weig.array() * kern.array()) / bw;   
    
    // X matrix
    for (int d = 0; d < deg; d++) {
      X.col(d+1) = (Lt.array() - newt(t)).pow(d+1);
    }
    
    // Huber regression
    fit = IRLScpp(Ly, X, W, 30, 0.0001, k);
    
    beta_hat = fit["beta"];
    mu_hat(t) = beta_hat(0);
  }
  
  // Rcpp::NumericVector T = Rcpp::wrap(Ly);
  // // Rcpp::as<NumericVector>(Ly);
  // double Y = median(T);
  
  return mu_hat;
}





/*** R
locpolysmooth(Lt = 1:5,
              Ly = 1:5,
              newt = 1:5)
# library(MASS)
# IRLScpp(Y = stackloss$stack.loss,
#         X = as.matrix(stackloss[, 1:3]),
#         weight_ = NULL)
# IRLS(Y = stackloss$stack.loss, 
#      X = stackloss[, 1:3],
#      method = "huber")
*/

