#include <RcppEigen.h>
//[[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

// Weight function I_W(x)
// [[Rcpp::export]]
double I_w(const Eigen::VectorXd& x, const Eigen::MatrixXd& S_inv) {
  return exp(-0.5 * x.transpose() * S_inv * x);
}

// [[Rcpp::export]]
double computeGamma(const Eigen::MatrixXd& X1, const Eigen::MatrixXd& X2, const Eigen::MatrixXd& S_inv) {
  int n1 = X1.rows();
  int n2 = X2.rows();
  double gamma = 0.0;
  
  // 1st term
  for (int i = 0; i < n1; ++i) {
    for (int j = 0; j < n1; ++j) {
      if (i == j) {
        gamma += I_w(X1.row(i) - X1.row(j), S_inv) / (n1 * n1);
      } else if (i > j) {
        gamma += 2. * I_w(X1.row(i) - X1.row(j), S_inv) / (n1 * n1);
      }
    }
    // 3rd term
    if (n1 >= n2) {
      for (int j = 0; j < n2; ++j) {
        if (i == j || ((i > j && i > n2-1))) {
          gamma -= I_w(X1.row(i) - X2.row(j), S_inv) * (2.0 / (n1 * n2));
        } else if (i > j && i <= n2-1) {
          gamma -= 2. * I_w(X1.row(i) - X2.row(j), S_inv) * (2.0 / (n1 * n2));
        }
      }
    }
  }
  
  // 2nd term
  for (int i = 0; i < n2; ++i) {
    for (int j = 0; j < n2; ++j) {
      if (i == j) {
        gamma += I_w(X2.row(i) - X2.row(j), S_inv) / (n2 * n2);
      } else if (i > j) {
        gamma += 2. * I_w(X2.row(i) - X2.row(j), S_inv) / (n2 * n2);
      }
    }
    // 3rd term
    if (n1 < n2) {
      for (int j = 0; j < n1; ++j) {
        if (i == j || ((i > j && i > n1-1))) {
          gamma -= I_w(X1.row(j) - X2.row(i), S_inv) * (2.0 / (n1 * n2));
        } else if (i > j && i <= n1-1) {
          gamma -= 2. * I_w(X1.row(j) - X2.row(i), S_inv) * (2.0 / (n1 * n2));
        }
      }
    }
  }
  
  return gamma;
}

