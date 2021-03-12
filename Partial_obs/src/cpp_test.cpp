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


// [[Rcpp::export]]
Rcpp::List get_positive_elements(Eigen::VectorXd Y,
                                 Eigen::MatrixXd X,
                                 Eigen::VectorXd W) {
  int n = W.size();
  
  Rcpp::NumericVector Y_pos(n);
  Rcpp::NumericMatrix X_pos(n, X.cols());
  Rcpp::NumericVector W_pos(n);
  Rcpp::NumericVector row_na(X.cols(), R_NaN);
  Rcpp::NumericVector row_x(X.cols());
  for (int i = 0; i < n; i++) {
    if (W(i) > 0) {
      Y_pos(i) = Y(i);
      row_x = Rcpp::wrap(X.row(i));
      X_pos.row(i) = row_x;
      W_pos(i) = W(i);
    } else {
      Y_pos(i) = R_NaN;
      X_pos.row(i) = row_na;
      W_pos(i) = R_NaN;
    }
  }
  
  Rcpp::NumericVector Y_rm_neg = na_omit(Y_pos);
  Rcpp::NumericVector W_rm_neg = na_omit(W_pos);
  
  int n_pos = Y_rm_neg.length();
  
  Rcpp::NumericMatrix X_rm_neg(n_pos, X.cols());
  Rcpp::NumericVector col_x(n_pos);
  for (int j = 0; j < X.cols(); j++) {
    col_x = na_omit(X_pos.column(j));
    X_rm_neg.column(j) = col_x;
  }
  
  return Rcpp::List::create(Rcpp::_("Y") = Y_rm_neg, 
                            Rcpp::_["X"] = X_rm_neg,
                            Rcpp::_["W"] = W_rm_neg,
                            Rcpp::_["n_pos"] = n_pos);
}


// [[Rcpp::export]]
Eigen::VectorXd test(Eigen::VectorXd Y) {
  Y.resize(3);
  
  return Y;
}


/*** R
i <- c(-1,1,-2,2,-3,3,-4,4,-5,5)
get_positive_elements(i,
                      cbind(i, i),
                      i)
# test(1:5)
*/
