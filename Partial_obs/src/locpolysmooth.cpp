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
inline Eigen::VectorXd locpolysmooth(NumericVector x) {
  return x * 2;
}



/*** R
library(MASS)
IRLScpp(Y = stackloss$stack.loss,
        X = as.matrix(stackloss[, 1:3]))
IRLS(Y = stackloss$stack.loss, 
     X = stackloss[, 1:3],
     method = "huber")
*/

