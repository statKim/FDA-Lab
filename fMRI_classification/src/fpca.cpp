#include <Rcpp.h>
//[[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "trapzRcpp.h"

using namespace Rcpp;
using namespace Eigen;


// [[Rcpp::export]]
// SVD for n x p matrix
List svd_cpp(Eigen::MatrixXd X) {
  X = X.transpose();  // p x n matrix
  Eigen::BDCSVD<Eigen::MatrixXd> svd(X, Eigen::ComputeFullU | Eigen::ComputeFullV);
  Eigen::MatrixXd U = svd.matrixU();
  Eigen::MatrixXd V = svd.matrixV();
  Eigen::VectorXd S = svd.singularValues();
  
  return List::create(Named("singular.value") = S,
                      Named("eig.ftn") = U);
}


// [[Rcpp::export]]
Rcpp::NumericMatrix normalize_phi(Rcpp::NumericMatrix phi, Rcpp::NumericVector work_grid) {
  Rcpp::NumericMatrix phi_normalized(phi.rows(), phi.cols());
  Rcpp::NumericVector phi_i(phi.rows());
  for (int i = 0; i < phi.cols(); i++) {
    phi_i = phi.column(i);
    // phi_i = phi_i / sqrt(grid_size);
    phi_i = phi_i / sqrt(trapzRcpp(work_grid, pow(phi_i, 2)));
    if (0 <= sum(phi_i * work_grid)) {
      phi_normalized.column(i) = phi_i;
    } else {
      phi_normalized.column(i) = -phi_i;
    }
  }
  // std::cout << work_grid;
  return phi_normalized;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix get_fpc_scores(Rcpp::NumericMatrix X, Rcpp::NumericMatrix phi, Rcpp::NumericVector work_grid) {
  int K = phi.cols();
  int n = X.rows();
  
  Rcpp::NumericMatrix xi(n, K);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < K; j++) {
      xi(i, j) = trapzRcpp(work_grid, X.row(i) * phi.column(j));
    }
  }
  
  return xi;
}

// // [[Rcpp::export]]
// Rcpp::NumericMatrix get_fpc_scores2(Eigen::MatrixXd X, Eigen::MatrixXd phi, Rcpp::NumericVector work_grid) {
//   int K = phi.cols();
//   int n = X.rows();
//   int m = X.cols();
// 
//   Rcpp::NumericMatrix xi(n, K);
//   Rcpp::NumericMatrix X_phi(m, n);
//   // Eigen::MatrixXd sub(X.cols(), K);
//   for (int i = 0; i < K; i++) {
//     // sub = X.transpose() * phi.col(i);
//     // X_phi = Rcpp::as<Rcpp::NumericVector>(sub);
//     // SEXP sub = Rcpp::wrap(X.transpose() * phi.col(i));
//     X_phi = Rcpp::wrap(X.transpose().array() * phi.col(i).array());
//     for (int j = 0; j < n; j++) {
//       xi(j, i) = trapzRcpp(work_grid, X_phi.column(j));
//     }
//   }
// 
//   return xi;
// }


/*** R
# system.time({
#   xi <- sapply(1:K, function(k) {
#     apply(t(X) * phi[, k], 2, function(x_i_phi_k) {
#       fdapace::trapzRcpp(work.grid, x_i_phi_k)
#     })
#   })
# })
# system.time({
#   xi2 <- get_fpc_scores(X, phi[, 1:K], work.grid)
# })
# all.equal(xi, xi2)

# # 약간 다름
# system.time({
#   phi <- svd.obj$u[, positive_ind]
#   phi <- phi / sqrt(grid_size)
#   phi <- apply(phi, 2, function(x) {
#     # x <- x / sqrt(grid_size)
#     x <- x / sqrt(fdapace::trapzRcpp(work.grid, x^2))
#     if ( 0 <= sum(x * work.grid) )
#       return(x)
#     else
#       return(-x)
#   })
# })
# system.time({
#   phi <- svd.obj$u[, positive_ind]
#   phi <- phi / sqrt(grid_size)
#   phi2 <- normalize_phi(phi, work.grid)
# })
# all.equal(phi, phi2)
# phi[1:5, 1:5]
# phi2[1:5, 1:5]

# # 많이 다름
# system.time({
#   svd.obj <- svd(t(X) / sqrt(n))
#   # positive_ind <- which(svd.obj$d > 0)   # indices of positive eigenvalues
#   # lambda <- svd.obj$d[positive_ind]^2   # eigenvalues
#   # phi <- svd.obj$u[, positive_ind]
# })
# 
# system.time({
#   svd.obj2 <- svd_cpp(X / sqrt(n))
#   # positive_ind <- which(svd.obj$d > 0)   # indices of positive eigenvalues
#   # lambda <- svd.obj$d[positive_ind]^2   # eigenvalues
#   # phi <- svd.obj$u[, positive_ind]
# })
# 
# svd.obj$d[1:10]
# svd.obj2$singular.val[1:10]
# 
# tail(svd.obj$d)
# tail(svd.obj2$singular.val)
# svd.obj$d %>% round(3)
# svd.obj2$singular.val %>% round(3)
*/
