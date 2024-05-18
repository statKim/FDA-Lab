#include <RcppEigen.h>
#include <Rcpp.h>
#include "trapzRcpp.h"
using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
List uFPCA_cpp(Eigen::MatrixXd X, Nullable<NumericVector> grid = R_NilValue, Nullable<int> K = R_NilValue, double FVE = 0.95, bool centered = false) {
  
  int n = X.rows();  // number of curves
  int m = X.cols();  // number of timepoints
  
  // Observed grid points
  NumericVector grid_r;
  if (grid.isNull()) {
    grid_r = seq(0, 1, length.out = m);
  } else {
    grid_r = grid;
  }
  
  // Centering the curves
  Eigen::VectorXd mu;
  if (!centered) {
    mu = X.colwise().mean();
    X = X.rowwise() - mu.transpose();
  } else {
    mu = Eigen::VectorXd::Zero(m);
  }
  
  // SVD
  BDCSVD<MatrixXd> svd(X.transpose() / sqrt(n), ComputeThinU | ComputeThinV);
  VectorXd d = svd.singularValues();
  MatrixXd u = svd.matrixU();
  
  // Eigen analysis
  // MatrixXd cov_mat = X * X.transpose() / n;
  // EigenSolver<MatrixXd> eig(cov_mat);
  // VectorXd d = eig.eigenvalues().real();
  // MatrixXd u = eig.eigenvectors().real();
  
  std::vector<int> positive_ind;
  for (int i = 0; i < d.size(); ++i) {
    if (d(i) > 0) {
      positive_ind.push_back(i);
    }
  }
  
  VectorXd lambda = d.segment(positive_ind[0], positive_ind.size()).array().square();
  MatrixXd phi = u.block(0, positive_ind[0], u.rows(), positive_ind.size());
  
  // Normalizing eigenvalues and eigenfunctions
  double grid_size = grid_r[1] - grid_r[0];
  NumericVector work_grid = seq(min(grid_r), max(grid_r), m);
  lambda = lambda * grid_size;
  for (int j = 0; j < phi.cols(); ++j) {
    phi.col(j) = phi.col(j) / sqrt(grid_size);
    // phi.col(j) = phi.col(j) / sqrt(phi.col(j).squaredNorm());  // normalize
    phi.col(j) = phi.col(j) / sqrt(trapzRcpp(work_grid, phi.col(j).array().square()));
    if (phi.col(j).dot(work_grid) < 0)
      phi.col(j) *= -1;
  }
  
  // Select the number of FPCs
  NumericVector cum_FVE = cumsum(lambda) / sum(lambda);  // cumulative FVE
  int K_cpp;
  if (K.isNull()) {
    K_cpp = std::distance(cum_FVE.begin(), std::find_if(cum_FVE.begin(), cum_FVE.end(), [FVE](double x){return x >= FVE;})) + 1;
  } else {
    K_cpp = as<int>(K);
  }
  if (K_cpp > lambda.size()) {
    stop("K is too large.");
  }
  
  // FPC scores
  lambda.conservativeResize(K_cpp);
  phi.conservativeResize(phi.rows(), K_cpp);
  MatrixXd xi(m, K_cpp);
  for (int k = 0; k < K_cpp; ++k) {
    for (int i = 0; i < n; ++i) {
      xi(i, k) = trapzRcpp(work_grid, X.row(i).array() * phi.col(k).array());
    }
  }
  
  return List::create(Named("eig.val") = lambda,
                      Named("eig.ftn") = phi,
                      Named("fpc.score") = xi,
                      Named("K") = K_cpp,
                      Named("FVE") = cum_FVE[K_cpp - 1],
                                            Named("work.grid") = work_grid,
                                            Named("mean.ftn") = mu);
}



/*** R
# timesTwo(42)
*/
