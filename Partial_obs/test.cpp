#include <RcppEigen.h>
#include <Rcpp.h>
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


// [[Rcpp::depends(RcppEigen)]]
// Get data with positive kernel weights for local regression
// - Y : observed data (unlist(Ly))
// - X : observed grid (unlist(Lt))
// - W : weight from kernel function
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
  
  return Rcpp::List::create(Rcpp::_["Y"] = Y_rm_neg,
                            Rcpp::_["X"] = X_rm_neg,
                            Rcpp::_["W"] = W_rm_neg,
                            Rcpp::_["n_pos"] = n_pos);
}


// Iteratively re-weighted least squares (IRLS) for robust regression (M-estimation)
// - Y : observed data (unlist(Ly))
// - X : observed grid (unlist(Lt))
// - weight_ : weight vector (Default = NULL)
// - maxit : maximun iteration numbers
// - weight : additional weight (for kernel regression)
// - tol : tolerence rate
// - k : delta for Huber function(or Tukey's biweight function)
// [[Rcpp::export]]
Rcpp::List IRLScpp(const Eigen::VectorXd Y,
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
  // Y_hat = X * beta.row(0);
  Rcpp::NumericVector resid = Rcpp::wrap(Y - Y_hat);
  Rcpp::NumericVector resid_med = abs(resid - median(resid));
  double s = median(resid_med) * 1.4826;  // re-scaled MAD by MAD*1.4826
  
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
    
    // update scale estimator
    resid = Rcpp::wrap(Y - Y_hat);
    resid_med = abs(resid - median(resid));
    s = median(resid_med) * 1.4826;  // re-scaled MAD by MAD*1.4826
    
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
    if ((beta.row(i+1) - beta.row(i)).cwiseAbs().sum() < tol) {
      iter = i+1;
      break;
    }
  }
  
  return Rcpp::List::create(Rcpp::_["beta"] = beta_hat,
                            Rcpp::_["iter"] = iter,
                            Rcpp::_["scale.est"] = s);
}


// Local polynomial kernel smoothing with huber loss (mean estimator)
// It is used for local_kern_smooth() in R function
// - Lt : a list of vectors or a vector containing time points for all curves
// - Ly : a list of vectors or a vector containing observations for all curves
// - newt : a vector containing time points to estimate
// - kernel : a kernel function for kernel smoothing ("epan", "gauss" are supported.)
// - bw : bandwidth
// - k : delta for Huber function (or Tukey's biweight function)
// - deg : degree of polynomial
// [[Rcpp::export]]
Eigen::VectorXd locpolysmooth(Eigen::VectorXd Lt,
                              Eigen::VectorXd Ly,
                              Eigen::VectorXd newt,
                              std::string kernel = "epanechnikov",
                              const double bw = 0.1,
                              const double k = 1.345,
                              const int deg = 1) {
  int n_newt = newt.size();   // number of grid which is predicted
  int n = Lt.size();   // number of Lt
  
  double weig = 1. / n;   // 1/length(Lt)
  
  Eigen::VectorXd tmp(n);
  Eigen::VectorXd kern(n);
  Rcpp::NumericVector W(n);
  // Eigen::VectorXd W(n);
  Eigen::MatrixXd X(n, deg+1);
  Eigen::VectorXd Y(n);
  Rcpp::List fit;
  Eigen::VectorXd beta_hat(deg+1);
  Eigen::VectorXd mu_hat(n_newt);
  Rcpp::List pos_idx_obj;
  int n_pos = 0;
  
  
  for (int t = 0; t < n_newt; t++) {
    tmp = (Lt.array() - newt(t)) / bw;
    
    if (kernel == "epanechnikov") {
      kern = (3./4.) * (1 - tmp.array().pow(2));   // Epanechnikov kernel
    } else if (kernel == "gauss") {
      // kern = (-1./2. * tmp.array().pow(2));
      // kern = kern.array().exp();
      kern = (1./sqrt(2.*M_PI)) * ((-0.5 * tmp.array().pow(2)).exp());   // gaussian kernel
    }
    
    Rcpp::Rcout << "sum =" << bw << "\t" << newt(t) << "\t" << kern << "\n";
    
    
    // X matrix
    X.resize(n, deg+1);
    X.setOnes();
    for (int d = 0; d < deg; d++) {
      X.col(d+1) = (Lt.array() - newt(t)).pow(d+1);
    }
    
    // obtain data which has non-negative weights
    pos_idx_obj = get_positive_elements(Ly, X, kern);
    n_pos = pos_idx_obj["n_pos"];
    X.resize(n_pos, deg+1);
    Y.resize(n_pos);
    kern.resize(n_pos);
    X = pos_idx_obj["X"];
    Y = pos_idx_obj["Y"];
    kern = pos_idx_obj["W"];
    
    // weight vector for w_i * K_h(x)
    W = (weig * kern.array()) / bw;
    
    
    
    // Huber regression
    fit = IRLScpp(Y, X, W, 30, 0.0001, k);
    
    beta_hat = fit["beta"];
    mu_hat(t) = beta_hat(0);
  }
  
  return mu_hat;
}




/*** R
library(mvtnorm)
library(robfpca)
source("R/sim_delaigle.R")
# source("R/sim_Lin_Wang(2020).R")
# source("load_source.R")

num_sim <- 30   # number of simulations
out_prop <- 0   # proportion of outliers
data_type <- "snippet"   # type of functional data
# kernel <- "epanechnikov"   # kernel function for local smoothing
kernel <- "gauss"   # kernel function for local smoothing
bw <- 0.1   # fixed bandwidth
n_cores <- 12   # number of threads for parallel computing

set.seed(1)
n <- 100
n.grid <- 51
# x.2 <- sim_kraus(n = n, out.prop = 0.2, out.type = 1, grid.length = n.grid)
x.2 <- sim_delaigle(n = n, 
                    model = 2,
                    type = data_type,
                    out.prop = out_prop, 
                    out.type = 1)
df <- data.frame(
  id = factor(unlist(sapply(1:length(x.2$Lt), 
                            function(id) { 
                              rep(id, length(x.2$Lt[[id]])) 
                            }) 
  )),
  y = unlist(x.2$Ly),
  t = unlist(x.2$Lt)
)

x <- list2matrix(x.2)
work.grid <- seq(0, 1, length.out = n.grid)



Lt <- x.2$Lt
Ly <- x.2$Ly
method <- "HUBER"
deg <- 1
K <- 5
bw_cand = 10^seq(-2, 0, length.out = 5)/3
cv_loss<-"HUBER"
delta <- 1.345


# get index for each folds
folds <- list()
n <- length(Lt)   # the number of curves
fold_num <- n %/% K   # the number of curves for each folds
fold_sort <- sample(1:n, n)
for (k in 1:K) {
  ind <- (fold_num*(k-1)+1):(fold_num*k)
  if (k == K) {
    ind <- (fold_num*(k-1)+1):n
  }
  folds[[k]] <- fold_sort[ind]
}
k <- 1


# for (k in 1:K) {
  print(k)
  Lt_train <- Lt[ -folds[[k]] ]
  Ly_train <- Ly[ -folds[[k]] ]
  Lt_test <- Lt[ folds[[k]] ]
  Ly_test <- Ly[ folds[[k]] ]
  
  for (i in 1:length(bw_cand)) {
    print(i)
    # y_hat <- local_kern_smooth(Lt = Lt_train,
    #                            Ly = Ly_train,
    #                            newt = Lt_test,
    #                            method = method,
    #                            bw = bw_cand[i],
    #                            kernel = kernel,
    #                            delta = 1.345)
    
    mu_hat <- locpolysmooth(Lt = unlist(Lt_train),
                            Ly = unlist(Ly_train),
                            newt = work.grid[30:31],
                            kernel = kernel,
                            bw = bw_cand[i],
                            k = delta,
                            deg = deg)
  }
# }

*/
