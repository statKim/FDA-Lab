#include <Rcpp.h>
#include <map>
#include <string>
#include <vector>
#include <algorithm>

using namespace Rcpp;

// expand.grid
// [[Rcpp::export]]
NumericMatrix expand_grid_cpp(NumericVector x,
                              NumericVector y) {
  int n = x.length();
  
  NumericMatrix res(pow(n, 2), 2);
  res.column(0) = rep_len(x, pow(n, 2));
  res.column(1) = rep_each(y, n);
  
  return res;
}



// Huber's location M-estimator
// [[Rcpp::export]]
double huber_cpp(NumericVector y,
                 double k = 1.5,
                 double tol = 1e-6) {
  y = y[!is_na(y)];
  int n = y.size();
  double mu = median(y);   // median
  
  if (n == 0) {
    Rcpp::stop("There is no available observation to compute Huber's location M-estimator.");
  }
  
  // MAD
  NumericVector s_ = abs(y - mu);
  double s = 1.4826 * median(s_);
  if (s == 0) {
    Rcpp::stop("Estimated MAD is 0 in Huber's location M-estimator.");
    // Rcout << "Estimated MAD is 0." << "\n";
    return 0;
  }
  
  // Obtain Huber's M-estimator
  NumericVector yy;
  double mu1;
  while(1) {   // infinite loop
    yy = pmin(pmax(mu - k * s, y), mu + k * s);
    mu1 = sum(yy) / n;
    
    if (abs(mu - mu1) < tol*s) {
      break;
    }
    
    mu = mu1;
  }
  
  return mu;
}


// // Local M-estimator
// // [[Rcpp::export]]
// NumericVector local_M_1D(NumericMatrix X,
//                          NumericVector obs_t,
//                          NumericVector new_t,
//                          double h) {
//   int p = new_t.length();
//   IntegerVector idx = seq(0, p-1);
//   NumericVector mu(p);
//   for (int i = 0; i < p; i++) {
//     IntegerVector i_neighbor = idx[obs_t < new_t[i] + h & obs_t > new_t[i] - h];
//     NumericVector X_i_neighbor;
//     mu[i] = huber_cpp(X_i_neighbor);
//   }
//   
//   return mu;
// }




// 같은 timepoint는 반복 안하도록 설정하기!!
// for loop of "get_raw_cov" in R function
// [[Rcpp::export]]
NumericMatrix get_raw_cov_cpp(NumericMatrix X,
                              NumericVector mu,
                              NumericVector gr,
                              bool diag = false) {
  int n = X.rows();
  int p = X.cols();
  IntegerVector ind_vec = seq(0, p-1);
  
  NumericMatrix raw_cov;
  for (int i = 0; i < n; i++) {
    IntegerVector ind = ind_vec[!is_na(X.row(i))];
    NumericVector gr_sub = gr[ind];
    NumericVector X_i = X.row(i);
    NumericVector mu_sub = mu[ind];
    X_i = X_i[ind];
    X_i = X_i - mu_sub;
    
    NumericMatrix exp_grid = expand_grid_cpp(gr_sub, gr_sub);
    NumericMatrix exp_X_i = expand_grid_cpp(X_i, X_i);
    NumericVector exp_grid_1 = exp_grid.column(0);
    NumericVector exp_grid_2 = exp_grid.column(1);
    
    if (diag == false) {
      // eliminate diagonal parts
      IntegerVector idx = seq(0, exp_X_i.nrow()-1);
      idx = idx[exp_grid_1 != exp_grid_2];
      NumericVector cov_i = exp_X_i.column(0) * exp_X_i.column(1);
      cov_i = cov_i[idx];
      NumericVector exp_grid_1_sub = exp_grid_1[idx];
      NumericVector exp_grid_2_sub = exp_grid_2[idx];
      
      // We use "cbind" in Rcpp, so we should transpose it.
      NumericMatrix raw_cov_i(3, cov_i.length());
      raw_cov_i(0, _) = cov_i;
      raw_cov_i(1, _) = exp_grid_1_sub;
      raw_cov_i(2, _) = exp_grid_2_sub;
      
      if (i == 0) {
        raw_cov = raw_cov_i;
      } else {
        raw_cov = cbind(raw_cov, raw_cov_i);
      }
    } else {
      NumericVector cov_i = exp_X_i.column(0) * exp_X_i.column(1);
      
      // We use "cbind" in Rcpp, so we should transpose it.
      NumericMatrix raw_cov_i(3, cov_i.length());
      raw_cov_i(0, _) = cov_i;
      raw_cov_i(1, _) = exp_grid_1;
      raw_cov_i(2, _) = exp_grid_2;
      
      if (i == 0) {
        raw_cov = raw_cov_i;
      } else {
        raw_cov = cbind(raw_cov, raw_cov_i);
      }
    }
  }
  
  // transpose because we cbind the raw_cov_i
  raw_cov = transpose(raw_cov);
  
  return raw_cov;
}


// for loop of "cov_local_M" in R function
// [[Rcpp::export]]
NumericMatrix cov_local_M_cpp(NumericVector raw_cov,
                              NumericVector s,
                              NumericVector t,
                              NumericVector gr,
                              double h = 0.02) {
  int p = gr.length();
  IntegerVector idx = seq(0, s.length()-1);
  
  // IntegerVector i_neighbor;
  // IntegerVector j_neighbor;
  // IntegerVector ind;
  // NumericVector raw_cov_sub;
  NumericMatrix cov_mat(p, p);
  for (int i = 0; i < p; i++) {
    IntegerVector i_neighbor = idx[s < gr[i] + h & s > gr[i] - h];
    
    for (int j = 0; j < p; j++) {
      if (i <= j) {
        IntegerVector j_neighbor = idx[t < gr[j] + h & t > gr[j] - h];

        IntegerVector ind = intersect(i_neighbor, j_neighbor);
        // ind = idx[t < gr[j] + h & t > gr[j] - h & s < gr[i] + h & s > gr[i] - h];
        NumericVector raw_cov_sub = raw_cov[ind];
        
        cov_mat(i, j) = huber_cpp(raw_cov_sub);
      } else {
        cov_mat(i, j) = cov_mat(j, i);
      }
    }
  }

  return cov_mat;
}




// 
// /*** R
// x <- matrix(c(1,NA,NA,10,2,3,4,NA), nrow = 2)
// set.seed(100)
// x <- matrix(rnorm(510), ncol = 51)
// gr <- seq(0, 1, length.out = 51)
// x.2 <- sim_delaigle(n = ,
//                     model = 2,
//                     type = data_type,
//                     out.prop = out_prop,
//                     out.type = out_type)
// x <- list2matrix(x.2)
// 
// # system.time({
// #   raw_cov <- get_raw_cov_cpp(x,
// #                              rep(0, 51),
// #                              gr, TRUE)
// #   a1 <- cov_local_M_cpp(raw_cov[, 1],
// #                         raw_cov[, 2],
// #                         raw_cov[, 3],
// #                         gr,
// #                         0.1)
// # })
// # system.time({
// #   a2 <- test(x, 0.1, gr)
// # })
// 
// all.equal(a1, a2)
// 
// # GA::persp3D(gr, gr, a1,
// #             theta = -70, phi = 30, expand = 1)
// # GA::persp3D(gr, gr, a2,
// #             theta = -70, phi = 30, expand = 1)
// */
