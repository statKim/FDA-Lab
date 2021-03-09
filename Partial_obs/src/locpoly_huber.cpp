#include <RcppEigen.h>
#include <map>          
#include <string>       
#include <vector>       
#include <algorithm>    

typedef Eigen::Map<Eigen::MatrixXd> MapMatd;
typedef Eigen::Map<Eigen::VectorXd> MapVecd;

using Rcpp::List;
using Rcpp::as;




inline double IRLS_huber(const Eigen::VectorXd& Y, 
                         const Eigen::MatrixXd& X, 
                         const Eigen::VectorXd& weight,  // weight for w_i*K_h
                         const int maxit,
                         const double tol,
                         const double k) {
  
  int n = Y.size();
  
  // initial value for beta (LSE estimator)
  Eigen::MatrixXd beta(X.rows(), X.cols());
  Eigen::LDLT<Eigen::MatrixXd> ldlt_XTX(X.transpose() * X);
  Eigen::VectorXd b = ldlt_XTX.solve(X.transpose() * Y);
  beta.rows(1) = b;   // LSE estimator
  
  Eigen::VectorXd resid = Y - X * beta.rows(1);
  double s = 1.4826 * abs(resid);  // median 빼야됨
  
  for (int i = 1; i <= maxit; i++) {
    Eigen::MatrixXd W(n, n);
    Eigen::VectorXd Y_hat = X * beta.rows(i);
    
    Eigen::VectorXd tmp = (Y - Y_hat) / s;
    
    Eigen::VectorXd psi = 
    
  }
  
  
  
  Eigen::VectorXd x_tmp((x.array()-xi)/h);
  
  Eigen::VectorXd temp = get_w(kcode,x_tmp,w);
  
  //Rcpp::Rcout << "temp.size()=" << temp.size() << std::endl;
  
  Eigen::MatrixXd X(z.size() , 1+d);
  X.setOnes();
  // if(d > 0)
  // {
  //     auto v = x.array()-xi;
  //     for(int j=1; j <= d; j++)
  //     {
  //         x.col(j) = v;
  //         if(j!=d) v = v * v;
  //     }
  // }
  
  for(int j=1; j <= d; j++) X.col(j) = (x.array()-xi).array().pow(j);
  
  //Rcpp::Rcout << "d=" << d << std::endl;
  
  //X.col(1) = x.array() - xi;
  Eigen::LDLT<Eigen::MatrixXd> ldlt_XTWX(X.transpose() * temp.asDiagonal() *X);
  Eigen::VectorXd b = ldlt_XTWX.solve(X.transpose() * temp.asDiagonal() * z);
  return b(0);
}






















const double SQRT2 = sqrt(2);
const double invSqrt2pi=  1./(sqrt(2.*M_PI));


inline Eigen::VectorXd get_w(const int kcode, const Eigen::VectorXd& x, const Eigen::VectorXd& w) {
  Eigen::VectorXd temp(w.size());
  
  switch (kcode){
    case 1: // Epan
    temp=  (1-x.array().pow(2)) * ((3./4)*w).array();
    break;
    case 2 : // Rect
    temp=(w.array())*.5 ;
    break;
    case 3:  // triangular
    temp = (1-x.array().abs()) * (w.array());
    break;
    case 4:  // quartic
    temp = (w.array()) * ((1.-x.array().pow(2)).array().pow(2)).array() * (15./16.);
    break;
    case 5: // triweight
    temp = (w.array()) *
      ((1.-x.array().pow(2)).array().pow(3)).array() * (35./32.);
    break;
    case 6: // tricube
    temp = (w.transpose().array()) * ((1.-x.array().pow(3)).array().pow(3)).array() * (70./81.);
    break;
    case 7: // cosine
    temp = (w.array()) * (((M_PI/2)*x.array()).cos()).array() * (M_PI/4);
    break;
    case 21: // gauss
    temp = ((-.5*(x.array().pow(2))).exp()) * invSqrt2pi *  (w.array());
    break;
    case 22: // logistic
    temp = w.array() / ((2+x.array().exp()+(-x.array()).exp()));
    break;
    case 23: // sigmoid
    temp = ((2./(M_PI)) * w.array()) / ((x.array().exp()+(-x.array()).exp()));
    break;
    case 24:  // silverman
    temp = (0.5 * w.array()) * ((x.array().abs() / (-SQRT2)).exp() * ((M_PI/4)+(x.array().abs() / SQRT2)).sin());
    break;
    default: Rcpp::stop("Unsupported kernel");
    
  }
  
  return temp;
}

inline double locpoly(const Eigen::VectorXd& x, 
                      const Eigen::VectorXd& z, 
                      const Eigen::VectorXd& w, 
                      const double h, 
                      const double xi,
                      const int d, // degree
                      const int kcode) {
  
  Eigen::VectorXd x_tmp((x.array()-xi)/h);
  
  Eigen::VectorXd temp = get_w(kcode,x_tmp,w);
  
  //Rcpp::Rcout << "temp.size()=" << temp.size() << std::endl;
  
  Eigen::MatrixXd X(z.size() , 1+d);
  X.setOnes();
  // if(d > 0)
    // {
      //     auto v = x.array()-xi;
      //     for(int j=1; j <= d; j++)
        //     {
          //         x.col(j) = v;
          //         if(j!=d) v = v * v;
          //     }
      // }
  
  for(int j=1; j <= d; j++) X.col(j) = (x.array()-xi).array().pow(j);
  
  //Rcpp::Rcout << "d=" << d << std::endl;
  
  //X.col(1) = x.array() - xi;
  Eigen::LDLT<Eigen::MatrixXd> ldlt_XTWX(X.transpose() * temp.asDiagonal() *X);
  Eigen::VectorXd b = ldlt_XTWX.solve(X.transpose() * temp.asDiagonal() * z);
  return b(0);
}


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]

Eigen::MatrixXd csmoothmean(const Eigen::Map<Eigen::VectorXd> & x, // sorted in increasing order
                            const Eigen::Map<Eigen::VectorXd> & z, // Y(x)
                            const Eigen::Map<Eigen::VectorXd> & w, // weight for each observation
                            const double h,
                            const std::string kernel,
                            const int d, // degree
                            Eigen::Map<Eigen::VectorXd> & newx // sorted in increasing order
)
{
  
  // kcode <= 20 is reserved for compactly supported kernel on [-1,1]
  std::map<std::string,int> supported_kernels;
  supported_kernels["epanechnikov"]    = 1;
  supported_kernels["rectangular"]    = 2;
  supported_kernels["triangular"]   = 3;
  supported_kernels["quartic"]    = 4;
  supported_kernels["triweight"]    = 5;
  supported_kernels["tricube"]   = 6;
  supported_kernels["cosine"] = 7;
  supported_kernels["gauss"]   = 21;
  supported_kernels["logistic"] = 22;
  supported_kernels["sigmoid"] = 23;
  supported_kernels["silverman"] = 24;
  
  int kcode = 1;
  if ( supported_kernels.count( kernel ) != 0){
    kcode = supported_kernels.find( kernel )->second;
  }
  else
  {
    Rcpp::stop("The specified kernel is not supported");
  }
  
  
  const unsigned int xn = newx.size();
  const unsigned int n = x.size();
  
  const double* xDat = x.data();
  
  //Rcpp::Rcout << "xn=" << xn << ", n=" << n << std::endl;
  
  
  Eigen::VectorXd muhat(xn);// = MatrixXd::Zero(xn,yn);
  muhat.setOnes();
  
  
  
  if(kcode > 20) // noncompact kernel
  {
    for (unsigned int j = 0; j < xn; ++j)
      muhat(j) = locpoly(x,z,w,h,newx(j),d,kcode);
  }
  else // compactly supported kernel
  {
    
    unsigned int idxl = 0;
    unsigned int idxu = 0;
    
    for (unsigned int j = 0; j < xn; ++j) {
      
      const double xl = newx(j) - h;
      const double xu = newx(j) + h;
      
      idxu = std::upper_bound(xDat+idxu, xDat + n, xu) - xDat;
      idxl = std::lower_bound(xDat+idxl, xDat + n, xl) - xDat;
      
      
      if(idxu - idxl < 5)
      {
        //Rcpp::stop("Less than 5 points in the local window.");
        muhat(j) = R_PosInf;
        continue;
      }
      
      unsigned int m = idxu - idxl;
      
      //Rcpp::Rcout << "m=" << m << std::endl;
      
      Eigen::VectorXd w_tmp(m);
      Eigen::VectorXd z_tmp(m);
      Eigen::VectorXd x_tmp(m);
      
      for (unsigned int k = 0; k < m; ++k){
        x_tmp(k) = x(idxl+k);
        w_tmp(k) = w(idxl+k);
        z_tmp(k) = z(idxl+k);
      }
      
      //Rcpp::Rcout << "x_tmp.size()=" << x_tmp.size() << std::endl;
      
      muhat(j) = locpoly(x_tmp,z_tmp,w_tmp,h,newx(j),d,kcode);
      
    }
  }
  
  
  return muhat;
}