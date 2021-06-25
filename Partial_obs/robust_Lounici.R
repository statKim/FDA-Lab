################################
### Lounichi (2014)
################################

### Frubenius norm
frobenius_norm <- function(A) {
  # svd_fit <- svd(A)
  # return( sqrt(sum(svd_fit$d^2)) )
  return( sqrt(sum(A^2)) )
}

### Nuclear norm
nuclear_norm <- function(A) {
  return( sum(diag(A)) )
}

### Sup-norm (Infinite norm)
sup_norm <- function(A) {
  svd_fit <- svd(A)
  return(svd_fit$d[1])
  
  # eig_fit <- eigen(A)
  # return(eig_fit$values[1])
}

### Covariance LASSO estimator via Singular Value Theresholding (SVT) algorithm
### Matrix LASSO
cov.lasso <- function(X, lambda = NULL, 
                      cov = NULL,
                      robust = FALSE, smooth = FALSE, noise.var = 0) {
  # proportion of observed entries
  delta <- 1 - sum(is.na(X)) / length(X)
  
  if (is.null(cov)) {
    if (robust == TRUE) {
      Sigma_delta <- var.rob.missfd(X, smooth = smooth, noise.var = noise.var)
    } else {
      X[is.na(X)] <- 0
      Sigma_delta <- cov(X)
      diag(Sigma_delta) <- diag(Sigma_delta) - noise.var
    }
  } else {
    Sigma_delta <- cov
    diag(Sigma_delta) <- diag(Sigma_delta) - noise.var
  }

  # equation (1.4)
  Sigma_tilde <- (1/delta - 1/delta^2)*diag(diag(Sigma_delta)) + 1/delta^2*Sigma_delta
  
  if (is.null(lambda)) {
    C <- 10
    lambda <- C*sqrt(sum(diag(Sigma_tilde))*sup_norm(Sigma_tilde))/delta * sqrt(log(2*ncol(X)) / nrow(X))
  }

  # # Singular value thresholding (SVT) algorithm
  # Sigma_hat <- filling::fill.SVT(Sigma_tilde,
  #                                maxiter = 1000,
  #                                lambda = lambda/2)$X

  # Optimal solution of equation (5.4) in page 1040 of paper
  eig_obj <- eigen(Sigma_tilde)
  ind_pos <- 1:max(which(eig_obj$values > 0))
  d <- eig_obj$values[ind_pos]
  u <- eig_obj$vectors[, ind_pos]
  
  Sigma_hat <- u %*% diag(ifelse(d > lambda/2, d-lambda/2, 0)) %*% t(u)
  
  return(Sigma_hat)
}


cov.lasso2 <- function(X, lambda = NULL, robust = FALSE, 
                       smooth = FALSE, noise.var = 0, C = 10) {
  delta <- 1 - sum(is.na(X)) / length(X)   # missingness
  
  if (robust == TRUE) {
    Sigma_delta <- var.rob.missfd(X, smooth = smooth, noise.var = noise.var)
  } else {
    X[is.na(X)] <- 0
    Sigma_delta <- cov(X)
  }
  
  Sigma_tilde <- (1/delta - 1/delta^2)*diag(diag(Sigma_delta)) + 1/delta^2*Sigma_delta
  
  if (is.null(lambda)) {
    lambda <- C*sqrt(sum(diag(Sigma_tilde))*sup_norm(Sigma_tilde))/delta * sqrt(log(2*ncol(X)) / nrow(X))
  }
  

  ### 박연주 교수님 코드
  Sigma.tilde <- Sigma_tilde
  # optimization solution based on page 1040 of the paper
  res=svd(Sigma.tilde)
  e.val = res$d
  e.vec = res$u
  ind = which(e.val < lambda/2)
  ind2 = which(e.val >= lambda/2)
  
  W = matrix(0,nrow=length(t), ncol=length(t))
  for(i in 1:length(ind)){
    tmp.ind = ind[i]
    W = W + 2*e.val[tmp.ind]*((e.vec[,tmp.ind])%*%t(e.vec[,tmp.ind]))/lambda
  }
  
  V.hat0 = matrix(0,nrow=length(t), ncol=length(t))
  for(i in 1:length(ind2)){
    tmp.ind = ind2[i]
    V.hat0 = V.hat0 +((e.vec[,tmp.ind])%*%t(e.vec[,tmp.ind]))
  }
  
  V.hat = V.hat0 + W 
  
  Sigma_hat = (2*Sigma.tilde - lambda*V.hat)/2
  
  return(Sigma_hat)
}



####################
### Test
####################

lambda <- 0.01
par(mfcol = c(2, 3))
# GA::persp3D(gr, gr, cov.true,
#             main = "True", xlab = "", ylab = "", zlab = "",
#             theta = -70, phi = 30, expand = 1)

cov.kraus <- var.missfd(x)
GA::persp3D(gr, gr, cov.kraus,
            main = "Kraus", xlab = "", ylab = "", zlab = "",
            theta = -70, phi = 30, expand = 1)

cov.lounici <- cov.lasso(x, lambda = lambda, robust = F) 
GA::persp3D(gr, gr, cov.lounici,
            main = "Lounici", xlab = "", ylab = "", zlab = "",
            theta = -70, phi = 30, expand = 1)

cov.Mest <- var.rob.missfd(x, smooth = F)
GA::persp3D(gr, gr, cov.Mest,
            main = "M-est", xlab = "", ylab = "", zlab = "",
            theta = -70, phi = 30, expand = 1)

cov.svt <- cov.lasso(x, lambda = lambda, 
                     cov = cov.Mest, robust = T) 
GA::persp3D(gr, gr, cov.svt,
            main = "M-Lounici", xlab = "", ylab = "", zlab = "",
            theta = -70, phi = 30, expand = 1)

cov.Mest.sm <- var.rob.missfd(x, smooth = T)
GA::persp3D(gr, gr, cov.Mest.sm,
            main = "M-est(smooth)", xlab = "", ylab = "", zlab = "",
            theta = -70, phi = 30, expand = 1)

cov.svt.sm <- cov.lasso(x, lambda = lambda, 
                        cov = cov.Mest.sm, robust = T, smooth = T) 
GA::persp3D(gr, gr, cov.svt.sm,
            main = "M-Lounici(smooth)", xlab = "", ylab = "", zlab = "",
            theta = -70, phi = 30, expand = 1)


