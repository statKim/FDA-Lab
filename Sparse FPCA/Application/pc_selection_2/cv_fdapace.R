#---------------------------
# Cross-validation for FPCA
#---------------------------
require(fdapace)

### generalized inverse approach (approximate correct PRESS)
# https://stats.stackexchange.com/questions/93845/how-to-perform-cross-validation-for-pca-to-determine-the-number-of-principal-com/115477#115477
cv.pca.sparse <- function(X, packages=c("fdapace"), plot=T) {
  fit <- FPCA(X$Ly, 
              X$Lt, 
              optns=list(dataType="Sparse",
                         methodXi="CE",
                         FVEthreshold=1,   # fit할 수 있는 최대의 k까지 fit
                         verbose=T))
  
  N <- length(X$Ly)   # number of curves
  
  # calculate AIC, BIC
  aic <- c()
  bic <- c()
  # AIC <- foreach(k=1:length(fit$lambda), .combine="c", .packages=packages) %dopar% {
  for (k in 1:length(fit$lambda)) {
    aic[k] <- GetLogLik(fit, k) + 2*k
    bic[k] <- GetLogLik(fit, k) + k*log(N)
  }
  
  # approximate PCA PRESS
  cv.error <- foreach(k=1:length(fit$lambda), .combine="c", .packages=packages) %dopar% {
    U <- fit$phi[, 1:k]   # eigenvector matrix without ith curve
    if (!is.matrix(U)) { U <- matrix(U, ncol=1) }  # transform to matrix
    uu.t <- U %*% t(U)
    w <- diag(1, nrow(uu.t)) - uu.t + diag(diag(uu.t))
    
    # calculate PRESS
    err <- mapply(function(t, y){
      time.ind <- sapply(t, function(x){ which.min( abs(x - fit$workGrid) ) })   # time point에 해당되는 grid 찾기
      sum( (w[, time.ind] %*% y)^2 )
    },
    X$Lt,
    X$Ly)
    
    return( sum(err) )
  }
  
  # 각 measure로 plot
  if (plot==T) {
    par(mfrow=c(1,3))
    plot(1:length(aic), 
         aic, 
         type="o",
         xlab="Nunber of PCs",
         ylab="AIC",
         main="AIC")
    points(which.min(aic),
           aic[which.min(aic)],
           col="red",
           pch=19,
           cex=2)    
    
    plot(1:length(bic), 
         bic, 
         type="o",
         xlab="Nunber of PCs",
         ylab="BIC",
         main="BIC")
    points(which.min(bic),
           bic[which.min(bic)],
           col="red",
           pch=19,
           cex=2)
    
    plot(1:length(cv.error), 
         cv.error, 
         type="o",
         xlab="Nunber of PCs",
         ylab="Cross-validation error",
         main="Leave-one-curve-out CV")
    points(which.min(cv.error),
           cv.error[which.min(cv.error)],
           col="red",
           pch=19,
           cex=2)
  }
  
  
  return( list(cv.error=cv.error,
               AIC=aic,
               BIC=bic,
               selected.K=data.frame(AIC=which.min(aic),
                                     BIC=which.min(bic),
                                     LOOCV=which.min(cv.error)),
               model=fit) )
  
  
  # cv.error <- foreach(i=1:length(X$Lt), .combine="rbind", .packages=packages) %dopar% {
  #   x.i <- X$Ly[[i]]   # ith curve
  #   t.i <- X$Lt[[i]]   # ith curve's time points
  #   
  #   fit <- FPCA(X$Ly[-i], 
  #               X$Lt[-i], 
  #               optns=list(dataType="Sparse",
  #                          methodXi="CE",
  #                          FVEthreshold=1,   # fit할 수 있는 최대의 k까지 fit
  #                          rho="cv",
  #                          verbose=T))
  #   
  #   time.grid <- fit$workGrid   # time grid
  #   err <- c()
  #   for (k in 1:length(fit$lambda)) {
  #     U <- fit$phi[, 1:k]   # eigenvector matrix without ith curve
  #     dec <- svd(U[-i, ])
  #     U.plus <- dec$v %*% diag(1/dec$d, ncol(dec$v), ncol(dec$u)) %*% t(dec$u)
  #     
  #     # 아마 time point 찾아서 계산해줘야 될듯
  #     
  #     time.ind <- time.grid[ sapply(t.i, function(x){ which.min( abs(x - time.grid) ) }) ]
  #     sub <- U %*% U.plus %*% x.i
  #     
  #     err[k] <- sum( (x.i - sub)^2 )
  #   }
  # }
  # 
  # cv.error <- data.frame(cv.error)
  # colnames(cv.error) <- 1:ncol(fit$phi)
  # return( colMeans(cv.error) )
}


### function from github of fdapace package
## https://github.com/functionaldata/tPACE/blob/master/R/GetLogLik.R
# K: input denoting number of components used
# returns -2 times log-likelihood
GetLogLik = function(fpcaObj, K, Ly = NULL, Lt = NULL){
  if(fpcaObj$optns$lean == TRUE && (is.null(Ly) || is.null(Lt))){
    stop("Option lean is TRUE, need input data Ly and measurement time list Lt to calculate log-likelihood.")
  }
  if(fpcaObj$optns$lean == FALSE){ # when input data is in fpcaObj
    Ly <- fpcaObj$inputData$Ly
    Lt <- fpcaObj$inputData$Lt
  }
  lambda = fpcaObj$lambda[1:K]
  sigma2 = fpcaObj$sigma2
  if(is.null(sigma2) && fpcaObj$optns$dataType == "Dense"){
    ymat = matrix(unlist(Ly),nrow=length(Ly), byrow=TRUE)
    sddiag = sqrt(diag(var(ymat)))
    sigma2 = sddiag*1e-4
    sigma2 = ConvertSupport(fromGrid = fpcaObj$obsGrid, toGrid = fpcaObj$workGrid, mu = sigma2)
  }
  logLik = 0
  phi = fpcaObj$phi[,1:K, drop=FALSE]
  
  if(fpcaObj$optns$dataType %in% c('Dense'
                                   #, 'DenseWithMV' # need extra imputation step
  )){
    if(K == 1){
      Sigma_y = phi %*% (lambda*diag(K)) %*% t(phi) + sigma2*diag(rep(1,nrow(phi)))
    } else {
      Sigma_y = phi %*% diag(lambda, length(lambda)) %*% t(phi) + sigma2*diag(rep(1,nrow(phi)))
    }
    detSigma_y = prod(c(lambda,rep(0,nrow(phi)-K))[1:length(lambda)]+sigma2)
    #detSigma_y = det(Sigma_y)
    if(detSigma_y == 0){
      logLik = NULL
      return(logLik)
    }
    # calculate loglikelihood via matrix multiplication
    ymatcenter = matrix(unlist(Ly)-fpcaObj$mu, nrow = length(Ly), byrow = TRUE)
    svd_Sigma_y = svd(Sigma_y)
    Sigma_y_inv = svd_Sigma_y$v %*% diag(1/svd_Sigma_y$d, length(svd_Sigma_y$d)) %*% t(svd_Sigma_y$u)
    logLik = sum(diag(t(Sigma_y_inv %*% t(ymatcenter)) %*% t(ymatcenter))) + length(Ly)*log(detSigma_y)
    return(logLik)
  } else { # Sparse case
    if(is.null(sigma2)){ sigma2 <- fpcaObj$rho }
    if(fpcaObj$optns$error == TRUE && sigma2 <= fpcaObj$rho){
      # especially for the case when sigma2 is estimated to be <=0 and set to 1e-6
      sigma2 <- fpcaObj$rho
    }
    for(i in 1:length(Ly)){
      if(length(Lt[[i]]) == 1){
        phi_i = t(as.matrix(ConvertSupport(fromGrid = fpcaObj$workGrid, toGrid = Lt[[i]],
                                           phi = phi)))        
      } else {
        phi_i = ConvertSupport(fromGrid = fpcaObj$workGrid, toGrid = Lt[[i]],
                               phi = phi)
      }
      mu_i = ConvertSupport(fromGrid = fpcaObj$workGrid, toGrid = Lt[[i]],
                            mu = fpcaObj$mu)
      if(K == 1){
        Sigma_yi = phi_i %*% (lambda*diag(K)) %*% t(phi_i) + sigma2 * diag(rep(1,length(mu_i)))
      } else{
        Sigma_yi = phi_i %*% diag(lambda, length(lambda)) %*% t(phi_i) + sigma2 * diag(rep(1,length(mu_i)))
      }
      detSigma_yi = det(Sigma_yi)
      if(detSigma_yi == 0){
        logLik = NULL
        return(logLik)
      }
      invtempi = solve(Sigma_yi, Ly[[i]] - mu_i)
      logLik = logLik + log(detSigma_yi) + invtempi %*% (Ly[[i]] - mu_i)
    }
    return(logLik)
  }
}