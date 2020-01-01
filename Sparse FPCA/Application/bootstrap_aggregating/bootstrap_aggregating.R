setwd("C:\\Users\\user\\Desktop\\KHS\\FDA-Lab\\Sparse FPCA\\Application")

library(tidyverse)
library(doParallel)   # parallel computing
library(fdapace)
library(e1071)

# load simulated data
load("simulated_data.RData")

# parallel computing setting
ncores <- detectCores() - 2
registerDoParallel(ncores)

packages <- c("fdapace","e1071")   # foreach에서 사용할 package 정의

s <- 10   # sparsity

system.time({   # running time check
for (s in 2:18) {
  print(paste("s=", s, sep=""))
  for (i in 1:100) {
    ## generate data
    # i <- 1
    set.seed(100)
    y <- factor(c(rep(0, 100), rep(1, 100)), level=c(0,1))
    ind <- sample(1:length(y), 100)   # train set index
    X.train <- as.matrix( X.curves[[i]][ind, ] )
    X.test <- as.matrix( X.curves[[i]][-ind, ] )
    y.train <- y[ind]
    y.test <- y[-ind]
    
    # sparsify train, test set
    train.data <- Sparsify(X.train, time, sparsity = s)
    test.data <- Sparsify(X.test, time, sparsity = s)
    
    ## Bootstrap aggregating for classification
    B <- 100
    N <- nrow(X.train)
    y.pred <- foreach(b=1:B, .packages=packages) %dopar% {
      # fit FPCA for bootstrapped data
      boot.ind <- sample(1:N, N, replace=T)
      fpca.fit <- FPCA(train.data$Ly[boot.ind], 
                       train.data$Lt[boot.ind], 
                       optns=list(dataType="Sparse",
                                  methodXi="CE",
                                  # methodSelectK="BIC",
                                  FVEthreshold=0.99,
                                  verbose=T))
      k <- fpca.fit$selectK   # optimal number of PCs when using BIC
      fpc.score <- fpca.fit$xiEst
      
      # train FPC scores
      train.fpc <- data.frame(y=y.train, 
                              x=fpc.score)
      colnames(train.fpc) <- c("y", paste("FPC", 1:k, sep=""))
      
      # fit the classification model using bootstrapped FPC scores
      fit.logit <- glm(y~., train.fpc, family=binomial)
      fit.svm.linear <- tune(svm, y ~ ., data=train.fpc, kernel="linear",
                             ranges=list(cost=10^(-1:2), gamma=c(.5,1,2)))$best.model
      fit.svm.radial <- tune(svm, y ~ ., data=train.fpc, kernel="radial",
                             ranges=list(cost=10^(-1:2), gamma=c(.5,1,2)))$best.model
      fit.svm.sigmoid <- tune(svm, y ~ ., data=train.fpc, kernel="sigmoid",
                              ranges=list(cost=10^(-1:2), gamma=c(.5,1,2)))$best.model
      fit.svm.poly <- tune(svm, y ~ ., data=train.fpc, kernel="poly",
                           ranges=list(cost=10^(-1:2), gamma=c(.5,1,2)))$best.model
      
      # test FPC scores
      fpc.score.test <- predict(fpca.fit, test.data$Ly, test.data$Lt, K=k, xiMethod="CE")
      test.fpc <- data.frame(y=y.test, 
                             x=fpc.score.test)
      colnames(test.fpc) <- c("y", paste("FPC", 1:k, sep=""))
      
      # predict
      pred.logit <- predict(fit.logit, test.fpc, type="response")
      pred.logit <- factor(ifelse(pred.logit > 0.5, 1, 0), levels=c(0, 1))
      pred.svm.linear <- predict(fit.svm.linear, test.fpc)
      pred.svm.radial <- predict(fit.svm.radial, test.fpc)
      pred.svm.sigmoid <- predict(fit.svm.sigmoid, test.fpc)
      pred.svm.poly <- predict(fit.svm.poly, test.fpc)
      
      return( list(logit=pred.logit,
                   svm.linear=pred.svm.linear,
                   svm.radial=pred.svm.radial,
                   svm.sigmoid=pred.svm.sigmoid,
                   svm.poly=pred.svm.poly) )
    }
    
    
    # majority voting functon
    majority_vote <- function(x) {
      key <- sort(unique(x))
      val <- key[which.max( table(x) )]
      return(val)
    }
    
    # save the accuracy
    cname <- c("logit","svm.linear","svm.radial","svm.sigmoid","svm.poly")
    acc <- c()
    for (j in 1:5) {
      pred <- apply(sapply(y.pred, function(x){ cbind( x[[j]] ) }), 
                    1, 
                    majority_vote)
      pred <- factor(pred-1, levels=c(0,1))
      
      acc[j] <- mean(pred == y.test)
    }
    
    if (i == 1) {
      acc.result <- acc
    } else {
      acc.result <- rbind(acc.result, acc)
    }
  }
  
  acc.result <- colMeans(acc.result)
    
  if (s == 2) {
    result <- data.frame(matrix(acc.result, nrow=1))
    colnames(result) <- cname
  } else {
    result <- rbind(result, acc.result)
  }
  
  if (s == 18) { rownames(result) <- paste("s=", 2:18, sep="") }
}
})   # running time check


table(Pred=pred, 
      True=y.test)
mean(pred == y.test)



