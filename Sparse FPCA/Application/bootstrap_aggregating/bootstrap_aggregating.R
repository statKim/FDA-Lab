setwd("C:\\Users\\user\\Desktop\\KHS\\FDA-Lab\\Sparse FPCA\\Application\\bootstrap_aggregating")

library(tidyverse)
library(doParallel)   # parallel computing
library(fdapace)
library(e1071)

# load simulated data
load("../simulated_data.RData")

# parallel computing setting
ncores <- detectCores() - 2
registerDoParallel(ncores)

packages <- c("fdapace","e1071","class","MASS")   # foreach에서 사용할 package 정의

# majority voting functon
majority_vote <- function(x) {
  key <- sort(unique(x))
  val <- key[which.max( table(x) )]
  return(val)
}


## Bagging to curves
system.time({   # running time check
for (s in 2:18) {
  print(paste("s=", s, sep=""))
  # train, test split
  i <- 1
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

  # Bootstrap aggregating
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
    k <- fpca.fit$selectK   # optimal number of PCs
    fpc.score <- fpca.fit$xiEst
    
    # train FPC scores
    train.fpc <- data.frame(y=y.train[boot.ind], 
                            x=fpc.score)
    colnames(train.fpc) <- c("y", paste("FPC", 1:k, sep=""))
    
    # classifiers with bootstrapped FPC scores
    tune.linear <- tune.svm(y ~ ., data=train.fpc, 
                            kernel="linear",
                            cost=c(10^(-3:1), 2^(5:10)))$best.model
    tune.radial <- tune.svm(y ~ ., data=train.fpc, 
                            kernel="radial",
                            cost=c(10^(-3:1), 2^(5:10)), 
                            gamma=c(10^(-3:1), 2^(5:10)) )$best.model
    tune.sigmoid <- tune.svm(y ~ ., data=train.fpc, 
                             kernel="sigmoid",
                             cost=c(10^(-3:1), 2^(5:10)), 
                             gamma=c(10^(-3:1), 2^(5:10)) )$best.model
    
    fit.logit <- glm(y~., train.fpc, family=binomial)
    fit.svm.linear <- svm(y ~ ., train.fpc, 
                          kernel="linear", 
                          cost=tune.linear$cost)
    fit.svm.radial <- svm(y ~ ., train.fpc, 
                          kernel="radial", 
                          cost=tune.radial$cost,
                          gamma=tune.radial$gamma)
    fit.svm.sigmoid <- svm(y ~ ., train.fpc, 
                           kernel="sigmoid", 
                           cost=tune.sigmoid$cost,
                           gamma=tune.sigmoid$gamma)
    # fit.svm.poly <- svm(y ~ ., train.fpc,
    #                     kernel="polynomial",
    #                     cost=tune.poly$cost,
    #                     gamma=tune.poly$gamma,
    #                     degree=tune.poly$degree)
    fit.knn <- tune.knn(x = train.fpc[, -1],
                        y = train.fpc[, 1],
                        k = 1:ceiling(N/2))$best.model
    fit.lda <- lda(y~., train.fpc)
    fit.qda <- qda(y~., train.fpc)
    fit.nb <- naiveBayes(y~., train.fpc)

    
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
    # pred.svm.poly <- predict(fit.svm.poly, test.fpc)
    pred.knn <- knn(train = cbind(y=fit.knn$cl,
                                  fit.knn$train), 
                    test = test.fpc, cl = fit.knn$cl, k = fit.knn$k)
    pred.lda <- predict(fit.lda, test.fpc)$class
    pred.qda <- predict(fit.qda, test.fpc)$class
    pred.nb <- predict(fit.nb, test.fpc, type="class")

    return( list(logit = pred.logit,
                 svm.linear = pred.svm.linear,
                 svm.radial = pred.svm.radial,
                 svm.sigmoid = pred.svm.sigmoid,
                 knn = pred.knn,
                 lda = pred.lda,
                 qda = pred.qda,
                 nb = pred.nb) )
  }
  
  # save the accuracy
  cname <- c("logit","svm.linear","svm.radial","svm.sigmoid","knn","lda","qda","naivebayes")
  acc <- c()
  for (j in 1:8) {
    pred <- apply(sapply(y.pred, function(x){ cbind( x[[j]] ) }),
                  1,
                  majority_vote)
    pred <- factor(pred-1, levels=c(0,1))
    acc[j] <- mean(pred == y.test)
  }

  # apply(sapply(y.pred, function(x){ cbind( x[[j]] ) })-1,
  #       2,
  #       function(x){mean(x==y.test)}
  # )
  
  if (s == 2) {
    result <- data.frame(matrix(acc, nrow=1))
    colnames(result) <- cname
  } else {
    result <- rbind(result, acc)
  }
  
  if (s == 18) { rownames(result) <- paste("s=", 2:18, sep="") }
}
})   # running time check

save(result, file="results/result_bag_raw.RData")




## Not bootstrap aggregating
system.time({   # running time check
for (s in 2:18) {
  print(paste("s=", s, sep=""))
  # train, test split
  acc <- foreach(i=1:100, .combine="rbind", .packages=packages) %dopar% {
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
    
    # fit FPCA
    fpca.fit <- FPCA(train.data$Ly, 
                     train.data$Lt, 
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
    
    # classifiers with FPC scores
    fit.logit <- glm(y~., train.fpc, family=binomial)
    fit.svm.linear <- tune(svm, y ~ ., data=train.fpc, kernel="linear",
                           ranges=list(cost=c(10^(-3:1), 2^(5:10))))$best.model
    fit.svm.radial <- tune(svm, y ~ ., data=train.fpc, kernel="radial",
                           ranges=list(cost=c(10^(-3:1), 2^(5:10)), 
                                       gamma=c(10^(-3:1), 2^(5:10))))$best.model
    fit.svm.sigmoid <- tune(svm, y ~ ., data=train.fpc, kernel="sigmoid",
                            ranges=list(cost=c(10^(-3:1), 2^(5:10)), 
                                        gamma=c(10^(-3:1), 2^(5:10))))$best.model
    fit.svm.poly <- tune(svm, y ~ ., data=train.fpc, kernel="poly",
                         ranges=list(cost=c(10^(-3:1), 2^(5:10)), 
                                     gamma=c(10^(-3:1), 2^(5:10))))$best.model
    
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
    
    y.pred <- list(logit=pred.logit,
                   svm.linear=pred.svm.linear,
                   svm.radial=pred.svm.radial,
                   svm.sigmoid=pred.svm.sigmoid,
                   svm.poly=pred.svm.poly)
    
    # save the accuracy
    cname <- c("logit","svm.linear","svm.radial","svm.sigmoid","svm.poly")
    acc <- c()
    for (j in 1:5) {
      pred <- factor(y.pred[[j]], levels=c(0,1))
      acc[j] <- mean(pred == y.test)
    }
    
    return( acc ) 
  }

  
  if (s == 2) {
    result <- data.frame(matrix(colMeans(acc), nrow=1))
    colnames(result) <- cname
  } else {
    result <- rbind(result, colMeans(acc))
  }
  
  if (s == 18) { rownames(result) <- paste("s=", 2:18, sep="") }
}
})   # running time check

save(result, file="results/result_fpc.RData")
