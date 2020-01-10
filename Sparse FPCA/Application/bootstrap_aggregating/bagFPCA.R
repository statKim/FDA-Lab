#---------------------------------------------------------------
# Bootstrap aggregated classification models using sparse FPCA
#---------------------------------------------------------------
require(fdapace)
require(e1071)
require(class)
require(MASS)
require(doParallel)

# The packages to use "foreach" statement
packages <- c("fdapace","e1071","class","MASS")

# # parallel computing setting
# ncores <- detectCores() - 2
# cl <- makeCluster(ncores)
# registerDoParallel(cl)
# stopCluster(cl)

### majority voting functon
majority_vote <- function(x) {
  key <- sort(unique(x))
  val <- key[which.max( table(x) )]
  return(val)
}

### The LSE-based weighting("Support Vector Machine Ensemble with Bagging"(2003), Kim et al.)
LSE_weight <- function(x) {
  A <- ifelse(x == 1, -1, 1)   # transform to -1 and 1
  w <- ginv(A) %*% ifelse(as.numeric(y.test) == 1, -1, 1)   # using generalized inverse(Not symmetric)
  val <- factor(ifelse(A %*% w < 0, 0, 1), levels=c(0,1))
  return(val)
}

### Fit the bagged classifier using sparse FPCA
# - Provide classification models
#   - Logistic regression("glm" function)
#   - SVM("svm" function in "e1071" package)
#   - KNN("knn" function in "class" package)
#   - LDA("lda" function in "e1071" package)
#   - qDA("qda" function in "e1071" package)
#   - Naive bayes("naiveBayes" function in "e1071" package)
# - Input parameters
#   - X.train : input data for "FPCA" function("fdapace" package)
#   - y.train : factor type vector of train set
#   - FPC.method : data type of FPCA("sparse", "dense), "sparse" is default. 
#   - selectK : the method of select the number of FPCs("FVE","AIC","BIC"), "FVE" is default.
#   - method : classification method("logit", "svm", "knn)
#   - kernel : (only method="svm")option of "svm" function("radial","linear","polynomial","sigmoid"), "radial" is default.
#   - B : the number of bootstrap replication
#   - ncores : the number of cores to use parallel computing
fit.fpca.bag <- function(X.train, y.train, FPC.method="Sparse", selectK="FVE", method="logit", kernel=NULL, B=10, ncores=NULL) {
  set.seed(777)
  
  # parallel computing setting
  ncores <- ceiling(detectCores() / 3)
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  # Bootstrap aggregating
  N <- length(y.train)
  y.pred <- foreach(b=1:B, .packages=packages) %dopar% {
    # fit FPCA for bootstrapped data
    boot.ind <- sample(1:N, N, replace=T)
    fpca.fit <- FPCA(X.train$Ly[boot.ind], 
                     X.train$Lt[boot.ind], 
                     optns=list(dataType=FPC.method,
                                methodXi="CE",
                                methodSelectK=selectK,
                                FVEthreshold=0.99,
                                verbose=T))
    k <- fpca.fit$selectK   # optimal number of PCs
    fpc.score <- fpca.fit$xiEst
    
    # train FPC scores
    train.fpc <- data.frame(y=y.train[boot.ind], 
                            x=fpc.score)
    colnames(train.fpc) <- c("y", paste("FPC", 1:k, sep=""))
    
    # fit the classification model
    if (method == "svm") {
      if (kernel == "linear") {
        tune.linear <- tune.svm(y ~ ., data=train.fpc, 
                                kernel="linear",
                                cost=c(10^(-3:1), 2^(5:10)) )$best.model
        fit.class <- svm(y ~ ., train.fpc, 
                         kernel="linear", 
                         probability=T,
                         cost=tune.linear$cost)
      } else if (kernel == "sigmoid") {
        tune.sigmoid <- tune.svm(y ~ ., data=train.fpc, 
                                 kernel="sigmoid",
                                 cost=c(10^(-3:1), 2^(5:10)), 
                                 gamma=c(10^(-3:1), 2^(5:10)) )$best.model
        fit.class <- svm(y ~ ., train.fpc, 
                         kernel="sigmoid", 
                         probability=T,
                         cost=tune.sigmoid$cost,
                         gamma=tune.sigmoid$gamma)
      } else if (kernel == "polynomial") {
        tune.poly <- tune.svm(y ~ ., data=train.fpc, 
                              kernel="polynomial",
                              cost=c(10^(-3:1), 2^(5:10)), 
                              gamma=c(10^(-3:1), 2^(5:10)) )$best.model
        fit.class <- svm(y ~ ., train.fpc,
                         kernel="polynomial",
                         probability=T,
                         cost=tune.poly$cost,
                         gamma=tune.poly$gamma,
                         degree=tune.poly$degree)
      } else {
        tune.radial <- tune.svm(y ~ ., data=train.fpc, 
                                kernel="radial",
                                cost=c(10^(-3:1), 2^(5:10)), 
                                gamma=c(10^(-3:1), 2^(5:10)) )$best.model
        fit.class <- svm(y ~ ., train.fpc, 
                         kernel="radial", 
                         probability=T,
                         cost=tune.radial$cost,
                         gamma=tune.radial$gamma)
      }
    } else if (method == "logit") {
      fit.class <- glm(y~., train.fpc, family=binomial)
    } else if (method == "knn") {
      fit.class <- tune.knn(x = train.fpc[, -1],
                            y = train.fpc[, 1],
                            k = 1:ceiling(N/2))$best.model
    } else if (method == "lda") {
      fit.class <- lda(y~., train.fpc)
    } else if (method == "qda") {
      fit.class <- qda(y~., train.fpc)
    } else if (method == "naivebayes") {
      fit.class <- naiveBayes(y~., train.fpc)
    }
    
    # OOB FPC scores
    fpc.score.oob <- predict(fpca.fit, 
                             X.train$Ly[-unique(boot.ind)], 
                             X.train$Lt[-unique(boot.ind)], 
                             K=k, 
                             xiMethod="CE")
    oob.fpc <- data.frame(y=y.train[-unique(boot.ind)], 
                          x=fpc.score.oob)
    colnames(oob.fpc) <- c("y", paste("FPC", 1:k, sep=""))
    
    # OOB error rate
    if (method == "svm") {
      pred.oob <- predict(fit.class, oob.fpc)
    } else if (method == "logit") {
      pred.oob <- predict(fit.class, oob.fpc, type="response")
      pred.oob <- factor(ifelse(pred.oob > 0.5, 1, 0), levels=c(0, 1))
    } else if (method == "knn") {
      pred.oob <- knn(train = train.fpc, test = oob.fpc, cl = train.fpc$y, k = fit.class$k)
    } else if (method %in% c("lda","qda")) {
      pred.oob <- predict(fit.class, oob.fpc)$class
    } else if (method == "naivebayes") {
      pred.oob <- predict(fit.class, oob.fpc, type='class')
    }
    oob.error <- mean(pred.oob != oob.fpc$y)
    
    # training error rate
    if (method == "svm") {
      pred.train <- predict(fit.class, train.fpc)
    } else if (method == "logit") {
      pred.train <- predict(fit.class, train.fpc, type="response")
      pred.train <- factor(ifelse(pred.train > 0.5, 1, 0), levels=c(0, 1))
    } else if (method == "knn") {
      pred.train <- knn(train = train.fpc, test = train.fpc, cl = train.fpc$y, k = fit.class$k)
    } else if (method %in% c("lda","qda")) {
      pred.train <- predict(fit.class, train.fpc)$class
    } else if (method == "naivebayes") {
      pred.train <- predict(fit.class, train.fpc, type='class')
    }
    train.error <- mean(pred.train != train.fpc$y)
    
    return( list(fpca.fit = fpca.fit,
                 K = k,
                 boot.index = boot.ind,
                 model = fit.class,
                 oob.error = oob.error,
                 train.error = train.error) )
  }
  
  # define method of SVM for kernel
  if (method == "svm") {
    if (is.null(kernel)) {
      method <- "svm.radial"
    } else {
      method <- paste(method, kernel, sep="")
    }
  }

  
  K <- sapply(y.pred, function(x) { x$K }) 
  fpca.fit <- lapply(y.pred, function(x) { x$fpca.fit })
  boot.index <- lapply(y.pred, function(x) { x$boot.index })
  model <- lapply(y.pred, function(x) { x$model })
  OOB.error <- sapply(y.pred, function(x) { x$oob.error }) 
  train.error <- sapply(y.pred, function(x) { x$train.error }) 
  
  result <- list(FPC.method = FPC.method,
                 method = method,
                 fpca.fit = fpca.fit,
                 K = K,
                 boot.index = boot.index,
                 model = model,
                 B = B,
                 OOB.error = OOB.error,
                 train.error = train.error)
  class(result) <- "fpca.bag.classifier"   # define class of object
  
  stopCluster(cl)    # End the parallel computing
  
  return( result )
}


### predict function of "fpca.bag.classifier" object(the bagged classifier using FPCA)
# - X.test : input data for "FPCA" function("fdapace" package)
predict.fpca.bag.classifier <- function(fit.fpca.bag.obj, X.test, y.test, probability=F, ncores=NULL) {
  # the number of bootstrap resamples
  B <- fit.fpca.bag.obj$B
  
  # parallel computing setting
  ncores <- ceiling(detectCores() / 3)
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  # predict for each model
  if (fit.fpca.bag.obj$method == "svm") {
    y.pred <- foreach(b=1:B, .combine="cbind", .packages=packages) %dopar% {
      fpca.fit <- fit.fpca.bag.obj$fpca.fit[[b]]
      K <- fit.fpca.bag.obj$K[b]
      model <- fit.fpca.bag.obj$model[[b]]
      
      # test FPC scores
      fpc.score.test <- predict(fpca.fit, X.test$Ly, X.test$Lt, K=K, xiMethod="CE")
      test.fpc <- data.frame(y=y.test, 
                             x=fpc.score.test)
      colnames(test.fpc) <- c("y", paste("FPC", 1:K, sep=""))
      
      pred.svm <- predict(model, test.fpc, probability=probability)
      if (probability == T) {
        pred.svm <- pred.svm[, 1]
      }
      return(pred.svm)
    }
  } else if (fit.fpca.bag.obj$method == "logit") {
    y.pred <- foreach(b=1:B, .combine="cbind", .packages=packages) %dopar% {
      fpca.fit <- fit.fpca.bag.obj$fpca.fit[[b]]
      K <- fit.fpca.bag.obj$K[b]
      model <- fit.fpca.bag.obj$model[[b]]
      
      # test FPC scores
      fpc.score.test <- predict(fpca.fit, X.test$Ly, X.test$Lt, K=K, xiMethod="CE")
      test.fpc <- data.frame(y=y.test, 
                             x=fpc.score.test)
      colnames(test.fpc) <- c("y", paste("FPC", 1:K, sep=""))
      
      pred.logit <- predict(model, test.fpc, type="response")
      if (probability == F) {
        pred.logit <- factor(ifelse(pred.logit > 0.5, 1, 0), levels=c(0, 1))
      }
      return(pred.logit)
    }
  } else if (fit.fpca.bag.obj$method == "knn") {
    y.pred <- foreach(b=1:B, .combine="cbind", .packages=packages) %dopar% {
      fpca.fit <- fit.fpca.bag.obj$fpca.fit[[b]]
      K <- fit.fpca.bag.obj$K[b]
      model <- fit.fpca.bag.obj$model[[b]]
      
      # test FPC scores
      fpc.score.test <- predict(fpca.fit, X.test$Ly, X.test$Lt, K=K, xiMethod="CE")
      test.fpc <- data.frame(y=y.test, 
                             x=fpc.score.test)
      colnames(test.fpc) <- c("y", paste("FPC", 1:K, sep=""))
        
      pred.knn <- knn(train = cbind(y=model$cl,
                                    model$train), test = test.fpc, cl = model$cl, k = model$k)
      return(pred.knn)
    }
  } else if (fit.fpca.bag.obj$method %in% c("lda","qda")) {
    y.pred <- foreach(b=1:B, .combine="cbind", .packages=packages) %dopar% {
      fpca.fit <- fit.fpca.bag.obj$fpca.fit[[b]]
      K <- fit.fpca.bag.obj$K[b]
      model <- fit.fpca.bag.obj$model[[b]]
      
      # test FPC scores
      fpc.score.test <- predict(fpca.fit, X.test$Ly, X.test$Lt, K=K, xiMethod="CE")
      test.fpc <- data.frame(y=y.test, 
                             x=fpc.score.test)
      colnames(test.fpc) <- c("y", paste("FPC", 1:K, sep=""))
      
      pred.lda.qda <- predict(model, test.fpc)$class
      return(pred.lda.qda)
    }
  } else if (fit.fpca.bag.obj$method == "naivebayes") {
    y.pred <- foreach(b=1:B, .combine="cbind", .packages=packages) %dopar% {
      fpca.fit <- fit.fpca.bag.obj$fpca.fit[[b]]
      K <- fit.fpca.bag.obj$K[b]
      model <- fit.fpca.bag.obj$model[[b]]
      
      # test FPC scores
      fpc.score.test <- predict(fpca.fit, X.test$Ly, X.test$Lt, K=K, xiMethod="CE")
      test.fpc <- data.frame(y=y.test, 
                             x=fpc.score.test)
      colnames(test.fpc) <- c("y", paste("FPC", 1:K, sep=""))
      
      pred.nb <- predict(model, test.fpc, type="class")
      return(pred.nb)
    }
  }
  
  if (probability == T) {
    pred <- rowMeans(y.pred)
  } else {
    # # marjority vote
    # pred <- apply(y.pred,
    #               1,
    #               majority_vote)
    # pred <- factor(pred-1, levels=c(0,1))
    
    # The LSE-based weighting
    pred <- LSE_weight( sapply(y.pred, function(x){ cbind( x[[j]] ) }) )
  }
  
  stopCluster(cl)   # End the parallel computing
  
  return(pred)
}


### print function of "fpca.bag.classifier" object
print.fpca.bag.classifier <- function(obj) {
  cat( "The bootstrap aggregating classifier using", obj$FPC.method, "FPCA \n",
       # "Type :", class(obj), "\n",
       obj$B, "Iterations", "\n", 
       "FPCA method :", obj$FPC.method, "\n" ,
       "Classification method :", obj$method, "\n",
       "OOB error rate :", mean(obj$OOB.error), "\n",
       "training error rate :", mean(obj$train.error) )
}


# Not bagged classifier
fit.fpca.classifier <- function(X.train, y.train, FPC.method="Sparse", selectK="FVE", method="logit", kernel=NULL, ncores=NULL) {
  set.seed(1000)
  
  # parallel computing setting
  ncores <- ceiling(detectCores() / 3)
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  N <- length(y.train)
  
  # fit FPCA for bootstrapped data
  boot.ind <- sample(1:N, N, replace=T)
  fpca.fit <- FPCA(X.train$Ly, 
                   X.train$Lt, 
                   optns=list(dataType=FPC.method,
                              methodXi="CE",
                              methodSelectK=selectK,
                              FVEthreshold=0.99,
                              verbose=T))
  k <- fpca.fit$selectK   # optimal number of PCs
  fpc.score <- fpca.fit$xiEst
  
  # train FPC scores
  train.fpc <- data.frame(y=y.train, 
                          x=fpc.score)
  colnames(train.fpc) <- c("y", paste("FPC", 1:k, sep=""))
  
  if (method == "svm") {
    if (kernel == "linear") {
      tune.linear <- tune.svm(y ~ ., data=train.fpc, 
                              kernel="linear",
                              cost=c(10^(-3:1), 2^(5:10)) )$best.model
      fit.class <- svm(y ~ ., train.fpc, 
                       kernel="linear", 
                       probability=T,
                       cost=tune.linear$cost)
    } else if (kernel == "sigmoid") {
      tune.sigmoid <- tune.svm(y ~ ., data=train.fpc, 
                               kernel="sigmoid",
                               cost=c(10^(-3:1), 2^(5:10)), 
                               gamma=c(10^(-3:1), 2^(5:10)) )$best.model
      fit.class <- svm(y ~ ., train.fpc, 
                       kernel="sigmoid", 
                       probability=T,
                       cost=tune.sigmoid$cost,
                       gamma=tune.sigmoid$gamma)
    } else if (kernel == "polynomial") {
      fit.svm.poly <- tune.svm(y ~ ., data=train.fpc, 
                               kernel="polynomial",
                               cost=c(10^(-3:1), 2^(5:10)), 
                               gamma=c(10^(-3:1), 2^(5:10)) )$best.model
      fit.class <- svm(y ~ ., train.fpc,
                       kernel="polynomial",
                       probability=T,
                       cost=tune.poly$cost,
                       gamma=tune.poly$gamma,
                       degree=tune.poly$degree)
    } else {
      tune.radial <- tune.svm(y ~ ., data=train.fpc, 
                              kernel="radial",
                              cost=c(10^(-3:1), 2^(5:10)), 
                              gamma=c(10^(-3:1), 2^(5:10)) )$best.model
      fit.class <- svm(y ~ ., train.fpc, 
                       kernel="radial", 
                       probability=T,
                       cost=tune.radial$cost,
                       gamma=tune.radial$gamma)
    }
  } else if (method == "logit") {
    fit.class <- glm(y~., train.fpc, family=binomial)
  } else if (method == "knn") {
    fit.class <- tune.knn(x = train.fpc[, -1],
                          y = train.fpc[, 1],
                          k = 1:ceiling(N/2))$best.model
  } else if (method == "lda") {
    fit.class <- lda(y~., train.fpc)
  } else if (method == "qda") {
    fit.class <- qda(y~., train.fpc)
  } else if (method == "naivebayes") {
    fit.class <- naiveBayes(y~., train.fpc)
  }
  
  result <- list(FPC.method=FPC.method,
                 method = method,
                 fpca.fit = fpca.fit,
                 K = k,
                 model = fit.class)
  class(result) <- "fpca.classifier"   # define class of object
  
  stopCluster(cl)    # End the parallel computing
  
  return( result )
}

### predict function of "fpca.classifier" object(the classifier using FPCA)
# - X.test : input data for "FPCA" function("fdapace" package)
predict.fpca.classifier <- function(fit.fpca.obj, X.test, y.test, probability=F, ncores=NULL) {
  # the number of bootstrap resamples
  B <- fit.fpca.obj$B
  
  # parallel computing setting
  ncores <- ceiling(detectCores() / 3)
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  fpca.fit <- fit.fpca.obj$fpca.fit
  K <- fit.fpca.obj$K
  model <- fit.fpca.obj$model
  
  # test FPC scores
  fpc.score.test <- predict(fpca.fit, X.test$Ly, X.test$Lt, K=K, xiMethod="CE")
  test.fpc <- data.frame(y=y.test, 
                         x=fpc.score.test)
  colnames(test.fpc) <- c("y", paste("FPC", 1:K, sep=""))
  
  
  # predict for model
  if (fit.fpca.obj$method == "svm") {
    pred <- predict(model, test.fpc, probability=probability)
    if (probability == T) {
      pred <- pred[, 1]
    }
  } else if (fit.fpca.obj$method == "logit") {
    pred <- predict(model, test.fpc, type="response")
    if (probability == F) {
      pred <- factor(ifelse(pred > 0.5, 1, 0), levels=c(0, 1))
    }
  } else if (fit.fpca.obj$method == "knn") {
    pred <- knn(train = cbind(y=model$cl,
                              model$train), test = test.fpc, cl = model$cl, k = model$k)
  } else if (fit.fpca.obj$method %in% c("lda","qda")) {
    pred <- predict(model, test.fpc)$class
  } else if (fit.fpca.obj$method == "naivebayes") {
    pred <- predict(model, test.fpc, type="class")
  }

  stopCluster(cl)   # End the parallel computing
  
  return(pred)
}




# -----------
# Example
# -----------
setwd("C:\\Users\\user\\Desktop\\KHS\\FDA-Lab\\Sparse FPCA\\Application\\bootstrap_aggregating")

library(doParallel)   # parallel computing
library(fdapace)
library(e1071)
library(class)
library(MASs)

# load simulated data
load("../simulated_data.RData")

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
s <- 10   # sparsity
train.data <- Sparsify(X.train, time, sparsity = s)
test.data <- Sparsify(X.test, time, sparsity = s)

ncores <- detectCores() - 2

# Bagged logistic regression
bag.logit <- fit.fpca.bag(train.data, y.train, method="logit", B=100, selectK="FVE", ncores=ncores)
pred <- predict(bag.logit, test.data, y.test)
# pred <- predict.fit.fpca.bag(bag.logit, test.data, y.test, probability=T)
# pred <- factor(ifelse(pred > 0.5, 1, 0), levels=c(0, 1))
acc.logit <- mean(pred == y.test)
acc.logit

# Bagged SVM
bag.svm <- fit.fpca.bag(train.data, y.train, method="svm", kernel="radial", B=100, selectK="BIC", ncores=ncores)
pred <- predict(bag.svm, test.data, y.test)
acc.svm <- mean(pred == y.test)
acc.svm

# Bagged KNN
bag.knn <- fit.fpca.bag(train.data, y.train, method="knn", B=100, selectK="FVE", ncores=ncores)
pred <- predict(bag.knn, test.data, y.test)
acc.knn <- mean(pred == y.test)
acc.knn

# Bagged LDA
bag.lda <- fit.fpca.bag(train.data, y.train, method="lda", B=100, selectK="FVE", ncores=ncores)
pred <- predict(bag.lda, test.data, y.test)
acc.lda <- mean(pred == y.test)
acc.lda

# Bagged QDA
bag.qda <- fit.fpca.bag(train.data, y.train, method="qda", B=100, selectK="FVE", ncores=ncores)
pred <- predict(bag.qda, test.data, y.test)
acc.qda <- mean(pred == y.test)
acc.qda

# Bagged Naive bayes
bag.nb <- fit.fpca.bag(train.data, y.train, method="naivebayes", B=100, selectK="FVE", ncores=ncores)
pred <- predict(bag.nb, test.data, y.test)
acc.nb <- mean(pred == y.test)
acc.nb



# train, test split
ncores <- detectCores() - 2

i <- 1
set.seed(100)
y <- factor(c(rep(0, 100), rep(1, 100)), level=c(0,1))
ind <- sample(1:length(y), 100)   # train set index
X.train <- as.matrix( X.curves[[i]][ind, ] )
X.test <- as.matrix( X.curves[[i]][-ind, ] )
y.train <- y[ind]
y.test <- y[-ind]

for (s in 2:18) {
  print(s)
  set.seed(100)
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
  # test FPC scores
  fpc.score.test <- predict(fpca.fit, test.data$Ly, test.data$Lt, K=k, xiMethod="CE")
  test.fpc <- data.frame(y=y.test, 
                         x=fpc.score.test)
  colnames(test.fpc) <- c("y", paste("FPC", 1:k, sep=""))
  
  
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
  # fit.svm.poly <- tune(svm, y ~ ., data=train.fpc, kernel="poly",
  #                      ranges=list(cost=c(10^(-3:1), 2^(5:10)), 
  #                                  gamma=c(10^(-3:1), 2^(5:10))))$best.model
  
  # predict
  pred.logit <- predict(fit.logit, test.fpc, type="response")
  pred.logit <- factor(ifelse(pred.logit > 0.5, 1, 0), levels=c(0, 1))
  pred.svm.linear <- predict(fit.svm.linear, test.fpc)
  pred.svm.radial <- predict(fit.svm.radial, test.fpc)
  pred.svm.sigmoid <- predict(fit.svm.sigmoid, test.fpc)
  # pred.svm.poly <- predict(fit.svm.poly, test.fpc)
  
  acc.logit <- mean(pred.logit == y.test)
  acc.svm.linear <- mean(pred.svm.linear == y.test)
  acc.svm.radial <- mean(pred.svm.radial == y.test) 
  acc.svm.sigmoid <- mean(pred.svm.sigmoid == y.test) 
  # acc.svm.polynomial <- mean(acc.svm.polynomial == y.test)
  
  # KNN
  fit.knn <- tune.knn(x = train.fpc[, -1],
                      y = train.fpc[, 1],
                      k = 1:ceiling(N/2))$best.model
  pred.knn <- knn(train = cbind(y=fit.knn$cl,
                                fit.knn$train), 
                  test = test.fpc, cl = fit.knn$cl, k = fit.knn$k)
  acc.knn <- mean(pred.knn == y.test)
  # LDA
  fit.lda <- lda(y~., train.fpc)
  pred.lda <- predict(fit.lda, test.fpc)$class
  acc.lda <- mean(pred.lda == y.test)
  # QDA
  fit.qda <- qda(y~., train.fpc)
  pred.qda <- predict(fit.qda, test.fpc)$class
  acc.qda <- mean(pred.qda == y.test)
  # Naive bayes
  fit.nb <- naiveBayes(y~., train.fpc)
  pred.nb <- predict(fit.nb, test.fpc, type="class")
  acc.nb <- mean(pred.nb == y.test)


  if (s == 2) {
    result <- c(acc.logit, acc.svm.linear, acc.svm.radial, acc.svm.sigmoid,  
                acc.knn, acc.lda, acc.qda, acc.nb)
  } else {
    result <- rbind(result,
                    c(acc.logit, acc.svm.linear, acc.svm.radial, acc.svm.sigmoid,  
                      acc.knn, acc.lda, acc.qda, acc.nb))
  }
}

save(result, file="results/result_fpc.RData")
