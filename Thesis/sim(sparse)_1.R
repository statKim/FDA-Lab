#########################################
### Simulation
###   Probability-enhanced effective dimension reduction
###   for classifying sparse functional data
###   Yao et al.
###   Link: https://link.springer.com/content/pdf/10.1007%2Fs11749-015-0470-2.pdf
#########################################
setwd("C:\\Users\\user\\Desktop\\KHS\\FDA-Lab\\Thesis")

library(tidyverse)    # dplyr, ggplot2, etc
library(doParallel)   # parallel computing
library(fdapace)      # functional PCA
library(e1071)        # SVM, Naive bayes
library(MASS)         # LDA, QDA
library(data.table)   # list rbind
library(gridExtra)    # subplot in ggplot2
library(reshape2)     # melt function

# parallel computing setting
ncores <- detectCores() - 2
registerDoParallel(ncores)

packages <- c("fdapace","e1071","MASS","tidyverse")   # foreach에서 사용할 package 정의

# function of calculating mode
Mode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# res.sim <- list()  # classification result of each bootstrap resamples
result <- list()   # error rate result
simm <- 0      # loop index
num.sim <- 0   # number of simulations
while (num.sim < 100) {
  simm <- simm + 1
  print(simm)
  ### generate the simulated dataset
  set.seed(simm)
  n <- 700   # number of observations
  j <- 1:50  # number of basis
  
  # parameters when generate class label
  b <- matrix(ifelse(j <= 2, 1, (j-2)^(-3)), ncol=1)
  
  ## generate curves
  data <- list()
  for (i in 1:n) {
    # random sparsify
    num.obs <- sample(10:20, 1)
    t <- sort(runif(num.obs, 0, 10))
    
    # 201~700 are test set
    if (i > 200) {
      range.train <- range(unlist(data$Lt[1:200]))
      while( (max(t) > range.train[2]) | (min(t) < range.train[1]) ) {
        t <- sort(runif(num.obs, 0, 10))
      }
    }
    
    # eigenfunctions
    phi <- sapply(j, function(x){
      if (x %% 2 == 0) {
        sqrt(1/5)*sin(pi*t*x/5)
      } else {
        sqrt(1/5)*cos(pi*t*x/5)
      }
    })
    
    # generate PC scores
    xi <- sapply(j, function(x){ rnorm(1, 0, sqrt( x^(-1.5) )) })
    
    # parameters when generate class label
    beta.1 <- phi %*% b
    beta.2 <- matrix(sqrt(3/10)*(t/5-1), ncol=1)
    
    # measurement error
    eps <- rnorm(num.obs, 0, sqrt(0.1))
    
    # generate the curve
    X <- xi %*% t(phi) + eps
    
    # generate class label
    eps <- rnorm(1, 0, sqrt(0.1))   # model error
    # fx <- sin(pi*(X %*% beta.1)/4)  # model 1
    fx <- exp((X %*% beta.1)/2)-1   # model 2
    # fx <- sin(pi*(X %*% beta.1)/3) + exp((X %*% beta.2)/3) - 1   # model 3
    # fx <- atan(pi*(X %*% beta.1)) + exp((X %*% beta.2)/3) - 1    # model 4
    y <- factor(ifelse(fx + eps < 0, 0, 1), levels=c(0, 1))
    
    data$id[[i]] <- rep(i, num.obs)
    data$y[[i]] <- rep(y, num.obs)
    data$Lt[[i]] <- t
    data$Ly[[i]] <- X
  }
  
  # # plot the generated curves
  # sapply(data, unlist) %>% 
  #   data.frame() %>% 
  #   mutate(y = ifelse(unlist(data$y) == 0, "G1", "G2")) %>% 
  #   ggplot(aes(x=Lt, y=Ly, group=id, color=y)) +
  #   geom_line() +
  #   xlab("Time") +
  #   ylab("") +
  #   theme_bw() +
  #   theme(legend.title = element_blank())
  
  
  ### train, test split => train: 1~200, test: 201~700
  id <- 1:n
  ind.train <- 1:200
  ind.test <- 201:700
  
  # split to train, test set
  y <- sapply(data$y, unique)
  X.train <- list(Lid = lapply(data$id[ind.train], unique),
                  Lt = data$Lt[ind.train],
                  Ly = data$Ly[ind.train])
  X.test <- list(Lid = lapply(data$id[ind.test], unique),
                 Lt = data$Lt[ind.test],
                 Ly = data$Ly[ind.test])
  y.train <- y[ind.train]
  y.test <- y[ind.test]
  
  
  ### Single classifier
  # fit functional PCA
  fpca.fit <- tryCatch({ 
    FPCA(X.train$Ly, 
         X.train$Lt, 
         optns=list(dataType="Sparse",
                    methodXi="CE",
                    # methodSelectK="BIC",
                    FVEthreshold=0.99,
                    verbose=F)) 
  }, error = function(e) { 
    print(e)
    # print("error on FPCA!!")
    return(FALSE) 
  })
  
  # fpca.fit == FALSE인 경우
  if ( isFALSE(fpca.fit) ) {
    next
    # return(NA)
  } else {
    num.sim <- num.sim + 1
  }
  
  k <- fpca.fit$selectK   # optimal number of PCs
  fpc.score <- fpca.fit$xiEst   # FPC scores
  
  # train FPC scores
  train.fpc <- data.frame(y = y.train, 
                          x = fpc.score)
  colnames(train.fpc) <- c("y", paste("FPC", 1:k, sep=""))
  
  # tune the hyperparmeters
  tune.linear <- tune.svm(y ~ .,
                          data = train.fpc,
                          kernel = "linear",
                          cost = c(10^(-3:1), 2^(5:10)))
  tune.radial <- tune.svm(y ~ .,
                          data = train.fpc,
                          kernel = "radial",
                          cost = c(10^(-3:1), 2^(5:10)),
                          gamma = c(10^(-3:1), 2^(5:10)))
  
  # fit classifiers
  fit.logit <- glm(y~., train.fpc, family=binomial)
  fit.svm.linear <- svm(y ~ .,
                        data = train.fpc, 
                        kernel = "linear", 
                        cost = tune.linear$best.model$cost)
  fit.svm.radial <- svm(y ~ .,
                        data = train.fpc, 
                        kernel = "radial", 
                        cost = tune.radial$best.model$cost,
                        gamma = tune.radial$best.model$gamma)
  fit.lda <- lda(y~., train.fpc)
  fit.qda <- qda(y~., train.fpc)
  fit.nb <- naiveBayes(y~., train.fpc)
  
  
  # test FPC scores
  fpc.score.test <- predict(fpca.fit, X.test$Ly, X.test$Lt, K=k, xiMethod="CE")
  test.fpc <- data.frame(y = y.test, 
                         x = fpc.score.test)
  colnames(test.fpc) <- c("y", paste("FPC", 1:k, sep=""))
  
  # predict
  pred <- list()
  pred.logit <- predict(fit.logit, test.fpc, type="response")
  pred[[1]] <- factor(ifelse(pred.logit > 0.5, 1, 0), levels=c(0, 1))
  pred[[2]] <- predict(fit.svm.linear, test.fpc)
  pred[[3]] <- predict(fit.svm.radial, test.fpc)
  pred[[4]] <- predict(fit.lda, test.fpc)$class
  pred[[5]] <- predict(fit.qda, test.fpc)$class
  pred[[6]] <- predict(fit.nb, test.fpc, type="class")
  
  # save the error rate
  err.single <- sapply(pred, function(x){ mean(x != y.test) })
  
  
  ### Bootstrap aggregating
  B <- 100
  N <- length(X.train$Ly)
  boot.mat <- matrix(sample(1:N, N*B, replace=T), N, B)
  y.pred <- foreach(b=1:B, .packages=packages) %dopar% {
    boot.ind <- boot.mat[, b]   # bootstraped sample index
    
    # bootstrap sample의 time range를 벗어나는 test set sample 제거
    id.test <- NULL
    if (min(unlist(X.test$Lt)) < min(unlist(X.train$Lt[boot.ind]))) {
      # train set에서의 index(Not id)
      over.ind <- which(sapply(X.test$Lt, 
                               function(x){ sum( x < min(unlist(X.train$Lt[boot.ind])) ) }) != 0)
      id.test <- unlist(X.test$Lid)[-over.ind]
    }
    if (max(unlist(X.test$Lt)) > max(unlist(X.train$Lt[boot.ind]))) {
      # train set에서의 index(Not id)
      over.ind <- which(sapply(X.test$Lt, 
                               function(x){ sum( x > max(unlist(X.train$Lt[boot.ind])) ) }) != 0)
      if (is.numeric(id.test)) {
        id.test <- intersect(id.test,
                             unlist(X.test$Lid)[-over.ind])
      } else {
        id.test <- unlist(X.test$Lid)[-over.ind]
      }
    }
    if (!is.numeric(id.test)) {
      id.test <- unlist(X.test$Lid)
    }
    
    # bootstrap sample의 time range를 벗어나는 OOB sample 제거
    id.oob <- NULL
    if (min(unlist(X.train$Lt[-unique(boot.ind)])) < min(unlist(X.train$Lt[boot.ind]))) {
      # train set에서의 index(Not id)
      over.ind <- which(sapply(X.train$Lt[-unique(boot.ind)], 
                               function(x){ sum( x < min(unlist(X.train$Lt[boot.ind])) ) }) != 0)
      id.oob <- unlist(X.train$Lid)[-unique(boot.ind)][-over.ind]
    }
    if (max(unlist(X.train$Lt[-unique(boot.ind)])) > max(unlist(X.train$Lt[boot.ind]))) {
      # train set에서의 index(Not id)
      over.ind <- which(sapply(X.train$Lt[-unique(boot.ind)], 
                               function(x){ sum( x > max(unlist(X.train$Lt[boot.ind])) ) }) != 0)
      if (is.numeric(id.oob)) {
        id.oob <- intersect(id.oob,
                            unlist(X.train$Lid)[-unique(boot.ind)][-over.ind])
      } else {
        id.oob <- unlist(X.train$Lid)[-unique(boot.ind)][-over.ind]
      }
    }
    if (!is.numeric(id.oob)) {
      id.oob <- unlist(X.train$Lid)[-unique(boot.ind)]
    }
    
    # fit functional PCA
    fpca.fit <- tryCatch({ 
      FPCA(X.train$Ly[boot.ind], 
           X.train$Lt[boot.ind], 
           optns=list(dataType="Sparse",
                      methodXi="CE",
                      # methodSelectK="BIC",
                      FVEthreshold=0.99,
                      verbose=F)) 
    }, error = function(e) { 
      print(e)
      # print("error on FPCA!!")
      return(FALSE) 
    })
    
    # fpca.fit == FALSE인 경우
    if ( isFALSE(fpca.fit) ) {
      # next
      return(NA)
    }
    
    k <- fpca.fit$selectK   # optimal number of PCs
    fpc.score <- fpca.fit$xiEst
    
    # train FPC scores
    train.fpc <- data.frame(y = y.train[boot.ind], 
                            x = fpc.score)
    colnames(train.fpc) <- c("y", paste("FPC", 1:k, sep=""))
    
    # fit classifiers
    fit.logit <- glm(y~., train.fpc, family=binomial)
    fit.svm.linear <- svm(y ~ ., 
                          data = train.fpc,
                          kernel = "linear",
                          cost = tune.linear$best.model$cost)
    fit.svm.radial <- svm(y ~ ., 
                          data = train.fpc,
                          kernel = "radial",
                          cost = tune.radial$best.model$cost,
                          gamma = tune.radial$best.model$gamma)
    fit.lda <- lda(y~., train.fpc)
    fit.qda <- qda(y~., train.fpc)
    fit.nb <- naiveBayes(y~., train.fpc)
    
    # test FPC scores
    ind <- which(unlist(X.test$Lid) %in% id.test)
    fpc.score.test <- predict(fpca.fit, X.test$Ly[ind], X.test$Lt[ind], K=k, xiMethod="CE")
    test.fpc <- data.frame(y = y.test[ind], 
                           x = fpc.score.test)
    colnames(test.fpc) <- c("y", paste("FPC", 1:k, sep=""))
    
    # predict
    pred.logit <- predict(fit.logit, test.fpc, type="response")
    pred.logit <- factor(ifelse(pred.logit > 0.5, 1, 0), levels=c(0, 1))
    pred.svm.linear <- predict(fit.svm.linear, test.fpc)
    pred.svm.radial <- predict(fit.svm.radial, test.fpc)
    pred.lda <- predict(fit.lda, test.fpc)$class
    pred.qda <- predict(fit.qda, test.fpc)$class
    pred.nb <- predict(fit.nb, test.fpc, type="class")
    
    # OOB FPC scores
    oob.ind <- which(X.train$Lid %in% id.oob)
    fpc.score.oob <- predict(fpca.fit,
                             X.train$Ly[oob.ind],
                             X.train$Lt[oob.ind],
                             K=k,
                             xiMethod="CE")
    oob.fpc <- data.frame(y = y.train[oob.ind],
                          x = fpc.score.oob)
    colnames(oob.fpc) <- c("y", paste("FPC", 1:k, sep=""))
    
    # OOB error rate
    oob.error <- c(mean(factor(ifelse(predict(fit.logit, oob.fpc, type="response") > 0.5, 1, 0), levels=c(0, 1)) != oob.fpc$y),
                   mean(predict(fit.svm.linear, oob.fpc) != oob.fpc$y),
                   mean(predict(fit.svm.radial, oob.fpc) != oob.fpc$y),
                   mean(predict(fit.lda, oob.fpc)$class != oob.fpc$y),
                   mean(predict(fit.qda, oob.fpc)$class != oob.fpc$y),
                   mean(predict(fit.nb, oob.fpc, type='class') != oob.fpc$y) )
    
    # train error rate
    train.error <- c(mean(factor(ifelse(predict(fit.logit, train.fpc, type="response") > 0.5, 1, 0), levels=c(0, 1)) != train.fpc$y),
                     mean(predict(fit.svm.linear, train.fpc) != train.fpc$y),
                     mean(predict(fit.svm.radial, train.fpc) != train.fpc$y),
                     mean(predict(fit.lda, train.fpc)$class != train.fpc$y),
                     mean(predict(fit.qda, train.fpc)$class != train.fpc$y),
                     mean(predict(fit.nb, train.fpc, type='class') != train.fpc$y) )
    
    return( list(oob.error = oob.error,
                 train.error = train.error,
                 boot = data.frame(id = id.test,
                                   y = y.test[ind],
                                   # predicted value
                                   logit = pred.logit,
                                   svm.linear = pred.svm.linear,
                                   svm.radial = pred.svm.radial,
                                   lda = pred.lda,
                                   qda = pred.qda,
                                   nb = pred.nb,
                                   # predicted value * OOB accuracy
                                   logit.oob = as.numeric(pred.logit) * (1 - oob.error[1]),
                                   svm.linear.oob = as.numeric(pred.svm.linear) * (1 - oob.error[2]),
                                   svm.radial.oob = as.numeric(pred.svm.radial) * (1 - oob.error[3]),
                                   lda.oob = as.numeric(pred.lda) * (1 - oob.error[4]),
                                   qda.oob = as.numeric(pred.qda) * (1 - oob.error[5]),
                                   nb.oob = as.numeric(pred.nb) * (1 - oob.error[6]))) )
  }
  
  
  ## save the accuracy
  res <- as.data.frame(rbindlist(lapply(y.pred, function(x){ x$boot })))
  # majority voting
  pred <- res %>% 
    group_by(id, y) %>% 
    summarise(logit = Mode(logit),
              svm.linear = Mode(svm.linear),
              svm.radial = Mode(svm.radial),
              lda = Mode(lda),
              qda = Mode(qda),
              nb = Mode(nb))
  err.majority <- 1 - apply(pred[, 3:8], 2, function(x){ mean(x == pred$y) })
  
  # oob error
  oob.acc <- colSums( 1 - t(sapply(y.pred, function(x){ x$oob.error })) )
  pred <- res %>% 
    group_by(id, y) %>% 
    summarise(logit = factor(ifelse(sum(logit.oob)/oob.acc[1] > 1.5, 1, 0), 
                             levels=c(0, 1)),
              svm.linear = factor(ifelse(sum(svm.linear.oob)/oob.acc[2] > 1.5, 1, 0), 
                                  levels=c(0, 1)),
              svm.radial = factor(ifelse(sum(svm.radial.oob)/oob.acc[3] > 1.5, 1, 0), 
                                  levels=c(0, 1)),
              lda = factor(ifelse(sum(lda.oob)/oob.acc[4] > 1.5, 1, 0), 
                           levels=c(0, 1)),
              qda = factor(ifelse(sum(qda.oob)/oob.acc[5] > 1.5, 1, 0), 
                           levels=c(0, 1)),
              nb = factor(ifelse(sum(nb.oob)/oob.acc[6] > 1.5, 1, 0), 
                          levels=c(0, 1)))
  err.oob <- 1 - apply(pred[, 3:8], 2, function(x){ mean(x == pred$y) })
  
  # save result
  res <- as.data.frame(rbind(err.single,
                             err.majority,
                             err.oob))
  colnames(res) <- c("Logit","SVM(Linear)","SVM(Gaussian)","LDA","QDA","NaiveBayes")
  result[[simm]] <- res
  
  print(paste("# of sim =", num.sim))
  # res.sim[[simm]] <- y.pred
  
  # 10개마다 RData 저장
  if (num.sim %% 10 == 0) {
    save(result, file="RData/sim(sparse)_1.RData")
  }
}

save(result, file="RData/sim(sparse)_1.RData")


## 결과 정리
result <- result[!sapply(result, is.null)]

res <- sapply(1:3, function(i){
  paste(lapply(result, function(x){ x[i, ]*100 }) %>% 
          rbindlist %>% 
          colMeans %>% 
          round(1),
        " (",
        apply(lapply(result[!sapply(result, is.null)], 
                     function(x){ x[i, ]*100 }) %>% 
                rbindlist, 2, sd) %>% 
          round(2),
        ")",
        sep="")
}) %>% t()

library(xtable)
xtable(res)



