setwd("C:\\Users\\user\\Desktop\\KHS\\FDA-Lab\\Sparse FPCA\\Application\\bootstrap_aggregating")

library(tidyverse)
library(doParallel)   # parallel computing
library(fdapace)
library(e1071)
library(class)
library(MASS)


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

# load data
data <- read.csv("../spnbmd.csv", header=T)

ind <- data %>%
  group_by(idnum) %>%
  summarise(n=n(),
            max.age=max(age)) %>%
  filter(n >= 2 & max.age < 18) %>%
  dplyr::select(idnum)
data <- data[which(data$idnum %in% ind$idnum), ]

### train, test split
range(data$age)   # range check => train set에 없는 time point 포함시 에러
data[which(data$age == 17.9), "idnum"]   # 54 68 99 222
data[which(data$age == 8.8), "idnum"]    # 274  408
N <- length(unique(data$idnum))   # 160

# 위의 timepoint가 train set에 반드시 들어가도록 split
set.seed(100)
ind <- sample(unique(data$idnum), 100)
while((sum(c(54, 68, 99, 222) %in% ind) == 0) | (sum(c(274, 408) %in% ind) == 0)) {
  ind <- sample(unique(data$idnum), 100)
}
ind.train <- sample(ind, 60)
while((sum(c(54, 68, 99, 222) %in% ind.train) == 0) | (sum(c(274, 408) %in% ind.train) == 0)) {
  ind.train <- sample(ind, 60)
}
  
# response class 0, 1 coding
data$sex <- factor(ifelse(data$sex == "fem", 1, 0), levels=c(0,1))

# transform to FPCA input
train <- data[which(data$idnum %in% ind.train), ]
valid <- data[which(data$idnum %in% ind & !(data$idnum %in% ind.train)), ]
test <- data[-which(data$idnum %in% ind), ]

X.train <- MakeFPCAInputs(IDs = train$idnum,
                          tVec = train$age,
                          yVec = train$spnbmd)
X.test <- MakeFPCAInputs(IDs = test$idnum,
                         tVec = test$age,
                         yVec = test$spnbmd)
X.valid <- MakeFPCAInputs(IDs = valid$idnum,
                          tVec = valid$age,
                          yVec = valid$spnbmd)
y.train <- unique(train[, c("idnum", "sex")])[, "sex"]
y.test <- unique(test[, c("idnum", "sex")])[, "sex"]
y.valid <- unique(valid[, c("idnum", "sex")])[, "sex"]


# plot
par(mfrow=c(1,2))
plot(train$age[which(train$idnum == unique(train$idnum)[1])], 
     train$spnbmd[which(train$idnum == unique(train$idnum)[1])], 
     col=as.numeric(train$sex)[which(train$idnum == unique(train$idnum)[1])[1]],
     type="o", xlim=c(8, 27), ylim=c(0.5, 1.7),
     xlab="Age", ylab="Bone mineral density",
     main="Train set")
for (i in 2:length(unique(train$idnum))) {
  lines(train$age[which(train$idnum == unique(train$idnum)[i])], 
        train$spnbmd[which(train$idnum == unique(train$idnum)[i])], 
        col=as.numeric(train$sex)[which(train$idnum == unique(train$idnum)[i])[1]],
        type="o")
}
plot(test$age[which(test$idnum == unique(test$idnum)[1])], 
     test$spnbmd[which(test$idnum == unique(test$idnum)[1])], 
     col=as.numeric(train$sex)[which(test$idnum == unique(test$idnum)[1])[1]],
     type="o", xlim=c(8, 27), ylim=c(0.5, 1.7),
     xlab="Age", ylab="Bone mineral density",
     main="Test set")
for (i in 2:length(unique(test$idnum))) {
  lines(test$age[which(test$idnum == unique(test$idnum)[i])], 
        test$spnbmd[which(test$idnum == unique(test$idnum)[i])], 
        col=as.numeric(test$sex)[which(test$idnum == unique(test$idnum)[i])[1]],
        type="o")
}


### Single classifier
set.seed(100)
fpca.fit <- FPCA(X.train$Ly, 
                 X.train$Lt, 
                 optns=list(dataType="Sparse",
                            methodXi="CE",
                            # methodSelectK="BIC",
                            FVEthreshold=0.99,
                            verbose=T))
k <- fpca.fit$selectK   # optimal number of PCs
fpc.score <- fpca.fit$xiEst

# train FPC scores
train.fpc <- data.frame(y=y.train, 
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
                      cost=tune.linear$cost,
                      probability=T)
fit.svm.radial <- svm(y ~ ., train.fpc, 
                      kernel="radial", 
                      cost=tune.radial$cost,
                      gamma=tune.radial$gamma,
                      probability=T)
fit.svm.sigmoid <- svm(y ~ ., train.fpc, 
                       kernel="sigmoid", 
                       cost=tune.sigmoid$cost,
                       gamma=tune.sigmoid$gamma,
                       probability=T)
fit.knn <- tune.knn(x = train.fpc[, -1],
                    y = train.fpc[, 1],
                    k = 1:ceiling(length(y.train)/2))$best.model
fit.lda <- lda(y~., train.fpc)
fit.qda <- qda(y~., train.fpc)
fit.nb <- naiveBayes(y~., train.fpc)


# test FPC scores
fpc.score.test <- predict(fpca.fit, X.test$Ly, X.test$Lt, K=k, xiMethod="CE")
test.fpc <- data.frame(y=y.test, 
                       x=fpc.score.test)
colnames(test.fpc) <- c("y", paste("FPC", 1:k, sep=""))

# predict
pred <- list()
pred.logit <- predict(fit.logit, test.fpc, type="response")
pred[[1]] <- factor(ifelse(pred.logit > 0.5, 1, 0), levels=c(0, 1))
pred[[2]] <- predict(fit.svm.linear, test.fpc)
pred[[3]] <- predict(fit.svm.radial, test.fpc)
pred[[4]] <- predict(fit.svm.sigmoid, test.fpc)
pred[[5]] <- knn(train = cbind(y=fit.knn$cl,
                               fit.knn$train),
                 test = test.fpc, cl = fit.knn$cl, k = fit.knn$k)
pred[[6]] <- predict(fit.lda, test.fpc)$class
pred[[7]] <- predict(fit.qda, test.fpc)$class
pred[[8]] <- predict(fit.nb, test.fpc, type="class")

# save the accuracy
cname <- c("logit","svm.linear","svm.radial","svm.sigmoid","knn","lda","qda","naivebayes")
acc <- sapply(pred, function(x){ mean(x == y.test) })
acc
# 0.7333333 0.7666667 0.7000000 0.5666667 1.0000000 0.7833333 0.6666667 0.7500000

# valid FPC scores
fpc.score.valid <- predict(fpca.fit, X.valid$Ly, X.valid$Lt, K=k, xiMethod="CE")
valid.fpc <- data.frame(y=y.valid, 
                        x=fpc.score.valid)
colnames(valid.fpc) <- c("y", paste("FPC", 1:k, sep=""))
valid.err <- c(mean(factor(ifelse(predict(fit.logit, valid.fpc, type="response") > 0.5, 1, 0), levels=c(0, 1)) != y.valid),
               mean(predict(fit.svm.linear, valid.fpc) != y.valid),
               mean(predict(fit.svm.radial, valid.fpc) != y.valid),
               mean(predict(fit.svm.sigmoid, valid.fpc) != y.valid),
               mean(knn(train = cbind(y=fit.knn$cl,
                                      fit.knn$train),
                        test = valid.fpc, cl = fit.knn$cl, k = fit.knn$k) != y.valid),
               mean(predict(fit.lda, valid.fpc)$class != y.valid),
               mean(predict(fit.qda, valid.fpc)$class != y.valid),
               mean(predict(fit.nb, valid.fpc, type="class") != y.valid) )
valid.err



## Robust Bagging
system.time({   # running time check
  # Bootstrap aggregating
  B <- 200
  N <- length(X.train$Ly)
  set.seed(100)
  y.pred <- foreach(b=1:B, .packages=packages) %dopar% {
    # fit FPCA for bootstrapped data
    boot.ind <- sample(1:N, N, replace=T) 
    while( !( (min(unlist(X.train$Lt[boot.ind])) <= min(unlist(X.test$Lt))) &
              (max(unlist(X.train$Lt[boot.ind])) >= max(unlist(X.test$Lt))) ) |
           !( (min(unlist(X.train$Lt[boot.ind])) <= min(unlist(X.train$Lt))) &
              (max(unlist(X.train$Lt[boot.ind])) >= max(unlist(X.train$Lt))) ) ) {
      boot.ind <- sample(1:N, N, replace=T)
    }
    
    fpca.fit <- FPCA(X.train$Ly[boot.ind], 
                     X.train$Lt[boot.ind], 
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
    
    
    # # classifiers with bootstrapped FPC scores
    # tune.linear <- tune.svm(y ~ ., data=train.fpc,
    #                         kernel="linear",
    #                         cost=c(10^(-3:1), 2^(5:10)))$best.model
    # tune.radial <- tune.svm(y ~ ., data=train.fpc,
    #                         kernel="radial",
    #                         cost=c(10^(-3:1), 2^(5:10)),
    #                         gamma=c(10^(-3:1), 2^(5:10)) )$best.model
    # tune.sigmoid <- tune.svm(y ~ ., data=train.fpc,
    #                          kernel="sigmoid",
    #                          cost=c(10^(-3:1), 2^(5:10)),
    #                          gamma=c(10^(-3:1), 2^(5:10)) )$best.model
    
    fit.logit <- glm(y~., train.fpc, family=binomial)
    fit.svm.linear <- svm(y ~ ., train.fpc,
                          kernel="linear",
                          cost=tune.linear$cost,
                          probability=T)
    fit.svm.radial <- svm(y ~ ., train.fpc,
                          kernel="radial",
                          cost=tune.radial$cost,
                          gamma=tune.radial$gamma,
                          probability=T)
    fit.svm.sigmoid <- svm(y ~ ., train.fpc,
                           kernel="sigmoid",
                           cost=tune.sigmoid$cost,
                           gamma=tune.sigmoid$gamma,
                           probability=T)
    # fit.knn <- tune.knn(x = train.fpc[, -1],
    #                     y = train.fpc[, 1],
    #                     k = 1:ceiling(length(y.train)/2))$best.model
    fit.lda <- lda(y~., train.fpc)
    fit.qda <- qda(y~., train.fpc)
    fit.nb <- naiveBayes(y~., train.fpc)
    

    # test FPC scores
    fpc.score.test <- predict(fpca.fit, X.test$Ly, X.test$Lt, K=k, xiMethod="CE")
    test.fpc <- data.frame(y=y.test, 
                           x=fpc.score.test)
    colnames(test.fpc) <- c("y", paste("FPC", 1:k, sep=""))
    
    # predict
    pred.logit <- predict(fit.logit, test.fpc, type="response")
    pred.logit <- factor(ifelse(pred.logit > 0.5, 1, 0), levels=c(0, 1))
    pred.svm.linear <- predict(fit.svm.linear, test.fpc)
    pred.svm.radial <- predict(fit.svm.radial, test.fpc)
    pred.svm.sigmoid <- predict(fit.svm.sigmoid, test.fpc)
    # pred.knn <- knn(train = train.fpc,
    #                 test = test.fpc,
    #                 cl = train.fpc$y, k = fit.knn$k)
    pred.knn <- knn(train = train.fpc,
                    test = test.fpc, cl = train.fpc$y, k = 1)
    pred.lda <- predict(fit.lda, test.fpc)$class
    pred.qda <- predict(fit.qda, test.fpc)$class
    pred.nb <- predict(fit.nb, test.fpc, type="class")
    
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
    oob.error <- c(mean(factor(ifelse(predict(fit.logit, oob.fpc, type="response") > 0.5, 1, 0), levels=c(0, 1)) != oob.fpc$y),
                   mean(predict(fit.svm.linear, oob.fpc) != oob.fpc$y),
                   mean(predict(fit.svm.radial, oob.fpc) != oob.fpc$y),
                   mean(predict(fit.svm.sigmoid, oob.fpc) != oob.fpc$y),
                   # mean(knn(train = train.fpc, test = oob.fpc, cl = train.fpc$y, k = fit.knn$k) != oob.fpc$y),
                   mean(knn(train = train.fpc, test = oob.fpc, cl = train.fpc$y, k = 1) != oob.fpc$y),
                   mean(predict(fit.lda, oob.fpc)$class != oob.fpc$y),
                   mean(predict(fit.qda, oob.fpc)$class != oob.fpc$y),
                   mean(predict(fit.nb, oob.fpc, type='class') != oob.fpc$y) )
    
    # # train error rate
    # train.error <- c(mean(predict(fit.svm.linear, train.fpc) != train.fpc$y),
    #                  mean(predict(fit.svm.radial, train.fpc) != train.fpc$y),
    #                  mean(predict(fit.svm.sigmoid, train.fpc) != train.fpc$y),
    #                  mean(factor(ifelse(predict(fit.logit, train.fpc, type="response") > 0.5, 1, 0), levels=c(0, 1)) != train.fpc$y),
    #                  mean(knn(train = train.fpc, test = train.fpc, cl = train.fpc$y, k = fit.knn$k) != train.fpc$y),
    #                  mean(predict(fit.lda, train.fpc)$class != train.fpc$y),
    #                  mean(predict(fit.qda, train.fpc)$class != train.fpc$y),
    #                  mean(predict(fit.nb, train.fpc, type='class') != train.fpc$y) )
    # 
    # # 전체 train set으로 pred
    # fpc.score.train <- predict(fpca.fit, X.train$Ly, X.train$Lt, K=k, xiMethod="CE")
    # train.fpc2 <- data.frame(y=y.train, 
    #                          x=fpc.score.train)
    # colnames(train.fpc2) <- c("y", paste("FPC", 1:k, sep=""))
    
    return( list(logit = pred.logit,
                 svm.linear = pred.svm.linear,
                 svm.radial = pred.svm.radial,
                 svm.sigmoid = pred.svm.sigmoid,
                 knn = pred.knn,
                 lda = pred.lda,
                 qda = pred.qda,
                 nb = pred.nb,
                 oob.error = oob.error,
                 # train.error = train.error,
                 # fitted.logit = factor(ifelse(predict(fit.logit, train.fpc2, type="response") > 0.5, 1, 0), levels=c(0, 1)),
                 # fitted.svm.linear = predict(fit.svm.linear, train.fpc2),
                 # fitted.svm.radial = predict(fit.svm.radial, train.fpc2),
                 # fitted.svm.sigmoid = predict(fit.svm.sigmoid, train.fpc2),
                 # fitted.knn = knn(train = train.fpc, test = train.fpc2, cl = train.fpc$y, k = fit.knn$k),
                 # fitted.lda = predict(fit.lda, train.fpc2)$class,
                 # fitted.qda = predict(fit.qda, train.fpc2)$class,
                 # fitted.nb = predict(fit.nb, train.fpc2, type='class'),
                 K = k) )
  }
})   # running time check



# save the accuracy
library(glmnet)
cname <- c("logit","svm.linear","svm.radial","svm.sigmoid","knn","lda","qda","naivebayes")
acc <- numeric(8)
for (j in 1:8) {
  # robust bagging
  pred <- data.frame(id = 1:B,
                     error = sapply(y.pred, function(x){ cbind( x$oob.error[j] ) })) %>%
    filter(error < valid.err[j]) %>%
    dplyr::select(id) %>%
    unlist()
  pred <- apply(sapply(y.pred, function(x){ cbind( x[[j]] ) })[, pred],
                1,
                median)
  pred <- factor(ifelse(pred > 1, 1, 0), levels=c(0,1))
  
  # # majority voting
  # pred <- apply(sapply(y.pred, function(x){ cbind( x[[j]] ) }),
  #               1,
  #               majority_vote)
  # pred <- factor(pred-1, levels=c(0,1))
  
  # # oob error
  # pred <- apply(sapply(y.pred, function(x){ cbind( x[[j]] ) }),
  #               1,
  #               function(x){
  #                 oob.error <- 1 - sapply(y.pred, function(x){ cbind( x$oob.error[j] ) })
  #                 # oob.error <- (oob.error - min(oob.error)) / (max(oob.error) - min(oob.error))
  #                 oob.error <- oob.error / sum(oob.error)
  #                 res <- factor(ifelse(x %*% oob.error > 1.5, 1, 0), levels=c(0, 1))
  #                 return(res)
  #               })
  
  acc[j] <- mean(pred == y.test)
}
acc

# Tuned
# 0.7333333 0.7666667 0.7000000 0.5666667 1.0000000 0.7833333 0.6666667 0.7500000  Single
# 0.7333333 0.7333333        NA 0.7500000        NA 0.7333333 0.7333333 0.7166667  robust bagging
# 0.8000000 0.7666667 0.6500000 0.7666667 1.0000000 0.7833333 0.7333333 0.7500000  majority vote
# 0.8000000 0.7666667 0.6500000 0.7666667 1.0000000 0.7833333 0.7333333 0.7500000  oob error weight

# Not tuned
# 0.7333333 0.7666667 0.7000000 0.5666667 1.0000000 0.7833333 0.6666667 0.7500000  Single
# 0.7500000 0.7500000 0.7333333 0.7833333        NA 0.7500000 0.7833333 0.7333333  robust bagging
# 0.8166667 0.7666667 0.7333333 0.7666667 1.0000000 0.7833333 0.7333333 0.7666667  majority vote
# 0.8166667 0.7500000 0.7333333 0.8000000 1.0000000 0.8166667 0.7500000 0.7500000  oob error weight


## iteration마다의 error rate 그래프
for (i in 2:B) {
  acc <- numeric(8)
  for (j in 1:8) {
    # majority voting
    pred <- apply(sapply(y.pred[1:i], function(x){ cbind( x[[j]] ) }),
                  1,
                  majority_vote)
    pred <- factor(pred-1, levels=c(0,1))
    
    acc[j] <- mean(pred == y.test)
  }
  if (i == 2) {
    res <- acc
  } else {
    res <- rbind(res, acc)
  }
}

single <- c(0.73, 0.77, 0.70, 0.57, 1.00, 0.78, 0.67, 0.75)
par(mfrow=c(2,4))
for (j in 1:8) {
  plot(1-res[,j], type="l", col=j+1, main=cname[j], ylim=c(0, 0.5), ylab="classification error rate")
  abline(h=1-single[j], lwd=2)
}
