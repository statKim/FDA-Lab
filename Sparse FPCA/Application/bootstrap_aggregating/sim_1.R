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

# ind <- data %>% 
#   group_by(idnum) %>% 
#   summarise(n=n(),
#             max.age=max(age)) %>% 
#   filter(n >= 2 & max.age < 18) %>% 
#   dplyr::select(idnum)
# data <- data[which(data$idnum %in% ind$idnum), ]
# 
# ### train, test split
# range(data$age)   # range check => train set에 없는 time point 포함시 에러
# data[which(data$age == 17.9), "idnum"]   # 54 68 99 222
# data[which(data$age == 8.8), "idnum"]    # 274  408
# N <- length(unique(data$idnum))   # 160
# 
# # 위의 timepoint가 train set에 반드시 들어가도록 split
# set.seed(100)
# ind <- sample(unique(data$idnum), 100)
# while((sum(c(54, 68, 99, 222) %in% ind) == 0) | (sum(c(274, 408) %in% ind) == 0)) {
#   ind <- sample(unique(data$idnum), 100)
# }

ind <- data %>%
  group_by(idnum) %>%
  summarise(n=n()) %>%
  filter(n >= 2) %>%
  dplyr::select(idnum)
data <- data[which(data$idnum %in% ind$idnum), ]

### train, test split
range(data$age)   # range check => train set에 없는 time point 포함시 에러
data[which(data$age == 26.2), "idnum"]   # 29
data[which(data$age == 8.8), "idnum"]    # 274  408
N <- length(unique(data$idnum))   # 280

# 위의 timepoint가 train set에 반드시 들어가도록 split
set.seed(100)
ind <- sample(unique(data$idnum), 180)
while((sum(c(29) %in% ind) == 0) | (sum(c(274, 408) %in% ind) == 0)) {
  ind <- sample(unique(data$idnum), 300)
}

# response class 0, 1 coding
data$sex <- factor(ifelse(data$sex == "fem", 1, 0), levels=c(0,1))

# transform to FPCA input
train <- data[which(data$idnum %in% ind), ]
test <- data[-which(data$idnum %in% ind), ]

X.train <- MakeFPCAInputs(IDs = train$idnum,
                          tVec = train$age,
                          yVec = train$spnbmd)
X.test <- MakeFPCAInputs(IDs = test$idnum,
                         tVec = test$age,
                         yVec = test$spnbmd)
y.train <- unique(train[, c("idnum", "sex")])[, "sex"]
y.test <- unique(test[, c("idnum", "sex")])[, "sex"]


### Single classifier
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
# 0.7666667 0.7166667 0.7333333 0.6500000 1.0000000 0.7666667 0.7500000 0.7666667  (Sub)
# 0.66 0.65 0.66 0.65 1.00 0.67 0.62 0.59  (Full - n >= 2)



## Bagging
system.time({   # running time check
# Bootstrap aggregating
B <- 200
N <- length(X.train$Ly)
set.seed(100)
y.pred <- foreach(b=1:B, .packages=packages) %dopar% {
  # fit FPCA for bootstrapped data
  boot.ind <- sample(1:N, N, replace=T) 
  while( !( (min(unlist(X.train$Lt[boot.ind])) <= min(unlist(X.test$Lt))) &
            (max(unlist(X.train$Lt[boot.ind])) >= max(unlist(X.test$Lt))) ) ) {
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
  pred.logit <- predict(fit.logit, test.fpc, type="response")
  pred.logit <- factor(ifelse(pred.logit > 0.5, 1, 0), levels=c(0, 1))
  pred.svm.linear <- predict(fit.svm.linear, test.fpc)
  pred.svm.radial <- predict(fit.svm.radial, test.fpc)
  pred.svm.sigmoid <- predict(fit.svm.sigmoid, test.fpc)
  pred.knn <- knn(train = cbind(y=fit.knn$cl,
                                fit.knn$train),
                  test = test.fpc, cl = fit.knn$cl, k = fit.knn$k)
  pred.lda <- predict(fit.lda, test.fpc)$class
  pred.qda <- predict(fit.qda, test.fpc)$class
  pred.nb <- predict(fit.nb, test.fpc, type="class")
  
  # # OOB FPC scores
  # fpc.score.oob <- predict(fpca.fit, 
  #                          X.train$Ly[-unique(boot.ind)], 
  #                          X.train$Lt[-unique(boot.ind)], 
  #                          K=k, 
  #                          xiMethod="CE")
  # oob.fpc <- data.frame(y=y.train[-unique(boot.ind)], 
  #                       x=fpc.score.oob)
  # colnames(oob.fpc) <- c("y", paste("FPC", 1:k, sep=""))
  # 
  # # OOB error rate
  # oob.error <- c(mean(predict(fit.svm.linear, oob.fpc) != oob.fpc$y),
  #                mean(predict(fit.svm.radial, oob.fpc) != oob.fpc$y),
  #                mean(predict(fit.svm.sigmoid, oob.fpc) != oob.fpc$y),
  #                mean(factor(ifelse(predict(fit.logit, oob.fpc, type="response") > 0.5, 1, 0), levels=c(0, 1)) != oob.fpc$y),
  #                mean(knn(train = train.fpc, test = oob.fpc, cl = train.fpc$y, k = fit.knn$k) != oob.fpc$y),
  #                mean(predict(fit.lda, oob.fpc)$class != oob.fpc$y),
  #                mean(predict(fit.qda, oob.fpc)$class != oob.fpc$y),
  #                mean(predict(fit.nb, oob.fpc, type='class') != oob.fpc$y) )
  
  # train error rate
  train.error <- c(mean(predict(fit.svm.linear, train.fpc) != train.fpc$y),
                   mean(predict(fit.svm.radial, train.fpc) != train.fpc$y),
                   mean(predict(fit.svm.sigmoid, train.fpc) != train.fpc$y),
                   mean(factor(ifelse(predict(fit.logit, train.fpc, type="response") > 0.5, 1, 0), levels=c(0, 1)) != train.fpc$y),
                   mean(knn(train = train.fpc, test = train.fpc, cl = train.fpc$y, k = fit.knn$k) != train.fpc$y),
                   mean(predict(fit.lda, train.fpc)$class != train.fpc$y),
                   mean(predict(fit.qda, train.fpc)$class != train.fpc$y),
                   mean(predict(fit.nb, train.fpc, type='class') != train.fpc$y) )
  
  return( list(logit = pred.logit,
               svm.linear = pred.svm.linear,
               svm.radial = pred.svm.radial,
               svm.sigmoid = pred.svm.sigmoid,
               knn = pred.knn,
               lda = pred.lda,
               qda = pred.qda,
               nb = pred.nb,
               oob.error = 0,
               train.error = train.error,
               fitted.logit = predict(fit.logit, train.fpc, type="response"),
               fitted.logit = factor(ifelse(predict(fit.logit, train.fpc, type="response") > 0.5, 1, 0), levels=c(0, 1)),
               fitted.svm.linear = predict(fit.svm.linear, train.fpc),
               fitted.svm.radial = predict(fit.svm.radial, train.fpc),
               fitted.svm.sigmoid = predict(fit.svm.sigmoid, train.fpc),
               fitted.knn = knn(train = train.fpc, test = train.fpc, cl = train.fpc$y, k = fit.knn$k),
               fitted.lda = predict(fit.lda, train.fpc)$class,
               fitted.qda = predict(fit.qda, train.fpc)$class,
               fitted.nb = predict(fit.nb, train.fpc, type='class')) )
}
})   # running time check

save(y.pred, file="sim_1.RData")

# save the accuracy
library(glmnet)
cname <- c("logit","svm.linear","svm.radial","svm.sigmoid","knn","lda","qda","naivebayes")
w <- list()
acc <- numeric(8)
for (j in 1:8) {
  # majority voting
  pred <- apply(sapply(y.pred, function(x){ cbind( x[[j]] ) }),
                1,
                majority_vote)
  pred <- factor(pred-1, levels=c(0,1))
  
  ### B < n
  # # regression
  # A <- ifelse(sapply(y.pred[1:100], function(x){ as.numeric( x[[j+10]] ) }) == 1, -1, 1)   # transform to -1 and 1
  # w[[j]] <- ginv(t(A) %*% A) %*% t(A) %*% matrix(ifelse(as.numeric(y.train) == 1, -1, 1), ncol=1)
  # A <- ifelse(sapply(y.pred[1:100], function(x){ as.numeric( x[[j]] ) }) == 1, -1, 1)   # transform to -1 and 1
  # pred <- factor(ifelse(A %*% matrix(w[[j]], ncol=1) < 0, 0, 1), levels=c(0, 1))
  
  # # logistic - numeric
  # A <- ifelse(sapply(y.pred[1:100], function(x){ as.numeric( x[[j+10]] ) }) == 1, -1, 1)   # transform to -1 and 1
  # data.w <- data.frame(y=factor(ifelse(as.numeric(y.train) == 1, -1, 1), levels=c(-1,1)),
  #                      x=A)
  # fit.w <- glm(y~., data.w, family=binomial)
  # w[[j]] <- fit.w
  # A <- ifelse(sapply(y.pred[1:100], function(x){ as.numeric( x[[j]] ) }) == 1, -1, 1)   # transform to -1 and 1
  # data.w <- data.frame(x=A)
  # pred <- predict(w[[j]], data.w, type="response")
  # pred <- factor(ifelse(pred > 0.5, 1, 0), levels=c(0,1))
  
  # # logistic - factor
  # A <- ifelse(sapply(y.pred[1:100], function(x){ as.numeric( x[[j+10]] ) }) == 1, -1, 1)   # transform to -1 and 1
  # A <- as.data.frame(A) %>% mutate_if(is.numeric, function(x){ factor(x, levels=c(-1,1)) })
  # data.w <- data.frame(y=factor(ifelse(as.numeric(y.train) == 1, -1, 1), levels=c(-1,1)),
  #                      x=A)
  # fit.w <- glm(y~., data.w, family=binomial)
  # w[[j]] <- fit.w
  # A <- ifelse(sapply(y.pred[1:100], function(x){ as.numeric( x[[j]] ) }) == 1, -1, 1)   # transform to -1 and 1
  # A <- as.data.frame(A) %>% mutate_if(is.numeric, function(x){ factor(x, levels=c(-1,1)) })
  # data.w <- data.frame(x=A)
  # pred <- predict(w[[j]], data.w, type="response")
  # pred <- factor(ifelse(pred > 0.5, 1, 0), levels=c(0,1))
  
  ## B > n : ridge
  # # regression
  # A <- ifelse(sapply(y.pred, function(x){ as.numeric( x[[j+10]] ) }) == 1, -1, 1)   # transform to -1 and 1
  # y.w <- ifelse(as.numeric(y.train) == 1, -1, 1)
  # cv.ridge <- cv.glmnet(A, y.w, alpha = 0)
  # fit.w <- glmnet(A, y.w, alpha = 0, lambda = cv.ridge$lambda.min)
  # w[[j]] <- fit.w
  # A <- ifelse(sapply(y.pred, function(x){ as.numeric( x[[j]] ) }) == 1, -1, 1)   # transform to -1 and 1
  # pred <- predict(w[[j]], A)
  # pred <- factor(ifelse(pred > 0, 1, 0), levels=c(0,1))
  
  # # logistic - numeric
  # A <- ifelse(sapply(y.pred, function(x){ as.numeric( x[[j+10]] ) }) == 1, -1, 1)   # transform to -1 and 1
  # y.w <- as.numeric(y.train) - 1
  # cv.ridge <- cv.glmnet(A, y.w, alpha = 0)
  # fit.w <- glmnet(A, y.w, alpha = 0, lambda = cv.ridge$lambda.min, family="binomial")
  # w[[j]] <- fit.w
  # A <- ifelse(sapply(y.pred, function(x){ as.numeric( x[[j]] ) }) == 1, -1, 1)   # transform to -1 and 1
  # pred <- predict(w[[j]], A, type="response")
  # pred <- factor(ifelse(pred > 0.5, 1, 0), levels=c(0,1))
  
  # # logistic - factor
  # A <- sapply(y.pred, function(x){ as.numeric( x[[j+10]] ) }) - 1
  # y.w <- as.numeric(y.train) - 1
  # cv.ridge <- cv.glmnet(A, y.w, alpha = 0)
  # fit.w <- glmnet(A, y.w, alpha = 0, lambda = cv.ridge$lambda.min)
  # w[[j]] <- fit.w
  # A <- sapply(y.pred, function(x){ as.numeric( x[[j]] ) }) - 1
  # pred <- predict(w[[j]], A, type="response")
  # pred <- factor(ifelse(pred > 0.5, 1, 0), levels=c(0,1))
  
  
  acc[j] <- mean(pred == y.test)
}
acc
# 0.7666667 0.7166667 0.7333333 0.6500000 1.0000000 0.7666667 0.7500000 0.7666667  Single
# 0.7500000 0.6833333 0.5666667 0.7166667 1.0000000 0.7500000 0.7166667 0.7166667  majority - N
# 0.7166667 0.7000000 0.5333333 0.5833333 0.0000000 0.2333333 0.4166667 0.3833333  reg
# 0.5166667 0.7500000 0.5666667 0.6833333 1.0000000 0.4333333 0.7333333 0.6666667  logit

# 0.66 0.65 0.66 0.65 1.00 0.67 0.62 0.59  Single
# 0.65 0.63 0.51 0.64 1.00 0.64 0.63 0.61  majority
# 0.65 0.70 0.49 0.67 1.00 0.71 0.57 0.58  reg
# 0.49 0.65 0.54 0.65 0.00 0.65 0.59 0.53  logit

## iteration마다의 error rate 그래프
for (i in 2:B) {
  acc <- numeric(8)
  for (j in 1:8) {
    # majority voting
    pred <- apply(sapply(y.pred[1:i], function(x){ cbind( x[[j]] ) }),
                  1,
                  majority_vote)
    pred <- factor(pred-1, levels=c(0,1))
    
    # # regression
    # A <- ifelse(sapply(y.pred[1:i], function(x){ as.numeric( x[[j+10]] ) }) == 1, -1, 1)   # transform to -1 and 1
    # w[[j]] <- ginv(t(A) %*% A) %*% t(A) %*% matrix(ifelse(as.numeric(y.train) == 1, -1, 1), ncol=1)
    # A <- ifelse(sapply(y.pred[1:i], function(x){ as.numeric( x[[j]] ) }) == 1, -1, 1)   # transform to -1 and 1
    # pred <- factor(ifelse(A %*% matrix(w[[j]], ncol=1) < 0, 0, 1), levels=c(0, 1))
    
    acc[j] <- mean(pred == y.test)
  }
  if (i == 2) {
    res <- acc
  } else {
    res <- rbind(res, acc)
  }
}

single <- c(0.7666667, 0.7166667, 0.7333333, 0.6500000, 1.0000000, 0.7666667, 0.7500000, 0.7666667)
par(mfrow=c(2,4))
for (j in 1:8) {
  plot(1-res[,j], type="l", col=j+1, main=cname[j], ylim=c(0, 0.5), ylab="classification error rate")
  abline(h=1-single[j], lwd=2)
}
