library(ddalpha)

names(dataf.growth())

dataf <- dataf.growth()

Ly <- lapply(dataf$dataf, function(x){ x$vals })
time <- lapply(dataf$dataf, function(x){ x$args })[[1]]

X <- t(sapply(Ly, cbind))

# train, test split
set.seed(100)
y <- factor(ifelse(unlist(dataf$labels) == "boy", 1, 0), levels=c(0, 1))
n <- length(y)
ind <- sample(1:n, 60)   # train set index

X.train <- as.matrix( X[ind, ] )
X.test <- as.matrix( X[-ind, ] )
y.train <- y[ind]
y.test <- y[-ind]

# sparsify train, test set
s <- 2:6
train.data <- Sparsify(X.train, time, sparsity = s)
test.data <- Sparsify(X.test, time, sparsity = s)


### Single classifier
fpca.fit <- FPCA(train.data$Ly, 
                 train.data$Lt, 
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
fpc.score.test <- predict(fpca.fit, test.data$Ly, test.data$Lt, K=k, xiMethod="CE")
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

# 0.7575758 0.7272727 0.7272727 0.7272727 0.6060606 0.7272727 0.6969697 0.7272727  Single
# 0.8181818 0.7878788 0.7575758 0.8181818 0.5757576 0.8181818 0.7575758 0.7575758  majority



## Bagging to curves
system.time({   # running time check
# Bootstrap aggregating
B <- 200
N <- nrow(X.train)
set.seed(100)
y.pred <- foreach(b=1:B, .packages=packages) %dopar% {
  # fit FPCA for bootstrapped data
  boot.ind <- sample(1:N, N, replace=T)
  while( !(min(unlist(train.data$Lt[boot.ind])) <= min(unlist(test.data$Lt))) &
         !(max(unlist(train.data$Lt[boot.ind])) >= max(unlist(test.data$Lt))) ) {
    boot.ind <- sample(1:N, N, replace=T)
  }
  
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
  pred.knn <- knn(train = cbind(y=fit.knn$cl,
                                fit.knn$train),
                  test = test.fpc, cl = fit.knn$cl, k = fit.knn$k)
  pred.lda <- predict(fit.lda, test.fpc)$class
  pred.qda <- predict(fit.qda, test.fpc)$class
  pred.nb <- predict(fit.nb, test.fpc, type="class")
  
  
  # OOB FPC scores
  fpc.score.oob <- predict(fpca.fit, 
                           train.data$Ly[-unique(boot.ind)], 
                           train.data$Lt[-unique(boot.ind)], 
                           K=k, 
                           xiMethod="CE")
  oob.fpc <- data.frame(y=y.train[-unique(boot.ind)], 
                        x=fpc.score.oob)
  colnames(oob.fpc) <- c("y", paste("FPC", 1:k, sep=""))
  
  # OOB error rate
  oob.error <- c(mean(predict(fit.svm.linear, oob.fpc) != oob.fpc$y),
                 mean(predict(fit.svm.radial, oob.fpc) != oob.fpc$y),
                 mean(predict(fit.svm.sigmoid, oob.fpc) != oob.fpc$y),
                 mean(factor(ifelse(predict(fit.logit, oob.fpc, type="response") > 0.5, 1, 0), levels=c(0, 1)) != oob.fpc$y),
                 mean(knn(train = train.fpc, test = oob.fpc, cl = train.fpc$y, k = fit.knn$k) != oob.fpc$y),
                 mean(predict(fit.lda, oob.fpc)$class != oob.fpc$y),
                 mean(predict(fit.qda, oob.fpc)$class != oob.fpc$y),
                 mean(predict(fit.nb, oob.fpc, type='class') != oob.fpc$y) )
  
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
               oob.error = oob.error,
               train.error = train.error,
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

save(y.pred, file="sim_2.RData")

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
  # A <- ifelse(sapply(y.pred, function(x){ as.numeric( x[[j+10]] ) }) == 1, -1, 1)   # transform to -1 and 1
  # w[[j]] <- solve(t(A) %*% A) %*% t(A) %*% matrix(ifelse(as.numeric(y.train) == 1, -1, 1), ncol=1)
  # A <- ifelse(sapply(y.pred, function(x){ as.numeric( x[[j]] ) }) == 1, -1, 1)   # transform to -1 and 1
  # pred <- factor(ifelse(A %*% matrix(w[[j]], ncol=1) < 0, 0, 1), levels=c(0, 1))
  
  # # logistic - numeric
  # A <- ifelse(sapply(y.pred, function(x){ as.numeric( x[[j+10]] ) }) == 1, -1, 1)   # transform to -1 and 1
  # data.w <- data.frame(y=factor(ifelse(as.numeric(y.train) == 1, -1, 1), levels=c(-1,1)),
  #                      x=A)
  # fit.w <- glm(y~., data.w, family=binomial)
  # w[[j]] <- fit.w
  # A <- ifelse(sapply(y.pred, function(x){ as.numeric( x[[j]] ) }) == 1, -1, 1)   # transform to -1 and 1
  # data.w <- data.frame(x=A)
  # pred <- predict(w[[j]], data.w, type="response")
  # pred <- factor(ifelse(pred > 0.5, 1, 0), levels=c(0,1))
  
  # # logistic - factor
  # A <- ifelse(sapply(y.pred, function(x){ as.numeric( x[[j+10]] ) }) == 1, -1, 1)   # transform to -1 and 1
  # A <- as.data.frame(A) %>% mutate_if(is.numeric, function(x){ factor(x, levels=c(-1,1)) })
  # data.w <- data.frame(y=factor(ifelse(as.numeric(y.train) == 1, -1, 1), levels=c(-1,1)),
  #                      x=A)
  # fit.w <- glm(y~., data.w, family=binomial)
  # w[[j]] <- fit.w
  # A <- ifelse(sapply(y.pred, function(x){ as.numeric( x[[j]] ) }) == 1, -1, 1)   # transform to -1 and 1
  # A <- as.data.frame(A) %>% mutate_if(is.numeric, function(x){ factor(x, levels=c(-1,1)) })
  # data.w <- data.frame(x=A)
  # pred <- predict(w[[j]], data.w, type="response")
  # pred <- factor(ifelse(pred > 0.5, 1, 0), levels=c(0,1))
  
  ### B > n : ridge
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
# 0.8181818 0.7878788 0.7575758 0.8181818 0.5757576 0.8181818 0.7575758 0.7575758  majority vote(30)
# 0.7878788 0.7878788 0.7878788 0.7878788 0.6363636 0.8181818 0.6969697 0.7575758  majority vote(100)
# 0.7878788 0.7575758 0.7575758 0.7575758 0.6363636 0.7878788 0.7878788 0.7878788  majority(200)

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

single <- c(0.7575758, 0.7272727, 0.7272727, 0.7272727, 0.6060606, 0.7272727, 0.6969697, 0.7272727)
par(mfrow=c(2,4))
# plot(res[, 1], type="o", ylim=c(0.55, 0.85))
for (j in 1:8) {
  plot(1-res[,j], type="l", col=j+1, main=cname[j], ylim=c(0.1, 0.5), ylab="classification error rate")
  abline(h=1-single[j], lwd=2)
}



