
s <- 6:12
## Bagging to curves
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


system.time({   # running time check
# Bootstrap aggregating
B <- 300
N <- nrow(X.train)
set.seed(100)
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
  
  # # predict
  # pred.logit <- predict(fit.logit, test.fpc, type="response")
  # pred.logit <- factor(ifelse(pred.logit > 0.5, 1, 0), levels=c(0, 1))
  # pred.svm.linear <- predict(fit.svm.linear, test.fpc)
  # pred.svm.radial <- predict(fit.svm.radial, test.fpc)
  # pred.svm.sigmoid <- predict(fit.svm.sigmoid, test.fpc)
  # # pred.svm.poly <- predict(fit.svm.poly, test.fpc)
  # pred.knn <- knn(train = cbind(y=fit.knn$cl,
  #                               fit.knn$train),
  #                 test = test.fpc, cl = fit.knn$cl, k = fit.knn$k)
  # pred.lda <- predict(fit.lda, test.fpc)$class
  # pred.qda <- predict(fit.qda, test.fpc)$class
  # pred.nb <- predict(fit.nb, test.fpc, type="class")
  
  # predict
  pred.logit <- predict(fit.logit, test.fpc, type="response")
  pred.svm.linear <- attr(predict(fit.svm.linear, test.fpc, probability = T), "prob")[, 1]
  pred.svm.radial <- attr(predict(fit.svm.radial, test.fpc, probability = T), "prob")[, 1]
  pred.svm.sigmoid <- attr(predict(fit.svm.sigmoid, test.fpc, probability = T), "prob")[, 1]
  # pred.svm.poly <- predict(fit.svm.poly, test.fpc)
  pred.knn <- knn(train = cbind(y=fit.knn$cl,
                                fit.knn$train),
                  test = test.fpc, cl = fit.knn$cl, k = fit.knn$k,
                  prob=T)
  pred.knn <- ifelse(pred.knn == 1, attr(pred.knn, "prob"), 1-attr(pred.knn, "prob"))
  pred.lda <- predict(fit.lda, test.fpc)$posterior[, 2]
  pred.qda <- predict(fit.qda, test.fpc)$posterior[, 2]
  pred.nb <- predict(fit.nb, test.fpc, type="raw")[, 2]
  
  
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
  
  
  
  fitted.knn <- knn(train = cbind(y=fit.knn$cl,
                                fit.knn$train),
                  test = test.fpc, cl = fit.knn$cl, k = fit.knn$k,
                  prob=T)
  
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
               fitted.logit = predict(fit.logit, train.fpc, type="response"),
               fitted.svm.linear = attr(predict(fit.svm.linear, train.fpc, probability = T), "prob")[, 1],
               fitted.svm.radial = attr(predict(fit.svm.radial, train.fpc, probability = T), "prob")[, 1],
               fitted.svm.sigmoid = attr(predict(fit.svm.sigmoid, train.fpc, probability = T), "prob")[, 1],
               fitted.knn = ifelse(fitted.knn == 1, attr(fitted.knn, "prob"), 1-attr(fitted.knn, "prob")),
               fitted.lda = predict(fit.lda, train.fpc)$posterior[, 2],
               fitted.qda = predict(fit.qda, train.fpc)$posterior[, 2],
               fitted.nb = predict(fit.nb, train.fpc, type="raw")[, 2],
               # fitted.logit = factor(ifelse(predict(fit.logit, train.fpc, type="response") > 0.5, 1, 0), levels=c(0, 1)),
               # fitted.svm.linear = predict(fit.svm.linear, train.fpc),
               # fitted.svm.radial = predict(fit.svm.radial, train.fpc),
               # fitted.svm.sigmoid = predict(fit.svm.sigmoid, train.fpc),
               # fitted.knn = knn(train = train.fpc, test = train.fpc, cl = train.fpc$y, k = fit.knn$k),
               # fitted.lda = predict(fit.lda, train.fpc)$class,
               # fitted.qda = predict(fit.qda, train.fpc)$class,
               # fitted.nb = predict(fit.nb, train.fpc, type='class'),
               prop = sum(as.numeric(y.train[boot.ind])-1) / N) )
}
})   # running time check


# save the accuracy
library(glmnet)
cname <- c("logit","svm.linear","svm.radial","svm.sigmoid","knn","lda","qda","naivebayes")
w <- list()
acc <- numeric(8)
for (j in 1:8) {
  # # majority voting
  # pred <- apply(sapply(y.pred, function(x){ cbind( x[[j]] ) }),
  #               1,
  #               majority_vote)
  # pred <- factor(pred-1, levels=c(0,1))
  
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

# 0.81 0.79 0.78 0.79 0.86 0.78 0.81 0.74  Single 
# 0.79 0.77 0.74 0.77 0.81 0.78 0.84 0.74  reg(100)
# 0.81 0.80 0.75 0.81 0.88 0.81 0.83 0.75  logit-numeric(100)
# 0.81 0.80 0.75 0.81 0.88 0.81 0.83 0.75  logit-factor(100)
# 0.78 0.78 0.74 0.76 0.82 0.79 0.84 0.74  ridge-reg(300)
# 0.78 0.78 0.77 0.75 0.82 0.78 0.84 0.74  ridge-logit-numeric(300)
# 0.78 0.78 0.74 0.76 0.82 0.78 0.83 0.74  ridge-logit-factor(300)


y.pred2 <- y.pred



for (j in 2:5) {
acc <- numeric(300)
for (i in 2:300) {
  y.pred <- y.pred2[1:i]
  
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
  # w <- fit.w
  # A <- ifelse(sapply(y.pred, function(x){ as.numeric( x[[j]] ) }) == 1, -1, 1)   # transform to -1 and 1
  # pred <- predict(w, A)
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

    
  acc[i] <- mean(pred == y.test)
}
cat("method :", cname[j],
    "\nB =", which.max(acc), 
    "\nAccuracy =", acc[which.max(acc)], "\n")
}

## majority vote
# method :  svm.linear 
# B = 4 
# Accuracy = 0.81 
# method :  svm.radial 
# B = 6 
# Accuracy = 0.83 
# method :  svm.sigmoid 
# B = 9 
# Accuracy = 0.8 
# method :  knn 
# B = 7 
# Accuracy = 0.89 

## reg
# method : svm.linear
# B = 28
# Accuracy = 0.81
# method : svm.radial
# B = 12
# Accuracy = 0.82
# method : svm.sigmoid
# B = 76
# Accuracy = 0.83
# method : knn
# B = 10
# Accuracy = 0.88

## logit-numeric
# method : svm.linear 
# B = 40 
# Accuracy = 0.83 
# method : svm.radial 
# B = 12 
# Accuracy = 0.82 
# method : svm.sigmoid 
# B = 42 
# Accuracy = 0.84 
# method : knn 
# B = 9 
# Accuracy = 0.88 

## logit-factor
# method : svm.linear 
# B = 40 
# Accuracy = 0.83 
# method : svm.radial 
# B = 12 
# Accuracy = 0.82 
# method : svm.sigmoid 
# B = 42 
# Accuracy = 0.84 
# method : knn 
# B = 9 
# Accuracy = 0.88 


## ridge-reg
# method : svm.linear 
# B = 4 
# Accuracy = 0.81 
# method : svm.radial 
# B = 12 
# Accuracy = 0.83 
# method : svm.sigmoid 
# B = 266 
# Accuracy = 0.83 
# method : knn 
# B = 12 
# Accuracy = 0.88 

## ridge-logit-numeric
# method : svm.linear 
# B = 4 
# Accuracy = 0.81 
# method : svm.radial 
# B = 3 
# Accuracy = 0.81 
# method : svm.sigmoid 
# B = 262 
# Accuracy = 0.81 
# method : knn 
# B = 19 
# Accuracy = 0.88 

## ridge-logit-factor
# method : svm.linear 
# B = 12 
# Accuracy = 0.82 
# method : svm.radial 
# B = 12 
# Accuracy = 0.83 
# method : svm.sigmoid 
# B = 265 
# Accuracy = 0.83 
# method : knn 
# B = 9 
# Accuracy = 0.88 






