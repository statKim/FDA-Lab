#########################################
### Simulation
###   Methods for Sparse Functional Data
###   by Edwin Kam Fai Lei
#########################################
library(tidyverse)    # dplyr, ggplot2, etc
library(doParallel)   # parallel computing
library(fdapace)      # functional PCA
library(e1071)        # SVM, Naive bayes
# library(class)        # kNN
library(MASS)         # LDA, QDA
library(data.table)   # list rbind
library(gridExtra)    # subplot in ggplot2
library(reshape2)     # melt function

# parallel computing setting
ncores <- detectCores() - 2
registerDoParallel(ncores)

packages <- c("fdapace","e1071","MASS","tidyverse","gridExtra")   # foreach에서 사용할 package 정의

# function of calculating mode
Mode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


### generate the simulated dataset
n <- 100   # number of observations
j <- 1:40  # number of eigenvalues
# t <- seq(0, 1, length.out = 101)   # time points

# parameters of 3 different simulation studies
A <- list(lambda = list(lambda.0 = j^(-3),
                        lambda.1 = 4*j^(-2)),
          mu = list(mu.0 = c(rep(1, 4), rep(0, 36)),
                    mu.1 = rep(0, 40)))
B <- list(lambda = list(lambda.0 = j^(-3),
                        lambda.1 = 4*j^(-2)),
          mu = list(mu.0 = rep(0, 40),
                    mu.1 = rep(0, 40)))
C <- list(lambda = list(lambda.0 = 3*j^(-2),
                        lambda.1 = 3*j^(-2)),
          mu = list(mu.0 = c(rep(1, 4), rep(0, 36)),
                    mu.1 = rep(0, 40)))

y.pred <- foreach(simm=1:30, .combine="rbind", .packages=packages) %dopar% {
print(simm)
set.seed(simm)
lambda <- A$lambda
mu <- A$mu

# # eigenfunctions
# phi <- sapply(j, function(x){ 
#   if (x %% 2 == 0) {
#     sin(pi*t*x/2) / sqrt(2)
#   } else {
#     cos(pi*t*x/2) / sqrt(2)
#   } 
# })

## generate curves
y <- factor(c(rep(0, n/2), rep(1, n/2)), levels=c(0, 1))
data <- list()
# random sparsify => generate curve
for (i in 1:n) {
  k <- ifelse(i > n/2, 2, 1)   # different mu and theta for each k
  
  # random parameter
  theta <- sapply(j, function(x){
    p <- runif(1)
    if (p > 1/2) {
      rnorm(1, -sqrt(lambda[[k]][x]/2), sqrt(lambda[[k]][x]/2))
    } else {
      rnorm(1, sqrt(lambda[[k]][x]/2), sqrt(lambda[[k]][x]/2))
    }
  })
  
  # random sparsify
  num.obs <- sample(5:14, 1)
  t <- sort(runif(num.obs, 0, 1))
  
  # eigenfunctions
  phi <- sapply(j, function(x){ 
    if (x %% 2 == 0) {
      sin(pi*t*x/2) / sqrt(2)
    } else {
      cos(pi*t*x/2) / sqrt(2)
    } 
  })
  
  # measurement error
  eps <- rnorm(num.obs, 0, 0.1)
  
  # generate the curve with measurement error
  X <- as.numeric( matrix((theta + mu[[k]]), nrow=1) %*% t(phi) ) + eps
  
  data$id[[i]] <- rep(i, num.obs)
  data$y[[i]] <- rep(y[i], num.obs)
  data$Lt[[i]] <- t
  data$Ly[[i]] <- X
}

# # random sparsify => generate curve
# for (i in 1:n) {
#   k <- ifelse(i > n/2, 2, 1)   # different mu and theta for each k
# 
#   # random parameter
#   theta <- sapply(j, function(x){
#     p <- runif(1)
#     if (p > 1/2) {
#       rnorm(1, -sqrt(lambda[[k]][x] / 2), sqrt(lambda[[k]][x] / 2))
#     } else {
#       rnorm(1, sqrt(lambda[[k]][x] / 2), sqrt(lambda[[k]][x] / 2))
#     }
#   })
# 
#   # random sparsify
#   num.obs <- sample(5:14, 1)
#   ind <- sort(sample(1:length(t), num.obs))
# 
#   # measurement error
#   eps <- rnorm(num.obs, 0, 0.1)
# 
#   # generate the curve with measurement error
#   X <- as.numeric( matrix((theta + mu[[k]]), nrow=1) %*% t(phi[ind, ]) ) + eps
# 
#   data$id[[i]] <- rep(i, num.obs)
#   data$y[[i]] <- rep(y[[i]], num.obs)
#   data$Lt[[i]] <- t[ind]
#   data$Ly[[i]] <- X
# }

# # generate curve => random sparsify
# for (i in 1:n) {
#   k <- ifelse(i > n/2, 2, 1)   # different mu and theta for each k
# 
#   # random parameter
#   theta <- sapply(j, function(x){
#     p <- runif(1)
#     if (p > 1/2) {
#       rnorm(1, -sqrt(lambda[[k]][x] / 2), sqrt(lambda[[k]][x] / 2))
#     } else {
#       rnorm(1, sqrt(lambda[[k]][x] / 2), sqrt(lambda[[k]][x] / 2))
#     }
#   })
# 
#   # measurement error
#   eps <- rnorm(length(t), 0, 0.1)
# 
#   # generate each curves with measurement error
#   X <- as.numeric( matrix((theta + mu[[k]]), nrow=1) %*% t(phi) ) + eps
# 
#   # random sparsify
#   num.obs <- sample(5:14, 1)
#   ind <- sort(sample(1:length(t), num.obs))
# 
#   data$id[[i]] <- rep(i, num.obs)
#   data$y[[i]] <- rep(y[[i]], num.obs)
#   data$Lt[[i]] <- t[ind]
#   data$Ly[[i]] <- X[ind]
# }

# plot the generated curves
sapply(data, unlist) %>%
  data.frame() %>%
  mutate(y = ifelse(unlist(data$y) == 0, "G1", "G2")) %>%
  ggplot(aes(x=Lt, y=Ly, group=id, color=y)) +
  geom_line() +
  xlab("Time") +
  ylab("") +
  theme_bw() +
  theme(legend.title = element_blank())



### train, test split
# => train set에 없는 time point 포함시 FPCA할 때, error 발생
# range(test set) > range(train set) => 넘는 부분에 해당하는 id 제거
id <- 1:n
ind.train <- sample(id, ceiling(n/3))
ind.test <- NULL
range.train <- range(unlist(data$Lt[ind.train]))
range.test <- range(unlist(data$Lt[-ind.train]))
if (range.test[1] < range.train[1]) {   # min(test) < min(train) 인 경우
  ind <- sapply(data$Lt[-ind.train], function(x){ sum(x < range.train[1]) == 0 })
  ind.test <- id[-ind.train][ind]
}
if (range.test[2] > range.train[2]) {   # max(test) > max(train) 인 경우
  ind <- sapply(data$Lt[-ind.train], function(x){ sum(x > range.train[2]) == 0 })
  if (is.numeric(ind.test)) {
    ind.test <- intersect(ind.test,
                          id[-ind.train][ind])
  } else {
    ind.test <- id[-ind.train][ind]
  }
}
if (!is.numeric(ind.test)) {   # range(test) < range(train) 인 경우
  ind.test <- id[-ind.train]
}

# split to train, test set
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
fpca.fit <- FPCA(X.train$Ly, 
                 X.train$Lt, 
                 optns=list(dataType="Sparse",
                            methodXi="CE",
                            # methodSelectK="BIC",
                            FVEthreshold=0.99,
                            verbose=F))
k <- fpca.fit$selectK   # optimal number of PCs
# if (k > 2) {
#   k <- 2
# }
fpc.score <- fpca.fit$xiEst[, 1:k]   # FPC scores

# train FPC scores
train.fpc <- data.frame(y = y.train, 
                        x = fpc.score)
colnames(train.fpc) <- c("y", paste("FPC", 1:k, sep=""))

# # tune the hyperparmeters
# tune.linear <- tune.svm(y ~ .,
#                         data = train.fpc,
#                         kernel = "linear",
#                         cost = c(10^(-3:1), 2^(5:10)))
# tune.radial <- tune.svm(y ~ .,
#                         data = train.fpc,
#                         kernel = "radial",
#                         cost = c(10^(-3:1), 2^(5:10)),
#                         gamma = c(10^(-3:1), 2^(5:10)))

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
                       x = fpc.score.test[, 1:k])
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

# save the accuracy
err.single <- sapply(pred, function(x){ mean(x != y.test) })
err.single
}
y.pred
colMeans(y.pred)
apply(y.pred, 2, sd)
round(100* (1-colMeans(y.pred)), 3)


# plot first 2 PCs
p1 <- train.fpc %>% 
  mutate(y = ifelse(y == 0, "G1", "G2")) %>% 
  ggplot(aes(x=FPC1, y=FPC2, color=y)) +
  geom_point() +
  ggtitle("Training set") +
  theme_bw() +
  theme(legend.title = element_blank())
p2 <- test.fpc %>% 
  mutate(y = ifelse(y == 0, "G1", "G2")) %>% 
  ggplot(aes(x=FPC1, y=FPC2, color=y)) +
  geom_point() +
  ggtitle("Test set") +
  theme_bw() +
  theme(legend.title = element_blank())
grid.arrange(p1, p2, ncol=1)



### Bootstrap aggregating
B <- 100
N <- length(X.train$Ly)
set.seed(100)
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
acc.majority <- apply(pred[, 3:8], 2, function(x){ mean(x == pred$y) })
acc.majority

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
acc.oob <- apply(pred[, 3:8], 2, function(x){ mean(x == pred$y) })
acc.oob



## iteration마다의 error rate 그래프
for (i in 2:B) {
  res <- as.data.frame(rbindlist(lapply(y.pred[1:i], function(x){ x$boot })))
  
  # majority voting
  pred <- res %>% 
    group_by(id, y) %>% 
    summarise(logit = Mode(logit),
              svm.linear = Mode(svm.linear),
              svm.radial = Mode(svm.radial),
              lda = Mode(lda),
              qda = Mode(qda),
              nb = Mode(nb))
  err <- apply(pred[, 3:8], 2, function(x){ mean(x != pred$y) })
  
  if (i == 2) {
    result <- err
  } else {
    result <- rbind(result, err)
  }
}
rownames(result) <- 2:B
colnames(result) <- c("Logit","SVM(Linear)","SVM(Gaussian)","LDA","QDA","NaiveBayes")

res2 <- reshape2::melt(result)
colnames(res2) <- c("iter","method","error")
ggplot(res2, aes(x=iter, y=error, color=method)) +
  geom_line(size=1) +
  ylab("Classification Error Rate") +
  xlab("Iterations") +
  theme_bw() +
  theme(legend.title = element_blank())



