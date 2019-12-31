setwd("C:\\Users\\user\\Desktop\\KHS\\FDA-Lab\\Sparse FPCA\\Application")

library(tidyverse)
library(doParallel)   # parallel computing
library(fdapace)
library(rlist)
library(MASS)

# load simulated data
load("simulated_data.RData")

# parallel computing setting
ncores <- detectCores() - 3
registerDoParallel(ncores)

packages <- c("fdapace","MASS")   # foreach에서 사용할 package 정의


set.seed(100)
s <- 10   # sparsity
# generate data
i <- 1
y <- factor(c(rep(0, 100), rep(1, 100)), level=c(0, 1))
data <- cbind(y, X.curves[[i]])
ind <- sample(1:nrow(data), 100)   # train set index
X.train <- as.matrix(data[ind, -1])
X.test <- as.matrix(data[-ind, -1])
y.train <- y[ind]
y.test <- y[-ind]

# sparsify train, test set
train.data <- Sparsify(X.train, time, sparsity = s)
test.data <- Sparsify(X.test, time, sparsity = s)

## Bootstrap aggregating
B <- 10
N <- nrow(train)
y.pred <- foreach(b=1:B, .combine="cbind", .packages=packages) %dopar% {
  boot.ind <- sample(1:N, N, replace=T)
  
  fit <- FPCA(train.data$Ly[boot.ind], 
              train.data$Lt[boot.ind], 
              optns=list(dataType="Sparse",
                         methodXi="CE",
                         methodSelectK="BIC",
                         verbose=T))
  K <- fit$selectK
  fpc.score <- fit$xiEst
  
  train <- data.frame(y=y.train, 
                      x=fpc.score)
  colnames(train) <- c("y", paste("FPC", 1:K, sep=""))
  
  logit <- glm(y~., train, family=binomial)
  
  # test data
  fpc.score.test <- predict(fit, test_data$Ly, test_data$Lt, K=K, xiMethod="CE")
  test <- data.frame(y=y.test, 
                     x=fpc.score.test)
  colnames(test) <- c("y", paste("FPC", 1:K, sep=""))
  pred <- predict(logit, test, type="response")
  pred <- ifelse(y.pred > 0.5, 1, 0)
  
  return(pred)
  # if (b==1) {
  #   y.pred <- pred
  # } else {
  #   y.pred <- cbind(y.pred, 
  #                   pred)
  # }
}

# 2 class majority voting functon
majority_vote <- function(x) {
  key <- unique(x)
  if (length(key) > 1) {
    val <- which.max( table(y.pred[1, ]) ) - 1
  } else {
    val <- key
  }
  
  return(val)
}

pred <- apply(y.pred, 1, majority_vote)

table(pred, y.test)
