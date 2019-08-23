setwd("C:\\Users\\user\\Desktop\\KHS\\Thesis\\seminar\\application\\logistic_reg")

#####################
### Data generating
#####################

library(fda.usc)
library(dplyr)
library(ggplot2)
library(data.table)

set.seed(100)
cell <- read.delim("mbc_9_12_3273__CDCDATA.txt")
gene <- cell[, c(1, 6:23)] # 0~120min, 18 points
gene <- na.omit(gene)

train <- gene[, -1] # remove the column of gene's names
train.fdata <- fdata(train)
train.fd <- fdata2fd(train.fdata, nbasis=5)
train.pca <- create.pc.basis(train.fdata, l=1:5, lambda=1)
train.pca$basis$data <- -train.pca$basis$data

# epsilon 1~5, group label을 성분으로 하는 100개의 datasets(trian 100, test 100) 생성
means <- seq(0.6, 0.2, -0.1)
variances <- c(2.6950, 0.8850, 0.1957, 0.1266, 0.1079)
random.sets <- list()
for (i in 1:100) {
  e1 <- append( rnorm(100, means[1], sqrt(variances[1])), rnorm(100, -means[1], sqrt(variances[1])) )
  e2 <- append( rnorm(100, means[2], sqrt(variances[2])), rnorm(100, -means[2], sqrt(variances[2])) )
  e3 <- append( rnorm(100, means[3], sqrt(variances[3])), rnorm(100, -means[3], sqrt(variances[3])) )
  e4 <- append( rnorm(100, means[4], sqrt(variances[4])), rnorm(100, -means[4], sqrt(variances[4])) )
  e5 <- append( rnorm(100, means[5], sqrt(variances[5])), rnorm(100, -means[5], sqrt(variances[5])) )
  group <- as.factor(c( rep("g1", 100), rep("g0", 100) ))
  random.set <- list(e1=e1, e2=e2, e3=e3, e4=e4, e5=e5, group=group)
  set.i <- paste("set", as.character(i), sep="")
  random.sets[[set.i]] <- random.set
}
random.sets$set1$e1[1:10]
append( rnorm(10, means[1], sqrt(variances[1])), rnorm(10, -means[1], sqrt(variances[1])) )
e1.new <- as.data.frame((random.sets$set1$e1[c(1:5)] %o% train.pca$basis$data[1,]))
e2.new <- as.data.frame(t(random.sets$set1$e2[c(1:5)] %o% train.pca$basis$data[2,]))

# mean curve에 mth eigenfunction_m과 epsilon_m을 곱한 값을 더한 X samples 생성
X.curves <- list()
mean.curve <- as.vector(train.pca$mean$data)
pca.coefs <- train.pca$basis$data
i <- 1
for (set.i in random.sets) {
  sum.e <- mean.curve # 1 datase합t당 PC X e들의 합
  for (j in 1:5) {
    ej <- paste("e", as.character(j), sep="")
    ej.new <- as.data.frame( (set.i[[ej]] %o% pca.coefs[j,]) )
    sum.e <- ej.new + sum.e
  }
  X.curves[[names(random.sets)[i]]] <- sum.e
  i <- i + 1
}

length(X.curves)
dim(X.curves[[1]])
random.set$group


#######################################
### functional logistic regression
#######################################
# time points 변수 만들기
time <- colnames(gene)[-1]
time <- as.integer(gsub(pattern="alpha|min",
                        replacement="",
                        x=time))

# reduced rank model에 맞도록 데이터 변환
source("../EM_reduced_rank.R")
library(dplyr)
library(reshape2)

system.time({
  for (k in 1:5) {
    # 100개의 generated data로 accuracy 계산
    k <- 5
    accuracy <- c()
    result <- list()
    for (i in 1:100) {
      set.seed(100)
      y <- factor(c(rep(0, 100), rep(1, 100)), level=c(0, 1))
      data <- cbind(y, X.curves[[i]])
      ind <- sample(1:nrow(data), 100)   # train set index
      train <- data[ind, ]
      test <- data[-ind, ]
      
      
      # PC score 계산(PACE, Yao et al.)
      # pc.score <- matrix(0, nrow=nrow(train), ncol=k)
      # y.center <- as.matrix(train[,-1]) - t( matrix(rep(colMeans(train[,-1]), 100), nrow=18) )
      # cov.fn <- t(y.center) %*% y.center / nrow(train)
      # eig <- eigen(cov.fn)
      # lambda <- eig$values
      # phi <- eig$vectors
      # for (j in 1:nrow(train)) {
      #   # j <- 1
      #   y.center <- as.matrix( train[j,-1] - t( colMeans(train[,-1]) ), ncol=1)
      #   # PC score 계산(PACE, Yao et al.)
      #   cov.fn <- t(y.center) %*% y.center
      #   
      #   pc.score[j, ] <- diag(lambda[1:k]) %*% t(phi[, 1:k]) %*% ginv(cov.fn) %*% t(y.center)
      # }
      # 
      # pred.curve <- colMeans(train[,-1]) + phi[, 1:k] %*% as.matrix(pc.score[1, ], ncol=1)
      # 
      # round(lambda / sum(lambda), 2)
      # 
      # plot(1:18, train[1, -1], type="l")
      # lines(1:18, pred.curve, col="2")
      
      # PC score 계산(PACE, Yao et al.)
      pc.score <- matrix(0, nrow=nrow(train), ncol=k)
      y.center <- as.matrix(train[,-1]) - t( matrix(rep(colMeans(train[,-1]), 100), nrow=18) )
      cov.fn <- t(y.center) %*% y.center / nrow(train)
      eig <- eigen(cov.fn)
      lambda.train <- eig$values
      phi.train <- eig$vectors
      for (j in 1:nrow(train)) {
        # j <- 1
        y.center <- as.matrix( train[j,-1] - t( colMeans(train[,-1]) ), ncol=1)
        # PC score 계산(PACE, Yao et al.)
        cov.fn <- t(y.center) %*% y.center

        pc.score[j, ] <- diag(lambda.train[1:k]) %*% t(phi.train[, 1:k]) %*% ginv(cov.fn) %*% t(y.center)
      }
      train.new <- cbind(train$y, as.data.frame(pc.score))
      colnames(train.new) <- c("y", paste("PC", 1:k))

      
      pc.score <- matrix(0, nrow=nrow(train), ncol=k)
      y.center <- as.matrix(test[,-1]) - t( matrix(rep(colMeans(test[,-1]), 100), nrow=18) )
      cov.fn <- t(y.center) %*% y.center / nrow(test)
      eig <- eigen(cov.fn)
      lambda.test <- eig$values
      phi.test <- eig$vectors
      
      # PC function의 부호 같게 바꾸기(각 PC function의 첫 번째 값의 부호에 따라)
      if (is.null(dim(phi.train))) {   # k=1인 경우
        if (sign(phi.train[1]) != sign(phi.test[1])) {
          phi.test <- -phi.test
        }
      } else {
        for (col in 1:k) {
          if (sign(phi.train[1, col]) != sign(phi.test[1, col])) {
            phi.test[, col] <- -phi.test[, col]
          }
        }
      }
      
      for (j in 1:nrow(test)) {
        # j <- 1
        y.center <- as.matrix( test[j,-1] - t( colMeans(test[,-1]) ), ncol=1)
        # PC score 계산(PACE, Yao et al.)
        cov.fn <- t(y.center) %*% y.center
        
        pc.score[j, ] <- diag(lambda.test[1:k]) %*% t(phi.test[, 1:k]) %*% ginv(cov.fn) %*% t(y.center)
      }
      test.new <- cbind(test$y, as.data.frame(pc.score))
      colnames(test.new) <- c("y", paste("PC", 1:k))    

      
      
      # fit logistic regression
      fit.logit <- glm(y ~ ., data=train.new, family = "binomial")
      # summary(fit.logit)
      
      # predict
      test.data <- test.new[, -1]
      if (!is.data.frame(test.data)) {   # k=1인 경우(vector로 변환되기 때문)
        test.data <- data.frame(test.data)
        colnames(test.data) <- "PC 1"
      }
      pred <- predict(fit.logit, test.data, type="response")
      # pred <- factor(ifelse(pred > 0.5, 1, 0), levels=c(0, 1))
      pred <- factor(round(pred), levels=c(0, 1))
      acc <- mean(pred == test.new$y)   # accuracy 계산
      
      
      # output 저장
      accuracy[i] <- acc   # accuracy
      res <- list("accuracy"=acc,
                  "cross.table"=table(pred, true=test.new$y))
      result[[i]] <- res
      
      print( paste(i, "th data's Accuracy :", acc) )
      # print( rbind(pred, true=test.new$y) )
    }
    # mean accuracy
    print( paste("Mean accuracy :", mean(accuracy)) )
    
    # RData로 저장
    # save(list=c("result", "accuracy"), file=paste("result_PC", k, ".RData", sep=""))
  }
})





{
  # train FPC plot
  par(mfrow=c(2,3))
  # PC function
  for (i in 1:5){
    plot(1:18, phi.train[,i], type="l", lwd=3,
         xlab="minutes", ylab="Princ. comp.", main=paste("PC", i))
  }
  
  # test plot
  par(mfrow=c(2,3))
  # PC function
  for (i in 1:5){
    plot(1:18, phi.test[,i], type="l", lwd=3,
         xlab="minutes", ylab="Princ. comp.", main=paste("PC", i))
  }
}
