setwd("C:\\Users\\user\\Desktop\\KHS\\FDA-Lab\\Sparse FPCA\\Application\\pc_selection_2")


# load simulated data
load("../simulated_data.RData")

# 10-fold CV
library(tidyverse)
library(doParallel)   # parallel computing
library(fdapace)
library(rlist)
library(MASS)

# parallel computing setting
ncores <- detectCores() - 2
registerDoParallel(ncores)


i <- 1
set.seed(100)
y <- factor(c(rep(0, 100), rep(1, 100)), level=c(0, 1))
data <- cbind(y, X.curves[[i]])
ind <- sample(1:nrow(data), 100)   # train set index
train <- as.matrix(data[ind, -1])
test <- as.matrix(data[-ind, -1])

# sparsify train, test set
train_data <- Sparsify(train, time, sparsity = c(2:10))
test_data <- Sparsify(test, time, sparsity = c(2:10))


packages <- c("fdapace","MASS")   # foreach에서 사용할 package 정의
k <- 3

X <- train_data
N <- length(X$Ly)

system.time(
# Eigenvector cross-validation
eig.cv <- foreach(i=1:N, .combine="rbind", .packages=packages) %dopar% {
  cv.error <- c()
  fit <- FPCA(X$Ly[-i], 
              X$Lt[-i], 
              optns=list(dataType="Sparse",
                         methodXi="CE",
                         FVEthreshold=1,   # fit할 수 있는 최대의 k까지 fit
                         verbose=T))
  
  time.grid <- fit$workGrid
  ind <- sapply(X$Lt[[i]], function(x){which.min( abs(x - time.grid) )})   # time point index
  x.center <- X$Ly[[i]] - fit$mu[ind]
  J <- length(X$Lt[[i]])
  
  for (k in 1:20) {
    if (k == 1) {
      eig.mat <- matrix(fit$phi[, 1:k], ncol=1)
    } else {
      eig.mat <- fit$phi[, 1:k]
    }
    
    loocv <- 0
    for (j in 1:J) {
      # predict FPC score of test set
      if (J > 2) {
        score.pred <- t(x.center[-j]) %*% eig.mat[ind[-j], ] %*% ginv( t(eig.mat[ind[-j], ]) %*% eig.mat[ind[-j], ] )
      } else {
        score.pred <- t(x.center[-j]) %*% eig.mat[ind[-j], ] %*% ginv( matrix(eig.mat[ind[-j], ], ncol=1) %*% matrix(eig.mat[ind[-j], ], nrow=1) )
      }
      
      x.hat <- score.pred %*% eig.mat[ind[j], ]
      loocv <- loocv + (x.center[j] - x.hat)^2
    }
    cv.error[k] <- loocv
  }
  return(cv.error)
}
)
cv.error <- colSums(eig.cv)
plot(1:length(cv.error), cv.error, type="o", xlab="No. of PCs", ylab="Cross-validation Error")
points(which.min(cv.error),
       cv.error[which.min(cv.error)],
       col="blue", pch=19, cex=1.5)
# score.pred <- predict(fit, 
#                       list(X$Ly[[i]][-j]), 
#                       list(X$Lt[[i]][-j]), 
#                       K=k)


### generalized inverse approach (approximate correct PRESS)
# https://stats.stackexchange.com/questions/93845/how-to-perform-cross-validation-for-pca-to-determine-the-number-of-principal-com/115477#115477
source("cv_fdapace.R")
fit <- cv.pca.sparse(train_data)
cv.error <- fit$cv.error
qplot(1:length(cv.error), 
      cv.error, 
      # type="o",
      xlab="Nunber of PCs",
      ylab="Cross-validation error",
      main="Leave-one-curve-out CV") + geom_line(colour="red")



### cross-validation => Overfitting
i <- 1
set.seed(100)
y <- factor(c(rep(0, 100), rep(1, 100)), level=c(0, 1))
ind <- sample(1:200, 100)
data <- cbind(y[ind], X.curves[[i]][ind, ])


start.time <- Sys.time()   # 연산속도 체크
# apply 10-fold cross validation
set.seed(777)
k <- 10  # CV parameter
samp <- sample(1:nrow(data), nrow(data))   # train set index
cv.error <- foreach(i = 1:k, .packages=packages) %dopar% {
  ind <- samp[(((i-1)*k)+1):(i*k)]
  # ind <- i   # LOOCV
  
  train <- as.matrix(data[-ind, -1])
  test <- as.matrix(data[ind, -1])
  
  # sparsify train, test set
  train_data <- Sparsify(train, time, sparsity=10)
  test_data <- Sparsify(test, time, sparsity=10)
  
  # fit FPCA
  fpca.train <- FPCA(train_data$Ly, 
                     train_data$Lt, 
                     optns=list(dataType="Sparse",
                                methodXi="CE",
                                FVEthreshold=1,   # fit할 수 있는 최대의 k까지 fit
                                rho="cv",
                                verbose=T))
  
  # pc 개수별로 mse 확인
  cv.error <- 0
  cv.error.train <- 0
  for (num.pc in 1:length(time)) {
    # predict FPC score of test set
    fpc.test <- predict(fpca.train, test_data$Ly, test_data$Lt, K=num.pc, xiMethod="CE")
    
    # plot predicted curve of test set
    pred.test <- matrix(rep(fpca.train$mu, nrow(test)), nrow(test), byrow=T) + fpc.test %*% t(fpca.train$phi[, 1:num.pc])
    pred.test <- as.list(as.data.frame(t(pred.test)))
    
    time.grid <- fpca.train$workGrid
    ind <- lapply(test_data$Lt, function(y){ sapply(y, function(x){which.min( abs(x - time.grid) )}) } )
    
    # calculate mse
    mse <- mapply(function(x, y, z){ mean( (z[x] - y)^2 ) }, 
                  ind,
                  test_data$Ly,
                  pred.test)
    cv.error[num.pc] <- mean(mse)
    
    
    # calculate traininig error
    fpc.train <- predict(fpca.train, train_data$Ly, train_data$Lt, K=num.pc, xiMethod="CE")
    
    # plot predicted curve of test set
    # pred.train <- matrix(rep(fpca.train$mu, nrow(train)), nrow(train), byrow=T) + fpc.train %*% t(fpca.train$phi[, 1:num.pc])
    pred.train <- fitted(fpca.train, 
                         ciOptns=list(alpha=0.05,
                                      cvgMethod = "interval"),                    
                         K=num.pc)
    pred.train <- as.list(as.data.frame(t(pred.train$fitted)))
    
    time.grid <- fpca.train$workGrid
    ind <- lapply(train_data$Lt, function(y){ sapply(y, function(x){which.min( abs(x - time.grid) )}) } )
    
    # calculate mse
    mse <- mapply(function(x, y, z){ mean( (z[x] - y)^2 ) }, 
                  ind,
                  train_data$Ly,
                  pred.train)
    cv.error.train[num.pc] <- mean(mse)
  }
  # # predict FPC score of test set
  # fpc.test <- predict(fpca.train, test_data$Ly, test_data$Lt, K=5)
  # 
  # # plot predicted curve of test set
  # pred.test <- matrix(rep(fpca.train$mu, nrow(test)), nrow(test), byrow=T) + fpc.test %*% t(fpca.train$phi[, 1:5])
  # pred.test <- as.list(as.data.frame(t(pred.test)))
  # 
  # time.grid <- fpca.train$workGrid
  # ind <- lapply(test_data$Lt, function(y){ sapply(y, function(x){which.min( abs(x - time.grid) )}) } )
  # 
  # # calculate mse
  # mse <- mean( mapply(function(x, y, z){ mean( (z[x] - y)^2 ) }, 
  #                     ind,
  #                     test_data$Ly,
  #                     pred.test) )
  return( list(test.error=cv.error,
               train.error=cv.error.train) )
}
end.time <- Sys.time()
end.time - start.time   # 연산 속도 측정


# CreatePathPlot( fpca.train, subset = 1:3, main = "GCV bandwidth", pch = 16)

test.error <- colMeans( list.rbind( lapply(cv.error, function(x){ x$test.error }) ) )
train.error <- colMeans( list.rbind( lapply(cv.error, function(x){ x$train.error }) ) )
which.min(test.error)
which.min(train.error)
train.error
test.error
plot(1:18, test.error, type="o", col="red", xlab="No. of PCs", ylab="CV Error")
lines(1:18, train.error, type="o")
points(which.min(test.error),
       test.error[which.min(test.error)],
       col="red", pch=19, cex=1.5)
points(which.min(train.error),
       train.error[which.min(train.error)],
       pch=19, cex=1.5)
legend("topright",
       c("Training CV error", "Test CV error"),
       col=c("black","red"),
       lty=c(1,1))

## Real data analysis => train set에 max time point가 없으면 test set에서 predict 못합
data <- read.csv("../spnbmd.csv", stringsAsFactors=F)
data <- data[which(data$sex=="fem" & data$ethnic=="White"), ]
dim(data)
library(dplyr)
source("GetLogLik.R")

ind <- data %>% 
  group_by(idnum) %>% 
  summarise(n=n())
ind <- ind[which(ind$n != 1), "idnum"]
data <- data[which(data$idnum %in% t(ind)), ]
dim(data)
data <- data[,c(1,3,5)]   # 필요없는 변수 제거
length(unique(data$idnum))   # 48


start.time <- Sys.time()   # 연산속도 체크
# apply 10-fold cross validation
set.seed(777)
id <- unique(data$idnum)
k <- length(id)
cv.error <- foreach(i = 1:k, .packages=packages) %dopar% {
  ind <- id[i]   # LOOCV
  
  # input data 형태로 변환
  train_data <- list(Lt=lapply(id[id != ind], function(x){ data$age[ which(data$idnum %in% x) ] }),
                     Ly=lapply(id[id != ind], function(x){ data$spnbmd[ which(data$idnum %in% x) ] }))
  test_data <- list(Lt=lapply(ind, function(x){ data$age[ which(data$idnum %in% x) ] }),
                    Ly=lapply(ind, function(x){ data$spnbmd[ which(data$idnum %in% x) ] }))
  
  # fit FPCA
  fpca.train <- FPCA(train_data$Ly, 
                     train_data$Lt, 
                     list(dataType="Sparse",
                          methodXi="CE",
                          FVEthreshold=1,   # fit할 수 있는 최대의 k까지 fit
                          rho="cv",
                          verbose=T))
  
  # pc 개수별로 mse 확인
  cv.error <- 0
  cv.error.train <- 0
  for (num.pc in 1:length(fpca.train$lambda)) {
    # predict FPC score of test set
    fpc.test <- predict(fpca.train, test_data$Ly, test_data$Lt, K=num.pc, xiMethod="CE")
    
    # plot predicted curve of test set
    pred.test <- matrix(rep(fpca.train$mu, length(test_data$Ly)), length(test_data$Ly), byrow=T) + fpc.test %*% t(fpca.train$phi[, 1:num.pc])
    pred.test <- as.list(as.data.frame(t(pred.test)))
    
    time.grid <- fpca.train$workGrid
    ind <- lapply(test_data$Lt, function(y){ sapply(y, function(x){which.min( abs(x - time.grid) )}) } )
    
    # calculate mse
    mse <- mapply(function(x, y, z){ mean( (z[x] - y)^2 ) }, 
                  ind,
                  test_data$Ly,
                  pred.test)
    cv.error[num.pc] <- mean(mse)
    
    
    # calculate traininig error
    fpc.train <- predict(fpca.train, train_data$Ly, train_data$Lt, K=num.pc, xiMethod="CE")
    
    # plot predicted curve of test set
    pred.train <- fitted(fpca.train, 
                         ciOptns=list(alpha=0.05,
                                      cvgMethod = "interval"),                    
                         K=num.pc)
    pred.train <- as.list(as.data.frame(t(pred.train$fitted)))
    
    time.grid <- fpca.train$workGrid
    ind <- lapply(train_data$Lt, function(y){ sapply(y, function(x){which.min( abs(x - time.grid) )}) } )
    
    # calculate mse
    mse <- mapply(function(x, y, z){ mean( (z[x] - y)^2 ) }, 
                  ind,
                  train_data$Ly,
                  pred.train)
    cv.error.train[num.pc] <- mean(mse)
  }
  return( list(test.error=cv.error,
               train.error=cv.error.train) )
}
end.time <- Sys.time()
end.time - start.time   # 연산 속도 측정


test.error <- colMeans( list.rbind( lapply(cv.error, function(x){ x$test.error }) ) )
train.error <- colMeans( list.rbind( lapply(cv.error, function(x){ x$train.error }) ) )
which.min(test.error)
which.min(train.error)
plot(1:18, test.error, type="o")
lines(1:18, train.error, type="o", col="red")













### Example
i <- 1
set.seed(100)
y <- factor(c(rep(0, 100), rep(1, 100)), level=c(0, 1))
data <- cbind(y, X.curves[[i]])
ind <- sample(1:nrow(data), 100)   # train set index
train <- as.matrix(data[ind, -1])
test <- as.matrix(data[-ind, -1])


# # create input shape data
# data <- MakeFPCAInputs(IDs = rep(1:100, each=length(X.curves[[i]])),
#                        tVec = rep(time, times = 100),
#                        yVec = t(X.curves[[i]]))

# sparsify train, test set
train_data <- Sparsify(train, time, sparsity = c(2:10))
test_data <- Sparsify(test, time, sparsity = c(2:10))

# fit FPCA
fpca.train <- FPCA(train_data$Ly, 
                   train_data$Lt, 
                   list(dataType="Sparse",
                        methodXi="CE",
                        rho="cv",
                        verbose=T))

# plot(fpca.train)

# CreateModeOfVarPlot(fpca.train, main="Confidence Interval of 1st FPC")
# CreatePathPlot(fpca.train, subset = 1:3, main = "GCV bandwidth", pch = 16)

pred.train <- fitted(fpca.train, ciOptns = list(alpha=0.05,
                                                cvgMethod = "interval"))

# # plot fitted curve and 95% CI of train set
# time.grid <- fpca.train$workGrid
# upper.limit <- pred.train$cvgUpper
# lower.limit <- pred.train$cvgLower
# par(mfrow=c(3,3))
# ind <- sample(1:n, 9)
# for (i in 1:9) {
#   j <- ind[i]
#   plot(time.grid, upper.limit[j,], type='l', col=4, lty=2,
#        ylim=c(min(cvgLower[j,]), max(cvgUpper[j,])),
#        xlab='minutes', ylab='X(t)', 
#        main=paste(j,'th predicted trajectory of curve', sep=''))
#   points(time.grid, lower.limit[j,], type='l', col=4, lty=2)
#   points(train_data$Lt[[j]], train_data$Ly[[j]])    # true points
#   points(time.grid, pred.train$fitted[j,], type='l', col="red")   # fitted curve
# }


# predict FPC score of test set
fpc.test <- predict(fpca.train, test_data$Ly, test_data$Lt, K=5)
plot(fpc.test[, 1], fpc.test[, 2])


# # plot predicted curve of test set
# pred.test <- matrix(rep(fpca.train$mu, nrow(test)), nrow(test), byrow=T) + fpc.test %*% t(fpca.train$phi[, 1:5])
# 
# par(mfrow=c(2,3))
# for (i in 1:6){
#   plot(time.grid, pred.test[i, ], type="l", col="red", xlab='minutes', ylab='X(t)')  
#   points(test_data$Lt[[i]], test_data$Ly[[i]])  
# }
