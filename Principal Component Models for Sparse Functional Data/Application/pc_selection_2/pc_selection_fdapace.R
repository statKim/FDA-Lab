library(fdapace)

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

# fit FPCA
fpca.train <- FPCA(train_data$Ly, 
                   train_data$Lt, 
                   list(dataType="Sparse",
                        methodXi="CE",
                        FVEthreshold=1,   # fit할 수 있는 최대의 k까지 fit
                        rho="cv",
                        verbose=T))

# predict FPC score of test set
fpc.test <- predict(fpca.train, test_data$Ly, test_data$Lt, K=5)

# plot predicted curve of test set
pred.test <- matrix(rep(fpca.train$mu, nrow(test)), nrow(test), byrow=T) + fpc.test %*% t(fpca.train$phi[, 1:5])
pred.test <- as.list(as.data.frame(t(pred.test)))

time.grid <- fpca.train$workGrid
ind <- lapply(test_data$Lt, function(y){ sapply(y, function(x){which.min( abs(x - time.grid) )}) } )

# calculate mse
mse <- mean( mapply(function(x, y, z){ mean( (z[x] - y)^2 ) }, 
                    ind,
                    test_data$Ly,
                    pred.test) )
mse



# 10-fold CV
library(tidyverse)
library(e1071)   # SVM fit하기 위함
library(doParallel)   # parallel computing
library(fdapace)

# parallel computing setting
ncores <- detectCores() - 3
registerDoParallel(ncores)

packages <- c("tidyverse","fdapace","e1071")   # foreach에서 사용할 package 정의


i <- 1
set.seed(100)
y <- factor(c(rep(0, 100), rep(1, 100)), level=c(0, 1))
data <- cbind(y, X.curves[[i]])

k <- 10  # CV parameter
samp <- sample(1:nrow(data), nrow(data))   # train set index

# apply cross validation
res <- foreach(i = 1:10, .combine="rbind", .packages=packages) %dopar% {
  ind <- samp[(((i-1)*k)+1):(i*k)]
  
  train <- as.matrix(data[-ind, -1])
  test <- as.matrix(data[ind, -1])
  
  # sparsify train, test set
  train_data <- Sparsify(train, time, sparsity = c(2:10))
  test_data <- Sparsify(test, time, sparsity = c(2:10))
  
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
  for (num.pc in 2:10) {
    # predict FPC score of test set
    fpc.test <- predict(fpca.train, test_data$Ly, test_data$Lt, K=num.pc)
    
    # plot predicted curve of test set
    pred.test <- matrix(rep(fpca.train$mu, nrow(test)), nrow(test), byrow=T) + fpc.test %*% t(fpca.train$phi[, 1:num.pc])
    pred.test <- as.list(as.data.frame(t(pred.test)))
    
    time.grid <- fpca.train$workGrid
    ind <- lapply(test_data$Lt, function(y){ sapply(y, function(x){which.min( abs(x - time.grid) )}) } )
    
    # calculate mse
    mse <- mean( mapply(function(x, y, z){ mean( (z[x] - y)^2 ) }, 
                        ind,
                        test_data$Ly,
                        pred.test) )
    cv.error[num.pc] <- mse
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
  return( cv.error )
}












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
