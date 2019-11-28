setwd("C:\\Users\\user\\Desktop\\KHS\\FDA-Lab\\Principal Component Models for Sparse Functional Data\\Application\\pc_selection_2")


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



### pseudo inverse approach (approximate correct PRESS)
# https://stats.stackexchange.com/questions/93845/how-to-perform-cross-validation-for-pca-to-determine-the-number-of-principal-com/115477#115477
cv.pca.sparse <- function(X, packages=c("fdapace"), plot=T) {
  fit <- FPCA(X$Ly, 
              X$Lt, 
              optns=list(dataType="Sparse",
                         methodXi="CE",
                         FVEthreshold=1,   # fit할 수 있는 최대의 k까지 fit
                         verbose=T))
  
  N <- length(X$Ly)   # number of curves
  
  # calculate AIC, BIC
  aic <- c()
  bic <- c()
  # AIC <- foreach(k=1:length(fit$lambda), .combine="c", .packages=packages) %dopar% {
  for (k in 1:length(fit$lambda)) {
    aic[k] <- GetLogLik(fit, k) + 2*k
    bic[k] <- GetLogLik(fit, k) + k*log(N)
  }
  
  # approximate PCA PRESS
  cv.error <- foreach(k=1:length(fit$lambda), .combine="c", .packages=packages) %dopar% {
    U <- fit$phi[, 1:k]   # eigenvector matrix without ith curve
    if (!is.matrix(U)) { U <- matrix(U, ncol=1) }  # transform to matrix
    uu.t <- U %*% t(U)
    w <- diag(1, nrow(uu.t)) - uu.t + diag(diag(uu.t))
    
    # calculate PRESS
    err <- mapply(function(t, y){
      time.ind <- sapply(t, function(x){ which.min( abs(x - fit$workGrid) ) })   # time point에 해당되는 grid 찾기
      sum( (w[, time.ind] %*% y)^2 )
    },
    X$Lt,
    X$Ly)
    
    return( sum(err) )
  }
  
  # 각 measure로 plot
  if (plot==T) {
    par(mfrow=c(1,3))
    plot(1:length(aic), 
         aic, 
         type="o",
         xlab="Nunber of PCs",
         ylab="AIC",
         main="AIC")
    points(which.min(aic),
           aic[which.min(aic)],
           col="red",
           pch=19,
           cex=2)    
    
    plot(1:length(bic), 
         bic, 
         type="o",
         xlab="Nunber of PCs",
         ylab="BIC",
         main="BIC")
    points(which.min(bic),
           bic[which.min(bic)],
           col="red",
           pch=19,
           cex=2)
    
    plot(1:length(cv.error), 
         cv.error, 
         type="o",
         xlab="Nunber of PCs",
         ylab="Cross-validation error",
         main="Leave-one-curve-out CV")
    points(which.min(cv.error),
           cv.error[which.min(cv.error)],
           col="red",
           pch=19,
           cex=2)
  }

  
  return( list(cv.error=cv.error,
               AIC=aic,
               BIC=bic,
               selected.K=data.frame(AIC=which.min(aic),
                                     BIC=which.min(bic),
                                     LOOCV=which.min(cv.error)),
               model=fit) )
  
  
  
  
  
  # cv.error <- foreach(i=1:length(X$Lt), .combine="rbind", .packages=packages) %dopar% {
  #   x.i <- X$Ly[[i]]   # ith curve
  #   t.i <- X$Lt[[i]]   # ith curve's time points
  #   
  #   fit <- FPCA(X$Ly[-i], 
  #               X$Lt[-i], 
  #               optns=list(dataType="Sparse",
  #                          methodXi="CE",
  #                          FVEthreshold=1,   # fit할 수 있는 최대의 k까지 fit
  #                          rho="cv",
  #                          verbose=T))
  #   
  #   time.grid <- fit$workGrid   # time grid
  #   err <- c()
  #   for (k in 1:length(fit$lambda)) {
  #     U <- fit$phi[, 1:k]   # eigenvector matrix without ith curve
  #     dec <- svd(U[-i, ])
  #     U.plus <- dec$v %*% diag(1/dec$d, ncol(dec$v), ncol(dec$u)) %*% t(dec$u)
  #     
  #     # 아마 time point 찾아서 계산해줘야 될듯
  #     
  #     time.ind <- time.grid[ sapply(t.i, function(x){ which.min( abs(x - time.grid) ) }) ]
  #     sub <- U %*% U.plus %*% x.i
  #     
  #     err[k] <- sum( (x.i - sub)^2 )
  #   }
  # }
  # 
  # cv.error <- data.frame(cv.error)
  # colnames(cv.error) <- 1:ncol(fit$phi)
  # return( colMeans(cv.error) )
}
fit <- cv.pca.sparse(train_data)

cv.error <- fit$cv.error

qplot(1:length(cv.error), 
      cv.error, 
      # type="o",
      xlab="Nunber of PCs",
      ylab="Cross-validation error",
      main="Leave-one-curve-out CV") + geom_line(colour="red")

# load simulated data
load("../simulated_data.RData")

# 10-fold CV
library(tidyverse)
library(doParallel)   # parallel computing
library(fdapace)
library(rlist)

# parallel computing setting
ncores <- detectCores() - 3
registerDoParallel(ncores)

packages <- c("tidyverse","fdapace")   # foreach에서 사용할 package 정의


i <- 5
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
plot(1:18, test.error, type="o")
lines(1:18, train.error, type="o", col="red")




## Real data analysis => train set에 max time point가 없으면 test set에서 predict 못합
data <- read.csv("../spnbmd.csv", stringsAsFactors=F)
data <- data[which(data$sex=="fem" & data$ethnic=="White"), ]
dim(data)
library(dplyr)
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
