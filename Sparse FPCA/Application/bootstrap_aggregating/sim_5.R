#########################################
### Simulation
###   Probability-enhanced effective dimension reduction
###   for classifying sparse functional data
###   Yao et al.
###   Link: https://link.springer.com/content/pdf/10.1007%2Fs11749-015-0470-2.pdf
#########################################
library(tidyverse)    # dplyr, ggplot2, etc
library(doParallel)   # parallel computing
library(fdapace)      # functional PCA
library(e1071)        # SVM, Naive bayes
library(MASS)         # LDA, QDA
library(data.table)   # list rbind
library(gridExtra)    # subplot in ggplot2
library(reshape2)     # melt function

# parallel computing setting
ncores <- detectCores() - 2
registerDoParallel(ncores)

packages <- c("fdapace","e1071","class","MASS")   # foreach에서 사용할 package 정의

# y.pred <- foreach(num=1:100, .combine="rbind", .packages=packages) %dopar% {
### generate the simulated dataset
set.seed(100*num)
n <- 700   # number of observations
j <- 1:50  # number of eigenvalues
t <- seq(0, 10, length.out = 101)   # time points

# eigenfunctions
phi <- sapply(j, function(x){
  if (x %% 2 == 0) {
    sqrt(1/5)*sin(pi*t*x/5)
  } else {
    sqrt(1/5)*cos(pi*t*x/5)
  }
})

# parameters when generate class label
b <- matrix(ifelse(j <= 2, 1, (j-2)^(-3)), ncol=1)
beta.1 <- phi %*% b
beta.2 <- matrix(sqrt(3/10)*(t/5-1), ncol=1)

## generate curves
data <- list()
# # random sparsify => generate curve
# for (i in 1:n) {
#   # generate PC scores
#   xi <- sapply(j, function(x){ rnorm(1, 0, sqrt( x^(-1.5) )) })
# 
#   # random sparsify
#   num.obs <- sample(10:20, 1)
#   ind <- sort(sample(1:length(t), num.obs))
# 
#   # measurement error
#   eps <- rnorm(num.obs, 0, sqrt(0.1))
# 
#   # generate the curve
#   X <- xi %*% t(phi[ind, ]) + eps
# 
#   # generate class label
#   eps <- rnorm(1, 0, sqrt(0.1))   # model error
#   fx <- sin(pi*(X %*% beta.1[ind, ])/4)  # model 1
#   # fx <- exp((X %*% beta.1[ind, ])/2)-1   # model 2
#   # fx <- sin(pi*(X %*% beta.1[ind, ])/3) + exp((X %*% beta.2[ind, ])/3) - 1   # model 3
#   # fx <- atan(pi*(X %*% beta.1[ind, ])) + exp((X %*% beta.2[ind, ])/3) - 1    # model 4
#   y <- factor(ifelse(fx + eps < 0, 0, 1), levels=c(0, 1))
# 
#   data$id[[i]] <- rep(i, num.obs)
#   data$y[[i]] <- rep(y, num.obs)
#   data$Lt[[i]] <- t[ind]
#   data$Ly[[i]] <- X
# }

# generate curve => random sparsify
for (i in 1:n) {
  # generate PC scores
  xi <- sapply(j, function(x){ rnorm(1, 0, sqrt( x^(-1.5) )) })

  # measurement error
  eps <- rnorm(length(t), 0, sqrt(0.1))

  # generate the curve
  X <- xi %*% t(phi) + eps

  # generate class label
  eps <- rnorm(1, 0, sqrt(0.1))   # model error
  fx <- sin(pi*(X %*% beta.1)/4)  # model 1
  # fx <- exp((X %*% beta.1)/2)-1   # model 2
  # fx <- sin(pi*(X %*% beta.1)/3) + exp((X %*% beta.2)/3) - 1   # model 3
  # fx <- atan(pi*(X %*% beta.1)) + exp((X %*% beta.2)/3) - 1    # model 4
  y <- factor(ifelse(fx + eps < 0, 0, 1), levels=c(0, 1))

  # random sparsify
  num.obs <- sample(10:20, 1)
  ind <- sort(sample(1:length(t), num.obs))

  data$id[[i]] <- rep(i, num.obs)
  data$y[[i]] <- rep(y, num.obs)
  data$Lt[[i]] <- t[ind]
  data$Ly[[i]] <- X[ind]
}

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



### train, test split => train set에 없는 time point 포함시 FPCA할 때, error 발생
# range(test set) > range(train set) => 넘는 부분에 해당하는 id 제거
id <- 1:n
ind.train <- 1:200
ind.test <- NULL
range.train <- range(unlist(data$Lt[ind.train]))
range.test <- range(unlist(data$Lt[-ind.train]))
if (range.test[1] < range.train[1]) {   # min(test) < min(train) 인 경우
  ind <- sapply(data$Lt[-ind.train], function(x){ sum(x < range.train[1]) == 0 })
  ind.test <- id[ind]
}
if (range.test[2] > range.train[2]) {   # max(test) > max(train) 인 경우
  ind <- sapply(data$Lt[-ind.train], function(x){ sum(x > range.train[2]) == 0 })
  if (is.numeric(ind.test)) {
    ind.test <- intersect(ind.test,
                          id[ind])
  } else {
    ind.test <- id[ind]
  }
}
if (!is.numeric(ind.test)) {   # range(test) < range(train) 인 경우
  ind.test <- id[-ind.train]
}

# split to train, test set
y <- sapply(data$y, unique)
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
fpc.score <- fpca.fit$xiEst   # FPC scores

# train FPC scores
train.fpc <- data.frame(y = y.train, 
                        x = fpc.score)
colnames(train.fpc) <- c("y", paste("FPC", 1:k, sep=""))

# tune the hyperparmeters
tune.linear <- tune.svm(y ~ ., 
                        data = train.fpc, 
                        kernel = "linear",
                        cost = c(10^(-3:1), 2^(5:10)))
tune.radial <- tune.svm(y ~ ., 
                        data = train.fpc, 
                        kernel = "radial",
                        cost = c(10^(-3:1), 2^(5:10)), 
                        gamma = c(10^(-3:1), 2^(5:10)))

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
                       x = fpc.score.test)
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
acc.single <- sapply(pred, function(x){ mean(x == y.test) })
acc.single
# }
# y.pred
# colMeans(y.pred)
  
  
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
  
  
