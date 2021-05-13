setwd("C:\\Users\\user\\Desktop\\KHS\\FDA-Lab\\Sparse FPCA\\Application\\bootstrap_aggregating")

library(tidyverse)
library(doParallel)   # parallel computing
library(fdapace)
library(e1071)
# library(class)
# library(MASS)
library(data.table)

# parallel computing setting
ncores <- detectCores() - 2
registerDoParallel(ncores)

packages <- c("fdapace","e1071")   # foreach에서 사용할 package 정의

# calculate mode
Mode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# load data
data <- read.csv("../spnbmd.csv", header=T)

# ind <- data %>%
#   group_by(idnum) %>%
#   summarise(n=n(),
#             max.age=max(age)) %>%
#   filter(n >= 2 & max.age < 18) %>%
#   dplyr::select(idnum)
ind <- data %>%
  group_by(idnum) %>%
  summarise(n=n()) %>%
  filter(n >= 2) %>%
  dplyr::select(idnum)
data <- data[which(data$idnum %in% ind$idnum), ]

# gender classification
data <- data %>% 
  mutate(y = factor(ifelse(sex == "fem", 1, 0), levels=c(0, 1))) %>% 
  dplyr::select(y, idnum, age, spnbmd)
length(unique(data$idnum))   # 280


acc.single <- numeric(30)
acc.majority <- numeric(30)
acc.oob <- numeric(30)

for (num in 1:30) {
  print(num)
  
  ### train, test split => train set에 없는 time point 포함시 에러
  # range(test set) > range(train set) => 넘는 부분에 해당하는 id 제거
  set.seed(100*num)
  id <- unique(data$idnum)
  id.train <- sample(id, floor(length(id) * 2/3))
  id.test <- NULL
  range.train <- range(data$age[which(data$idnum %in% id.train)])
  range.test <- range(data$age[-which(data$idnum %in% id.train)])
  if (range.test[1] < range.train[1]) {
    over.ind <- which(data$age[-which(data$idnum %in% id.train)] < range.train[1] )
    over.ind <- data$idnum[-which(data$idnum %in% id.train)][over.ind]
    id.test <- id[-which(id %in% c(id.train, over.ind))]
  }
  if (range.test[2] > range.train[2]) {
    over.ind <- which(data$age[-which(data$idnum %in% id.train)] > range.train[2] )
    over.ind <- data$idnum[-which(data$idnum %in% id.train)][over.ind]
    if (is.numeric(id.test)) {
      id.test <- intersect(id.test,
                           id[-which(unique(data$idnum) %in% c(id.train, over.ind))])
    } else {
      id.test <- id[-which(id %in% c(id.train, over.ind))]
    }
  }
  if (!is.numeric(id.test)) {
    id.test <- id[-which(id %in% id.train)]
  }
  
  # transform to FPCA input
  train <- data[which(data$idnum %in% id.train), ]
  test <- data[which(data$idnum %in% id.test), ]
  
  X.train <- MakeFPCAInputs(IDs = train$idnum,
                            tVec = train$age,
                            yVec = train$spnbmd)
  X.test <- MakeFPCAInputs(IDs = test$idnum,
                           tVec = test$age,
                           yVec = test$spnbmd)
  y.train <- unique(train[, c("idnum", "y")])[, "y"]
  y.test <- unique(test[, c("idnum", "y")])[, "y"]
  
  
  ### Single classifier
  fpca.fit <- FPCA(X.train$Ly, 
                   X.train$Lt, 
                   optns=list(dataType="Sparse",
                              methodXi="CE",
                              # methodSelectK="BIC",
                              FVEthreshold=0.99,
                              verbose=F))
  k <- fpca.fit$selectK   # optimal number of PCs
  fpc.score <- fpca.fit$xiEst
  
  # train FPC scores
  train.fpc <- data.frame(y=y.train, 
                          x=fpc.score)
  colnames(train.fpc) <- c("y", paste("FPC", 1:k, sep=""))
  
  # classifiers with bootstrapped FPC scores
  # fit.logit <- glm(y~., train.fpc, family=binomial)
  tuned.svm <- tune.svm(y ~ .,
                        data = train.fpc,
                        # kernel = "linear",
                        # kernel = "radial",
                        # kernel = "sigmoid",
                        kernel = "polynomial",
                        gamma = c(10^(-3:1), 2^(5:10)),
                        # coef0 = ,
                        degree = 2,
                        cost = c(10^(-3:1), 2^(5:10)) )
  fit.svm <- svm(y ~ .,
                 data = train.fpc,
                 # kernel = "linear",
                 # kernel = "radial",
                 # kernel = "sigmoid",
                 kernel = "polynomial",
                 gamma = tuned.svm$best.model$gamma,
                 # coef0 = tuned.svm$best.model$coef0,
                 degree = tuned.svm$best.model$degree,
                 cost = tuned.svm$best.model$cost)
  
  # test FPC scores
  fpc.score.test <- predict(fpca.fit, X.test$Ly, X.test$Lt, K=k, xiMethod="CE")
  test.fpc <- data.frame(y = y.test, 
                         x = fpc.score.test)
  colnames(test.fpc) <- c("y", paste("FPC", 1:k, sep=""))
  
  # predict
  # pred <- predict(fit.logit, test.fpc, type="response")
  # pred <- factor(ifelse(pred > 0.5, 1, 0), levels=c(0, 1))
  pred <- predict(fit.svm, test.fpc)
  
  # save the accuracy
  acc.single[num] <- mean(pred == y.test)
  
  
  ## Bagging
  # Bootstrap aggregating
  B <- 100
  N <- length(X.train$Ly)
  set.seed(100*num)
  boot.mat <- matrix(sample(1:N, N*B, replace=T), N, B)
  y.pred <- foreach(b=1:B, .packages=packages) %dopar% {
    boot.ind <- boot.mat[, b]
    
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
    
    # fit FPCA for bootstrapped data
    fpca.fit <- FPCA(X.train$Ly[boot.ind], 
                     X.train$Lt[boot.ind], 
                     optns=list(dataType="Sparse",
                                methodXi="CE",
                                # methodSelectK="BIC",
                                FVEthreshold=0.99,
                                verbose=F))
    k <- fpca.fit$selectK   # optimal number of PCs
    fpc.score <- fpca.fit$xiEst
    
    # train FPC scores
    train.fpc <- data.frame(y = y.train[boot.ind], 
                            x = fpc.score)
    colnames(train.fpc) <- c("y", paste("FPC", 1:k, sep=""))
    
    # classifiers with bootstrapped FPC scores
    # fit.logit <- glm(y~., train.fpc, family=binomial)
    # tuned.svm <- tune.svm(y ~ .,
    #                       data = train.fpc,
    #                       # kernel = "linear",
    #                       kernel = "radial",
    #                       gamma = c(10^(-3:1), 2^(5:10)),
    #                       cost = c(10^(-3:1), 2^(5:10)) )
    fit.svm <- svm(y ~ .,
                   data = train.fpc,
                   # kernel = "linear",
                   # kernel = "radial",
                   # kernel = "sigmoid",
                   kernel = "polynomial",
                   gamma = tuned.svm$best.model$gamma,
                   # coef0 = tuned.svm$best.model$coef0,
                   degree = tuned.svm$best.model$degree,
                   cost = tuned.svm$best.model$cost)
    
    # test FPC scores
    ind <- which(unlist(X.test$Lid) %in% id.test)
    fpc.score.test <- predict(fpca.fit, X.test$Ly[ind], X.test$Lt[ind], K=k, xiMethod="CE")
    test.fpc <- data.frame(y=y.test[ind], 
                           x=fpc.score.test)
    colnames(test.fpc) <- c("y", paste("FPC", 1:k, sep=""))
    
    # predict
    # pred <- predict(fit.logit, test.fpc, type="response")
    # pred <- factor(ifelse(pred > 0.5, 1, 0), levels=c(0, 1))
    pred <- predict(fit.svm, test.fpc)
    
    # OOB FPC scores
    oob.ind <- which(X.train$Lid %in% id.oob)
    fpc.score.oob <- predict(fpca.fit,
                             X.train$Ly[oob.ind],
                             X.train$Lt[oob.ind],
                             K=k,
                             xiMethod="CE")
    oob.fpc <- data.frame(y=y.train[oob.ind],
                          x=fpc.score.oob)
    colnames(oob.fpc) <- c("y", paste("FPC", 1:k, sep=""))
    
    # OOB error rate
    # oob.error <- mean(factor(ifelse(predict(fit.logit, oob.fpc, type="response") > 0.5, 1, 0), levels=c(0, 1)) != oob.fpc$y)
    oob.error <- mean(predict(fit.svm, oob.fpc) != oob.fpc$y)
    return( list(oob.error = oob.error,
                 boot = data.frame(id = id.test,
                                   y = y.test[ind],
                                   pred = pred,
                                   pred.oob = as.numeric(pred) * (1 - oob.error))) )
  }
  
  
  ## save the accuracy
  res <- as.data.frame(rbindlist(lapply(y.pred, function(x){ x$boot })))
  # majority voting
  pred <- res %>% 
    group_by(id, y) %>% 
    summarise(pred = Mode(pred))
  acc.majority[num] <- mean(pred$y == pred$pred)
  
  # oob error
  oob.acc <- 1 - sapply(y.pred, function(x){ x$oob.error })
  pred <- cbind(res[, 1:3],
                pred.oob = res$pred.oob / sum(oob.acc) ) %>% 
    group_by(id, y) %>% 
    summarise(pred2 = factor(ifelse(sum(pred.oob) > 1.5, 1, 0), levels=c(0, 1)))
  acc.oob[num] <- mean(pred$y == pred$pred2)
}

acc.single
acc.majority
acc.oob
c(mean(acc.single),
  mean(acc.majority),
  mean(acc.oob))

plot(1:length(acc.single), acc.single, 
     type="l", lwd=2, col=2, xlab="", ylab="Accuracy", ylim=c(0.45, 0.9))
lines(acc.majority, type="l", col=3, lwd=2)
lines(acc.oob, type="l", col=4, lwd=2)
legend("topleft", 
       paste(c("Single","Majority","OOB"), "(", 
             round(c(mean(acc.single),
                     mean(acc.majority),
                     mean(acc.oob)), 3), ")", sep=""), 
       lty=c(1,1,1), 
       lwd=c(2,2,2),
       col=c(2,3,4))

data.frame(Iteration = 1:30,
           Single = acc.single,
           Majority = acc.majority,
           OOB = acc.oob) %>% 
  reshape2::melt(id="Iteration", value.name="Accuracy") %>% 
  ggplot(aes(Iteration, Accuracy, group=variable, colour=variable)) +
  geom_line(size=1) +
  ylim(0.45, 0.9) +
  xlab("") +
  ylab("Accuracy") +
  theme_bw() +
  theme(legend.title = element_blank())
  
res.logit <- data.frame(seed = 100*1:30,
                        Single = acc.single,
                        Majority = acc.majority,
                        OOB = acc.oob)
res.svm <- data.frame(seed = 100*1:30,
                      Single = acc.single,
                      Majority = acc.majority,
                      OOB = acc.oob)
res.svm.tune <- data.frame(seed = 100*1:30,
                           Single = acc.single,
                           Majority = acc.majority,
                           OOB = acc.oob)
res.svm.linear <- data.frame(seed = 100*1:30,
                             Single = acc.single,
                             Majority = acc.majority,
                             OOB = acc.oob)
res.svm.sigmoid <- data.frame(seed = 100*1:30,
                              Single = acc.single,
                              Majority = acc.majority,
                              OOB = acc.oob)
res.svm.ploy <- data.frame(seed = 100*1:30,
                           Single = acc.single,
                           Majority = acc.majority,
                           OOB = acc.oob)
res.svm.ploy.2 <- data.frame(seed = 100*1:30,
                             Single = acc.single,
                             Majority = acc.majority,
                             OOB = acc.oob)
# result <- list(logit = res.logit,
#                svm = res.svm,
#                svm.tune = res.svm.tune,
#                svm.linear = res.svm.linear,
#                svm.sigmoid = res.svm.sigmoid,
#                svm.poly = res.svm.ploy)
result$svm.poly.degree2 <- res.svm.ploy.2
save(result, file="result.RData")

lapply(result, function(x){colMeans(x[,-1])})
