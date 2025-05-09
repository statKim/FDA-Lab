setwd("C:\\Users\\user\\Desktop\\KHS\\Thesis\\Principal Component Models for Sparse Functional Data\\Application\\logistic_reg")

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
source("../sparseFPCA.R")
library(dplyr)
library(reshape2)

system.time({
for (k in 1:5) {
  # 100개의 generated data로 accuracy 계산
  # k <- 5
  kn <- 7
  accuracy <- c()
  result <- list()
  for (i in 1:100) {
    set.seed(100)
    y <- factor(c(rep(0, 100), rep(1, 100)), level=c(0, 1))
    data <- cbind(y, X.curves[[i]])
    ind <- sample(1:nrow(data), 100)   # train set index
    train <- data[ind, ]
    test <- data[-ind, ]
    
    # functional PC 계산하기 위해 데이터 변환
    train.fpc.data <- melt(cbind(id=1:nrow(train), train[, -1]),
                           "id", variable.name="time", value.name="gene_exp") %>% 
                      arrange(id, time)
    train.fpc.data$time <- time
    test.fpc.data <- melt(cbind(id=1:nrow(test), test[, -1]),
                          "id", variable.name="time", value.name="gene_exp") %>% 
                     arrange(id, time)
    test.fpc.data$time <- time
    
    # functional PC fitting
    fpca.train <- fpca.fit(train.fpc.data, iter=100, init_value=c(.1, .1, .1, .1, .1), num_knots=kn, num_pc=k, grid_sep=1)
    fpca.test <- fpca.fit(test.fpc.data, iter=100, init_value=c(.1, .1, .1, .1, .1), num_knots=kn, num_pc=k, grid_sep=1)
    
    # FPC로 train, test set 만들기
    pc_func_train <- fpca.train$PCfunction[which(fpca.train$Timegrid %in% time), ]
    pc_func_test <- fpca.test$PCfunction[which(fpca.test$Timegrid %in% time), ]
    
    
    # 만약 PC function이 complex인 경우, 넘어가기
    if (is.complex(pc_func_train) | is.complex(pc_func_test)) {
      next
    }
    # PC function의 부호 같게 바꾸기(각 PC function의 첫 번째 값의 부호에 따라)
    if (is.null(dim(pc_func_train))) {   # k=1인 경우
      L2.same <- sum( (pc_func_train - pc_func_test)^2 )
      L2.diff <- sum( (pc_func_train + pc_func_test)^2 )
      # L2 norm 비교
      if (L2.same > L2.diff) {
        pc_func_test <- -pc_func_test
        fpca.test$PCfunction <- -fpca.test$PCfunction
      }
    } else {
      for (col in 1:k) {
        L2.same <- sum( (pc_func_train[, col] - pc_func_test[, col])^2 )
        L2.diff <- sum( (pc_func_train[, col] + pc_func_test[, col])^2 )
        # L2 norm 비교
        if (L2.same > L2.diff) {
          pc_func_test[, col] <- -pc_func_test[, col]
          fpca.test$PCfunction[, col] <- -fpca.test$PCfunction[, col]
        }
      }
    }

    # PC score 계산
    train.new <- as.matrix(train[,-1]) %*% pc_func_train  # PC score
    # train.new <- (as.matrix(train[,-1]) - fpca.train$MeanFunction[which(fpca.train$Timegrid %in% time),]) %*% fpca.train$PCfunction[which(fpca.train$Timegrid %in% time),]  # PC score
    train.new <- cbind(train$y, as.data.frame(train.new))
    colnames(train.new) <- c("y", paste("PC", 1:k))
    
    
    test.new <- as.matrix(test[,-1]) %*% pc_func_test
    # test.new <- (as.matrix(test[,-1]) - fpca.test$MeanFunction[which(fpca.test$Timegrid %in% time),]) %*% fpca.test$PCfunction[which(fpca.test$Timegrid %in% time),]
    test.new <- cbind(test$y, as.data.frame(test.new))
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
  save(list=c("result", "accuracy"), file=paste("dense_result/result_PC", k, ".RData", sep=""))
}
})

# error rate 계산
load("dense_result/result_PC5.RData")
# G1 error rate
mean( sapply(result, function(x){x$cross.table[2,1] / sum(x$cross.table[,1])}))
sd( sapply(result, function(x){x$cross.table[2,1] / sum(x$cross.table[,1])}) )
# G2 error rate
mean( sapply(result, function(x){x$cross.table[1,2] / sum(x$cross.table[,2])}) )
sd( sapply(result, function(x){x$cross.table[1,2] / sum(x$cross.table[,2])}) )
# overall
mean(1- accuracy)
sd(1- accuracy)


# misclassify된 curve를 구분지어서 그리기
par(mfrow=c(1,2))
for (j in c(0,1)) {
  ind_X <- which(test.new$y != pred & test.new$y == j)   # misclassify된 test set의 curve
  ind_O <- which(test.new$y == pred & test.new$y == j)
  sub <- test.fpc.data[which(test.fpc.data[,1]==ind_O[1]), ]
  if (j == 0) main.name <- "Group 1"
  else main.name <- "Group 2"
  plot(sub[,2], sub[,3],
       type="o",
       xlim=c(range(test.fpc.data[,2])[1], range(test.fpc.data[,2])[2]),
       ylim=c(range(test.fpc.data[,3])[1], range(test.fpc.data[,3])[2]),
       xlab="minutes",
       ylab="gene expression level",
       main=main.name)
  for (i in ind_O) {
    sub <- test.fpc.data[which(test.fpc.data[,1]==i), ]
    lines(sub[,2], sub[,3], type="l", col="gray", lwd=2)
    # if (test.new$y[i] == 0) {
    #   lines(sub[,2], sub[,3], type="l", col="gray", lwd=2)
    # } else {
    #   lines(sub[,2], sub[,3], type="l", lwd=1)
    # }
  }
  for (i in ind_X) {
    sub <- test.fpc.data[which(test.fpc.data[,1]==i), ]
    lines(sub[,2], sub[,3], type="o", lwd=2)
  }
  legend("topright", 
         c("Correct", "Misclassification"),
         col=c("gray","black"),
         lty=1)
}


# scatter plot of pairwise FPC scores(training set)
par(mfrow=c(1,2))
plot(train.new[,2], train.new[,3], col=as.integer(train.new$y), pch=as.integer(train.new$y), cex=1.5,
     xlab="First FPC score", ylab="Second FPC score")
abline(h=0, v=0)
legend("topright", 
       c("Group 1", "Group 2"),
       col=c("black","red"),
       pch=c(1,2))
plot(train.new[,3], train.new[,4], col=as.integer(train.new$y), pch=as.integer(train.new$y), cex=1.5,
     xlab="Second FPC score", ylab="Third FPC score")
abline(h=0, v=0)
legend("topright", 
       c("Group 1", "Group 2"),
       col=c("black","red"),
       pch=c(1,2))

# PVE
pve <- apply(test.new[, -1], 2, var)
cumsum(pve)

# train FPC plot
{
  par(mfrow=c(2,3))
  fit <- fpca.train
  
  # curves and mean function
  # par(mfrow=c(1,2))
  ind <- unique(train.fpc.data[,1])
  sub <- train.fpc.data[which(train.fpc.data[,1]==ind[1]), ]
  plot(sub[,2], sub[,3],
       type="o",
       xlim=c(range(train.fpc.data[,2])[1], range(train.fpc.data[,2])[2]),
       ylim=c(range(train.fpc.data[,3])[1], range(train.fpc.data[,3])[2]),
       xlab="minutes",
       ylab="gene expression level")
  for (i in ind) {
    sub <- train.fpc.data[which(train.fpc.data[,1]==i), ]
    lines(sub[,2], sub[,3], type="o")
  }
  lines(fit$Timegrid, fit$MeanFunction, col="red", lwd=3, type="l")
  
  # PC function
  for (i in 1:k){
    plot(fit$Timegrid, fit$PCfunction[,i], type="l", lwd=3, ylim=c(-.4, .2),
         xlab="minutes", ylab="Princ. comp.", main=paste("FPC", i, "(", round(pve[i]*100, 1), "%)"))
  }
}

# test FPC plot
{
  par(mfrow=c(2,3))
  fit <- fpca.test
  
  # curves and mean function
  # par(mfrow=c(1,2))
  ind <- unique(test.fpc.data[,1])
  sub <- test.fpc.data[which(test.fpc.data[,1]==ind[1]), ]
  plot(sub[,2], sub[,3],
       type="o",
       xlim=c(range(test.fpc.data[,2])[1], range(test.fpc.data[,2])[2]),
       ylim=c(range(test.fpc.data[,3])[1], range(test.fpc.data[,3])[2]),
       xlab="minutes",
       ylab="gene expression level")
  for (i in ind) {
    sub <- test.fpc.data[which(test.fpc.data[,1]==i), ]
    lines(sub[,2], sub[,3], type="o")
  }
  lines(fit$Timegrid, fit$MeanFunction, col="red", lwd=3, type="l")
  
  # PC function
  for (i in 1:k){
    plot(fit$Timegrid, fit$PCfunction[,i], type="l", lwd=3, ylim=c(-.4, .2),
         xlab="minutes", ylab="Princ. comp.", main=paste("PC", i))
  }
}


# predicted curve
fit <- fpca.train
par(mfrow=c(3,3))
for (i in 1:length(unique(train.fpc.data[, 1]))){
  ind <- which(fit$Timegrid %in% train.fpc.data[which(train.fpc.data[, 1]==unique(train.fpc.data[, 1])[i]), 2])
  time_point <- fit$Timegrid[ind]
  t_range <- range(ind)
  y_hat <- fit$MeanFunction[t_range[1]:t_range[2]] + fit$PCfunction[t_range[1]:t_range[2],]%*%fit$alpha[i,]
  
  y <- matrix(train.fpc.data[which(train.fpc.data[, 1]==unique(train.fpc.data[, 1])[i]), 3], ncol=1)
  
  plot(time_point, y,
       xlim=range(train.fpc.data[, 2]), ylim=range(train.fpc.data[, 3]),
       xlab="minutes", ylab="gene expression level",
       main=paste("predicted trajectory of curve", unique(train.fpc.data[, 1])[i]))
  points(fit$Timegrid[t_range[1]:t_range[2]], y_hat, col=3, type='l')
}



# PVE
pve <- apply(test.new[, -1], 2, var)
cumsum(pve)


x <- diag(fpca.train$D)
cumsum(x) / sum(x)
x/ sum(x)
