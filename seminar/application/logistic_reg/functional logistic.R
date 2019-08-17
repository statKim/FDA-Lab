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
# 100개의 generated data로 accuracy 계산
k <- 5
accuracy <- c()
result <- list()
for (i in 1:100) {
  set.seed(777)
  y <- factor(c(rep(0, 100), rep(1, 100)), level=c(0, 1))
  data <- cbind(y, X.curves[[i]])
  ind <- sample(1:nrow(data), 100)   # train set index
  train <- data[ind, ]
  test <- data[-ind, ]
  
  train.fpc.data <- melt(cbind(id=1:nrow(train), train[, -1]), "id", variable.name="time", value.name="gene_exp") %>% 
    arrange(id, time)
  train.fpc.data$time <- time
  test.fpc.data <- melt(cbind(id=1:nrow(test), test[, -1]), "id", variable.name="time", value.name="gene_exp") %>% 
    arrange(id, time)
  test.fpc.data$time <- time
  
  # functional PC 계산
  fpca.train <- fpca.fit(train.fpc.data, iter=100, init_value=c(.1, .1, .1, .1, .1), num_knots=9, num_pc=k, grid_sep=1)
  fpca.test <- fpca.fit(test.fpc.data, iter=100, init_value=c(.1, .1, .1, .1, .1), num_knots=9, num_pc=k, grid_sep=1)
  
  # FPC로 trian, test set 만들기
  # train.new <- as.matrix(train[,-1]) %*% fpca.train$PCfunction[which(fpca.train$Timegrid %in% time),]  # PC score
  train.new <- (as.matrix(train[,-1]) - fpca.train$MeanFunction[which(fpca.train$Timegrid %in% time),]) %*% fpca.train$PCfunction[which(fpca.train$Timegrid %in% time),]  # PC score
  train.new <- cbind(train$y, as.data.frame(train.new))
  colnames(train.new) <- c("y", paste("PC", 1:k))
  # train.new$y <- factor(train.new$y, levels=c(0, 1))
  # train.new$y <- as.factor(train.new$y)
  # head(train.new)
  # str(train.new)
  
  # test.new <- as.matrix(test[,-1]) %*% fpca.test$PCfunction[which(fpca.test$Timegrid %in% time),] 
  test.new <- (as.matrix(test[,-1]) - fpca.test$MeanFunction[which(fpca.test$Timegrid %in% time),]) %*% fpca.test$PCfunction[which(fpca.test$Timegrid %in% time),]
  test.new <- cbind(test$y, as.data.frame(test.new))
  colnames(test.new) <- c("y", paste("PC", 1:k))
  # test.new$y <- factor(test.new$y, levels=c(0, 1))
  # test.new$y <- as.factor(test.new$y)
  
  # fit logistic regression
  fit.logit <- glm(y ~ ., data=train.new, family = "binomial")
  # summary(fit.logit)
  
  # predict
  pred <- predict(fit.logit, test.new[,-1], type="response")
  # pred <- factor(ifelse(pred > 0.5, 1, 0), levels=c(0, 1))
  pred <- factor(round(pred), levels=c(0, 1))
  acc <- mean(pred == test.new$y)
  
  if (acc < 0.5) {   # PC curve가 eigenvector때문에 부호가 바뀌는 경우
    pred <- predict(fit.logit, -test.new[,-1], type="response")
    pred <- factor(ifelse(pred > 0.5, 1, 0), levels=c(0, 1))
    acc <- mean(pred == test.new$y)
  }
  result[[i]] <- table(pred, true=test.new$y)
  # rbind(pred, true=test.new$y)
  accuracy[i] <- acc   # accuracy
  
  print( paste(i, "th data's Accuracy :", accuracy[i]) )
  # print( rbind(pred, true=test.new$y) )
}
})
# mean accuracy
mean(accuracy)


# PVE
x <- apply(test.new[, -1], MARGIN=2, FUN=var)
cumsum(x) / var(c( as.matrix(train[, -1]) ))

par(mfrow=c(2,1))
plot(train.new[,2], train.new[,3], col=as.integer(train.new$y), pch=as.integer(train.new$y), cex=1.5,
     xlab="PC1", ylab="PC2")
plot(test.new[,2], test.new[,3], col=as.integer(test.new$y), pch=as.integer(test.new$y), cex=1.5,
     xlab="PC2", ylab="PC3")



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
  for (i in 1:5){
    plot(fit$Timegrid, fit$PCfunction[,i], type="l", lwd=3, ylim=c(-.4, .2),
         xlab="minutes", ylab="Princ. comp.", main=paste("PC", i))
  }
  
  # loglikelihood 비교
  print( paste("Loglikelihood : ", round(fit$Loglik, 2), sep="") )
  
  # PVE
  print(fit$PVE)
}


x <- diag(fpca.train$D)
cumsum(x) / sum(x)
x/ sum(x)
