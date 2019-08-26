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


#######################################
### functional logistic regression
#######################################
source("../sparseFPCA.R")
library(dplyr)
library(reshape2)
library(e1071)

# time points 변수 만들기
time <- colnames(gene)[-1]
time <- as.integer(gsub(pattern="alpha|min",
                        replacement="",
                        x=time))
system.time({
# 100개의 generated data로 accuracy 계산
k <- 2
kn <- 8
# mean.acc.logit <- c()
# mean.acc.svm <- c()
# for (kn in 4:12) {
# print(paste("=====", kn, "====="))
  
accuracy.logit <- c()
accuracy.svm <- c()
# result <- list()
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
  
  ## 데이터 sparse하게 변환 => curve의 연속된 index random하게 구하기
  # train set
  samp <- sample(2:6, length(unique(train.fpc.data$id)), replace=T)   # curve의 obs 개수 선택
  ind.sparse <- c()
  for (j in unique(train.fpc.data$id)) {   # 각 curve 별로 개수별로 obs 뽑기
    ind <- which(train.fpc.data$id == j)
    m <- sort(sample(ind, samp[j]))
    ind.sparse <- append(ind.sparse,
                         m)
  }
  train.fpc.data <- train.fpc.data[ind.sparse, ]

  # test set
  samp <- sample(2:6, length(unique(test.fpc.data$id)), replace=T)   # curve의 obs 개수 선택
  ind.sparse <- c()
  for (j in unique(test.fpc.data$id)) {   # 각 curve 별로 개수별로 obs 뽑기
    ind <- which(test.fpc.data$id == j)
    m <- sort(sample(ind, samp[j]))
    ind.sparse <- append(ind.sparse,
                         m)
  }
  test.fpc.data <- test.fpc.data[ind.sparse, ]
  
  
  ### fpca package 사용한 경우 - PACE 방법으로 PC score 계산
  # column 순서 변경
  train.fpc.data <- as.matrix( train.fpc.data[, c(1,3,2)] )
  test.fpc.data <- as.matrix( test.fpc.data[, c(1,3,2)] )
  
  # parameter 정의
  ini.method <- "EM"        # optimization method
  basis.method <- "bs"      # spline basis - "ns"의 경우 현재 오류 발생
  sl.v <- rep(0.5, 10)      # Neuton-Raphson
  max.step <- 50            # EM iteration 횟수
  grid.l <- seq(0,1,0.01)   # EM에서는 사용 X
  grids <- seq(0,1,0.002)   # EM의 grid
  
  # fit fpca.mle
  fpca.train <- fpca.mle(train.fpc.data, kn, k, ini.method, basis.method, sl.v, max.step, grid.l, grids)
  fpca.test <- fpca.mle(test.fpc.data, kn, k, ini.method, basis.method, sl.v, max.step, grid.l, grids)
  
  # FPC로 train, test set 만들기
  pc_func_train <- fpca.train$eigenfunctions
  pc_func_test <- fpca.test$eigenfunctions
  for (col in 1:k) {   # 부호 비교해서 PC function 부호 통일
    # if (sign(pc_func_train[col, 1]) != sign(pc_func_test[col, 1])) {
    if ( (which.max(table(sign(pc_func_train[col, ]))) != which.max(table(sign(pc_func_test[col, ]))) &
         (table(sign(pc_func_train[col, ]))[which.max(table(sign(pc_func_train[col, ])))] > 334 &
          table(sign(pc_func_test[col, ]))[which.max(table(sign(pc_func_test[col, ])))] > 334) ) |
    # if (((table(sign(pc_func_train[col, ]))[1] / length(fpca.train$grid) > 2/3) &
    #      !(table(sign(pc_func_test[col, ]))[1] / length(fpca.test$grid) > 2/3)) |
    #     ((table(sign(pc_func_train[col, ]))[2] / length(fpca.train$grid) > 2/3) &
    #      !(table(sign(pc_func_test[col, ]))[2] / length(fpca.test$grid) > 2/3)) |
        sum( sign(pc_func_train[col, c(1, ncol(pc_func_train))]) == sign(pc_func_test[col, c(1, ncol(pc_func_test))]) ) == 0) {
    # if ( length( setdiff(which.peaks(pc_func_train[col,]), c(1, ncol(pc_func_train))) ) !=
    #         length( setdiff(which.peaks(pc_func_test[col,]), c(1, ncol(pc_func_test))) ) &
    #      length( setdiff(which.peaks(-pc_func_train[col,]), c(1, ncol(pc_func_train))) ) !=
    #      length( setdiff(which.peaks(-pc_func_test[col,]), c(1, ncol(pc_func_test))) ) ) {
       pc_func_test[col, ] <- -pc_func_test[col, ]
    }
  }
  fpca.train$eigenfunctions <- pc_func_train
  fpca.test$eigenfunctions <- pc_func_test
  
  # PC scores - 마지막 원소(grid_sep)는 timepoints 간격이 1, 0.1, 0.01, ...를 지정해줘야함
  train.new <- fpca.score(train.fpc.data, fpca.train$grid, fpca.train$fitted_mean, fpca.train$eigenvalues, 
                          fpca.train$eigenfunctions, fpca.train$error_var, k, 1)   
  train.new <- cbind(train$y, as.data.frame(train.new))
  colnames(train.new) <- c("y", paste("PC", 1:k))
  
  test.new <- fpca.score(test.fpc.data, fpca.test$grid, fpca.test$fitted_mean, fpca.test$eigenvalues, 
                         fpca.test$eigenfunctions, fpca.test$error_var, k, 1)
  test.new <- cbind(test$y, as.data.frame(test.new))
  colnames(test.new) <- c("y", paste("PC", 1:k))

  
  
  # ### 내가 만든 function으로 한 것
  # # functional PC fitting
  # fpca.train <- fpca.fit(train.fpc.data, iter=50, init_value=c(.1, .1, .1, .1, .1), num_knots=kn, num_pc=k, grid_sep=1)
  # fpca.test <- fpca.fit(test.fpc.data, iter=50, init_value=c(.1, .1, .1, .1, .1), num_knots=kn, num_pc=k, grid_sep=1)
  # 
  # # FPC로 train, test set 만들기
  # pc_func_train <- fpca.train$PCfunction[which(fpca.train$Timegrid %in% time),]
  # pc_func_test <- fpca.test$PCfunction[which(fpca.test$Timegrid %in% time),]
  # 
  # # 만약 PC function이 complex인 경우, 넘어가기
  # if (is.complex(pc_func_train) | is.complex(pc_func_test)) {
  #   next
  # }
  # # PC function의 부호 같게 바꾸기(각 PC function의 첫 번째 값의 부호에 따라)
  # if (is.null(dim(pc_func_train))) {   # k=1인 경우
  #   if (sign(pc_func_train[1]) != sign(pc_func_test[1])) {
  #     pc_func_test <- -pc_func_test
  #   }
  # } else {
  #   for (col in 1:k) {
  #     if (sign(pc_func_train[1, col]) != sign(pc_func_test[1, col])) {
  #       pc_func_test[, col] <- -pc_func_test[, col]
  #     }
  #   }
  # }
  # 
  # # train, test set 각각 PC score 계산
  # train.new <- matrix(NA, nrow(train), k)
  # for (j in unique(train.fpc.data$id)) {
  #   sub <- train.fpc.data[which(train.fpc.data$id == j), ]
  #   if (k == 1) {
  #     train.new[j, ] <- matrix(sub[, 3] - fpca.train$MeanFunction[which(fpca.train$Timegrid %in% sub$time)], nrow=1) %*% 
  #       pc_func_train[which(time %in% sub$time)]   # PC score
  #   } else {
  #     train.new[j, ] <- matrix(sub[, 3] - fpca.train$MeanFunction[which(fpca.train$Timegrid %in% sub$time)], nrow=1) %*% 
  #       pc_func_train[which(time %in% sub$time), ]   # PC score
  #   }
  # }
  # train.new <- cbind(train$y, as.data.frame(train.new))
  # colnames(train.new) <- c("y", paste("PC", 1:k))
  # 
  # test.new <- matrix(NA, nrow(test), k)
  # for (j in unique(test.fpc.data$id)) {
  #   sub <- test.fpc.data[which(test.fpc.data$id == j), ]
  #   if (k == 1) {
  #     test.new[j, ] <- matrix(sub[, 3] - fpca.train$MeanFunction[which(fpca.train$Timegrid %in% sub$time)], nrow=1) %*% 
  #       pc_func_test[which(time %in% sub$time)]   # PC score
  #   } else {
  #     test.new[j, ] <- matrix(sub[, 3] - fpca.train$MeanFunction[which(fpca.train$Timegrid %in% sub$time)], nrow=1) %*% 
  #       pc_func_test[which(time %in% sub$time), ]   # PC score
  #   }
  # }
  # test.new <- cbind(test$y, as.data.frame(test.new))
  # colnames(test.new) <- c("y", paste("PC", 1:k))
  #
  # FPC.plot(fpca.train, fpca.test)
  
  # fit logistic regression
  fit.logit <- glm(y ~ ., data=train.new, family = "binomial")
  
  # fit svm
  fit.svm <- svm(y ~ ., data=train.new)
  tune.svm <- tune(svm, train.x=train.new[, -1], train.y=train.new[, 1], kernel="radial", 
                   ranges=list(cost=10^(-1:2), gamma=c(.5,1,2)))
  fit.svm <- svm(y ~ ., data=train.new, kernel="radial",
                 cost=tune.svm$best.parameters[1], gamma=tune.svm$best.parameters[2])
  
  # predict
  test.data <- test.new[, -1]
  if (!is.data.frame(test.data)) {   # k=1인 경우(vector로 변환되기 때문)
    test.data <- data.frame(test.data)
    colnames(test.data) <- "PC 1"
  }
  pred.logit <- predict(fit.logit, test.data, type="response")
  pred.logit <- factor(ifelse(pred.logit > 0.5, 1, 0), levels=c(0, 1))
  # pred.logit <- factor(round(pred.logit), levels=c(0, 1))
  acc.logit <- mean(pred.logit == test.new$y)   # accuracy 계산
  
  pred.svm <- predict(fit.svm, test.data)
  acc.svm <- mean(pred.svm == test.new$y)   # accuracy 계산
  
  # output 저장
  accuracy.logit[i] <- acc.logit   # accuracy
  accuracy.svm[i] <- acc.svm
  # res <- list("accuracy"=acc.logit,
  #             "cross.table"=table(pred, true=test.new$y))
  # result[[i]] <- res
  
  print( paste(i, "th data's Accuracy(logit) :", acc.logit) )
  print( paste(i, "th data's Accuracy(svm) :", acc.svm) )
  
  par(mfrow=c(2, k+1))
  # train set curves and mean function, PC function
  fpc.plot(train.fpc.data, fpca.train, xlab="minute", ylab="gene expression level")
  # test set curves and mean function, PC function
  fpc.plot(test.fpc.data, fpca.test, xlab="minute", ylab="gene expression level")
  mtext(paste(i, "th :", acc.logit),
               side = 3, # which margin to place text. 1=bottom, 2=left, 3=top, 4=right
               line = 3, # to indicate the line in the margin starting with 0 and moving out
               adj = 1, # adj=0 for left/bottom alignment or adj=1 for top/right alignment
               cex = 2, # font size
               outer = F)
}
# mean accuracy
print( paste("Mean accuracy(logit) :", mean(accuracy.logit)) )
print( paste("Mean accuracy(svm) :", mean(accuracy.svm)) )

# mean.acc.logit <- append(mean.acc.logit, mean(accuracy.logit))
# mean.acc.svm <- append(mean.acc.svm, mean(accuracy.svm))
# }
# RData로 저장
# save(list=c("result", "accuracy"), file=paste("sparse_result/result_PC", k, ".RData", sep=""))
})

acc.knots <- setNames(data.frame(rbind(mean.acc.logit, mean.acc.svm), row.names = c("logit","svm")),4:12)
save(acc.knots, file="sparse_result/accuracy_knots.RData")

# find local maxima index from https://gist.github.com/jamiefolson/5831746
which.peaks <- function(x,partial=TRUE,decreasing=FALSE){
  if (decreasing){
    if (partial){
      which(diff(c(TRUE,diff(x)<=0,FALSE))>0)
    }else {
      which(diff(diff(x)<=0)>0)
    }    
  }else {
    if (partial){
      which(diff(c(TRUE,diff(x)>=0,FALSE))<0)
    }else {
      which(diff(diff(x)>=0)<0)
    }
  }
}

# predicted curve
pred.plot(train.fpc.data, fpca.train, xlab="minutes", ylab="gene expression level")
pred.plot(test.fpc.data, fpca.test, xlab="minutes", ylab="gene expression level")

# error rate 계산
load("result_PC5.RData")
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
FPC.plot <- function(fpca.train, fpca.test){
  par(mfrow=c(2,k+1))
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
    plot(time, pc_func_train[,i], type="l", lwd=3, ylim=c(-.4, .2),
         xlab="minutes", ylab="Princ. comp.", main=paste("FPC", i))
  }
 
  # test FPC plot
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
    plot(time, pc_func_test[,i], type="l", lwd=3, ylim=c(-.4, .2),
         xlab="minutes", ylab="Princ. comp.", main=paste("FPC", i))
  }
}










##############################
### bone mineral density data
##############################
setwd("C:\\Users\\user\\Desktop\\KHS\\Thesis\\Principal Component Models for Sparse Functional Data\\Application")

# White이면서 obs 개수가 1개인 경우 제외하면 딱 48개 curve가 됨!!
data <- read.csv("spnbmd.csv", stringsAsFactors=F)
data <- data[which(data$ethnic %in% c("Asian","Black")), ]
dim(data)
library(dplyr)
ind <- data %>% 
  group_by(idnum) %>% 
  summarise(n=n())
ind <- ind[which(ind$n != 1), "idnum"]

data <- data[which(data$idnum %in% t(ind)), ]
dim(data)
data <- data[,c(2,1,3,5)]   # 필요없는 변수 제거
length(unique(data$idnum))   # 90
head(data)
# data$sex <- factor( ifelse(data$sex == "fem", 1, 0), levels=c(0,1) )
data$ethnic <- factor( ifelse(data$ethnic == "Asian", 1, 0), levels=c(0,1) )

source("sparseFPCA.R")
kn <- 4
k <- 5

# training set
set.seed(100)
samp <- sample(unique(data$idnum), round(length(unique(data$idnum))/2))   # curve의 obs 개수 선택
train.fpc.data <- data[which(data$id %in% samp), ]
test.fpc.data <- data[-which(data$id %in% samp), ]

# functional PC fitting
fpca.train <- fpca.fit(train.fpc.data[, -1], iter=100, init_value=c(.1, .1, .1, .1, .1), num_knots=kn, num_pc=k, grid_sep=.1)
fpca.test <- fpca.fit(test.fpc.data[, -1], iter=100, init_value=c(.1, .1, .1, .1, .1), num_knots=kn, num_pc=k, grid_sep=.1)

# FPC로 train, test set 만들기
pc_func_train <- fpca.train$PCfunction
pc_func_test <- fpca.test$PCfunction

# 만약 PC function이 complex인 경우, 넘어가기
if (is.complex(pc_func_train) | is.complex(pc_func_test)) {
  next
}
# PC function의 부호 같게 바꾸기(각 PC function의 첫 번째 값의 부호에 따라)
if (is.null(dim(pc_func_train))) {   # k=1인 경우
  if (sign(pc_func_train[1]) != sign(pc_func_test[1])) {
    pc_func_test <- -pc_func_test
  }
} else {
  for (col in 1:k) {
    if (sign(pc_func_train[1, col]) != sign(pc_func_test[1, col])) {
      pc_func_test[, col] <- -pc_func_test[, col]
    }
  }
}

# PC score 계산
for (j in unique(train.fpc.data$id)) {
  sub <- train.fpc.data[which(train.fpc.data$id == j), ]
  pc.score <- matrix(sub[, 3] - fpca.train$MeanFunction[which(fpca.train$Timegrid %in% sub$age)], nrow=1) %*% 
              pc_func_train[which(fpca.train$Timegrid %in% sub$age), ]   # PC score
  if (j == unique(train.fpc.data$id)[1]) {
    train.new <- pc.score
  } else {
    train.new <- rbind(train.new, pc.score)
  }
}
train.new <- cbind(unique(train.fpc.data[, 1:2])[, 1], as.data.frame(train.new))
colnames(train.new) <- c("y", paste("PC", 1:k))

for (j in unique(test.fpc.data$id)) {
  sub <- test.fpc.data[which(test.fpc.data$id == j), ]
  pc.score <- matrix(sub[, 3] - fpca.test$MeanFunction[which(fpca.test$Timegrid %in% sub$age)], nrow=1) %*% 
    pc_func_test[which(fpca.test$Timegrid %in% sub$age), ]   # PC score
  if (j == unique(test.fpc.data$id)[1]) {
    test.new <- pc.score
  } else {
    test.new <- rbind(test.new, pc.score)
  }
}
test.new <- cbind(unique(test.fpc.data[, 1:2])[, 1], as.data.frame(test.new))
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
pred <- factor(ifelse(pred > 0.5, 1, 0), levels=c(0, 1))
acc <- mean(pred == train.new$y)   # accuracy 계산
acc


{ # train FPC plot
  par(mfrow=c(2,k+1))
  # curves and mean function
  ind <- unique(train.fpc.data[,2])
  sub <- train.fpc.data[which(train.fpc.data[,2]==ind[1]), ]
  plot(sub[,3], sub[,4],
       type="o",
       xlim=c(range(train.fpc.data[,3])[1], range(train.fpc.data[,3])[2]),
       ylim=c(range(train.fpc.data[,4])[1], range(train.fpc.data[,4])[2]),
       xlab="minutes",
       ylab="gene expression level")
  for (i in ind) {
    sub <- train.fpc.data[which(train.fpc.data[,2]==i), ]
    lines(sub[,3], sub[,4], type="o")
  }
  lines(fpca.train$Timegrid, fpca.train$MeanFunction, col="red", lwd=3, type="l")
  
  # PC function
  for (i in 1:k){
    plot(fpca.train$Timegrid, pc_func_train[,i], type="l", lwd=3, ylim=c(-.2, .2),
         xlab="minutes", ylab="Princ. comp.", main=paste("FPC", i))
  }
}

{ # test FPC plot
  # curves and mean function
  ind <- unique(test.fpc.data[,2])
  sub <- test.fpc.data[which(test.fpc.data[,2]==ind[1]), ]
  plot(sub[,3], sub[,4],
       type="o",
       xlim=c(range(test.fpc.data[,3])[1], range(test.fpc.data[,3])[2]),
       ylim=c(range(test.fpc.data[,4])[1], range(test.fpc.data[,4])[2]),
       xlab="minutes",
       ylab="gene expression level")
  for (i in ind) {
    sub <- test.fpc.data[which(test.fpc.data[,2]==i), ]
    lines(sub[,3], sub[,4], type="o")
  }
  lines(fpca.test$Timegrid, fpca.test$MeanFunction, col="red", lwd=3, type="l")
  
  # PC function
  for (i in 1:k){
    plot(fpca.test$Timegrid, pc_func_test[,i], type="l", lwd=3, ylim=c(-.2, .2),
         xlab="minutes", ylab="Princ. comp.", main=paste("FPC", i))
  }
}
