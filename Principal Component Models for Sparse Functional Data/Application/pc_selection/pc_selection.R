setwd("C:\\Users\\user\\Desktop\\KHS\\Thesis\\Principal Component Models for Sparse Functional Data\\Application\\pc_selection")

#####################
### Data generating
#####################
library(fda.usc)
library(dplyr)
library(ggplot2)
library(data.table)

set.seed(100)
cell <- read.delim("../logistic_reg/mbc_9_12_3273__CDCDATA.txt")
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

# time points 변수 만들기
time <- colnames(gene)[-1]
time <- as.integer(gsub(pattern="alpha|min",
                        replacement="",
                        x=time))


###########################################
### PC selection (PVE, LOOCV)
###########################################
source("../sparseFPCA.R")
library(dplyr)
library(reshape2)
library(e1071)   # SVM fit하기 위함
library(doParallel)   # parallel computing
# 기존 fpca package의 함수를 재정의
path <- "C:/Users/user/Desktop/KHS/fpca/R/"
source( paste(path, "fpca_func.R", sep="") )
source( paste(path, "fpca_pred.R", sep="") )
source( paste(path, "fpca_score.R", sep="") )
source( paste(path, "functions_EM.R", sep="") )
source( paste(path, "functions_GenData.R", sep="") )
source( paste(path, "functions_LocLin.R", sep="") )
source( paste(path, "functions_Optimization.R", sep="") )


# parallel computing setting
ncores <- detectCores() - 3
registerDoParallel(ncores)

##############################################
## 1. PVE (Proportion of Variance Explained)
##############################################
packages <- c("splines","sm","dplyr","reshape2","e1071")   # foreach에서 사용할 package 정의
result <- list()   # obs 개수 별로 output 저장
start.time <- Sys.time()   # 연산속도 체크
for (num.obs in 2:18) {
  res <- foreach(i = 1:100, .combine="rbind", .packages=packages) %dopar% {
    print(i)
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
    samp <- rep(num.obs, length(unique(train.fpc.data$id)))   # curve의 obs 개수 선택(모두 같은 개수로 sparse하게 만듬)
    # train set
    # samp <- sample(2:6, length(unique(train.fpc.data$id)), replace=T)   # curve의 obs 개수 선택
    ind.sparse <- c()
    for (j in unique(train.fpc.data$id)) {   # 각 curve 별로 개수별로 obs 뽑기
      ind <- which(train.fpc.data$id == j)
      m <- sort(sample(ind, samp[j]))
      ind.sparse <- append(ind.sparse,
                           m)
    }
    train.fpc.data <- train.fpc.data[ind.sparse, ]
    
    # test set
    # samp <- sample(2:6, length(unique(test.fpc.data$id)), replace=T)   # curve의 obs 개수 선택
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
    max.step <- 100           # EM iteration 횟수
    grid.l <- seq(0,1,0.01)   # EM에서는 사용 X
    grids <- seq(0,1,0.002)   # EM의 grid
    kn <- 8                   # knots 개수 - 1
    
    # initial value
    PVE <- 0
    k <- 1
    
    # PVE > 0.95 일 때의 PC 개수 사용
    while (PVE < 0.95) {
      k <- k + 1
      
      if (k > kn) {
        break
      }
      
      # fit fpca.mle
      fpca.train <- fpca.mle(train.fpc.data, kn, k, ini.method, basis.method, sl.v, max.step, grid.l, grids)
      
      # PVE > 0.95인 경우만 고르기
      PVE[k] <- fpca.train$PVE[2, k]
      
      # PC 1개를 추가할 때, cumulative PVE가 0.01 이하만 증가하면 개수 추가 X
      if (PVE[k] - PVE[k-1] < 0.01) {
        k <- k - 1
        break
      }
    }
    
    # select된 k를 사용해서 fpca fit
    fpca.train <- fpca.mle(train.fpc.data, kn, k, ini.method, basis.method, sl.v, max.step, grid.l, grids)
    
    # PC scores - 마지막 원소(grid_sep)는 timepoints 간격이 1, 0.1, 0.01, ...를 지정해줘야함
    train.new <- tryCatch(fpca.score(train.fpc.data, fpca.train$grid, fpca.train$fitted_mean, fpca.train$eigenvalues, 
                                     fpca.train$eigenfunctions, fpca.train$error_var, k, 1),
                          error = function(e) {
                            print(paste(i, "th data's train set fpca error"))
                            return(TRUE)
                          })
    test.new <- tryCatch(fpca.score(test.fpc.data, fpca.train$grid, fpca.train$fitted_mean, fpca.train$eigenvalues, 
                           fpca.train$eigenfunctions, fpca.train$error_var, k, 1),
                         error = function(e) {
                           print(paste(i, "th data's train set fpca error"))
                           return(TRUE)
                         })

    
    # score 계산 중 에러가 발생하지 않을 때, fpca.train의 Newton-iteration이 수렴하지 못했을 때(error_var>>, cv_score 계산 X)
    if (!is.numeric(train.new) & !is.numeric(test.new) & fpca.train$converge < 1) {
      train.new <- cbind(train$y, as.data.frame(train.new))
      colnames(train.new) <- c("y", paste("PC", 1:k))
      test.new <- cbind(test$y, as.data.frame(test.new))
      colnames(test.new) <- c("y", paste("PC", 1:k))

      # fit logistic regression
      fit.logit <- glm(y ~ ., data=train.new, family = "binomial")

      # fit svm
      fit.svm.linear <- tune(svm, y ~ ., data=train.new, kernel="linear",
                             ranges=list(cost=10^(-1:2), gamma=c(.5,1,2)))$best.model
      fit.svm.radial <- tune(svm, y ~ ., data=train.new, kernel="radial",
                             ranges=list(cost=10^(-1:2), gamma=c(.5,1,2)))$best.model
      fit.svm.sigmoid <- tune(svm, y ~ ., data=train.new, kernel="sigmoid",
                              ranges=list(cost=10^(-1:2), gamma=c(.5,1,2)))$best.model
      fit.svm.poly <- tune(svm, y ~ ., data=train.new, kernel="poly",
                           ranges=list(cost=10^(-1:2), gamma=c(.5,1,2)))$best.model

      # predict
      test.data <- test.new[, -1]
      if (!is.data.frame(test.data)) {   # k=1인 경우(vector로 변환되기 때문)
        test.data <- data.frame(test.data)
        colnames(test.data) <- "PC 1"
      }

      # logit predict
      pred.logit <- predict(fit.logit, test.data, type="response")
      pred.logit <- factor(ifelse(pred.logit > 0.5, 1, 0), levels=c(0, 1))
      acc.logit <- mean(pred.logit == test.new$y)   # accuracy 계산

      # svm predict
      pred.svm.linear <- predict(fit.svm.linear, test.data)
      acc.svm.linear <- mean(pred.svm.linear == test.new$y)   # accuracy 계산
      pred.svm.radial <- predict(fit.svm.radial, test.data)
      acc.svm.radial <- mean(pred.svm.radial == test.new$y)   # accuracy 계산
      pred.svm.sigmoid <- predict(fit.svm.sigmoid, test.data)
      acc.svm.sigmoid <- mean(pred.svm.sigmoid == test.new$y)   # accuracy 계산
      pred.svm.poly <- predict(fit.svm.poly, test.data)
      acc.svm.poly <- mean(pred.svm.poly == test.new$y)   # accuracy 계산

      # output 저장
      return( data.frame(logit=acc.logit, 
                         svm.linear=acc.svm.linear,
                         svm.radial=acc.svm.radial,
                         svm.sigmoid=acc.svm.sigmoid,
                         svm.poly=acc.svm.poly,
                         K=k,
                         PVE=round(PVE[k], 3)) )
    } else {
      return( data.frame(logit=NA, 
                         svm.linear=NA,
                         svm.radial=NA,
                         svm.sigmoid=NA,
                         svm.poly=NA,
                         K=k,
                         PVE=round(PVE[k], 3)) )
    }
  }

  result[[paste(num.obs, "obs", sep="")]] <- res
}  
end.time <- Sys.time()
end.time - start.time   # 연산 속도 측정


save(result, file="results/pve_result.RData")

load("results/pve_result.RData")


sapply(result, function(x){colMeans(x, na.rm=T)})
lapply(result, function(x){table(x[, 6])})
sapply(result, function(x){sum(is.na(x[, 1]))})

library(xtable)
xtable(t(sapply(result, function(x){colMeans(x, na.rm=T)})), digits=c(rep(3, 6), 2,2))


###########################################
## 2. Leave-one-curve-out cross validation
###########################################
# parallel computing
library(doParallel)
ncores <- detectCores() - 4
registerDoParallel(ncores)

i <- 1
num.obs <- 10

# result <- list()
# result.sse <- list()
# result.kl <- list()
start.time <- Sys.time()   # 연산속도 체크
for (num.obs in 2:18) {
  # for (i in 1:100) {
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
    samp <- rep(num.obs, length(unique(train.fpc.data$id)))   # curve의 obs 개수 선택(모두 같은 개수로 sparse하게 만듬)
    # train set
    # samp <- sample(2:10, length(unique(train.fpc.data$id)), replace=T)   # curve의 obs 개수 선택
    ind.sparse <- c()
    for (j in unique(train.fpc.data$id)) {   # 각 curve 별로 개수별로 obs 뽑기
      ind <- which(train.fpc.data$id == j)
      m <- sort(sample(ind, samp[j]))
      ind.sparse <- append(ind.sparse,
                           m)
    }
    train.fpc.data <- train.fpc.data[ind.sparse, ]
    
    # test set
    # samp <- sample(2:10, length(unique(test.fpc.data$id)), replace=T)   # curve의 obs 개수 선택
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
    max.step <- 100           # EM iteration 횟수
    grid.l <- seq(0,1,0.01)   # EM에서는 사용 X
    grids <- seq(0,1,0.002)   # EM의 grid
    kn <- 8                   # knots 개수 - 1
    
    
    cv.loocv <- c()
    cv.fpca <- c()
    for (k in 2:8) {
      res <- foreach(i = 1:100, .combine="c", .packages=packages) %dopar% {
        ind <- which(train.fpc.data[, 1] == i)
        train.set <- train.fpc.data[-ind, ]
        valid.set <- train.fpc.data[ind, ]
        
        # fit fpca.mle
        fpca.train <- tryCatch(fpca.mle(train.set, kn, k, ini.method, basis.method, sl.v, max.step, grid.l, grids),
                               error = function(e) {
                                 print(paste(i, "th data's train set fpca error"))
                                 return(TRUE)
                               }) 
        
        # FPC score
        train.new <- tryCatch(fpca.score(valid.set, fpca.train$grid, fpca.train$fitted_mean, fpca.train$eigenvalues, 
                                         fpca.train$eigenfunctions, fpca.train$error_var, k, 1),
                              error = function(e) {
                                print(paste(i, "th data's train set fpca error"))
                                return(TRUE)
                              }) 
        
        if ( is.list(fpca.train) ) {   # fpca.mle()가 error가 나오지 않은 경우
          if (is.numeric(train.new) & fpca.train$converge < 1) {
            # Karhunen-Leave expansion
            curve.hat <- fpca.pred(train.new, fpca.train$fitted_mean, fpca.train$eigenfunctions)
            
            time.point <- valid.set[, 3]
            ind <- sapply(time.point, function(x){ceiling(median(which( round(fpca.train$grid) %in% x )))} )
            
            return( mean( (valid.set[, 2] - curve.hat[ind])^2 ) )
          } else {
            return(NA)
          }
        } else {
          return(NA)
        }
      }
      cv.loocv[k] <- mean(res, na.rm=T)
      
      fpca.train <- tryCatch(fpca.mle(train.fpc.data, kn, k, ini.method, basis.method, sl.v, max.step, grid.l, grids),
                             error = function(e) {
                               print(paste(i, "th data's train set fpca error"))
                               return(TRUE)
                             })
      if ( is.list(fpca.train) ) {
        if (fpca.train$converge < 1) {
          cv.fpca[k] <- fpca.train$cv_scores
        } else {
          cv.fpca[k] <- NA
        }
      } else {
        cv.fpca[k] <- NA
      }
    }
    
  
    # # selected k
    # k <- which.min(cv.loocv)
    # 
    # # select된 k를 사용해서 fpca fit
    # fpca.train <- fpca.mle(train.fpc.data, kn, k, ini.method, basis.method, sl.v, max.step, grid.l, grids)
    # 
    # # PC scores - 마지막 원소(grid_sep)는 timepoints 간격이 1, 0.1, 0.01, ...를 지정해줘야함
    # train.new <- tryCatch(fpca.score(train.fpc.data, fpca.train$grid, fpca.train$fitted_mean, fpca.train$eigenvalues, 
    #                                  fpca.train$eigenfunctions, fpca.train$error_var, k, 1),
    #                       error = function(e) {
    #                         print(paste(i, "th data's train set fpca error"))
    #                         return(TRUE)
    #                       })
    # test.new <- tryCatch(fpca.score(test.fpc.data, fpca.train$grid, fpca.train$fitted_mean, fpca.train$eigenvalues, 
    #                                 fpca.train$eigenfunctions, fpca.train$error_var, k, 1),
    #                      error = function(e) {
    #                        print(paste(i, "th data's train set fpca error"))
    #                        return(TRUE)
    #                      })
    # 
    # 
    # # score 계산 중 에러가 발생하지 않을 때, fpca.train의 Newton-iteration이 수렴하지 못했을 때(error_var>>, cv_score 계산 X)
    # if (!is.numeric(train.new) & !is.numeric(test.new) & fpca.train$converge < 1) {
    #   train.new <- cbind(train$y, as.data.frame(train.new))
    #   colnames(train.new) <- c("y", paste("PC", 1:k))
    #   test.new <- cbind(test$y, as.data.frame(test.new))
    #   colnames(test.new) <- c("y", paste("PC", 1:k))
    #   
    #   # fit logistic regression
    #   fit.logit <- glm(y ~ ., data=train.new, family = "binomial")
    #   
    #   # fit svm
    #   fit.svm.linear <- tune(svm, y ~ ., data=train.new, kernel="linear",
    #                          ranges=list(cost=10^(-1:2), gamma=c(.5,1,2)))$best.model
    #   fit.svm.radial <- tune(svm, y ~ ., data=train.new, kernel="radial",
    #                          ranges=list(cost=10^(-1:2), gamma=c(.5,1,2)))$best.model
    #   fit.svm.sigmoid <- tune(svm, y ~ ., data=train.new, kernel="sigmoid",
    #                           ranges=list(cost=10^(-1:2), gamma=c(.5,1,2)))$best.model
    #   fit.svm.poly <- tune(svm, y ~ ., data=train.new, kernel="poly",
    #                        ranges=list(cost=10^(-1:2), gamma=c(.5,1,2)))$best.model
    #   
    #   # predict
    #   test.data <- test.new[, -1]
    #   if (!is.data.frame(test.data)) {   # k=1인 경우(vector로 변환되기 때문)
    #     test.data <- data.frame(test.data)
    #     colnames(test.data) <- "PC 1"
    #   }
    #   
    #   # logit predict
    #   pred.logit <- predict(fit.logit, test.data, type="response")
    #   pred.logit <- factor(ifelse(pred.logit > 0.5, 1, 0), levels=c(0, 1))
    #   acc.logit <- mean(pred.logit == test.new$y)   # accuracy 계산
    #   
    #   # svm predict
    #   pred.svm.linear <- predict(fit.svm.linear, test.data)
    #   acc.svm.linear <- mean(pred.svm.linear == test.new$y)   # accuracy 계산
    #   pred.svm.radial <- predict(fit.svm.radial, test.data)
    #   acc.svm.radial <- mean(pred.svm.radial == test.new$y)   # accuracy 계산
    #   pred.svm.sigmoid <- predict(fit.svm.sigmoid, test.data)
    #   acc.svm.sigmoid <- mean(pred.svm.sigmoid == test.new$y)   # accuracy 계산
    #   pred.svm.poly <- predict(fit.svm.poly, test.data)
    #   acc.svm.poly <- mean(pred.svm.poly == test.new$y)   # accuracy 계산
    #   
    #   # output 저장
    #   output <- data.frame(logit=acc.logit, 
    #                        svm.linear=acc.svm.linear,
    #                        svm.radial=acc.svm.radial,
    #                        svm.sigmoid=acc.svm.sigmoid,
    #                        svm.poly=acc.svm.poly,
    #                        K=k,
    #                        PVE=round(fpca.train$PVE[2, k], 3))
    # } else {
    #   output <- data.frame(logit=NA, 
    #                        svm.linear=NA,
    #                        svm.radial=NA,
    #                        svm.sigmoid=NA,
    #                        svm.poly=NA,
    #                        K=k,
    #                        PVE=round(fpca.train$PVE[2, k], 3))
    # }
    
    # if (i > 1) {
    if (num.obs > 2) {
      cv.sse <- rbind(cv.sse, cv.loocv)
      cv.kl <- rbind(cv.kl, cv.fpca)
      # res <- rbind(res, output)
    } else {
      cv.sse <- cv.loocv
      cv.kl <- cv.fpca
      # res <- output
    }
  # }
  
#   result[[paste(num.obs, "obs", sep="")]] <- res
#   result.sse[[paste(num.obs, "obs", sep="")]] <- cv.sse
#   result.kl[[paste(num.obs, "obs", sep="")]] <- cv.kl
}
end.time <- Sys.time()
end.time - start.time   # 연산 속도 측정

rownames(cv.sse) <- 2:18
rownames(cv.kl) <- 2:18

save(list=c("cv.sse", "cv.kl"), file="results/cv_result.RData")

par(mfrow=c(3,2))
for (i in 2:18) {
  plot(2:8, cv.sse[i-1, -1], type="o", xlab="No. of FPC", ylab="MSE", main=paste(i, "obs"))
  points(which.min(cv.sse[i-1, -1])+1, 
         cv.sse[i-1, -1][which.min(cv.sse[i-1, -1])], 
         col="red", pch=19, cex=1.5)
  plot(2:8, cv.kl[i-1, -1], type="o", xlab="No. of FPC", ylab="Kullback–Leibler Loss", main=paste(i, "obs"))
  points(which.min(cv.kl[i-1, -1])+1,
         cv.kl[i-1, -1][which.min(cv.kl[i-1, -1])],
         col="red", pch=19, cex=1.5)
}


