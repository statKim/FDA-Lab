setwd("C:\\Users\\user\\Desktop\\KHS\\Thesis\\Principal Component Models for Sparse Functional Data\\Application")

## EM algorithm
# White이면서 obs 개수가 1개인 경우 제외하면 딱 48개 curve가 됨!!
data <- read.csv("spnbmd.csv", stringsAsFactors=F)
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
data <- data[, c(1,3,2)]   # column 순서 변경

### fpca package 이용해서 functional PC 계산
library(fpca)
# source("fpca/fpca_func.R")
# source("fpca/fpca_pred.R")
# source("fpca/fpca_score.R")
# source("fpca/functions_EM.R")
# source("fpca/functions_GenData.R")
# source("fpca/functions_LocLin.R")
# source("fpca/functions_Optimization.R")
# fpca 내의 잘못 정의된 함수 수정
fpca.score <- function(data.m, grids.u, muhat, eigenvals, eigenfuncs, sig2hat, K){
  ##estimated conditional principal component scores (BLUPs): \hat{E(\xi_k|Y)}, k=1,...,K
  ##Name:FPcScore
  ##para:
  ##     data.m -- data matrix; same as input for fpca.mle
  ##     grids.u -- grid of time points used in evaluating the mean and eigenfunctions (on the original scale); (returned by fpca. mle)
  ##     muhat,eigenvals, eigenfuncs, sig2hat -- (estimated) mean, eigenvalues, eigenfunctions and noise variance; (returned by fpca.mle)
  ##     K -- number of eigenfunctions used in the model, i.e., (estimated) dimension of the process
  ##return: first K conditional PC scores (the BLUP estimates): n by K
  temp <- table(data.m[,1])
  n <- length(temp)             # number of curves;
  m.l <- as.vector(temp)        # m.l -- number of time points per curve
  result <- matrix(0, n, K)       # First K FPC scores for each subject
  
  N <- length(grids.u)        # number of time points on the grid
  evalmat <- diag(eigenvals[1:K])  # diagonal matrix of the first K (estimated) eigenvalues
  current <- 0  # current index
  eigenfuncs.u <- t(eigenfuncs)   # dimmension: grid_length by K
  
  # 정확한 grid의 index 구하기 위해 실제 time points와 1:1 매칭
  tt <- sort(unique(data.m[,3]))
  tt <- cbind(tt, seq(0.0001, 1, length.out=length(tt)))
  data.u <- matrix(as.numeric(as.vector(data.m[,-1])),   # convert obs matrix to be numierc
                   nrow=nrow(data.m[,-1]),
                   ncol=ncol(data.m[,-1]))    
  
  for (i in 1:n){
    Y <- as.vector(data.u[(current+1):(current+m.l[i]),1])  # observed  measurements of ith curve
    meastime <- data.u[(current+1):(current+m.l[i]),2]      # measurement times of the ith curve
    meastime2 <- tt[which(tt[,1] %in% meastime), 2]
    # round한 gridpoint와 timepoints를 비교해서(같은 것 여러개 나옴) 같은 것의 2번째 값 사용
    gridtime <- apply(matrix(meastime), 1, function(x){ which( round(grids.u, 1) %in% x )[2] })   
    muy <- muhat[gridtime]
    Phiy  <- matrix(eigenfuncs.u[gridtime,1:K], ncol=K)
    Sigy <- Phiy %*% evalmat %*% t(Phiy) + sig2hat * diag(m.l[i])
    temp.y <- matrix(Y-muy)
    result[i,] <- evalmat %*% t(Phiy) %*% solve(Sigy, temp.y)
    current <- current + m.l[i]
    
    # # gridpoint와 실제 timepoint가 비슷한지 비교하기 위함
    # print(i)
    # print(rbind(meastime, grids.u[gridtime]))
  }
  return(result)
}



# spline knots 개수, PC 개수 정의
kn <- c(14)   # knots 개수
k <- c(2)   # PC 개수

# parameters for fpca.mle
ini.method <- "EM"        # optimization method
basis.method <- "bs"      # spline basis - "ns"의 경우 현재 오류 발생
sl.v <- rep(0.5, 10)      # Neuton-Raphson
max.step <- 50            # EM iteration 횟수
grid.l <- seq(0,1,0.01)   # EM에서는 사용 X
grids <- seq(0,1,0.002)   # EM의 grid

# fit candidate models by fpca.mle
fit <- fpca.mle(data, kn, k, ini.method, basis.method, sl.v, max.step, grid.l, grids)
summary(fit)

# fitted parameters
kn <- fit$selected_model[1]   # selected knots 개수
k <- fit$selected_model[2]    # selected PC 개수
grids.new <- fit$grid   # rescaled grid
evalest <- fit$eigenvalues   # estimated eigenvalues
sig2est <- fit$error_var   # estimated error variance
eigenfest <- fit$eigenfunctions   # estimated eigenfunctions(PC function)
muest <- fit$fitted_mean   # estimated mean function


# mean function
par(mfrow=c(1, k+1))
ind <- unique(data[, 1])
sub <- data[which(data[, 1]==ind[1]), ]
plot(sub[, 3], sub[, 2],
     type="o",
     xlim=c(range(data[, 3])[1], range(data[, 3])[2]),
     ylim=c(range(data[, 2])[1], range(data[, 2])[2]),
     xlab="Age (years)",
     ylab="Spinal bone density")
for (i in ind) {
  sub <- data[which(data[, 1]==i), ]
  lines(sub[, 3], sub[, 2], type="o")
}
lines(grids.new, fit$fitted_mean, col="red", lwd=3, type="l")

# PC function
for(i in 1:k){
  plot(grids.new, -eigenfest[i,], type="l", xlab="Age (years)", ylab="Princ. comp.", main=paste("PC", i))
}


##look at the CV scores and convergence for each model: note that model (M=5, r=4) does not converge.
result$cv_scores ##CV
result$converge ##convergence

# PC scores
fpcs <- fpca.score(data, grids.new, muest, evalest, eigenfest, sig2est, k)

# predicted curve
pred <- fpca.pred(fpcs, muest, eigenfest)

# true curve와 predicted curve 비교
tt <- sort(unique(data[,3]))
tt <- cbind(tt, seq(0.01, 1, length.out = length(tt)))
N <- length(grids.new)
par(mfrow=c(3,3))
for (i in 1:length(unique(data[,1]))){
  id <- unique(data[,1])[i]      # for curve i
  t.c <- data[data[,1] == id, 3]  # measurement points
  meastime2 <- tt[which(tt[,1] %in% t.c), 2]
  t.proj <- ceiling(N*meastime2)   # measurement points projected on the grid
  t.proj <- apply(matrix(t.c), 1, function(x){ which( round(grids.new, 1) %in% x )[2] })
  y.c <- data[data[,1]==id,2]   # obs
  y.pred.proj <- pred[t.proj, i]   # predicted obs on the measurement points
  
  # plots
  plot(t.c, y.c, ylim=c(0.3, 1.5), xlab="time", ylab="obs", main=paste("predicted trajectory of curve", id))
  points(grids.new, pred[, i], col=3, type='l')
}
par(mfrow=c(1,1))

