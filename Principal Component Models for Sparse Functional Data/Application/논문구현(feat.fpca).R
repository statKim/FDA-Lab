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
data <- as.matrix( data[, c(1,3,2)] )   # column 순서 변경

### fpca package 이용해서 functional PC 계산
source("sparseFPCA.R")
# 기존 fpca package의 함수를 재정의
path <- "C:/Users/user/Desktop/KHS/fpca/R/"
source( paste(path, "fpca_func.R", sep="") )
source( paste(path, "fpca_pred.R", sep="") )
source( paste(path, "fpca_score.R", sep="") )
source( paste(path, "functions_EM.R", sep="") )
source( paste(path, "functions_GenData.R", sep="") )
source( paste(path, "functions_LocLin.R", sep="") )
source( paste(path, "functions_Optimization.R", sep="") )


# spline knots 개수, PC 개수 정의
kn <- c(4, 9, 14)   # knots 개수
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
fit$selected_model[1]   # best knots 개수
fit$selected_model[2]   # best PC 개수
fit$grid                # rescaled grid
fit$error_var           # estimated error variance
fit$eigenvalues         # estimated eigenvalues (PC function)
fit$eigenfunctions      # estimated eigenfunctions 
fit$fitted_mean         # estimated mean function
fit$PVE                 # PVE
fit$cv_scores           # CV
fit$converge            # convergence


# mean function, PC functions
par(mfrow=c(1, k+1))
fpc.plot(data, fit, xlab="Age (years)", ylab="Spinal bone density")

# PC scores
fpcs <- fpca.score(data, fit$grid, fit$fitted_mean, fit$eigenvalues, 
                   fit$eigenfunctions, fit$error_var, k, .1)   

# predicted curve
pred.plot(data, fit, fpcs, xlab="Age (years)", ylab="Spinal bone density", xrange=T)
pred.plot(data, fit, fpcs, xlab="Age (years)", ylab="Spinal bone density")





