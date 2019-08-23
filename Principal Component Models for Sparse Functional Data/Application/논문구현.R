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


source("sparseFPCA.R")
library(splines)
library(far)    # spline basis orthonormal하게 변환하기 위함
par(mfrow=c(3,2))
for (kn in c(6)) {
  # values
  iter <- 100
  N <- length(unique(data$idnum))
  num_knots <- kn   # knots개수 = num_knots-2
  age_range <- range(data$age)
  knots <- seq(from=age_range[1], to=age_range[2], length=num_knots)[-c(1,num_knots)]
  Tgrid <- unique( round( seq(age_range[1], age_range[2], 0.05), 1) )
  B <- ns(Tgrid, knots=knots, intercept=T)
  q <- num_knots   # basis matrix dimension
  k <- 2   # k=q => mixed effects model
  
  # B <- orthonormalization(B, basis=F, norm=T)   # Gram-Schmidt orthogonalization
  QR <- qr(B)
  R <- qr.R(QR)
  g <- nrow(B)
  L <- age_range[2] - age_range[1]
  T_mat <- sqrt(g/L)*t(solve(R))
  L/g * t(B%*%t(T_mat))%*%B%*%t(T_mat)   # identity matrix => orthonormal
  # T_mat <- sqrt(nrow(B)/nrow(data))*t(solve(R))
  B <- B%*%t(T_mat) * sqrt(L/g)
  
  
  # inital value
  init_value <- .1
  
  # initial parameters
  sigma <- rep(init_value, iter)
  D <- matrix(init_value, iter, k)  # diagonal term  q x q
  alpha <- list(matrix(init_value, N, k))  # matrix list => ncol(B) 대신 pc 개수 k로 해야됨!!
  theta0 <- matrix(init_value, iter, q)
  Theta <- list(matrix(init_value, q, k))  # matrix list
  alpha_outer <- rep(list(matrix(init_value, k, k)), N)
  
  # apply EM-algorithm
  for (t in 1:(iter-1)) {
    # M-step
    sigma[t+1] <- proc_sigma(data, B, theta0[t,], Theta[[t]], alpha[[t]], D[t,], sigma[t], Tgrid)
    D[t+1,] <- proc_D(data, B, Theta[[t]], alpha[[t]], D[t,], sigma[t+1], Tgrid)
    
    theta0[t+1,] <- sol_theta0(data, B, Theta[[t]], alpha[[t]], Tgrid)
    Theta[[t+1]] <- sol_Theta(data, B, theta0[t,], Theta[[t]], D[t+1,], alpha[[t]], sigma[t+1], alpha_outer, Tgrid)
    
    # E-step (predict alpha)
    alpha[[t+1]] <- sol_alpha(data, B, theta0[t+1,], Theta[[t+1]], D[t+1,], sigma[t+1], Tgrid)
    alpha_outer <- sol_alpha_outprod(data, B, Theta[[t+1]], alpha[[t+1]], D[t+1,], sigma[t+1], Tgrid)
  }
  
  # estimate PC functions and mean function
  orth_Theta <- orthog_Theta(Theta[[t]], D[t,])
  fpc <- B%*%orth_Theta
  if (sum(fpc[,1]) < 0) {   # 1st PC의 부호를 양수로 만들기 위해
    orth_Theta <- -orthog_Theta(Theta[[t]], D[t,])
    fpc <- B%*%orth_Theta
  }
  mean_fn <- B%*%theta0[t,]
  
  
  # likelihood
  likelihood <- c()
  for (i in 1:N) {
    B_i <- B[which(Tgrid %in% data$age[which(data$idnum==unique(data$idnum)[i])]), ]
    y_i <- matrix(data$spnbmd[which(data$idnum==unique(data$idnum)[i])], ncol=1)
    n_i <- length(y_i)
    likelihood[i] <- exp( t(y_i-B_i%*%theta0[t,])%*%
                            ginv(diag(rep(sigma[t], n_i))+B_i%*%Theta[[t]]%*%diag(D[t,])%*%
                            t(Theta[[t]])%*%t(B_i))%*%(y_i-B_i%*%theta0[t,])/(-2) ) /
                     ( (2*pi)^(n_i/2) * sqrt( det( diag(rep(sigma[t], n_i))+B_i%*%Theta[[t]]%*%diag(D[t,])%*%t(Theta[[t]])%*%t(B_i) ) ) )
  }
  likelihood <- prod(likelihood)
  
  # 48 curves and mean function
  # par(mfrow=c(1,2))
  ind <- unique(data$idnum)
  sub <- data[which(data$idnum==ind[1]), ]
  plot(sub$age, sub$spnbmd,
       type="o",
       xlim=c(age_range[1], age_range[2]),
       ylim=c(range(data$spnbmd)[1], range(data$spnbmd)[2]),
       xlab="Age (years)",
       ylab="Spinal bone density")
  for (i in ind) {
    sub <- data[which(data$idnum==i), ]
    lines(sub$age, sub$spnbmd, type="o")
  }
  lines(Tgrid, mean_fn, col="red", lwd=3, type="l")
  
  # 1st PC function
  plot(Tgrid, fpc[,1], type="l", lwd=3, xlab="Age (years)", ylab="Princ. comp.")
  
  # # compute residuals
  # res <- c()
  # for (i in unique(data$idnum)){
  #   res[which(data$idnum==i)] <- data$spnbmd[which(data$idnum==i)]-mean_fn[which(Tgrid %in% data$age[which(data$idnum==i)])]
  # }
  # # residual plot
  # plot(data$age[which(data$idnum==ind[1])],
  #      res[which(data$idnum==ind[1])], type="o",
  #      xlim=c(age_range[1], age_range[2]), ylim=c(range(res)[1], range(res)[2]),
  #      xlab="Age (years)", ylab="Residuals")
  # for (i in unique(data$idnum)){
  #   lines(data$age[which(data$idnum==i)],
  #         res[which(data$idnum==i)], type="o")
  # }
  # print(sum(res^2))
}

source("sparseFPCA.R")
par(mfrow=c(3,2))
# for (start_point in seq(-5, 5, 0.5)) {
# (sigma,D,alpha,theta0,Theta)   
init_value <- c(.1, .1, .1, .1, .1)
for (kn in c(4)) {
  fit <- fpca.fit(data, iter=100, init_value=init_value, num_knots=kn, num_pc=2)
  
  # 48 curves and mean function
  # par(mfrow=c(1,2))
  ind <- unique(data$idnum)
  sub <- data[which(data$idnum==ind[1]), ]
  plot(sub$age, sub$spnbmd,
       type="o",
       xlim=c(range(data$age)[1], range(data$age)[2]),
       ylim=c(range(data$spnbmd)[1], range(data$spnbmd)[2]),
       xlab="Age (years)",
       ylab="Spinal bone density")
  for (i in ind) {
    sub <- data[which(data$idnum==i), ]
    lines(sub$age, sub$spnbmd, type="o")
  }
  lines(fit$Timegrid, fit$MeanFunction, col="red", lwd=3, type="l")
  
  # 1st PC function
  plot(fit$Timegrid, fit$PCfunction[,1], type="l", lwd=3, xlab="Age (years)", ylab="Princ. comp.", ylim=c(0, 0.12))
  
  # # mixed effects model과 비교
  # fit2 <- fpca.fit(data, iter=100, init_value=init_value, num_knots=kn, mixed.model=T)
  # lines(fit2$Timegrid, fit2$PCfunction[,1], lwd=3, lty=2, type="l")

  # loglikelihood 비교
  print( paste("Reduced rank(", kn, " knots) : ", round(fit$Loglik, 2), sep="") )
  # print( paste("mixed effects(", kn, " knots) : ", fit2$Loglik, sep="") )
}
# print("===================================================")
# }

# predicted curve
fit <- fpca.fit(data, iter=100, init_value=init_value, num_knots=4, num_pc=2)
par(mfrow=c(5,5))
for (i in 1:25){
  ind <- which(fit$Timegrid %in% data$age[which(data$idnum==unique(data$idnum)[i])])
  time_point <- fit$Timegrid[ind]
  t_range <- range(ind)
  y_hat <- fit$MeanFunction[t_range[1]:t_range[2]] + fit$PCfunction[t_range[1]:t_range[2],]%*%fit$alpha[i,]
  
  y <- matrix(data$spnbmd[which(data$idnum==unique(data$idnum)[i])], ncol=1)
  
  plot(time_point, y,
       xlim=range(data$age), ylim=range(data$spnbmd),
       xlab="Age (years)", ylab="Spinal bone density",
       main=paste("predicted trajectory of curve", unique(data$idnum)[i]))
  points(fit$Timegrid[t_range[1]:t_range[2]], y_hat, col=3, type='l')
}

# 1st PC function
plot(fit$Timegrid, fit$PCfunction[,1], type="l", lwd=3, xlab="Age (years)", ylab="Princ. comp.", ylim=c(0, 0.12))

