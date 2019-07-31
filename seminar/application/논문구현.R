setwd("C:\\Users\\user\\Desktop\\KHS\\Thesis\\seminar\\application")

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


par(mfrow=c(3,3))
library(splines)
library(MASS)   # ginv 함수 사용하기위함
library(far)    # spline basis orthonormal하게 변환하기 위함
# 각 방법별로 연산속도 측정
system.time(
for (kn in c(6,11,16)) {
  # values
  iter <- 100
  N <- length(unique(data$idnum))
  num_knots <- kn   # knots개수 = num_knots-2
  age_range <- range(data$age)
  knots <- seq(from=age_range[1], to=age_range[2], length=num_knots)[-c(1,num_knots)]
  # Tgrid <- sort( union( sort(unique(data$age)), seq(age_range[1], age_range[2], 0.1) ) )
  Tgrid <- unique( round( seq(age_range[1], age_range[2], 0.05), 1) )
  B <- ns(Tgrid, knots=knots, intercept=T)
  q <- num_knots   # basis matrix dimension
  k <- 2   # k=q => mixed effects model
  
  B <- orthonormalization(B, basis=F, norm=T)   # Gram-Schmidt orthogonalization
  
  # inital value
  init_value <- .1
  
  # initial parameters
  sigma <- rep(init_value, iter)
  D <- matrix(init_value, iter, k)  # diagonal term  q x q
  alpha <- list(matrix(init_value, N, k))  # matrix list => ncol(B) 대신 pc 개수 k로 해야됨!!
  theta0 <- matrix(init_value, iter, q)
  Theta <- list(matrix(init_value, q, k))  # matrix list
  
  
  # apply EM-algorithm
  for (t in 1:(iter-1)) {
    # M-step

    theta0[t+1,] <- sol_theta0(data, B, Theta[[t]], alpha[[t]])
    Theta[[t+1]] <- sol_Theta(data, B, theta0[t,], Theta[[t]], D[t,], alpha[[t]], sigma[t])
  
    # E-step (predict alpha)
    alpha[[t+1]] <- sol_alpha(data, B, theta0[t+1,], Theta[[t+1]], D[t,], sigma[t]) 
    sigma[t+1] <- proc_sigma(data, B, theta0[t+1,], Theta[[t+1]], alpha[[t]], D[t,], sigma[t])
    D[t+1,] <- proc_D(data, B, theta0[t+1,], Theta[[t+1]], alpha[[t]], D[t,], sigma[t])
  }
  
  # estimate PC functions and mean function
  orth_Theta <- orthog_Theta(Theta[[t]], D[t,])
  fpc <- B%*%orth_Theta
  mean_fn <- B%*%theta0[t,]
  
  # 48 curves and mean function
  # par(mfrow=c(1,2))
  ind <- unique(data$idnum)
  sub <- data[which(data$idnum==ind[1]), ]
  plot(sub$age, sub$spnbmd, type="o", xlim=c(age_range[1], age_range[2]),
       ylim=c(range(data$spnbmd)[1], range(data$spnbmd)[2]), xlab="Age (years)", ylab="Spinal bone density")
  for (i in ind) {
    sub <- data[which(data$idnum==i), ]
    lines(sub$age, sub$spnbmd, type="o")
  }
  lines(Tgrid, mean_fn, col="red", lwd=3, type="l")
  
  # 1st PC function
  plot(Tgrid, fpc[,1], type="l", lwd=3, xlab="Age (years)", ylab="Princ. comp.")
  
  # compute residuals
  res <- c()
  for (i in unique(data$idnum)){
    res[which(data$idnum==i)] <- data$spnbmd[which(data$idnum==i)]-mean_fn[which(Tgrid %in% data$age[which(data$idnum==i)])]
  }
  # residual plot
  plot(data$age[which(data$idnum==ind[1])],
       res[which(data$idnum==ind[1])], type="o",
       xlim=c(age_range[1], age_range[2]), ylim=c(range(res)[1], range(res)[2]),
       xlab="Age (years)", ylab="Residuals")
  for (i in unique(data$idnum)){
    lines(data$age[which(data$idnum==i)],
          res[which(data$idnum==i)], type="o")
  }
  print(sum(res^2))
}
)


# function 정의
proc_sigma <- function(data, B, theta0, Theta, alpha, D, sigma) {
  N <- length(unique(data$idnum))
  theta0 <- matrix(theta0, ncol=1)
  comp <- c()
  for (i in 1:N) {
    B_i <- B[which(Tgrid %in% data$age[which(data$idnum==unique(data$idnum)[i])]), ]
    y <- matrix(data$spnbmd[which(data$idnum==unique(data$idnum)[i])], ncol = 1)
    alpha_i <- matrix(alpha[i,], ncol=1)
    comp[i] <- t(y - B_i%*%theta0 - B_i%*%Theta%*%alpha_i) %*% (y - B_i%*%theta0 - B_i%*%Theta%*%alpha_i) +
      sum(diag( B_i%*%Theta%*%ginv( diag(1/D)+(t(Theta)%*%t(B_i)%*%B_i%*%Theta)/sigma )%*%t(Theta)%*%t(B_i) ))
  }
  result <- (1/nrow(data)) * sum(comp)
  return(result)
}

proc_D <- function(data, B, theta0, Theta, alpha, D, sigma) {
  N <- length(unique(data$idnum))
  K <- ncol(Theta)
  result <- c()
  for (j in 1:K) {
    comp <- c()
    for (i in 1:N) {
      B_i <- B[which(Tgrid %in% data$age[which(data$idnum==unique(data$idnum)[i])]), ]
      comp[i] <- alpha[i,j]^2 + (ginv( diag(1/D)+(t(Theta)%*%t(B_i)%*%B_i%*%Theta)/sigma ))[j,j]
    }
    result[j] <- sum(comp)
  }
  result <- result / N
  return(result)
}


sol_theta0 <- function(data, B, Theta, alpha) {
  # left_s <- solve( Reduce("+", lapply(B, FUN = function(x){ t(x)%*%x })) )   # matrix sum
  N <- length(unique(data$idnum))
  for (i in 1:N) {
    B_i <- B[which(Tgrid %in% data$age[which(data$idnum==unique(data$idnum)[i])]), ]
    y <- matrix(data$spnbmd[which(data$idnum==unique(data$idnum)[i])], ncol=1)
    alpha_i <- matrix(alpha[i,], ncol=1)
    if (i == 1) {
      left_s <- t(B_i)%*%B_i
      right_s <- t(B_i)%*%(y - B_i%*%Theta%*%alpha_i)
    } else {
      left_s <- left_s + t(B_i)%*%B_i
      right_s <- right_s + t(B_i)%*%(y - B_i%*%Theta%*%alpha_i)
    }
  }
  result <- ginv(left_s) %*% right_s
  return(result)
}

orthog_Theta <- function(Theta, D) {
  K <- ncol(Theta)   # number of PCs
  Gamma <- Theta %*% diag(D) %*% t(Theta)
  # tryCatch(
  #   { 
  #     r <- eigen(Gamma)
  #     result <- r$vectors 
  #   },
  #   error = function(e) {
  #     print(paste(t,"번째에서 에러다"))
  #     r <- svd(Gamma)
  #     result <- r$u
  #   }
  # )
  r <- eigen(Gamma)
  result <- r$vectors[,1:K]
  return(result)
}

sol_Theta <- function(data, B, theta0, Theta, D, alpha, sigma) {
  N <- length(unique(data$idnum))
  K <- ncol(Theta)
  q <- nrow(Theta)
  result <- matrix(0, q, K)
  for (j in 1:K) {
    for (i in 1:N) {
      B_i <- B[which(Tgrid %in% data$age[which(data$idnum==unique(data$idnum)[i])]), ]
      y <- matrix(data$spnbmd[which(data$idnum==unique(data$idnum)[i])], ncol=1)
      alp_outer <- sol_alpha_outprod(B_i, Theta, alpha[i,], D, sigma)
      if (i == 1) {
        # left_s <- alpha[i,j]^2 * t(B_i)%*%B_i
        left_s <- alp_outer[j,j] * t(B_i)%*%B_i
        right_s <- t(B_i)%*%( alpha[i,j] * (y - B_i%*%theta0) -
                                  sum(B_i%*%Theta%*%alp_outer[j, ]) +   # dimension??? => n_i
                                  alp_outer[j,j]*B_i%*%Theta[,j] )        
      } else {
        # left_s <- left_s + alpha[i,j]^2 * t(B_i)%*%B_i
        left_s <- left_s + alp_outer[j,j] * t(B_i)%*%B_i
        right_s <- right_s + t(B_i)%*%(alpha[i,j] * (y - B_i%*%theta0) -
                                            sum(B_i%*%Theta%*%alp_outer[j, ]) +
                                            alp_outer[j,j]*B_i%*%Theta[,j] )
      }
    }
    result[,j] <- ginv(left_s) %*% right_s
  }
  return(result)
}

# alpha_i의 outerproduct
sol_alpha_outprod <- function(B_i, Theta, alpha_i, D, sigma) {
  alpha_i <- matrix(alpha_i, ncol = 1)
  result <- alpha_i %*% t(alpha_i) + ginv( diag(1/D) + (t(Theta)%*%t(B_i)%*%B_i%*%Theta)/sigma )
  return(result)
}


sol_alpha <- function(data, B, theta0, Theta, D, sigma) {
  N <- length(unique(data$idnum))
  K <- ncol(Theta)
  result <- matrix(0, N, K)
  for (i in 1:N) {
    B_i <- B[which(Tgrid %in% data$age[which(data$idnum==unique(data$idnum)[i])]), ]
    y <- matrix(data$spnbmd[which(data$idnum==unique(data$idnum)[i])], ncol=1)
    alpha_i <- ginv( sigma*diag(1/D) + t(Theta)%*%t(B_i)%*%B_i%*%Theta ) %*%
               t(Theta)%*%t(B_i)%*%(y-B_i%*%theta0)
    result[i,] <- alpha_i
  }
  return(result)
}
