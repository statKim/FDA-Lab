########################################
### EM algorithm for reduced rank model
########################################

library(splines)
library(MASS)   # ginv 함수 사용하기위함
library(far)    # spline basis orthonormal하게 변환하기 위함

# function 정의
proc_sigma <- function(data, B, theta0, Theta, alpha, D, sigma) {
  N <- length(unique(data$idnum))
  theta0 <- matrix(theta0, ncol=1)
  comp <- c()
  for (i in 1:N) {
    B_i <- B[which(Tgrid %in% data$age[which(data$idnum==unique(data$idnum)[i])]), ]
    y <- matrix(data$spnbmd[which(data$idnum==unique(data$idnum)[i])], ncol=1)
    alpha_i <- matrix(alpha[i,], ncol=1)
    comp[i] <- t(y - B_i%*%theta0 - B_i%*%Theta%*%alpha_i) %*% (y - B_i%*%theta0 - B_i%*%Theta%*%alpha_i) +
      sum(diag( B_i%*%Theta%*%ginv( diag(1/D)+(t(Theta)%*%t(B_i)%*%B_i%*%Theta)/sigma )%*%t(Theta)%*%t(B_i) ))
  }
  result <- (1/nrow(data)) * sum(comp)
  return(result)
}

proc_D <- function(data, B, Theta, alpha, D, sigma) {
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
  # left_s <- list()
  # right_s <- list()
  for (i in 1:N) {
    B_i <- B[which(Tgrid %in% data$age[which(data$idnum==unique(data$idnum)[i])]), ]
    y <- matrix(data$spnbmd[which(data$idnum==unique(data$idnum)[i])], ncol=1)
    alpha_i <- matrix(alpha[i,], ncol=1)
    # left_s[i] <- t(B_i)%*%B_i
    # right_s[i] <- t(B_i)%*%(y - B_i%*%Theta%*%alpha_i)
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


sol_Theta <- function(data, B, theta0, Theta, D, alpha, sigma, alpha_outprod) {
  N <- length(unique(data$idnum))
  K <- ncol(Theta)
  q <- nrow(Theta)
  result <- matrix(0, q, K)
  for (j in 1:K) {
    for (i in 1:N) {
      B_i <- B[which(Tgrid %in% data$age[which(data$idnum==unique(data$idnum)[i])]), ]
      y <- matrix(data$spnbmd[which(data$idnum==unique(data$idnum)[i])], ncol=1)
      # alp_outer <- sol_alpha_outprod(B_i, Theta, alpha[i,], D, sigma)
      alp_outer <- alpha_outprod[[i]]
      if (i == 1) {
        # left_s <- alpha[i,j]^2 * t(B_i)%*%B_i
        left_s <- alp_outer[j,j] * t(B_i)%*%B_i
        right_s <- t(B_i)%*%( alpha[i,j] * (y - B_i%*%theta0) -
                                B_i%*%Theta%*%alp_outer[j, ] +
                                alp_outer[j,j]*B_i%*%Theta[,j] )        
      } else {
        # left_s <- left_s + alpha[i,j]^2 * t(B_i)%*%B_i
        left_s <- left_s + alp_outer[j,j] * t(B_i)%*%B_i
        right_s <- right_s + t(B_i)%*%(alpha[i,j] * (y - B_i%*%theta0) -
                                         B_i%*%Theta%*%alp_outer[j, ] +
                                         alp_outer[j,j]*B_i%*%Theta[,j] )
      }
    }
    result[,j] <- ginv(left_s) %*% right_s
  }
  return(result)
}

# alpha_i의 outerproduct
# sol_alpha_outprod <- function(B_i, Theta, alpha_i, D, sigma) {
#   alpha_i <- matrix(alpha_i, ncol = 1)
#   result <- alpha_i %*% t(alpha_i) + ginv( diag(1/D) + (t(Theta)%*%t(B_i)%*%B_i%*%Theta)/sigma )
#   return(result)
# }
sol_alpha_outprod <- function(data, B, Theta, alpha, D, sigma) {
  N <- length(unique(data$idnum))
  K <- ncol(Theta)
  result <- list()
  for (i in 1:N) {
    B_i <- B[which(Tgrid %in% data$age[which(data$idnum==unique(data$idnum)[i])]), ]
    alpha_i <- matrix(alpha[i,], ncol = 1)
    result[[i]] <- alpha_i %*% t(alpha_i) + ginv( diag(1/D) + (t(Theta)%*%t(B_i)%*%B_i%*%Theta)/sigma )
  }
  return(result)
}


sol_alpha <- function(data, B, theta0, Theta, D, sigma) {
  N <- length(unique(data$idnum))
  K <- ncol(Theta)
  result <- matrix(0, N, K)
  for (i in 1:N) {
    B_i <- B[which(Tgrid %in% data$age[which(data$idnum==unique(data$idnum)[i])]), ]
    y <- matrix(data$spnbmd[which(data$idnum==unique(data$idnum)[i])], ncol=1)
    alpha_i <- ginv( sigma*diag(1/D) + t(Theta)%*%t(B_i)%*%B_i%*%Theta ) %*% t(Theta)%*%t(B_i)%*%(y-B_i%*%theta0)
    result[i,] <- alpha_i
  }
  return(result)
}


orthog_Theta <- function(Theta, D) {
  K <- ncol(Theta)   # number of PCs
  Gamma <- Theta %*% diag(D) %*% t(Theta)
  r <- eigen(Gamma)
  result <- r$vectors[,1:K]
  return(result)
}



# fit reduced rank model for fpca
fpca.fit <- function(data, iter=100, init_value=.1, num_knots=4, num_pc=2, mixed.model=F) {
  # values
  iter <- 100
  N <- length(unique(data$idnum))
  num_knots <- num_knots + 2   # knots개수 = num_knots-2 => basis intercept때문
  age_range <- range(data$age)
  knots <- seq(from=age_range[1], to=age_range[2], length=num_knots)[-c(1,num_knots)]
  Tgrid <- unique( round( seq(age_range[1], age_range[2], 0.05), 1) )
  B <- ns(Tgrid, knots=knots, intercept=T)
  q <- num_knots   # basis matrix dimension
  k <- num_pc   # k=q => mixed effects model
  
  # mixed effects model fitting
  if (mixed.model == T) {
    k <- q
  }
  
  # orthnormalization for spline basis
  B <- orthonormalization(B, basis=F, norm=T)   # Gram-Schmidt orthogonalization
  
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
    theta0[t+1,] <- sol_theta0(data, B, Theta[[t]], alpha[[t]])
    Theta[[t+1]] <- sol_Theta(data, B, theta0[t+1,], Theta[[t]], D[t,], alpha[[t]], sigma[t], alpha_outer)
    
    # E-step (predict alpha)
    alpha[[t+1]] <- sol_alpha(data, B, theta0[t+1,], Theta[[t+1]], D[t,], sigma[t])
    alpha_outer <- sol_alpha_outprod(data, B, Theta[[t+1]], alpha[[t+1]], D[t,], sigma[t])
    
    sigma[t+1] <- proc_sigma(data, B, theta0[t+1,], Theta[[t+1]], alpha[[t+1]], D[t,], sigma[t])
    D[t+1,] <- proc_D(data, B, Theta[[t+1]], alpha[[t+1]], D[t,], sigma[t+1])
  }
  
  # estimate PC functions and mean function
  orth_Theta <- orthog_Theta(Theta[[t]], D[t,])
  fpc <- B%*%orth_Theta
  if (sum(fpc[,1]) < 0) {   # 1st PC의 부호를 양수로 만들기 위해
    orth_Theta <- -orthog_Theta(Theta[[t]], D[t,])
    fpc <- B%*%orth_Theta
  }
  mean_fn <- B%*%theta0[t,]
  
  # output
  result <- list()
  result[["basis"]] <- B
  result[["Timegrid"]] <- Tgrid
  result[["sigma^2"]] <- sigma[[t]]
  result[["alpha"]] <- alpha[[t]]
  result[["theta0"]] <- theta0[t,]
  result[["D"]] <- diag(D[t,])
  result[["eigenfunction"]] <- orth_Theta
  result[["MeanFunction"]] <- mean_fn
  result[["FPCscore"]] <- fpc
  
  return(result)
}



