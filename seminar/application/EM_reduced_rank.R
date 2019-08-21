########################################
### EM algorithm for reduced rank model
########################################
# 데이터를 정의해줘야함
# input data : 1st column(id), 2nd column(time points), 3rd column(y)


library(splines)
library(MASS)   # ginv 함수 사용하기위함
library(far)    # spline basis orthonormal하게 변환하기 위함

# function 정의
proc_sigma <- function(data, B, theta0, Theta, alpha, D, sigma, Tgrid) {
  N <- length(unique(data[,1]))
  K <- ncol(Theta)
  theta0 <- matrix(theta0, ncol=1)
  comp <- c()
  for (i in 1:N) {
    B_i <- B[which(Tgrid %in% data[,2][which(data[,1]==unique(data[,1])[i])]), ]
    y <- matrix(data[,3][which(data[,1]==unique(data[,1])[i])], ncol=1)
    alpha_i <- matrix(alpha[i,], ncol=1)
    if (K > 1) {   # PC 개수가 1인 경우의 조건
      comp[i] <- t(y - B_i%*%theta0 - B_i%*%Theta%*%alpha_i) %*% (y - B_i%*%theta0 - B_i%*%Theta%*%alpha_i) +
        sum(diag( B_i%*%Theta%*%ginv( diag(1/D)+(t(Theta)%*%t(B_i)%*%B_i%*%Theta)/sigma )%*%t(Theta)%*%t(B_i) ))
    } else {
      comp[i] <- t(y - B_i%*%theta0 - B_i%*%Theta%*%alpha_i) %*% (y - B_i%*%theta0 - B_i%*%Theta%*%alpha_i) +
        sum(diag( B_i%*%Theta%*%ginv( (1/D)+(t(Theta)%*%t(B_i)%*%B_i%*%Theta)/sigma )%*%t(Theta)%*%t(B_i) ))
    }
  }
  result <- (1/nrow(data)) * sum(comp)
  return(result)
}

proc_D <- function(data, B, Theta, alpha, D, sigma, Tgrid) {
  N <- length(unique(data[,1]))
  K <- ncol(Theta)
  result <- c()
  for (j in 1:K) {
    comp <- c()
    for (i in 1:N) {
      B_i <- B[which(Tgrid %in% data[,2][which(data[,1]==unique(data[,1])[i])]), ]
      if (K > 1) {   # PC 개수가 1인 경우의 조건
        comp[i] <- alpha[i,j]^2 + (ginv( diag(1/D)+(t(Theta)%*%t(B_i)%*%B_i%*%Theta)/sigma ))[j,j]
      } else {
        comp[i] <- alpha[i,j]^2 + (ginv( (1/D)+(t(Theta)%*%t(B_i)%*%B_i%*%Theta)/sigma ))[j,j]
      }
    }
    result[j] <- sum(comp)
  }
  result <- result / N
  return(result)
}


sol_theta0 <- function(data, B, Theta, alpha, Tgrid) {
  N <- length(unique(data[,1]))
  for (i in 1:N) {
    B_i <- B[which(Tgrid %in% data[,2][which(data[,1]==unique(data[,1])[i])]), ]
    y <- matrix(data[,3][which(data[,1]==unique(data[,1])[i])], ncol=1)
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


sol_Theta <- function(data, B, theta0, Theta, D, alpha, sigma, alpha_outprod, Tgrid) {
  N <- length(unique(data[,1]))
  K <- ncol(Theta)
  q <- nrow(Theta)
  result <- matrix(0, q, K)
  for (j in 1:K) {
    for (i in 1:N) {
      B_i <- B[which(Tgrid %in% data[,2][which(data[,1]==unique(data[,1])[i])]), ]
      y <- matrix(data[,3][which(data[,1]==unique(data[,1])[i])], ncol=1)
      alp_outer <- alpha_outprod[[i]]
      if (i == 1) {
        left_s <- alp_outer[j,j] * t(B_i)%*%B_i
        right_s <- t(B_i)%*%( alpha[i,j] * (y - B_i%*%theta0) -
                                B_i%*%Theta%*%alp_outer[j, ] +
                                alp_outer[j,j]*B_i%*%Theta[,j] )        
      } else {
        left_s <- left_s + alp_outer[j,j] * t(B_i)%*%B_i
        right_s <- right_s + t(B_i)%*%( alpha[i,j] * (y - B_i%*%theta0) -
                                          B_i%*%Theta%*%alp_outer[j, ] +
                                          alp_outer[j,j]*B_i%*%Theta[,j] )
      }
    }
    result[,j] <- ginv(left_s) %*% right_s
  }
  return(result)
}

# alpha_i의 outerproduct
sol_alpha_outprod <- function(data, B, Theta, alpha, D, sigma, Tgrid) {
  N <- length(unique(data[,1]))
  K <- ncol(Theta)
  result <- list()
  for (i in 1:N) {
    B_i <- B[which(Tgrid %in% data[,2][which(data[,1]==unique(data[,1])[i])]), ]
    alpha_i <- matrix(alpha[i,], ncol = 1)
    if (K > 1) {   # PC 개수가 1인 경우의 조건
      result[[i]] <- alpha_i %*% t(alpha_i) + ginv( diag(1/D) + (t(Theta)%*%t(B_i)%*%B_i%*%Theta)/sigma )
    } else {
      result[[i]] <- alpha_i %*% t(alpha_i) + ginv( (1/D) + (t(Theta)%*%t(B_i)%*%B_i%*%Theta)/sigma )
    }
    
  }
  return(result)
}


sol_alpha <- function(data, B, theta0, Theta, D, sigma, Tgrid) {
  N <- length(unique(data[,1]))
  K <- ncol(Theta)
  result <- matrix(0, N, K)
  for (i in 1:N) {
    B_i <- B[which(Tgrid %in% data[,2][which(data[,1]==unique(data[,1])[i])]), ]
    y <- matrix(data[,3][which(data[,1]==unique(data[,1])[i])], ncol=1)
    if (K > 1) {   # PC 개수가 1인 경우의 조건
      alpha_i <- ginv( sigma*diag(1/D) + t(Theta)%*%t(B_i)%*%B_i%*%Theta ) %*% t(Theta)%*%t(B_i)%*%(y-B_i%*%theta0)
    } else {
      alpha_i <- ginv( sigma*(1/D) + t(Theta)%*%t(B_i)%*%B_i%*%Theta ) %*% t(Theta)%*%t(B_i)%*%(y-B_i%*%theta0)
    }
    
    result[i,] <- alpha_i
  }
  return(result)
}


orthog_Theta <- function(Theta, D) {
  K <- ncol(Theta)   # number of PCs
  if (K > 1) {   # PC 개수가 1인 경우의 조건
    Gamma <- Theta %*% diag(D) %*% t(Theta)
  } else {
    Gamma <- D * Theta %*% t(Theta)
  }
  # set.seed(5)   # 같은 부호의 eigenvector 뽑기 위함
  r <- eigen(Gamma)
  eigenvector <- r$vectors[,1:K]
  eigenvalues <- r$values[1:K]
  return(list(eigenvector=eigenvector, eigenvalues=eigenvalues))
}



# fit reduced rank model for fpca
fpca.fit <- function(data, iter=100, init_value=c(.1, .1, .1, .1, .1), num_knots=4, num_pc=2, mixed.model=F, grid_sep=0.1) {
  # 연산속도 높이기 위해 matrix로 변환
  data <- as.matrix(data)
  
  # 필요한 변수 정의
  N <- length(unique(data[,1]))
  num_knots <- num_knots + 2   # knots개수 = num_knots-2 => basis intercept때문
  age_range <- range(data[,2])
  knots <- seq(from=age_range[1], to=age_range[2], length=num_knots)[-c(1,num_knots)]
  if (grid_sep >= 1) {
    interval <- 0
  } else if (grid_sep >= 0.1) {
    interval <- 1
  } else if (grid_sep >= 0.01) {
    interval <- 2
  } else if (grid_sep >= 0.001) {
    interval <- 3
  } else {
    interval <- 3
  }
  Tgrid <- unique( round( seq(age_range[1], age_range[2], grid_sep/2), interval) )
  # Tgrid <- unique( round( seq(age_range[1], age_range[2], 0.05), 1) )
  B <- ns(Tgrid, knots=knots, intercept=T)
  q <- num_knots   # basis matrix dimension
  k <- num_pc   # k=q => mixed effects model
  
  # mixed effects model fitting
  if (mixed.model == T) {
    k <- q
  }
  
  # orthnormalization for spline basis
  # B <- orthonormalization(B, basis=F, norm=T)   # Gram-Schmidt orthogonalization
  QR <- qr(B)
  R <- qr.R(QR)
  g <- nrow(B)
  L <- age_range[2] - age_range[1]
  T_mat <- sqrt(g/L)*t(solve(R))
  # L/g * t(B%*%t(T_mat))%*%B%*%t(T_mat)   # identity matrix => orthonormal
  B <- B%*%t(T_mat) * sqrt(L/g)
  
  # initial parameters
  sigma <- rep(init_value[1], iter)
  D <- matrix(init_value[2], iter, k)  # diagonal term  q x q
  alpha <- list(matrix(init_value[3], N, k))  # matrix list => ncol(B) 대신 pc 개수 k로 해야됨!!
  alpha_outer <- rep(list(matrix(init_value[3], k, k)), N)
  theta0 <- matrix(init_value[4], iter, q)
  Theta <- list(matrix(init_value[5], q, k))  # matrix list
  
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
  
  # estimate PC functions
  orth_Theta <- orthog_Theta(Theta[[t]], D[t,])
  pc_func <- B%*%orth_Theta$eigenvector   # PC function
  # if (sum(pc_func[,1]) < 0) {   # 1st PC의 부호를 양수로 만들기 위해
  #   orth_Theta$eigenvector <- -orth_Theta$eigenvector
  #   pc_func <- B%*%orth_Theta$eigenvector
  # }
  
  # estimate mean function
  mean_fn <- B%*%theta0[t,]
  
  # estimate PC score
  # fpc <- 
  
  # PVE 계산
  PVE <- t( data.frame(individual=orth_Theta$eigenvalues / sum(orth_Theta$eigenvalues),
                       cumulative=cumsum(orth_Theta$eigenvalues) / sum(orth_Theta$eigenvalues)) )
  colnames(PVE) <- paste("PC", 1:k, sep="")
  
  # log likelihood
  loglik <- c()
  # sloglik <- c()
  for (i in 1:N) {
    B_i <- B[which(Tgrid %in% data[,2][which(data[,1]==unique(data[,1])[i])]), ]
    y_i <- matrix(data[,3][which(data[,1]==unique(data[,1])[i])], ncol=1)
    n_i <- length(y_i)
    alpha_i <- matrix(alpha[[t]][i,], ncol=1)
    if (k > 1) {   # PC 개수가 1인 경우의 조건
      loglik[i] <- ( t(y_i-B_i%*%theta0[t,])%*%ginv(diag(rep(sigma[t], n_i))+B_i%*%Theta[[t]]%*%diag(D[t,])%*%
                                                      t(Theta[[t]])%*%t(B_i))%*%(y_i-B_i%*%theta0[t,])/(-2) ) -
        log( (2*pi)^(n_i/2) * sqrt( det( diag(rep(sigma[t], n_i))+B_i%*%Theta[[t]]%*%diag(D[t,])%*%t(Theta[[t]])%*%t(B_i) ) ) )
    } else {
      loglik[i] <- ( t(y_i-B_i%*%theta0[t,])%*%ginv(diag(rep(sigma[t], n_i))+D[t,]*B_i%*%Theta[[t]]%*%
                                                      t(Theta[[t]])%*%t(B_i))%*%(y_i-B_i%*%theta0[t,])/(-2) ) -
        log( (2*pi)^(n_i/2) * sqrt( det( diag(rep(sigma[t], n_i))+D[t,]*B_i%*%Theta[[t]]%*%t(Theta[[t]])%*%t(B_i) ) ) )
      
    }
  }
  loglik <- sum(loglik)
  
  # predicted curve
  pred.curve <- list()
  for (i in 1:N){
    ind <- which(Tgrid %in% data[,2][which(data[,1]==unique(data[,1])[i])])
    time_point <- Tgrid[ind]
    t_range <- range(ind)
    alpha_i <- matrix(alpha[[t]][i,], ncol=1)
    pc_curve <- pc_func[t_range[1]:t_range[2],]
    y_hat <- mean_fn[t_range[1]:t_range[2]] + pc_curve%*%alpha_i
    
    pred.curve[[ as.character( unique(data[,1])[i] ) ]] <- data.frame(time=Tgrid[t_range[1]:t_range[2]],
                                                                      y.hat=y_hat)
  }
  
  # output
  result <- list()
  result[["basis"]] <- B
  result[["knots"]] <- knots  
  result[["Timegrid"]] <- Tgrid
  result[["sigma2"]] <- sigma[[t]]
  result[["alpha"]] <- alpha[[t]]
  result[["theta0"]] <- theta0[t,]
  result[["D"]] <- diag(D[t,])
  result[["Theta"]] <- orth_Theta$eigenvector
  result[["MeanFunction"]] <- mean_fn
  result[["PCfunction"]] <- pc_func
  # result[["FPCscore"]] <- fpc
  result[["Loglik"]] <- loglik
  result[["predcurve"]] <- pred.curve
  result[["PVE"]] <- PVE
  
  return(result)
}



