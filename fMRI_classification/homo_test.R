Rcpp::sourceCpp("src/homo_test.cpp")
library(fda.usc)

# Calculate test statistic from permutation test
perm_stat <- function(X1, X2) {
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  n <- n1 + n2
  
  # Pooled sample covariance matrix
  S1 <- cov(X1)
  S2 <- cov(X2)
  S <- ((n1-1)*S1 + (n2-1)*S2) / (n-2)
  
  # Inverse of pooled sample covariance matrix
  eig.obj <- eigen(S)
  pos.eig.idx <- which(eig.obj$values > 1e-10)
  S_inv <- eig.obj$vectors[, pos.eig.idx] %*% diag(1/eig.obj$values[pos.eig.idx]) %*% t(eig.obj$vectors[, pos.eig.idx])
  
  # Test statistic using C++
  gamma <- computeGamma(X1, X2, S_inv)
  
  # # Weight function I_W(x)
  # I_w <- function(x, S_inv){ exp( -( as.numeric(t(x) %*% (S_inv %*% x)) )/2 ) }
  # 
  # # Test statistic
  # gamma <- 0
  # # 1st term
  # for (i in 1:n1) {
  #   for (j in 1:n1) {
  #     if (i == j) {
  #       # diagonal term
  #       gamma <- gamma + I_w(X1[i, ] - X1[j, ], S_inv) / n1^2
  #     } else if (i > j) {
  #       # off-diagonal term
  #       gamma <- gamma + 2*I_w(X1[i, ] - X1[j, ], S_inv) / n1^2
  #     }
  #   }
  #   # 3rd term
  #   if (n1 >= n2) {
  #     for (j in 1:n2) {
  #       if (i == j | ((i > j & i > n2))) {
  #         # diagonal term
  #         gamma <- gamma - I_w(X1[i, ] - X2[j, ], S_inv) * 2/(n1*n2)
  #       } else if (i > j & i <= n2) {
  #         # off-diagonal term
  #         gamma <- gamma - 2*I_w(X1[i, ] - X2[j, ], S_inv) * 2/(n1*n2)
  #       }
  #       # gamma <- gamma - I_w(X1[i, ] - X2[j, ], S_inv) * 2/(n1*n2)
  #     }
  #   }
  # }
  # 
  # # 2nd term
  # for (i in 1:n2) {
  #   for (j in 1:n2) {
  #     if (i == j) {
  #       # diagonal term
  #       gamma <- gamma + I_w(X2[i, ] - X2[j, ], S_inv) / n2^2
  #     } else if (i > j) {
  #       # off-diagonal term
  #       gamma <- gamma + 2*I_w(X2[i, ] - X2[j, ], S_inv) / n2^2
  #     }
  #   }
  #   # 3rd term
  #   if (n1 < n2) {
  #     for (j in 1:n1) {
  #       if (i == j | ((i > j & i > n1))) {
  #         # diagonal term
  #         gamma <- gamma - I_w(X1[j, ] - X2[i, ], S_inv) * 2/(n1*n2)
  #       } else if (i > j & i <= n1) {
  #         # off-diagonal term
  #         gamma <- gamma - 2*I_w(X1[j, ] - X2[i, ], S_inv) * 2/(n1*n2)
  #       }
  #       # gamma <- gamma - I_w(X1[i, ] - X2[j, ], S_inv) * 2/(n1*n2)
  #     }
  #   }
  # }
  # # # 3rd term
  # # for (i in 1:n1) {
  # #   for (j in 1:n2) {
  # #     if (i == j) {
  # #       # diagonal term
  # #       gamma <- gamma - I_w(X1[i, ] - X2[j, ]) / (2/n1*n2)
  # #     } else if (i > j) {
  # #       # off-diagonal term
  # #       gamma <- gamma - 2*I_w(X1[i, ] - X2[j, ]) / (2/n1*n2)
  # #     }
  # #     # gamma <- gamma - I_w(X1[i, ] - X2[j, ]) / (2/n1*n2)
  # #   }
  # # }
  
  return(gamma)
}


# Test homogeneity of 2-sample functional data
fun_homo_test <- function(X1, X2, method = "ECF", depth = depth.FM, trim = 0.1, alpha = 0.05, B = 500) {
  
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  n <- n1 + n2
  m <- ncol(X1)
  
  X <- rbind(X1, X2)
  group <- c(rep(0, n1), rep(1, n2))
  idx_g1 <- 1:n1
  idx_g2 <- (n1+1):n
  
  if (ncol(X2) != m) {
    stop("X1 and X2 have different timepoints!")
  }
  
  if (method == "ECF") {
    # Jiang, Q., Hušková, M., Meintanis, S. G., & Zhu, L. (2019). Asymptotics, finite-sample comparisons and applications for two-sample tests with functional data. Journal of Multivariate Analysis, 170, 202-220.
    
    # Test statistic
    test_stat <- perm_stat(X1, X2)
    
    # Permutation test statistics
    gamma <- rep(0, B)
    for (b in 1:B) {
      # Permutation
      perm <- sample(1:n, n)
      X1_perm <- X[perm[idx_g1], ]
      X2_perm <- X[perm[idx_g2], ]
      
      # Calculate test statstic from bth permutation
      gamma[b] <- perm_stat(X1_perm, X2_perm)
    }
    
    # Compute p-value
    # critical_value <- quantile(gamma, 1-alpha)
    # test_stat > critical_value
    p_value <- mean(gamma > test_stat)
    
    # hist(gamma)
    # abline(v = test_stat, col = 2, lwd = 3)
    
    res <- list(
      test_stat = test_stat,
      p_value = p_value,
      test_stat_perm = gamma
    )
  } else if (method == "ddplot") {
    # Calle-Saldarriaga, A., Laniado, H., Zuluaga, F., & Leiva, V. (2021). Homogeneity tests for functional data based on depth-depth plots with chemical applications. Chemometrics and Intelligent Laboratory Systems, 219, 104420.
    
    # Compute depth
    fdataobj <- fdata(rbind(X1, X2), gr)
    fdataobj_X1 <- fdata(X1, gr)
    fdataobj_X2 <- fdata(X2, gr)
    
    depth.obj <- depth(fdataobj, fdataobj_X1, trim = trim)
    depth_x <- depth.obj$dep
    depth.obj <- depth(fdataobj, fdataobj_X2, trim = trim)
    depth_y <- depth.obj$dep
    
    # # DD-plot
    # plot(depth_x, depth_y)
    # abline(a = 0, b = 1, col = 2, lwd = 2)
    
    # Least square estimate
    fit <- lm(depth_x ~ depth_y)
    fit_summary <- summary(fit)
    
    sigma_beta0 <- sqrt( sum(residuals(fit)^2) * sum(depth_y) / ((n-2) * sum((depth_y - mean(depth_y))^2)) )
    # sigma_beta0 <- sqrt( sum(residuals(fit)^2) * mean(depth_y)^2 / ((n-2) * sum((depth_y - mean(depth_y))^2)) )
    sigma_beta1 <- sqrt( sum(residuals(fit)^2) / ((n-2) * sum((depth_y - mean(depth_y))^2)) )
    
    T0 <- fit_summary$coefficients[1, 1] / sigma_beta0
    T1 <- fit_summary$coefficients[2, 1] / sigma_beta1
    
    # Bootstrap under H_0: X =_d Y
    set.seed(1000)
    T_star <- matrix(0, B, 2)
    for (b in 1:B) {
      # Bootstrap
      X1_boot <- X1[sample(1:n1, n1, replace = T), ]
      X2_boot <- X2[sample(1:n2, n2, replace = T), ]
      
      # Compute depth
      # fdataobj <- fdata(rbind(X1_boot, X2), gr)
      fdataobj <- fdata(rbind(X1_boot, X2_boot), gr)
      fdataobj_X1 <- fdata(X1_boot, gr)
      fdataobj_X2 <- fdata(X2_boot, gr)
      
      depth.obj <- depth(fdataobj, fdataobj_X1, trim = trim)
      depth_x <- depth.obj$dep
      depth.obj <- depth(fdataobj, fdataobj_X2, trim = trim)
      depth_y <- depth.obj$dep
      
      # Least square estimate
      fit <- lm(depth_x ~ depth_y)
      fit_summary <- summary(fit)
      
      T_star[b, ] <- (fit_summary$coefficients[, 1] - c(0, 1)) / fit_summary$coefficients[, 2]
    }
    
    p_value <- rbind(
      colMeans(T_star > c(T0, T1)),
      colMeans(T_star < c(T0, T1))
    )
    p_value <- 2 * apply(p_value, 2, min)
    
    adj_p_value <- min(2*min(p_value), max(p_value))
    
    res <- list(
      p_value = adj_p_value,
      p_value_list = p_value,
      test_stat_boot = T_star,
      test_stat = c(T0, T1)
    )
  }
  
  return(res)
}

# library(ddalpha)

##################################################
### Load fMRI data
##################################################
library(tidyverse)

# Class label of 191 subjects
y <- read.csv("./fMRI_Classification/class_label.csv", header = T)
y <- y$y

# Functional covariates from 82 retions
file_list <- list.files("./fMRI_Classification/AlignedSubject/")
for (i in 1:length(file_list)) {
  df <- read.csv(paste0("./fMRI_Classification/AlignedSubject/", file_list[i]), header = T)
  df <- df[, -1] %>% as.matrix()
  
  if (i == 1) {
    X <- array(0, c(length(file_list), ncol(df), nrow(df)))
  }
  X[i, , ] <- t(df)
}
dim(X)

n <- dim(X)[1]   # number of curves
m <- dim(X)[2]   # number of timepoints
p <- dim(X)[3]  # number of functional variables
gr <- seq(0, 1, length.out = m)


##################################################
### Homogeneity test for each region
##################################################

test_obj <- list()
for (i in 1:p) {
  print(paste(i, "th region"))
  
  set.seed(i)
  X1 <- X[y == 0, , i]
  X2 <- X[y == 1, , i]
  
  start_time <- Sys.time()
  # test_obj[[i]] <- fun_homo_test(X1, X2, method = "ECF", B = 10)
  test_obj[[i]] <- fun_homo_test(X1, X2, method = "ddplot", depth = depth.FM, B = 500)
  end_time <- Sys.time()
  
  # test_obj[[i]]
  
  print(test_obj[[i]]$p_value)
  print(end_time - start_time)
}
# save(test_obj, file = "RData/homo_test_ddplot.RData")

p_value <- sapply(test_obj, function(x){ x$p_value })
which(p_value < 0.05)

# Transform using symmetric distribution
hist(qnorm(p_value), breaks = 15)   # probit transform


i <- 2
X1 <- X[y == 0, , i]
X2 <- X[y == 1, , i]
perm_stat(X1, X2)






