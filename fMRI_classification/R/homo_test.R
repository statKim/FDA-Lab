Rcpp::sourceCpp("src/homo_test.cpp")
library(fda.usc)
library(ddalpha)
library(doParallel)
library(foreach)

# Remove the parallel backends parameter
unregister <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}


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
fun_homo_test <- function(X1, X2, method = "ddplot", 
                          depth = c("FM","RPD","RT","ABD","hM","fd2"), trim = 0.2, 
                          alpha = 0.05, B = 500) {
  
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
    gr <- seq(0, 1, length.out = m)
    fdataobj <- fdata(rbind(X1, X2), gr)
    fdataobj_X1 <- fdata(X1, gr)
    fdataobj_X2 <- fdata(X2, gr)
    
    # Function for computing depth
    get_depth <- function(fdataobj, fdataobj_X1, fdataobj_X2, 
                          depth, trim, 
                          group, idx_g1, idx_g2) {
      if (depth %in% c("FM","RPD","RT")) {
        if (depth == "FM") {
          depth <- depth.FM
        } else if (depth == "RPD") {
          depth <- depth.RPD
        } else if (depth == "RT") {
          depth <- depth.RT
        }
        
        depth.obj <- depth(fdataobj, fdataobj_X1, trim = trim)
        depth_x <- depth.obj$dep
        depth.obj <- depth(fdataobj, fdataobj_X2, trim = trim)
        depth_y <- depth.obj$dep
      } else if (depth %in% c("ABD","hM","fd2")) {
        datafobj <- dataf(fdataobj, group)$dataf
        datafobj_X1 <- dataf(fdataobj_X1, group[idx_g1])$dataf
        datafobj_X2 <- dataf(fdataobj_X2, group[idx_g2])$dataf
        
        if (depth == "fd2") {
          datafobj <- derivatives.est(datafobj, deriv = c(0,1))
          datafobj_X1 <- derivatives.est(datafobj_X1, deriv = c(0,1))
          datafobj_X2 <- derivatives.est(datafobj_X2, deriv = c(0,1))
          depth_x <- depthf.fd2(datafobj, datafobj_X1)$Half_FD
          depth_y <- depthf.fd2(datafobj, datafobj_X2)$Half_FD
        } else if (depth == "hM") {
          depth_x <- depthf.hM(datafobj, datafobj_X1, norm = "L2", q = trim)
          depth_y <- depthf.hM(datafobj, datafobj_X2, norm = "L2", q = trim)
        } else if (depth == "ABD") {
          depth_x <- depthf.ABD(datafobj, datafobj_X1, norm = "L2")
          depth_y <- depthf.ABD(datafobj, datafobj_X2, norm = "L2")
        }
      }
      
      return(list(depth_x = depth_x,
                  depth_y = depth_y))
    }
    
    # Compute depth
    depth_obj <- get_depth(fdataobj, fdataobj_X1, fdataobj_X2, 
                           depth, trim, 
                           group, idx_g1, idx_g2)
    depth_x <- depth_obj$depth_x
    depth_y <- depth_obj$depth_y
    
    # # DD-plot
    # plot(depth_x, depth_y)
    # abline(a = 0, b = 1, col = 2, lwd = 2)
    
    # Least square estimate
    fit <- lm(depth_x ~ depth_y)
    fit_summary <- summary(fit)
    
    # Test statistics
    T_value <- (fit_summary$coefficients[, 1] - c(0, 1)) / fit_summary$coefficients[, 2]
    
    # sigma_beta0 <- sqrt( sum(residuals(fit)^2) * sum(depth_y)^2 / ((n-2) * sum((depth_y - mean(depth_y))^2)) )
    # sigma_beta1 <- sqrt( sum(residuals(fit)^2) / ((n-2) * sum((depth_y - mean(depth_y))^2)) )
    # 
    # T0 <- fit_summary$coefficients[1, 1] / sigma_beta0
    # T1 <- (fit_summary$coefficients[2, 1] - 1) / sigma_beta1
    
    # Bootstrap under H_0: X =_d Y
    idx_boot_mat <- sample(1:n, n*B, replace = T)
    idx_boot_mat <- matrix(idx_boot_mat, n, B)
    T_star <- foreach(i=1:B, 
                      .packages=c("ddalpha","fda.usc"), 
                      .export=c("get_depth"), 
                      .combine = rbind) %dopar% {
       # Bootstrap
       idx_boot <- idx_boot_mat[, i]
       # idx_boot <- sample(1:n, n, replace = T)
       group_boot <- group[idx_boot]
       idx_g1_boot <- idx_boot[idx_g1]
       idx_g2_boot <- idx_boot[idx_g2]
       X1_boot <- X[idx_g1_boot, ]
       X2_boot <- X[idx_g2_boot, ]
       
       # Compute depth
       fdataobj <- fdata(rbind(X1_boot, X2_boot), gr)
       fdataobj_X1 <- fdata(X1_boot, gr)
       fdataobj_X2 <- fdata(X2_boot, gr)
       
       # Compute depth
       depth_obj <- get_depth(fdataobj, fdataobj_X1, fdataobj_X2, 
                              depth, trim, 
                              group_boot, idx_g1_boot, idx_g2_boot)
       depth_x <- depth_obj$depth_x
       depth_y <- depth_obj$depth_y
       
       # Least square estimates
       fit <- lm(depth_x ~ depth_y)
       fit_summary <- summary(fit)
       
       # Test statistics
       test_stat <- (fit_summary$coefficients[, 1] - c(0, 1)) / fit_summary$coefficients[, 2]
       
       return(test_stat)
     }
    
    # T_star <- matrix(0, B, 2)
    # for (b in 1:B) {
    #   # Bootstrap
    #   idx_boot <- sample(1:n, n, replace = T)
    #   group_boot <- group[idx_boot]
    #   idx_g1_boot <- idx_boot[idx_g1]
    #   idx_g2_boot <- idx_boot[idx_g2]
    #   X1_boot <- X[idx_g1_boot, ]
    #   X2_boot <- X[idx_g2_boot, ]
    #   
    #   # X1_boot <- X1[sample(1:n1, n1, replace = T), ]
    #   # X2_boot <- X2[sample(1:n2, n2, replace = T), ]
    #   
    #   # Compute depth
    #   fdataobj <- fdata(rbind(X1_boot, X2_boot), gr)
    #   fdataobj_X1 <- fdata(X1_boot, gr)
    #   fdataobj_X2 <- fdata(X2_boot, gr)
    # 
    #   # depth.obj <- depth(fdataobj, fdataobj_X1, trim = trim)
    #   # depth_x <- depth.obj$dep
    #   # depth.obj <- depth(fdataobj, fdataobj_X2, trim = trim)
    #   # depth_y <- depth.obj$dep
    #   
    #   # Compute depth
    #   depth_obj <- get_depth(fdataobj, fdataobj_X1, fdataobj_X2, 
    #                          depth, trim, 
    #                          group_boot, idx_g1_boot, idx_g2_boot)
    #   depth_x <- depth_obj$depth_x
    #   depth_y <- depth_obj$depth_y
    #   
    #   # Least square estimates
    #   fit <- lm(depth_x ~ depth_y)
    #   fit_summary <- summary(fit)
    #   
    #   # Test statistics
    #   T_star[b, ] <- (fit_summary$coefficients[, 1] - c(0, 1)) / fit_summary$coefficients[, 2]
    # }
    
    p_value <- rbind(
      rowMeans(t(T_star) > T_value),
      rowMeans(t(T_star) < T_value)
    )
    p_value <- 2 * apply(p_value, 2, min)
    
    adj_p_value <- min(2*min(p_value), max(p_value))
    
    res <- list(
      p_value = adj_p_value,
      p_value_list = p_value,
      test_stat_boot = T_star,
      test_stat = T_value
    )
  }
  
  return(res)
}


# Find the cutoff of mirror statistics controlling the FDR
mirror_stat_cutoff <- function(mirror_stat, q = 0.1, est = FALSE) {
  # Cutoff candidates
  cand <- sort(abs(mirror_stat), decreasing = T)   # sorted abs(mirror stat)
  cand <- (cand[-1] + cand[-(length(cand))]) / 2   # mean values of paired 'cand'
  fdp <- rep(1, length(cand))
  for (i in 1:length(cand)){
    if (is.infinite(cand[i])) {
      next
    }
    
    # Compute FDP
    fdp[i] <- sum(mirror_stat < -cand[i]) / max(1, sum(mirror_stat > cand[i]))
    # fdp[i] <- sum(mirror_stat <= -cand[i]) / max(1, sum(mirror_stat >= cand[i]))
    
    # # Check FDP < fdr
    # if (fdp[i] <= fdr) {
    # # if (fdp[i] <= fdr & fdp[i] > 0) {
    #   cutoff <- mean(cand[i], cand[i-1])
    #   if (cutoff < 0) {
    #     cutoff <- -cutoff
    #   }
    #   break
    # }
  }
  
  # print(fdp)
  
  # Find cutoff
  cutoff_idx_list <- which(fdp <= q)
  # print(cutoff_idx_list)
  if (length(cutoff_idx_list) == 0) {
    warning("Cannot control the FDR using the input 'mirror_stat' and 'fdr' level!")
    cutoff <- Inf  # Refered by "Knockoffs" paper (2015)
  } else {
    cutoff_idx <- max(cutoff_idx_list)   # index of the minimum value
    # cutoff <- mean(c(cand[cutoff_idx], cand[cutoff_idx-1]))
    cutoff <- cand[cutoff_idx]
  }
  
  # if (sum(fdp <= fdr) == 0) {
  #   stop("Cannot control the FDR using the input 'mirror_stat' and 'fdr' level!")
  # }
  
  if (est == TRUE) {
    return(list(cutoff = cutoff,
                fdp = fdp))
  } else {
    return(cutoff)  
  }
}

