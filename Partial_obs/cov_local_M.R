

get_kernel <- function(t_vec, t, bw, method = "epan") {
  u <- (t_vec - t) / bw
  kern <- (1-u^2) * (3/4)
  
  return(kern)
}


get_raw_cov <- function(x, 
                        mu = NULL, 
                        gr,
                        diag = FALSE,
                        engine = "C++") {
  # gr <- seq(0, 1, length.out = 51)
  n <- nrow(x)
  p <- ncol(x)
  
  if (is.null(mu)) {
    mu <- mean_Mest(x)
  }
  
  if (engine == "C++") {
    raw_cov <- get_raw_cov_cpp(x, mu, gr, diag)
    raw_cov <- as.data.frame(raw_cov)
  } else if (engine == "R") {
    for (i in 1:n) {
      ind <- which(!is.na(x[i, ]))
      gr_sub <- gr[ind]
      x_sub <- x[i, ind] - mu[ind]
      
      exp_grid <- expand.grid(gr_sub, gr_sub)
      exp_x <- expand.grid(x_sub, x_sub)
      
      if (diag == FALSE) {
        # eliminate diagonal parts
        idx <- which(exp_grid[, 1] == exp_grid[, 2])
        cov_i <- exp_x[-idx, 1] * exp_x[-idx, 2]
        
        if (i == 1) {
          raw_cov <- cbind(cov_i,
                           exp_grid[-idx, ])
        } else {
          raw_cov <- rbind(raw_cov,
                           cbind(cov_i,
                                 exp_grid[-idx, ]))
        }
      } else {
        cov_i <- exp_x[, 1] * exp_x[, 2]
        
        if (i == 1) {
          raw_cov <- cbind(cov_i,
                           exp_grid)
        } else {
          raw_cov <- rbind(raw_cov,
                           cbind(cov_i,
                                 exp_grid))
        }
      }
    }
  }

  colnames(raw_cov) <- c("cov","s","t")
  
  return(raw_cov)  
}


### Local covariance
cov_local_M <- function(x, 
                        h = 0.02,
                        grid = NULL,
                        diag = FALSE,
                        engine = "R") {
  if (is.null(grid)) {
    gr <- seq(0, 1, length.out = ncol(x))
  } else {
    gr <- grid
  }
  p <- length(gr)
  
  # raw covariance
  raw_cov <- get_raw_cov(x, 
                         mu = NULL, 
                         gr = gr,
                         diag = diag,
                         engine = engine)
  s <- raw_cov$s
  t <- raw_cov$t
  raw_cov <- raw_cov$cov
  
  if (engine == "C++") {
    cov_mat <- cov_local_M_cpp(raw_cov = raw_cov,
                               s = s,
                               t = t,
                               gr = gr,
                               h = h)
  } else if (engine == "R") {
    cov_mat <- matrix(NA, p, p)
    for (i in 1:p) {
      i_neighbor <- which(s < gr[i] + h & s > gr[i] - h)
      # kern_s <- get_kernel(s, gr[i], h)
      # i_neighbor <- which(kern_s > 0)
      for (j in 1:p) {
        if (i <= j) {
          j_neighbor <- which(t < gr[j] + h & t > gr[j] - h)
          
          ind <- intersect(i_neighbor, j_neighbor)
          
          cov_mat[i, j] <- MASS::huber(raw_cov[ind])$mu
        } else {
          cov_mat[i, j] <- cov_mat[j, i]  
        }
        
        # kern_t <- get_kernel(t, gr[j], h)
        # j_neighbor <- which(kern_t > 0)
        # 
        # ind <- intersect(i_neighbor, j_neighbor)
        # # kern_weig <- kern_s[ind] * kern_t[ind]
        # # kern_weig <- kern_weig / max(kern_weig)
        # # 
        # # cov_weig <- kern_weig * raw_cov$cov[ind]
        # cov_weig <- raw_cov$cov[ind]
        # 
        # cov_mat[i, j] <- MASS::huber(cov_weig)$mu
      }
    }
  }

  # # make positive definite
  # eig <- get_eigen(cov_mat, gr)
  # lambda <- eig$lambda
  # idx_positive <- which(lambda > 0)
  # lambda <- lambda[idx_positive]
  # phi <- eig$phi[, idx_positive]
  # cov_mat <- phi %*% diag(lambda, ncol = length(idx_positive)) %*% t(phi)
  
  return(cov_mat)  
}


# system.time({
#   cc <- cov_local_M(x, h = 0.05)  
# })
# dim(cc)
# 
# 
# gr <- work.grid
# GA::persp3D(gr, gr, cc,
#             theta = -70, phi = 30, expand = 1)
# 
# par(mfrow = c(2, 2))
# bw <- c(0.03, 0.05, 0.1, 0.2)
# for (i in 1:4) {
#   cc <- cov_local_M(x, h = bw[i])  
#   GA::persp3D(gr, gr, cc,
#               theta = -70, phi = 30, expand = 1)  
# }
# 
# 
# cc <- cov_local_M(x, h = 0.05)  
# 
# eig.true <- get_delaigle_eigen(work.grid, model = 2)
# # eig.true <- get_eigen(cov.true, gr)$phi[, 1:4]
# eig.huber <- get_eigen(cov.huber, gr)$phi[, 1:4]
# eig.Mest <- get_eigen(cov.Mest, gr)$phi[, 1:4]
# eig.Mest.sm <- get_eigen(cov.Mest.sm, gr)$phi[, 1:4]
# eig.cc.sm <- get_eigen(cc, gr)$phi[, 1:4]
# par(mfrow = c(2, 2))
# for (k in 1:4) {
#   matplot(work.grid,
#           cbind(
#             eig.true[, k],
#             check_eigen_sign(eig.huber, eig.true)[, k],
#             check_eigen_sign(eig.Mest, eig.true)[, k],
#             check_eigen_sign(eig.Mest.sm, eig.true)[, k],
#             check_eigen_sign(eig.cc.sm, eig.true)[, k]
#           ),
#           type = "l",
#           col = 1:5,
#           lty = 1:5,
#           main = paste("Eigenfunction", k),
#           xlab = "", ylab = "",
#           lwd = rep(2, 4))
#   if (k == 1) {
#     legend("bottomleft",
#            c("True","Huber","M-est","M-est(smooth)","cc_sm"),
#            col = 1:7,
#            lty = 1:7,
#            lwd = rep(2, 7))
#   }
# }







cv.cov_local_M <- function(x,
                           bw_cand = NULL,
                           ncores = 1) {
  if (is.list(x)) {
    gr <- sort(unique(unlist(x$Lt)))
    x <- list2matrix(x)
  } else {
    gr <- seq(0, 1, length.out = ncol(x))
  }
  
  n <- nrow(x)
  p <- ncol(x)
  
  mu <- mean_Mest(x)
  
  # bandwidth candidates
  if (is.null(bw_cand)) {
    a <- min(gr)
    b <- max(gr)
    bw_cand <- seq(min(diff(gr))*1.5, (b-a)/3, length.out = 5)
    # bw_cand <- 10^seq(-2, 0, length.out = 10) * (b - a)/3
  }
  
  
  # Leave-one-curve-out cross-validation
  cv_error <- matrix(0, n, length(bw_cand))
  for (i in 1:n) {
    print(i)
    # obtain the raw covariance of ith curve
    curve_i <- matrix(x[i, ], nrow = 1)
    ind_not_na <- which(!is.na(x[i, ]))
    cov_raw <- get_raw_cov(curve_i, mu, diag = TRUE)
    cov_raw <- matrix(cov_raw$cov, ncol = length(ind_not_na))
    diag(cov_raw) <- 0
    
    for (j in 1:length(bw_cand)) {
      # smoothed covariance without ith curve
      cov_sm <- cov_local_M(x[-i, ], 
                            h = bw_cand[j],
                            grid = gr[ind_not_na])
      diag(cov_sm) <- 0
      cv_error[i, j] <- mean((cov_sm - cov_raw)^2)
    }
  }
  
  # trimmed mean where cv_error <= 0.75 quantile
  cv_error <- apply(cv_error, 2, function(cv_error_bw) {
    cutoff_trimmed <- quantile(cv_error_bw, 0.75)
    return( mean(cv_error_bw[cv_error_bw < cutoff_trimmed]) )
  })
  
  bw <- list(selected_bw = bw_cand[ which.min(cv_error) ],
             cv.error = data.frame(bw = bw_cand,
                                   error = cv_error))
  
  return(bw)
  
  
  # # element-wise covariances
  # st <- expand.grid(gr, gr)
  # cov_st <- as.numeric(cov_hat)
  # 
  # # remove diagonal parts from raw covariance (See Yao et al.(2005))
  # ind <- which(st[, 1] == st[, 2], arr.ind = T)
  # st <- st[-ind, ]
  # cov_st <- cov_st[-ind]
  # 
  # # get index for each folds
  # n <- length(st)
  # folds <- list()
  # fold_num <- n %/% K   # the number of curves for each folds
  # fold_sort <- sample(1:n, n)
  # for (k in 1:K) {
  #   ind <- (fold_num*(k-1)+1):(fold_num*k)
  #   if (k == K) {
  #     ind <- (fold_num*(k-1)+1):n
  #   }
  #   folds[[k]] <- fold_sort[ind]
  # }
  # 
  # # K-fold cross validation
  # if (ncores > 1) {
  #   # Parallel computing setting
  #   if (ncores > parallel::detectCores()) {
  #     ncores <- parallel::detectCores() - 1
  #     warning(paste0("ncores is too large. We now use ", ncores, " cores."))
  #   }
  #   # ncores <- detectCores() - 3
  #   cl <- parallel::makeCluster(ncores)
  #   doParallel::registerDoParallel(cl)
  #   
  #   # matrix of bw_cand and fold
  #   bw_fold_mat <- data.frame(bw_cand = rep(bw_cand, each = K),
  #                             fold = rep(1:K, length(bw_cand)))
  #   
  #   cv_error <- foreach::foreach(i = 1:nrow(bw_fold_mat),
  #                                .combine = "c",
  #                                # .export = c("local_kern_smooth"),
  #                                .packages = c("robfpca"),
  #                                .errorhandling = "pass") %dopar% {
  #                                  
  #                                  bw <- bw_fold_mat$bw_cand[i]   # bandwidth candidate
  #                                  k <- bw_fold_mat$fold[i]   # fold for K-fold CV
  #                                  
  #                                  # data of kth fold
  #                                  st_train <- st[-folds[[k]], ]
  #                                  st_test <- st[folds[[k]], ]
  #                                  cov_train <- cov_st[-folds[[k]]]
  #                                  cov_test <- cov_st[folds[[k]]]
  #                                  
  #                                  # Bivariate Nadaraya-Watson smoothing
  #                                  cov_hat_sm <- fields::smooth.2d(cov_train,
  #                                                                  x = st_train, 
  #                                                                  surface = F,
  #                                                                  theta = bw, 
  #                                                                  nrow = p, 
  #                                                                  ncol = p)
  #                                  cov_hat_sm <- as.numeric(cov_hat_sm)[folds[[k]]]
  #                                  err <- sum((cov_test - cov_hat_sm)^2)   # squared errors
  #                                  
  #                                  return(err)
  #                                }
  #   parallel::stopCluster(cl)
  #   
  #   bw_fold_mat$cv_error <- cv_error
  #   cv_obj <- bw_fold_mat %>%
  #     dplyr::group_by(bw_cand) %>%
  #     dplyr::summarise(cv_error = sum(cv_error))
  #   
  #   bw <- list(selected_bw = cv_obj$bw_cand[ which.min(cv_obj$cv_error) ],
  #              cv.error = as.data.frame(cv_obj))
  # } else {
  #   cv_error <- rep(0, length(bw_cand))
  #   for (k in 1:K) {
  #     # data of kth fold
  #     st_train <- st[-folds[[k]], ]
  #     st_test <- st[folds[[k]], ]
  #     cov_train <- cov_st[-folds[[k]]]
  #     cov_test <- cov_st[folds[[k]]]
  #     
  #     for (i in 1:length(bw_cand)) {
  #       # Bivariate Nadaraya-Watson smoothing
  #       cov_hat_sm <- fields::smooth.2d(cov_train,
  #                                       x = st_train, 
  #                                       surface = F,
  #                                       theta = bw_cand[i], 
  #                                       nrow = p, 
  #                                       ncol = p)
  #       cov_hat_sm <- as.numeric(cov_hat_sm)[folds[[k]]]
  #       cv_error[i] <- cv_error[i] + sum((cov_test - cov_hat_sm)^2)   # squared errors
  #     }
  #   }
  #   
  #   bw <- list(selected_bw = bw_cand[ which.min(cv_error) ],
  #              cv.error = data.frame(bw = bw_cand,
  #                                    error = cv_error))
  # }
  # 
  # return(bw)
}


# cv.cov_local_M(x)
# 
# system.time({
#   a1 <- cov_local_M(x, 0.05, diag = F, engine = "C++")
# })
# system.time({
#   a2 <- cov_local_M(x, 0.05, diag = F, engine = "R")
# })
# all.equal(a1, a2)
# 
# par(mfrow = c(1, 2))
# GA::persp3D(gr, gr, a1,
#             theta = -70, phi = 30, expand = 1)
# GA::persp3D(gr, gr, a2,
#             theta = -70, phi = 30, expand = 1)
