# Rcpp::sourceCpp("../../robfpca/src/cov_Mest.cpp")
# source("../../robfpca/R/cov_Mest.R")

eig.obj <- get_eigen(cov.Mest, work.grid)
eig.obj.sm <- get_eigen(cov.Mest.sm.noise, work.grid)

phi_sm <- apply(eig.obj$phi, 2, function(phi){
  # pspline::sm.spline(work.grid, phi)$ysmth
  stats::smooth.spline(work.grid, phi)$y
})
cov_2 <- phi_sm %*% diag(eig.obj$lambda) %*% t(phi_sm)
eig.obj2 <- get_eigen(cov_2, work.grid)
phi_sm <- eig.obj2$phi[, 1:4]


par(mfrow = c(2, 2))
matplot(work.grid, get_delaigle_eigen(work.grid), type = "l")
matplot(work.grid, eig.obj$phi[, 1:4], type = "l")
matplot(work.grid, eig.obj.sm$phi[, 1:4], type = "l")
matplot(work.grid, phi_sm, type = "l")

all.equal(eig.obj$phi[, 1:4], phi_sm)


bw.cov_Mest <- function(x,
                        bw_cand = NULL,
                        K = 5,
                        ncores = 1,
                        noise.var = 0) {
  
  if (is.list(x)) {
    gr <- sort(unique(unlist(x$Lt)))
    x <- list2matrix(x)
  } else {
    gr <- seq(0, 1, length.out = ncol(x))
  }
  
  n <- nrow(x)
  p <- ncol(x)
  
  # bandwidth candidates
  if (is.null(bw_cand)) {
    a <- min(gr)
    b <- max(gr)
    bw_cand <- 10^seq(-2, 0, length.out = 10) * (b - a)/3
  }
  
  # get index for each folds
  folds <- list()
  fold_num <- n %/% K   # the number of curves for each folds
  fold_sort <- sample(1:n, n)
  for (k in 1:K) {
    ind <- (fold_num*(k-1)+1):(fold_num*k)
    if (k == K) {
      ind <- (fold_num*(k-1)+1):n
    }
    folds[[k]] <- fold_sort[ind]
  }
  
  # K-fold cross validation
  if (ncores > 1) {
    # Parallel computing setting
    if (ncores > parallel::detectCores()) {
      ncores <- parallel::detectCores() - 1
      warning(paste0("ncores is too large. We now use ", ncores, " cores."))
    }
    # ncores <- detectCores() - 3
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    
    # matrix of bw_cand and fold
    bw_fold_mat <- data.frame(bw_cand = rep(bw_cand, each = K),
                              fold = rep(1:K, length(bw_cand)))
    
    cv_error <- foreach::foreach(i = 1:nrow(bw_fold_mat),
                                 .combine = "c",
                                 # .export = c("local_kern_smooth"),
                                 .packages = c("robfpca"),
                                 .errorhandling = "pass") %dopar% {
                                   
      bw <- bw_fold_mat$bw_cand[i]   # bandwidth candidate
      k <- bw_fold_mat$fold[i]   # fold for K-fold CV
      
      # data of kth fold
      x_train <- x[-folds[[k]], ]
      x_test <- x[folds[[k]], ]
      
      # obtain the covariance based on marignal M-estimator via C++ code
      cov_hat = cov_Mest_cpp(x)
      
      # subtract noise variance
      diag(cov_hat) <- diag(cov_hat) - noise.var
      
      # Bivariate Nadaraya-Watson smoothing
      cov_hat_sm <- fields::smooth.2d(as.numeric(cov_hat),
                                      x = expand.grid(gr, gr), 
                                      surface = F,
                                      theta = bw, 
                                      nrow = p, 
                                      ncol = p)
      err <- sum((cov_hat - cov_hat_sm)^2)   # squared errors
      
      return(err)
    }
    parallel::stopCluster(cl)
    
    bw_fold_mat$cv_error <- cv_error
    cv_obj <- bw_fold_mat %>%
      dplyr::group_by(bw_cand) %>%
      dplyr::summarise(cv_error = sum(cv_error))
    
    bw <- list(selected_bw = cv_obj$bw_cand[ which.min(cv_obj$cv_error) ],
               cv.error = as.data.frame(cv_obj))
  } else {
    cv_error <- rep(0, length(bw_cand))
    for (k in 1:K) {
      # data of kth fold
      x_train <- x[-folds[[k]], ]
      x_test <- x[folds[[k]], ]
      
      for (i in 1:length(bw_cand)) {
        # obtain the covariance based on marignal M-estimator via C++ code
        cov_hat = cov_Mest_cpp(x)
        
        # subtract noise variance
        diag(cov_hat) <- diag(cov_hat) - noise.var
        
        # Bivariate Nadaraya-Watson smoothing
        cov_hat_sm <- fields::smooth.2d(as.numeric(cov_hat),
                                        x = expand.grid(gr, gr), 
                                        surface = F,
                                        theta = bw_cand[i], 
                                        nrow = p, 
                                        ncol = p)
        cv_error[i] <- cv_error[i] + sum((cov_hat - cov_hat_sm)^2)   # squared errors
      }
    }
    
    bw <- list(selected_bw = bw_cand[ which.min(cv_error) ],
               cv.error = data.frame(bw = bw_cand,
                                     error = cv_error))
  }
  
  return(bw)
}



### Example
# system.time({
#   set.seed(1000)
#   cv.obj <- bw.cov_Mest(x)
# })
# cv.obj

