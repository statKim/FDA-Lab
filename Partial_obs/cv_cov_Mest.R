

### K-fold cross-validation for Bivariate Nadaraya-Watson smoothing for Covariance
### - Not exactly observation-wise cross-validation
### - It is conducted for element-wise covariance
cv.cov_Mest <- function(x,
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
  
  # n <- nrow(x)
  p <- ncol(x)
  
  # bandwidth candidates
  if (is.null(bw_cand)) {
    a <- min(gr)
    b <- max(gr)
    bw_cand <- seq(min(diff(gr)), (b-a)/3, length.out = 10)
    # bw_cand <- 10^seq(-2, 0, length.out = 10) * (b - a)/3
  }
  
  # obtain the raw covariance (Not smoothed)
  cov_hat <- cov_Mest(x,
                      smooth = FALSE,
                      make.pos.semidef = FALSE)
  
  # subtract noise variance
  diag(cov_hat) <- diag(cov_hat) - noise.var
  
  # element-wise covariances
  st <- expand.grid(gr, gr)
  cov_st <- as.numeric(cov_hat)
  
  # remove diagonal parts from raw covariance (See Yao et al.(2005))
  ind <- which(st[, 1] == st[, 2], arr.ind = T)
  st <- st[-ind, ]
  cov_st <- cov_st[-ind]
  
  # get index for each folds
  n <- length(st)
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
      st_train <- st[-folds[[k]], ]
      st_test <- st[folds[[k]], ]
      cov_train <- cov_st[-folds[[k]]]
      cov_test <- cov_st[folds[[k]]]
      
      # Bivariate Nadaraya-Watson smoothing
      cov_hat_sm <- fields::smooth.2d(cov_train,
                                      x = st_train, 
                                      surface = F,
                                      theta = bw, 
                                      nrow = p, 
                                      ncol = p)
      cov_hat_sm <- as.numeric(cov_hat_sm)[folds[[k]]]
      err <- sum((cov_test - cov_hat_sm)^2)   # squared errors
      
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
      st_train <- st[-folds[[k]], ]
      st_test <- st[folds[[k]], ]
      cov_train <- cov_st[-folds[[k]]]
      cov_test <- cov_st[folds[[k]]]
      
      for (i in 1:length(bw_cand)) {
        # Bivariate Nadaraya-Watson smoothing
        cov_hat_sm <- fields::smooth.2d(cov_train,
                                        x = st_train, 
                                        surface = F,
                                        theta = bw_cand[i], 
                                        nrow = p, 
                                        ncol = p)
        cov_hat_sm <- as.numeric(cov_hat_sm)[folds[[k]]]
        cv_error[i] <- cv_error[i] + sum((cov_test - cov_hat_sm)^2)   # squared errors
      }
    }
    
    bw <- list(selected_bw = bw_cand[ which.min(cv_error) ],
               cv.error = data.frame(bw = bw_cand,
                                     error = cv_error))
  }
  
  return(bw)
}

# parallel 안한 것이 더 빠름...
set.seed(100)
system.time({
  obj <- cv.cov_Mest(x, ncores = 1, noise.var = noise_var)
})
obj








par(mfrow = c(3, 3))
for (i in 1:9) {
  plot(cov.Mest[i*5, (i*5):51], type = "o")  
}

cc <- matrix(0, 51, 51)
for (i in 1:51) {
  if (i > 51-3) {
    cc[i, i:51] <- cov.Mest[i, i:51]
  } else {
    cc[i, i:51] <- smooth.spline(gr[i:51], cov.Mest[i, i:51])$y
  }
}

for (j in 4:51) {
  cc[1:j, j] <- smooth.spline(gr[1:j], cc[1:j, j])$y
}
cc <- (cc + t(cc)) / 2

par(mfrow = c(2, 2))
GA::persp3D(gr, gr, cov.true,
            theta = -70, phi = 30, expand = 1)
GA::persp3D(gr, gr, cov.Mest,
            theta = -70, phi = 30, expand = 1)
GA::persp3D(gr, gr, cov.Mest.sm,
            theta = -70, phi = 30, expand = 1)
GA::persp3D(gr, gr, cc,
            theta = -70, phi = 30, expand = 1)

eig <- get_eigen(cc - diag(noise_var, 51), work.grid)
eig$PVE


round(eigen(cov.Mest)$values, 5)
cc <- cov.Mest
diag(cc) <- diag(cc) + 0.0001
round(eigen(cc)$values, 5)


cc <- cov.Mest.sm
diag(cc) <- smooth.spline(gr, diag(cov.Mest))$y
par(mfrow = c(1, 2))
GA::persp3D(gr, gr, cov.Mest.sm,
            theta = -70, phi = 30, expand = 1)
GA::persp3D(gr, gr, cc,
            theta = -70, phi = 30, expand = 1)
round(eigen(cov.Mest.sm)$values, 3)
round(eigen(cc)$values, 3)

get_eigen(cov.Mest.sm, work.grid)$PVE
get_eigen(cc, work.grid)$PVE


mean((diag(cov.Mest.sm) - diag(cc))^2)

lambda.true <- get_eigen(cov.true, work.grid)$lambda
round(lambda.true, 3)
lambda <- get_eigen(cov.Mest, work.grid)$lambda
round(lambda, 3)

beta <- 0.9
lambda.shirink <- beta*lambda + (1-beta)*mean(lambda)
cumsum(lambda.shirink) / sum(lambda.shirink)


par(mfrow = c(2, 2))
eig.true <- get_delaigle_eigen(work.grid)
# cc <- robfpca::cov_Mest(x, smooth = F, make.pos.semidef = T)
# GA::persp3D(gr, gr, cc,
#             theta = -70, phi = 30, expand = 1)
# eig <- get_eigen(cc, work.grid)
# print(eig$PVE)
bw <- c(0.02, 0.05, 0.1, 0.2)
# plot(diag(cc), type = "o")
for (i in 1:4) {
  cc <- robfpca::cov_Mest(x, smooth = T, bw = bw[i], make.pos.semidef = F)
  GA::persp3D(gr, gr, cc,
              theta = -70, phi = 30, expand = 1)
  eig <- get_eigen(cc, work.grid)
  print(bw[i])
  print(eig$PVE)
  print(mean((eig.true - eig$phi[, 1:4])^2))
  # lines(diag(cc), type = "l", col = i+1)
}




x.3 <- x.2
for (i in 1:100) {
  x.3$Ly[[i]] <- pspline_curve(x.2$Lt[[i]], x.2$Ly[[i]])  
}
xx <- list2matrix(x.3)
matplot(t(x[81:100, ]), type = "l")
matplot(t(xx[81:100, ]), type = "l", ylim = c(-5, 5))
matplot(t(xx[1:80, ]), type = "l", ylim = c(-5, 5))


par(mfrow = c(1, 2))
plot(gr, diag(cov.Mest))
lines(gr, diag(cov.Mest))
lines(gr, smooth.spline(gr, diag(cov.Mest))$y, col = 2, lwd = 2)
plot(gr, diag(cov.Mest.sm))
lines(gr, diag(cov.Mest.sm))
lines(gr, smooth.spline(gr, diag(cov.Mest.sm))$y, col = 2, lwd = 2)

cc <- cov.Mest
cc.sm <- smooth.spline(gr, diag(cov.Mest))$y
diag(cc)

for (i in 1:51) {
  cc[i, ] <- cc[i, ] * sqrt(cc.sm[i]^2) / sqrt(cov.Mest[i, i]^2)
}
for (j in 1:51) {
  cc[, j] <- cc[, j] * sqrt(cc.sm[j]^2) / sqrt(cov.Mest[j, j]^2)
}
diag(cc)
cc.sm

par(mfrow = c(1, 2))
GA::persp3D(gr, gr, cov.Mest,
            theta = -70, phi = 30, expand = 1)
GA::persp3D(gr, gr, cc,
            theta = -70, phi = 30, expand = 1)

eig <- get_eigen(cc, work.grid)
eig$PVE





gr <- work.grid
par(mfrow = c(2, 3))
cov.true <- get_cov_fragm(gr)
GA::persp3D(gr, gr, cov_Mest(x,
                             make.pos.semidef = FALSE),
            theta = -70, phi = 30, expand = 1)
GA::persp3D(gr, gr, cov_Mest(x,
                             bw = 0.003,
                             smooth = T,
                             make.pos.semidef = FALSE),
            theta = -70, phi = 30, expand = 1)
GA::persp3D(gr, gr, cov_Mest(x,
                             bw = 0.01,
                             smooth = T,
                             make.pos.semidef = FALSE),
            theta = -70, phi = 30, expand = 1)
GA::persp3D(gr, gr, cov_Mest(x,
                             bw = 0.1,
                             smooth = T,
                             make.pos.semidef = FALSE),
            theta = -70, phi = 30, expand = 1)
GA::persp3D(gr, gr, cov_Mest(x,
                             bw = 0.3,
                             smooth = T,
                             make.pos.semidef = FALSE),
            theta = -70, phi = 30, expand = 1)
GA::persp3D(gr, gr, cov_Mest(x,
                             bw = obj$selected_bw,
                             smooth = T,
                             make.pos.semidef = FALSE),
            theta = -70, phi = 30, expand = 1)
### Example
# system.time({
#   set.seed(1000)
#   cv.obj <- bw.cov_Mest(x)
# })
# cv.obj

cov.Mest <- cov_Mest(x, make.pos.semidef = F)
cc <- as.numeric(cov.Mest)
xx <- expand.grid(gr, gr)
ind <- which(xx[, 1] == xx[, 2], arr.ind = T)
cov_hat_sm <- fields::smooth.2d(cc[-ind],
                                x = xx[-ind, ], 
                                surface = F,
                                theta = 0.04, 
                                nrow = 51, 
                                ncol = 51)
diag(cov_hat_sm)

par(mfrow = c(2, 2))
GA::persp3D(gr, gr, cov.Mest,
            theta = -70, phi = 30, expand = 1)
GA::persp3D(gr, gr, cov_hat_sm,
            theta = -70, phi = 30, expand = 1)
diag(cov_hat_sm) <- diag(cov.Mest)
GA::persp3D(gr, gr, cov_hat_sm,
            theta = -70, phi = 30, expand = 1)

rob.var <- cov_hat_sm
diag(rob.var) <- diag(rob.var) - noise_var
eig <- eigen(rob.var)
k <- which(eig$values > 0)
lambda <- eig$values[k]
phi <- matrix(eig$vectors[, k],
              ncol = length(k))
rob.var <- phi %*% diag(lambda, ncol = length(k)) %*% t(phi)
GA::persp3D(gr, gr, rob.var,
            theta = -70, phi = 30, expand = 1)




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




plot(mean_Mest(x), type = "l")
lines(mean_Mest(x, smooth = TRUE), col = 2)
lines(smooth.spline(mean_Mest(x))$y, col = 2)

gr <- work.grid
yy <- gr^2 + rnorm(51, 0, 0.1)
plot(gr, yy)

fit <- smooth.spline(gr, yy)
fit$y
lines(gr, fit$y, lwd = 2)



eig.obj <- eigen(cov.Mest)
lambda <- eig.obj$values
lambda[lambda <= 0] <- 0
pve <- cumsum(lambda) / sum(lambda)
# k <- sum(lambda > 0)
k <- which(pve > 0.9)[1]
phi <- eig.obj$vectors
phi_sm <- phi
phi_sm[, 1:k] <- apply(phi[, 1:k], 2, function(x){ smooth.spline(gr, x)$y })
# phi_sm <- pracma::gramSchmidt(phi_sm)$Q
phi_sm[, 1:k] <- pracma::gramSchmidt(phi_sm[, 1:k])$Q

matplot(gr, phi[, 1:k], type = "l")
matplot(gr, phi_sm[, 1:k], type = "l")

diag(phi %*% t(phi))
diag(phi_sm %*% t(phi_sm))

cc <- phi %*% diag(lambda) %*% t(phi)
cc_sm <- phi_sm %*% diag(lambda) %*% t(phi_sm)

par(mfrow = c(2, 2))
GA::persp3D(gr, gr, cov.Mest,
            theta = -70, phi = 30, expand = 1)
GA::persp3D(gr, gr, cov.Mest.sm,
            theta = -70, phi = 30, expand = 1)
GA::persp3D(gr, gr, cc_sm,
            theta = -70, phi = 30, expand = 1)

get_eigen(cc_sm, gr)$PVE

# cov.Mest.sm <- cov_Mest(x, smooth = T, bw = 0.3, make.pos.semidef = F)

### Eigen function trajectories
eig.true <- get_delaigle_eigen(work.grid, model = 2)
# eig.true <- get_eigen(cov.true, gr)$phi[, 1:4]
eig.huber <- get_eigen(cov.huber, gr)$phi[, 1:4]
eig.Mest <- get_eigen(cov.Mest, gr)$phi[, 1:4]
eig.Mest.sm <- get_eigen(cov.Mest.sm, gr)$phi[, 1:4]
eig.cc.sm <- get_eigen(cc_sm, gr)$phi[, 1:4]
par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(work.grid,
          cbind(
            eig.true[, k],
            check_eigen_sign(eig.huber, eig.true)[, k],
            check_eigen_sign(eig.Mest, eig.true)[, k],
            check_eigen_sign(eig.Mest.sm, eig.true)[, k],
            check_eigen_sign(eig.cc.sm, eig.true)[, k]
          ),
          type = "l",
          col = 1:5,
          lty = 1:5,
          main = paste("Eigenfunction", k),
          xlab = "", ylab = "",
          lwd = rep(2, 4))
  if (k == 1) {
    legend("bottomleft",
           c("True","Huber","M-est","M-est(smooth)","cc_sm"),
           col = 1:7,
           lty = 1:7,
           lwd = rep(2, 7))
  }
}







matplot(gr, yy, pch = 1)
matlines(gr, fit2, type = "l", lwd = 2)

fit3 <- fit2 %*% diag( sqrt(diag(t(yy) %*% yy)) )
matlines(gr, fit3, type = "l", lwd = 2)



eig.obj <- get_eigen(cov.Mest, work.grid)
eig.obj$lambda
eig.obj$phi <- apply(eig.obj$phi, 2, function(x){ smooth.spline(work.grid, x)$y })
yy <- eig.obj$phi %*% diag(eig.obj$lambda) %*% t(eig.obj$phi)

diag(cov.Mest)
diag(yy)

par(mfrow = c(1, 2))
GA::persp3D(gr, gr, cov.Mest,
            theta = -70, phi = 30, expand = 1)
GA::persp3D(gr, gr, yy,
            theta = -70, phi = 30, expand = 1)


rob.var <- cov_Mest_cpp(x)
eig <- eigen(rob.var)
k <- which(eig$values > 0)
lambda <- eig$values[k]
phi <- matrix(eig$vectors[, k],
              ncol = length(k))

dim(phi)
# smoothing spline is performed for each eigenfunction
if (smooth == TRUE) {
  gr <- seq(0, 1, length.out = p)
  phi <- apply(phi, 2, function(u_i){
    stats::smooth.spline(gr, u_i)$y
  })
  # phi <- pracma::householder(phi)$Q
  phi <- far::orthonormalization(phi, basis = FALSE)
}
rob.var <- phi %*% diag(lambda, ncol = length(k)) %*% t(phi)



