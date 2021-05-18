################################################
### Simulation for covariance estimation
### - Lin setting (Low missingness)
### 1. PACE
### 2. Lin & Wang (2020)
### 3. Robust method (Huber loss)
### 4. WRM
### - For all methods, 5-fold CV are performed.
################################################
library(GA)   # persp plot
library(mvtnorm)
library(fdapace)   # 1, 2
library(mcfda)   # 7
library(synfd)   # 7
library(doParallel)   # parallel computing
library(doRNG)   # set.seed for foreach
library(MASS)   # huber, rlm
library(latex2exp)
library(tidyverse)
library(robfilter)
# source("R/functions.R")
# source("R/utills.R")
library(robfpca)
source("R/sim_Delaigle(2020).R")
source("R/sim_Lin_Wang(2020).R")
source("R/sim_kraus.R")
source("Kraus(2015)/pred.missfd.R")
source("Kraus(2015)/simul.missfd.R")


#####################################
### Delaigle (2020) setting
#####################################
num_sim <- 100
mise_reconstr <- matrix(0, num_sim, 3)
mse_reconstr <- matrix(0, num_sim, 3)
mise_completion <- matrix(0, num_sim, 5)
mse_completion <- matrix(0, num_sim, 5)

colnames(mise_reconstr) <- c("Yao","Lin","Robust")
colnames(mse_reconstr) <- c("Yao","Lin","Robust")
colnames(mise_completion) <- c("Yao","Lin","Kraus","Robust","Robust+align")
colnames(mse_completion) <- c("Yao","Lin","Kraus","Robust","Robust+align")

# simulation result
pca.est <- list()
num.sim <- 0   # number of simulations
seed <- 0   # current seed
sim.seed <- rep(NA, 100)   # collection of seed with no error occurs

while (num.sim < 100) {
  seed <- seed + 1
  set.seed(seed)
  print(paste0("Seed: ", seed))
# for (sim in 1:num_sim) {
#   print(paste(sim, "th simulation"))
#   # sim <- 1
#   set.seed(sim)
  
  #############################
  ### Data generation
  #############################
  n <- 100
  n.grid <- 51
  x.2 <- sim.kraus(n = 100, out.prop = 0.2, out.type = 4, grid.length = n.grid)
  df <- data.frame(
    id = factor(unlist(sapply(1:length(x.2$Lt), 
                              function(id) { 
                                rep(id, length(x.2$Lt[[id]])) 
                              }) 
    )),
    y = unlist(x.2$Ly),
    t = unlist(x.2$Lt)
  )
  # ggplot(df, aes(t, y, color = id)) +
  #   geom_line() +
  #   theme_bw() +
  #   # ylim(-10, 10) +
  #   theme(legend.position = "none")
  
  # spread data
  x <- df %>% 
    spread(key = "t", value = "y")
  x <- x[, -1] %>% 
    as.matrix
  # matplot(t(x), type = "l")
  
  
  #############################
  ### Covariance estimation
  #############################
  skip_sim <- FALSE   # if skip_sim == TRUE, pass this seed
  bw <- 0.2
  kernel <- "epanechnikov"
  work.grid <- seq(0, 1, length.out = n.grid)
  
  ### Yao, Müller, and Wang (2005)
  start_time <- Sys.time()
  registerDoRNG(seed)
  kern <- ifelse(kernel == "epanechnikov", "epan", kernel)
  optns <- list(methodXi = "CE", dataType = "Sparse", kernel = kern, verbose = FALSE,
                userBwMu = bw, userBwCov = bw)
  # optns <- list(methodXi = "CE", dataType = "Sparse", kernel = kern, verbose = FALSE,
  #               kFoldMuCov = 5, methodBwMu = "CV", methodBwCov = "CV", useBinnedCov = FALSE)
  tryCatch({
    mu.yao.obj <- GetMeanCurve(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)
    cov.yao.obj <- GetCovSurface(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)
  }, error = function(e) { 
    print("Yao cov error")
    print(e)
    skip_sim <<- TRUE
  })
  if (skip_sim == TRUE) {
    next
  }
  mu.yao <- mu.yao.obj$mu
  cov.yao <- cov.yao.obj$cov
  if (length(work.grid) != 51) {
    mu.yao <- ConvertSupport(fromGrid = cov.yao.obj$workGrid, toGrid = work.grid,
                             mu = mu.yao.obj$mu)
    cov.yao <- ConvertSupport(fromGrid = cov.yao.obj$workGrid, toGrid = work.grid,
                              Cov = cov.yao.obj$cov)
  }
  end_time <- Sys.time()
  print(paste0("Yao et al. : ", 
               round(difftime(end_time, start_time, units = "secs"), 3),
               " secs"))
  
  
  ### Lin & Wang (2020)
  start_time <- Sys.time()
  registerDoRNG(seed)
  tryCatch({
    # # 5-fold CV (It took very long time when we use CV option in mcfda package.)
    # cv.mu.lin.obj <- meanfunc.rob(x$Lt, x$Ly, method = "L2", kernel = kernel,
    #                               cv_bw_loss = "L2", ncores = ncores,
    #                               bw = NULL)
    # cv.var.lin.obj <- varfunc.rob(x$Lt, x$Ly, mu = cv.mu.lin.obj, kernel = kernel,
    #                               method = "L2",  cv_bw_loss = "L2", ncores = ncores,
    #                               bw = NULL)
    # # estimate mean, variance, covariance
    # mu.lin.obj <- meanfunc(x$Lt, x$Ly, method = "PACE", kernel = kernel,
    #                        bw = cv.mu.lin.obj$bw)   # It occurs error or very slow.
    # var.lin.obj <- varfunc(x$Lt, x$Ly, method = "PACE", kernel = kernel,
    #                        mu = mu.lin.obj, bw = cv.var.lin.obj$obj$bw)
    # cov.lin.obj <- covfunc(x$Lt, x$Ly, method = "SP",
    #                        mu = mu.lin.obj, sig2x = var.lin.obj)
    # cov.lin <- predict(cov.lin.obj, work.grid)
    mu.lin.obj <- meanfunc(x.2$Lt, x.2$Ly, method = "PACE", kernel = kernel,
                           bw = bw)   # It occurs error or very slow.
    var.lin.obj <- varfunc(x.2$Lt, x.2$Ly, method = "PACE", kernel = kernel,
                           mu = mu.lin.obj, bw = bw)
    cov.lin.obj <- covfunc(x.2$Lt, x.2$Ly, method = "SP",
                           mu = mu.lin.obj, sig2x = var.lin.obj)
  }, error = function(e) { 
    print("Lin cov error")
    print(e)
    skip_sim <<- TRUE
  })
  if (skip_sim == TRUE) {
    next
  }
  mu.lin <- predict(mu.lin.obj, work.grid)
  cov.lin <- predict(cov.lin.obj, work.grid)
  end_time <- Sys.time()
  print(paste0("Lin & Wang : ", 
               round(difftime(end_time, start_time, units = "secs"), 3),
               " secs"))
  
  
  ### Huber loss
  start_time <- Sys.time()
  registerDoRNG(seed)
  tryCatch({
    # # For delta in Huber function and bandwidth are selected from 5-fold CV
    # mu.huber.obj <- meanfunc.rob(x$Lt, x$Ly, method = "huber", kernel = kernel, 
    #                              bw = NULL, delta = NULL)
    # # bandwidth are selected from 5-fold CV (almost 3 minutes)
    # cov.huber.obj <- covfunc.rob(x$Lt, x$Ly, method = "huber", kernel = kernel, 
    #                              mu = mu.huber.obj, 
    #                              bw = NULL, delta = NULL)
    # For delta in Huber function and bandwidth are selected from 5-fold CV
    mu.huber.obj <- meanfunc.rob(x.2$Lt, x.2$Ly, method = "huber", kernel = kernel, 
                                 bw = bw, delta = 1.345)
    # bandwidth are selected from 5-fold CV (almost 3 minutes)
    cov.huber.obj <- covfunc.rob(x.2$Lt, x.2$Ly, method = "huber", kernel = kernel, 
                                 mu = mu.huber.obj, 
                                 bw = bw, delta = 1.345)
  }, error = function(e) { 
    print("Huber cov error")
    print(e)
    skip_sim <<- TRUE
  })
  if (skip_sim == TRUE) {
    next
  }
  mu.huber <- predict(mu.huber.obj, work.grid)
  cov.huber <- predict(cov.huber.obj, work.grid)
  end_time <- Sys.time()
  print(paste0("Robust (Huber loss) : ", 
               round(difftime(end_time, start_time, units = "secs"), 3),
               " secs"))
  
  
  # if some covariances is a not finite value
  if (!is.finite(sum(cov.yao)) | !is.finite(sum(cov.lin)) | 
      !is.finite(sum(cov.huber))) {
    cat("Estimated covariances do not have finite values. \n")
    next
  }
  
  # if all covariances are 0
  if ((sum(cov.yao) == 0) | (sum(cov.lin) == 0) | 
      (sum(cov.huber) == 0)) {
    cat("Estimated covariance have all 0 values. \n")
    next
  }
  
  
  ### Principal component analysis
  pve <- 0.99   # Not used if K is given
  K <- 5   # fixed number of PCs
  
  # Yao
  pca.yao.obj <- funPCA(x.2$Lt, x.2$Ly, 
                        mu.yao, cov.yao, sig2 = cov.yao.obj$sigma2, 
                        work.grid, PVE = pve, K = K)
  # Lin
  pca.lin.obj <- funPCA(x.2$Lt, x.2$Ly, 
                        mu.lin, cov.lin, sig2 = cov.lin.obj$sig2e, 
                        work.grid, PVE = pve, K = K)
  # Huber
  pca.huber.obj <- funPCA(x.2$Lt, x.2$Ly, 
                          mu.huber, cov.huber, sig2 = cov.huber.obj$sig2e, 
                          work.grid, PVE = pve, K = K)
  
  # # WRM - 535.63 secs(guass) / 75.33 (epan)
  # system.time({
  #   mu.wrm <- predict(mu.wrm.obj, work.grid)
  #   cov.wrm <- predict(cov.wrm.obj, work.grid)
  #   pca.wrm.obj <- PCA_CE(x.2$Lt, x.2$Ly, mu.wrm, cov.wrm, 
  #                         sig2 = cov.wrm.obj$sig2e, work.grid, K = NULL)
  # })
  
  
  ### Curve reconstruction via PCA
  # index of non-outlier curves having missing values
  cand <- which(apply(x, 1, function(x){ sum(is.na(x)) }) > 0)
  cand <- cand[cand <= 80]   # exclude outlier curves
  
  # reconstructed curves
  pred_yao_mat <- predict(pca.yao.obj, K = NULL)
  pred_lin_mat <- predict(pca.lin.obj, K = NULL)
  pred_huber_mat <- predict(pca.huber.obj, K = NULL)
  
  ise_reconstr <- matrix(NA, length(cand), 3)
  sse_reconstr <- matrix(NA, length(cand), 3)
  ise_completion <- matrix(NA, length(cand), 5)
  sse_completion <- matrix(NA, length(cand), 5)
  
  for (i in 1:length(cand)) {
    ind <- cand[i]
    
    pred_yao <- pred_yao_mat[ind, ]
    pred_lin <- pred_lin_mat[ind, ]
    pred_huber <- pred_huber_mat[ind, ]
    pred_kraus <- pred.missfd(x[ind, ], x)
    
    # ISE for reconstruction of overall interval
    df <- cbind(pred_yao,
                pred_lin,
                pred_huber)
    ise_reconstr[i, ] <- apply(df, 2, function(pred) { 
      get_ise(x.2$x.full[ind, ], pred, work.grid) 
    })
    sse_reconstr[i, ] <- apply(df, 2, function(pred) { 
      sum((x.2$x.full[ind, ] - pred)^2)
    })
    
    # ISE for completion
    NA_ind <- which(is.na(x[ind, ]))
    df <- cbind(pred_missing_curve(x[ind, ], pred_yao, conti = FALSE),
                pred_missing_curve(x[ind, ], pred_lin, conti = FALSE),
                pred_kraus,
                pred_missing_curve(x[ind, ], pred_huber, conti = FALSE),
                pred_missing_curve(x[ind, ], pred_huber, conti = FALSE, align = TRUE))
    df <- df[NA_ind, ]
    if (length(NA_ind) == 1) {
      df <- matrix(df, nrow = 1)
    }
    ise_completion[i, ] <- apply(df, 2, function(pred) { 
      get_ise(x.2$x.full[ind, NA_ind], pred, work.grid[NA_ind]) 
    })
    sse_completion[i, ] <- apply(df, 2, function(pred) { 
      sum((x.2$x.full[ind, NA_ind] - pred)^2)
    })
  }
  
  # update number of simulations and save seed which does not occur errors
  num.sim <- num.sim + 1
  sim.seed[num.sim] <- seed
  print(paste0("Total # of simulations: ", num.sim))
  
  mise_reconstr[num.sim, ] <- colMeans(ise_reconstr)
  mse_reconstr[num.sim, ] <- colMeans(sse_reconstr)
  mise_completion[num.sim, ] <- colMeans(ise_completion)
  mse_completion[num.sim, ] <- colMeans(sse_completion)
  
  pca.est[[num.sim]] <- list(seed = seed,
                             work.grid = work.grid,
                             mu.obj = list(yao = mu.yao.obj,
                                           lin = mu.lin.obj,
                                           huber = mu.huber.obj),
                             cov.obj = list(yao = cov.yao.obj,
                                            lin = cov.lin.obj,
                                            huber = cov.huber.obj),
                             cov = list(yao = cov.yao,
                                        lin = cov.lin,
                                        huber = cov.huber),
                             pca.obj = list(yao = pca.yao.obj,
                                            lin = pca.lin.obj,
                                            huber = pca.huber.obj))
}
save(list = c("pca.est","mise_reconstr","mse_reconstr","mise_completion","mse_completion"),
     file = "RData/20210421_completion_fixed.RData")

colMeans(mise_reconstr)
colMeans(mse_reconstr)
colMeans(mise_completion)
colMeans(mse_completion)

apply(mise_reconstr, 2, sd)
apply(mse_reconstr, 2, sd)
apply(mise_completion, 2, sd)
apply(mse_completion, 2, sd)


df <- cbind(
  MISE = paste0(
    round(colMeans(mise_reconstr), 2),
    " (",
    round(apply(mise_reconstr, 2, sd), 2),
    ")"
  ), 
  MSE = paste0(
    round(colMeans(mse_reconstr), 2),
    " (",
    round(apply(mse_reconstr, 2, sd), 2),
    ")"
  )
) %>% 
  as.data.frame %>% 
  rownames_to_column("method")
df <- cbind(
  MISE = paste0(
    round(colMeans(mise_completion[, c(1,2,4,3,5)]), 2),
    " (",
    round(apply(mise_completion[, c(1,2,4,3,5)], 2, sd), 2),
    ")"
  ), 
  MSE = paste0(
    round(colMeans(mse_completion[, c(1,2,4,3,5)]), 2),
    " (",
    round(apply(mse_completion[, c(1,2,4,3,5)], 2, sd), 2),
    ")"
  )
) %>% 
  as.data.frame %>% 
  rownames_to_column("method") %>% 
  left_join(df, by = "method")
df


##################################
### Test version
##################################
### data generation
set.seed(100)
n <- 100
n.grid <- 51
x.2 <- sim.kraus(n = n, out.prop = 0.2, out.type = 4, grid.length = n.grid)
df <- data.frame(
  id = factor(unlist(sapply(1:length(x.2$Lt), 
                            function(id) { 
                              rep(id, length(x.2$Lt[[id]])) 
                            }) 
  )),
  y = unlist(x.2$Ly),
  t = unlist(x.2$Lt)
)
# ggplot(df, aes(t, y, color = id)) +
#   geom_line() +
#   theme_bw() +
#   # ylim(-10, 10) +
#   theme(legend.position = "none")

# spread data
x <- df %>% 
  spread(key = "t", value = "y")
x <- x[, -1] %>% 
  as.matrix
matplot(t(x), type = "l")

# test for fixed parameters
system.time({
  # bw <- 0.1
  # kern <- "gauss"
  bw <- 0.2
  kern <- "epan"
  optns <- list(methodXi = "CE", dataType = "Sparse", kernel = kern, verbose = FALSE,
                userBwMu = bw, userBwCov = bw)
  mu.yao.obj <- GetMeanCurve(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)
  cov.yao.obj <- GetCovSurface(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)

  work.grid <- seq(0, 1, length.out = n.grid)
  mu.yao <- ConvertSupport(fromGrid = cov.yao.obj$workGrid, toGrid = work.grid,
                           mu = mu.yao.obj$mu)
  cov.yao <- ConvertSupport(fromGrid = cov.yao.obj$workGrid, toGrid = work.grid,
                            Cov = cov.yao.obj$cov)
})
# user  system elapsed
# 0.3     0.0     0.3

system.time({
  # kernel <- "gauss"
  kernel <- "epanechnikov"
  # estimate mean, variance, covariance
  mu.lin.obj <- meanfunc(x.2$Lt, x.2$Ly, method = "PACE", kernel = kernel,
                         bw = bw)   # It occurs error or very slow.
  # print("mean finish")
  var.lin.obj <- varfunc(x.2$Lt, x.2$Ly, method = "PACE", kernel = kernel,
                         mu = mu.lin.obj, bw = bw)
  # print("var finish")
  cov.lin.obj <- covfunc(x.2$Lt, x.2$Ly, method = "SP",
                         mu = mu.lin.obj, sig2x = var.lin.obj)
  # print("cov finish")
  mu.lin <- predict(mu.lin.obj, work.grid)
  cov.lin <- predict(cov.lin.obj, work.grid)
})
# user  system elapsed
# 20.65    0.11   20.75

system.time({
  # For delta in Huber function and bandwidth are selected from 5-fold CV
  mu.huber.obj <- meanfunc.rob(x.2$Lt, x.2$Ly, method = "huber", kernel = kernel,
                               bw = bw, delta = 1.345, ncores = 10)
  # bandwidth are selected from 5-fold CV (almost 3 minutes)
  cov.huber.obj <- covfunc.rob(x.2$Lt, x.2$Ly, method = "huber", kernel = kernel,
                               mu = mu.huber.obj,
                               bw = bw, delta = 1.345, ncores = 10)
  mu.huber <- predict(mu.huber.obj, work.grid)
  cov.huber <- predict(cov.huber.obj, work.grid)
})
# user  system elapsed
# 1.44    0.00    1.44
# 59.39    0.00   59.39    # 5-fold CV

# system.time({
#   # For delta in Huber function and bandwidth are selected from 5-fold CV
#   mu.wrm.obj <- meanfunc.rob(x.2$Lt, x.2$Ly, method = "wrm", kernel = kernel,
#                              bw = bw, delta = 1.345)
#   # bandwidth are selected from 5-fold CV (almost 3 minutes)
#   cov.wrm.obj <- covfunc.rob(x.2$Lt, x.2$Ly, method = "wrm", kernel = kernel,
#                              mu = mu.huber.obj,
#                              bw = bw, delta = 1.345)
#   mu.wrm <- predict(mu.wrm.obj, work.grid)
#   cov.wrm <- predict(cov.wrm.obj, work.grid)
# })
# # user  system elapsed
# # 286.64    0.01  286.82 gauss
# # user  system elapsed
# # 41.69    0.14   41.82 epan

### Principal component analysis
pve <- 0.99   # Not used if K is given
K <- 5   # fixed number of PCs

# Yao
pca.yao.obj <- funPCA(x.2$Lt, x.2$Ly, 
                      mu.yao, cov.yao, sig2 = cov.yao.obj$sigma2, 
                      work.grid, PVE = pve, K = K)
# Lin
pca.lin.obj <- funPCA(x.2$Lt, x.2$Ly, 
                      mu.lin, cov.lin, sig2 = cov.lin.obj$sig2e, 
                      work.grid, PVE = pve, K = K)
# Huber
pca.huber.obj <- funPCA(x.2$Lt, x.2$Ly, 
                        mu.huber, cov.huber, sig2 = cov.huber.obj$sig2e, 
                        work.grid, PVE = pve, K = K)
# # WRM - 535.63 secs(guass) / 75.33 (epan)
# system.time({
#   pca.wrm.obj <- PCA_CE(x.2$Lt, x.2$Ly, mu.wrm, cov.wrm,
#                         sig2 = cov.wrm.obj$sig2e, work.grid, K = NULL)
# })


### Covariance surfaces
cov.true <- get_cov_fragm(work.grid, model = 2)   # true covariance
par(mfrow = c(2, 2))
persp3D(work.grid, work.grid, cov.true,
        main = "True",
        theta = -70, phi = 30, expand = 1)
persp3D(work.grid, work.grid, cov.yao, 
        zlim = range(cov.true), main = "Yao",
        theta = -70, phi = 30, expand = 1)
persp3D(work.grid, work.grid, cov.lin, 
        zlim = range(cov.true), main = "Lin",
        theta = -70, phi = 30, expand = 1)
persp3D(work.grid, work.grid, cov.huber, 
        zlim = range(cov.true), main = "Huber",
        theta = -70, phi = 30, expand = 1)
par(mfrow = c(1, 1))

matplot(work.grid,
        cbind(diag(cov.true),
              diag(cov.yao),
              diag(cov.lin),
              diag(cov.huber)),
        type = "l", ylim = range(cov.true),
        xlab = "", ylab = "")
matplot(work.grid,
        cbind(mu.yao,
              mu.lin,
              mu.huber),
        type = "l", xlab = "", ylab = "")


### Completion
par(mfrow = c(3, 3))
cand <- which(apply(x, 1, function(x){ sum(is.na(x)) }) > 0)
cand <- cand[cand <= 80]   # exclude outlier curves
# par(mfrow = c(1, 3))
# cand <- c(25, 70, 80)
for (ind in cand) {
  pred_yao <- predict(pca.yao.obj, K = NULL)[ind, ]
  pred_lin <- predict(pca.lin.obj, K = NULL)[ind, ]
  pred_huber <- predict(pca.huber.obj, K = NULL)[ind, ]
  # pred_wrm <- mu.wrm + matrix(pca.wrm.obj$pc.score[1, ], nrow = 1) %*% t(pca.wrm.obj$eig.fun)
  pred_kraus <- pred.missfd(x[ind, ], x)
  
  is_snippets <- (max( diff( which(!is.na(x[ind, ])) ) ) == 1)
  if (is_snippets) {
    obs_range <- range(which(!is.na(x[ind, ])))   # index range of observed periods
    
    if ((obs_range[1] > 1) & (obs_range[2] < n.grid)) {
      # start and end
      obs_range <- obs_range
    } else if ((obs_range[1] > 1) | (obs_range[2] < n.grid)) {
      if (obs_range[1] > 1) {
        # start periods
        obs_range <- obs_range[1]
      } else if (obs_range[2] < n.grid) {
        # end periods
        obs_range <- obs_range[2]
      }
    }
  } else {
    # missing is in the middle.
    obs_range <- range(which(is.na(x[ind, ])))
    # include last observed point
    obs_range <- c(obs_range[1] - 1,
                   obs_range[2] + 1)
  }
  
  df <- cbind(x.2$x.full[ind, ],
              pred_missing_curve(x[ind, ], pred_yao),
              pred_missing_curve(x[ind, ], pred_lin),
              pred_missing_curve(x[ind, ], pred_huber),
              pred_kraus,
              pred_missing_curve(x[ind, ], pred_huber, align = TRUE))
  matplot(work.grid, df, type = "o",
          pch = rep(1, 6), lty = rep(1, 6), 
          # lwd = c(1,1,1,2,1,2),
          xlab = "", ylab = "", main = paste0(ind, "th trajectory"))
  # points(work.grid, x.2$x.full[ind, ])
  abline(v = work.grid[obs_range],
         lty = 2, lwd = 2)
  grid()
  if (ind %in% cand[(0:6)*9 + 1]) {
    legend("topleft",
           c("True","Yao","Lin","Huber","Kraus","Huber-aligned"),
           col = 1:6,
           lty = rep(1, 6))
  }
}
par(mfrow = c(1, 1))




### reconstruction for missing parts
which(apply(x, 1, function(x){ sum(is.na(x)) }) > 10)
which(apply(x, 1, function(x){ sum(is.na(x)) }) > 10 &
        apply(x, 1, function(x){ is.na(x[51]) | is.na(x[1]) }) > 0)

ind <- 37

pred_yao <- predict(pca.yao.obj, K = NULL)[ind, ]
pred_lin <- predict(pca.lin.obj, K = NULL)[ind, ]
pred_huber <- predict(pca.huber.obj, K = NULL)[ind, ]
# pred_wrm <- mu.wrm + matrix(pca.wrm.obj$pc.score[1, ], nrow = 1) %*% t(pca.wrm.obj$eig.fun)
pred_kraus <- pred.missfd(x[ind, ], x)

par(mfrow = c(1, 2))
df <- cbind(x.2$x.full[ind, ],
            pred_missing_curve(x[ind, ], pred_yao),
            pred_missing_curve(x[ind, ], pred_lin),
            pred_missing_curve(x[ind, ], pred_huber),
            pred_kraus,
            pred_missing_curve(x[ind, ], pred_huber, align = TRUE),
            pred_missing_curve_2(x[ind, ], pred_huber, align = TRUE))
matplot(work.grid, df, type = "l",
        xlab = "", ylab = "", main = "Completion for missing parts")
abline(v = work.grid[ range(which(!is.na(x[ind, ]))) ],
       lty = 2, lwd = 2)
grid()
legend("topleft",
       c("True","Yao","Lin","Huber","Kraus","Huber-algined"),
       col = 1:6,
       lty = rep(1, 7))

df <- cbind(
  x.2$x.full[ind, ],
  pred_yao,
  pred_lin,
  pred_huber,
  pred_kraus
)
matplot(work.grid, df, type = "l",
        xlab = "", ylab = "", main = "Reconstruction")
abline(v = work.grid[ range(which(!is.na(x[ind, ]))) ],
       lty = 2, lwd = 2)
grid()
legend("topleft", 
       c("True","Yao","Lin","Huber","Kraus"),
       col = 1:5,
       lty = rep(1, 5))


pred <- pred_missing_curve_2(x[ind, ], pred_huber, align = TRUE)
ind_missing <- which(is.na(x[ind, ]))
z <- x[ind, ]
z[ind_missing] <- pred[ind_missing]

i <- (ind_missing[1]-4):(ind_missing[1]-1)
j <- ind_missing[1]:(ind_missing[1]+2)

get_deriv(z, work.grid)[c(i, j)]
get_deriv(get_deriv(z, work.grid),
          work.grid)[c(i, j)]


get_deriv(c(x[ind, i], pred[j]),
          work.grid[c(i, j)])
get_deriv(get_deriv(c(x[ind, i], pred[j]),
                    work.grid[c(i, j)]),
          work.grid[c(i, j)])
pred_missing_curve(x[ind, ], pred_huber, align = TRUE)
pred
lines(work.grid, pred, col = 2, lwd = 2)








# eigen analysis
cov.true <-  get_cov_fragm(work.grid, model = model.cov)
eig.true <- get_eigen(cov = cov.true, grid = work.grid)
eig.yao <- get_eigen(cov = cov.yao, grid = work.grid)
eig.lin <- get_eigen(cov = cov.lin, grid = work.grid)
eig.huber <- get_eigen(cov = cov.huber, grid = work.grid)
# eig.wrm <- get_eigen(cov = cov.wrm, grid = work.grid)

# change eigen direction(sign) for first K eigenvectors
K <- 3
eig.yao$phi[, 1:K] <- check_eigen_sign(eig.yao$phi[, 1:K], eig.true$phi[, 1:K])
eig.lin$phi[, 1:K] <- check_eigen_sign(eig.lin$phi[, 1:K], eig.true$phi[, 1:K])
eig.huber$phi[, 1:K] <- check_eigen_sign(eig.huber$phi[, 1:K], eig.true$phi[, 1:K])
# eig.wrm$phi[, 1:K] <- check_eigen_sign(eig.wrm$phi[, 1:K], eig.true$phi[, 1:K])

# fitst 3 eigenfunctions
fig <- list()
for (i in 1:K) {
  fig.data <- data.frame(
    work.grid = rep(work.grid, 4),
    phi = c(eig.true$phi[, i],
            eig.yao$phi[, i],
            eig.lin$phi[, i],
            eig.huber$phi[, i]),
    method = factor(rep(c("True","Yao","Lin","Huber"),
                        each = length(work.grid)),
                    levels = c("True","Yao","Lin","Huber"))
  )
  fig[[i]] <- ggplot(data = fig.data, 
                     mapping = aes(work.grid, phi, color = method, linetype = method)) +
    geom_line(size = 1) +
    labs(x = TeX("$t$"), y = TeX(paste0("$\\phi_", i, "(t)$"))) +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.position = "bottom")
}
gridExtra::grid.arrange(grobs = fig, nrow = 1)




#####################################
### Lin & Wang (2020) setting
#####################################
set.seed(1000)
x <- reg.fd(
  n = 10, 
  X = kl.process(domain = c(0, 1),
                 eigen.values = 1/(3^(1:25)),
                 eigen.functions = c("FOURIER"),
                 distribution = c("GAUSSIAN"))
)
matplot(t(x$y), type = "l")

mu <- function(s) sin(2*pi*s)
D <- irreg.fd(mu=mu, X=synfd::gaussian.process(), n=100, m=5)
plot(D$t[[1]], D$y[[1]], type = "o", xlim = c(0, 1), ylim = range(unlist(D$y)))
for (i in 2:100) {
  lines(D$t[[i]], D$y[[i]], col = i, type = "o")
}


n.grid <- 51
x <- sim.kraus(n = 100, out.prop = 0, out.type = 4, grid.length = n.grid)
plot(x$Lt[[i]], x$Ly[[i]], type = "l",
     xlim = c(0, 1), ylim = range(unlist(x$Ly)))
for (i in 2:10) {
  lines(x$Lt[[i]], x$Ly[[i]], col = i, lty = i)
}
length(which(cov(x$x.full) < 0))


### simulation data test
set.seed(100)
n <- 100
n.grid <- 51
x.2 <- sim_lin_wang(n = n, out.prop = 0, out.type = 4,
                    process = kl.process(domain = c(0, 1),
                                         eigen.values = 1/(2^(1:3)),
                                         eigen.functions = c("COS"),
                                         distribution = c("EXPONENTIAL")),
                    regular.grid = TRUE, grid.length = 51)
df <- data.frame(
  id = factor(unlist(sapply(1:length(x.2$Lt), 
                            function(id) { 
                              rep(id, length(x.2$Lt[[id]])) 
                            }) 
  )),
  y = unlist(x.2$Ly),
  t = unlist(x.2$Lt)
)

# spread data
x <- df %>% 
  spread(key = "t", value = "y")
x <- x[, -1] %>% 
  as.matrix
matplot(t(x), type = "l")


# test for fixed parameters
system.time({
  # bw <- 0.1
  # kern <- "gauss"
  bw <- 0.2
  kern <- "epan"
  optns <- list(methodXi = "CE", dataType = "Sparse", kernel = kern, verbose = FALSE,
                userBwMu = bw, userBwCov = bw)
  mu.yao.obj <- GetMeanCurve(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)
  cov.yao.obj <- GetCovSurface(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)
  
  work.grid <- seq(0, 1, length.out = n.grid)
  mu.yao <- ConvertSupport(fromGrid = cov.yao.obj$workGrid, toGrid = work.grid,
                           mu = mu.yao.obj$mu)
  cov.yao <- ConvertSupport(fromGrid = cov.yao.obj$workGrid, toGrid = work.grid,
                            Cov = cov.yao.obj$cov)
})

system.time({
  # kernel <- "gauss"
  kernel <- "epanechnikov"
  # estimate mean, variance, covariance
  mu.lin.obj <- meanfunc(x.2$Lt, x.2$Ly, method = "PACE", kernel = kernel,
                         bw = bw)   # It occurs error or very slow.
  # print("mean finish")
  var.lin.obj <- varfunc(x.2$Lt, x.2$Ly, method = "PACE", kernel = kernel,
                         mu = mu.lin.obj, bw = bw)
  # print("var finish")
  cov.lin.obj <- covfunc(x.2$Lt, x.2$Ly, method = "SP",
                         mu = mu.lin.obj, sig2x = var.lin.obj)
  # print("cov finish")
  mu.lin <- predict(mu.lin.obj, work.grid)
  cov.lin <- predict(cov.lin.obj, work.grid)
  
  if (sum(cov.lin) == 0) {
    stop("All of estimated covariances are 0.")
  }
})

system.time({
  # For delta in Huber function and bandwidth are selected from 5-fold CV
  mu.huber.obj <- meanfunc.rob(x.2$Lt, x.2$Ly, method = "huber", kernel = kernel,
                               bw = bw, delta = 1.345)
  # bandwidth are selected from 5-fold CV (almost 3 minutes)
  cov.huber.obj <- covfunc.rob(x.2$Lt, x.2$Ly, method = "huber", kernel = kernel,
                               mu = mu.huber.obj,
                               bw = bw, delta = 1.345)
  mu.huber <- predict(mu.huber.obj, work.grid)
  cov.huber <- predict(cov.huber.obj, work.grid)
  
  if (sum(cov.huber) == 0) {
    stop("All of estimated covariances are 0.")
  }
})

### Covariance surfaces
# cov.true <- cov(x)   # true covariance
par(mfrow = c(2, 2))
# persp3D(work.grid, work.grid, cov.true,
#         main = "True",
#         theta = -70, phi = 30, expand = 1)
persp3D(work.grid, work.grid, cov.yao, 
        # zlim = range(cov.true), main = "Yao",
        theta = -70, phi = 30, expand = 1)
persp3D(work.grid, work.grid, cov.lin, 
        # zlim = range(cov.true), main = "Lin",
        theta = -70, phi = 30, expand = 1)
persp3D(work.grid, work.grid, cov.huber, 
        # zlim = range(cov.true), main = "Huber",
        theta = -70, phi = 30, expand = 1)
par(mfrow = c(1, 1))


### Principal component analysis
pve <- 0.99   # Not used if K is given
K <- 2   # fixed number of PCs

# Yao
pca.yao.obj <- funPCA(x.2$Lt, x.2$Ly, 
                      mu.yao, cov.yao, sig2 = cov.yao.obj$sigma2, 
                      work.grid, PVE = pve, K = K)
# Lin
pca.lin.obj <- funPCA(x.2$Lt, x.2$Ly, 
                      mu.lin, cov.lin, sig2 = cov.lin.obj$sig2e, 
                      work.grid, PVE = pve, K = K)
# Huber
pca.huber.obj <- funPCA(x.2$Lt, x.2$Ly, 
                        mu.huber, cov.huber, sig2 = cov.huber.obj$sig2e, 
                        work.grid, PVE = pve, K = K)

### Completion
par(mfrow = c(3, 3))
cand <- which(apply(x, 1, function(x){ sum(is.na(x)) }) > 0)
cand <- cand[cand <= 80]   # exclude outlier curves
# par(mfrow = c(1, 3))
# cand <- c(25, 70, 80)
for (ind in cand) {
  pred_yao <- predict(pca.yao.obj, K = NULL)[ind, ]
  pred_lin <- predict(pca.lin.obj, K = NULL)[ind, ]
  pred_huber <- predict(pca.huber.obj, K = NULL)[ind, ]
  # pred_wrm <- mu.wrm + matrix(pca.wrm.obj$pc.score[1, ], nrow = 1) %*% t(pca.wrm.obj$eig.fun)
  pred_kraus <- pred.missfd(x[ind, ], x)
  
  is_snippets <- (max( diff( which(!is.na(x[ind, ])) ) ) == 1)
  if (is_snippets) {
    obs_range <- range(which(!is.na(x[ind, ])))   # index range of observed periods
    
    if ((obs_range[1] > 1) & (obs_range[2] < n.grid)) {
      # start and end
      obs_range <- obs_range
    } else if ((obs_range[1] > 1) | (obs_range[2] < n.grid)) {
      if (obs_range[1] > 1) {
        # start periods
        obs_range <- obs_range[1]
      } else if (obs_range[2] < n.grid) {
        # end periods
        obs_range <- obs_range[2]
      }
    }
  } else {
    # missing is in the middle.
    obs_range <- range(which(is.na(x[ind, ])))
    # include last observed point
    obs_range <- c(obs_range[1] - 1,
                   obs_range[2] + 1)
  }
  
  df <- cbind(x.2$x.full[ind, ],
              pred_missing_curve(x[ind, ], pred_yao),
              pred_missing_curve(x[ind, ], pred_lin),
              pred_missing_curve(x[ind, ], pred_huber),
              pred_kraus,
              pred_missing_curve(x[ind, ], pred_huber, align = TRUE))
  matplot(work.grid, df, type = "l",
          pch = rep(1, 6), lty = rep(1, 6), 
          # lwd = c(1,1,1,2,1,2),
          xlab = "", ylab = "", main = paste0(ind, "th trajectory"))
  points(work.grid, x.2$x.full[ind, ])
  abline(v = work.grid[obs_range],
         lty = 2, lwd = 2)
  grid()
  if (ind %in% cand[(0:6)*9 + 1]) {
    legend("topleft",
           c("True","Yao","Lin","Huber","Kraus","Huber-aligned"),
           col = 1:6,
           lty = rep(1, 6))
  }
}
par(mfrow = c(1, 1))

