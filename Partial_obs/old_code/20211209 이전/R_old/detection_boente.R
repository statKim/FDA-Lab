################################################
### Outlier Detection - Boente(2015) setting
### - Partially observed case
################################################
library(GA)   # persp plot
library(mvtnorm)
library(fdapace)   # 1, 2
library(mcfda)   # 7
library(synfd)   # 7
library(doParallel)   # parallel computing
library(doRNG)   # set.seed for foreach
library(MASS)   # huber, rlm
library(tidyverse)
library(latex2exp)
library(xtable)
library(robfpca)
source("R/sim_boente.R")
source("Kraus(2015)/pred.missfd.R")
source("Kraus(2015)/simul.missfd.R")
source("robust_Kraus.R")
source("Boente_cov.R")


# score distance
score_dist <- function(funPCA.obj) {
  K <- funPCA.obj$K
  score <- funPCA.obj$pc.score
  lambda <- funPCA.obj$lambda
  
  SD <- apply(score, 1, function(row){
    sqrt(sum(row^2 / lambda))
  })
  
  return(SD)
}

# orthogonal distance
orthogonal_dist <- function(funPCA.obj, X) {
  K <- funPCA.obj$K
  X_hat <- predict(funPCA.obj, K = K)
  
  OD <- (X - X_hat)^2
  OD <- apply(OD, 1, function(row){
    sqrt(sum(row, na.rm = T))
  })
  
  return(OD)
}




#####################################
### Simulation Parameters
#####################################
num_sim <- 30   # number of simulations
out_prop <- 0.2   # proportion of outliers
model <- 4   # type of outliers
data_type <- "partial"   # type of functional data
kernel <- "epanechnikov"   # kernel function for local smoothing
# kernel <- "gauss"   # kernel function for local smoothing
bw_boente <- 0.1   # bandwidth for Boente(2020) - Error occurs for small bw
n_cores <- 12   # number of threads for parallel computing
pve <- 0.95   # Not used if K is given
# fixed number of PCs (If NULL, it is selected by PVE)
if (model == 3) {
  K <- 5
} else {
  K <- 2
}

sens <- matrix(NA, num_sim, 6)
spec <- matrix(NA, num_sim, 6)
colnames(sens) <- c("Yao","Boente",
                    "M-est","M-est-noise",
                    "M-est(smooth)","M-est(smooth)-noise")
colnames(spec) <- c("Yao","Boente",
                    "M-est","M-est-noise",
                    "M-est(smooth)","M-est(smooth)-noise")
pca.est <- list()

# simulation result
num.sim <- 0   # number of simulations
seed <- 0   # current seed
sim.seed <- rep(NA, num_sim)   # collection of seed with no error occurs
while (num.sim < num_sim) {
  seed <- seed + 1
  set.seed(seed)
  print(paste0("Seed: ", seed))
  
  #############################
  ### Data generation
  #############################
  n <- 100
  n.grid <- 51
  x.2 <- sim_boente(n = n, 
                    model = model,
                    type = data_type,
                    eps1 = out_prop,
                    outlier = TRUE)
  
  x <- list2matrix(x.2)
  # matplot(t(x), type = "l")
  
  #############################
  ### Covariance estimation
  #############################
  skip_sim <- FALSE   # if skip_sim == TRUE, pass this seed
  work.grid <- seq(0, 1, length.out = n.grid)
  
  ### M-estimator
  start_time <- Sys.time()
  registerDoRNG(seed)
  tryCatch({
    mu.Mest <- mean_Mest(x)
    mu.Mest.sm <- mean_Mest(x, smooth = TRUE)
    # noise_var <- sigma2.rob(x.2$Lt, x.2$Ly)   # robust noise variance estimator
    
    # Not adjust noise
    cov.Mest <- cov_Mest(x)
    cov.Mest.sm <- cov_Mest(x, smooth = T)
    
    # adjust noise
    # noise_var <- noise_var_M(cov.Mest, cov.Mest.sm, work.grid)
    noise_var <- 1
    cov.Mest.noise <- cov_Mest(x, noise.var = noise_var)
    cov.Mest.sm.noise <- cov_Mest(x, smooth = T,
                                  noise.var = noise_var)
  }, error = function(e) { 
    print("M-est cov error")
    print(e)
    skip_sim <<- TRUE
  })
  if (skip_sim == TRUE) {
    next
  }
  end_time <- Sys.time()
  print(paste0("M-est : ", 
               round(difftime(end_time, start_time, units = "secs"), 3),
               " secs"))
  
  
  ### Yao, Müller, and Wang (2005)
  ## 30 secs
  start_time <- Sys.time()
  registerDoRNG(seed)
  tryCatch({
    kern <- ifelse(kernel == "epanechnikov", "epan", kernel)
    optns <- list(methodXi = "CE", dataType = "Sparse", kernel = kern, verbose = FALSE,
                  # userBwMu = bw, userBwCov = bw)
                  kFoldMuCov = 5, methodBwMu = "CV", methodBwCov = "CV", useBinnedCov = FALSE)
    mu.yao.obj <- GetMeanCurve(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)
    cov.yao.obj <- GetCovSurface(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)
    mu.yao <- mu.yao.obj$mu
    cov.yao <- cov.yao.obj$cov
    if (length(work.grid) != 51) {
      mu.yao <- ConvertSupport(fromGrid = cov.yao.obj$workGrid, toGrid = work.grid,
                               mu = mu.yao.obj$mu)
      cov.yao <- ConvertSupport(fromGrid = cov.yao.obj$workGrid, toGrid = work.grid,
                                Cov = cov.yao.obj$cov)
    }
  }, error = function(e) { 
    print("Yao cov error")
    print(e)
    skip_sim <<- TRUE
  })
  if (skip_sim == TRUE) {
    next
  }
  end_time <- Sys.time()
  print(paste0("Yao et al. : ", 
               round(difftime(end_time, start_time, units = "secs"), 3),
               " secs"))
  
  
  ### Boente et al. (2020)
  start_time <- Sys.time()
  registerDoRNG(seed)
  tryCatch({
    cov.boente.obj <- cov_boente(x.2, bw.mu = bw_boente, bw.cov = bw_boente)
    mu.boente <- cov.boente.obj$mu
    cov.boente <- cov.boente.obj$cov
    # noise var from source code of sparseFPCA package
    noise_boente <- eigen(cov.boente)$values[1] / (1e3 - 1)
  }, error = function(e) {
    print("Boente (2020) cov error")
    print(e)
    skip_sim <<- TRUE
  })
  if (skip_sim == TRUE) {
    next
  }
  end_time <- Sys.time()
  print(paste0("Boente (2020) : ", 
               round(difftime(end_time, start_time, units = "secs"), 3),
               " secs"))
  
  # if some covariances is a not finite value
  if (!is.finite(sum(cov.yao)) | !is.finite(sum(cov.boente)) | 
      !is.finite(sum(cov.Mest)) | !is.finite(sum(cov.Mest.sm)) |
      !is.finite(sum(cov.Mest.noise)) | !is.finite(sum(cov.Mest.sm.noise))) {
    cat("Estimated covariances do not have finite values. \n")
    next
  }
  
  # if all covariances are 0
  if ((sum(cov.yao) == 0) | (sum(cov.boente) == 0) |
      (sum(cov.Mest) == 0) | (sum(cov.Mest.sm) == 0) |
      (sum(cov.Mest.noise) == 0) | (sum(cov.Mest.sm.noise) == 0)) {
    cat("Estimated covariance have all 0 values. \n")
    next
  }
  
  
  ### Principal component analysis
  # Yao
  pca.yao.obj <- funPCA(x.2$Lt, x.2$Ly, 
                        mu.yao, cov.yao, sig2 = cov.yao.obj$sigma2, 
                        work.grid, PVE = pve, K = K)
  # Boente
  pca.boente.obj <- funPCA(x.2$Lt, x.2$Ly, 
                           mu.boente, cov.boente, sig2 = noise_boente, 
                           work.grid, PVE = pve, K = K)
  # M-est
  pca.Mest.obj <- funPCA(x.2$Lt, x.2$Ly,
                         mu.Mest, cov.Mest, sig2 = 0,
                         work.grid, PVE = pve, K = K)
  pca.Mest.sm.obj <- funPCA(x.2$Lt, x.2$Ly,
                            mu.Mest.sm, cov.Mest.sm, sig2 = 0,
                            work.grid, PVE = pve, K = K)
  # consider noise var
  pca.Mest.noise.obj <- funPCA(x.2$Lt, x.2$Ly,
                               mu.Mest, cov.Mest.noise, sig2 = noise_var,
                               work.grid, PVE = pve, K = K)
  pca.Mest.sm.noise.obj <- funPCA(x.2$Lt, x.2$Ly,
                                  mu.Mest.sm, cov.Mest.sm.noise, sig2 = noise_var,
                                  work.grid, PVE = pve, K = K)
  
  ### Outlier detection
  SD <- cbind(
    score_dist(pca.yao.obj),
    score_dist(pca.boente.obj),
    score_dist(pca.Mest.obj),
    score_dist(pca.Mest.noise.obj),
    score_dist(pca.Mest.sm.obj),
    score_dist(pca.Mest.sm.noise.obj)
  )
  for (i in 1:6) {
    cut_off <- boxplot(SD[, i], plot = F)$stats[5, 1]
    
    y_hat <- rep(0, n)
    y_hat[which(SD[, i] > cut_off)] <- 1
    
    sens[num.sim+1, i] <- caret::sensitivity(factor(y_hat, levels = c(1, 0)), 
                                             factor(x.2$y, levels = c(1, 0)))
    spec[num.sim+1, i] <- caret::specificity(factor(y_hat, levels = c(1, 0)), 
                                             factor(x.2$y, levels = c(1, 0)))
  }
  
 
  num.sim <- num.sim + 1 
  print(paste0("Total # of simulations: ", num.sim))
  print(colMeans(sens, na.rm = T)) 
  
  pca.est[[num.sim]] <- list(seed = seed,
                             x.2 = x.2,
                             work.grid = work.grid,
                             pca.obj = list(pca.yao.obj,
                                            pca.boente.obj,
                                            pca.Mest.obj,
                                            pca.Mest.noise.obj,
                                            pca.Mest.sm.obj,
                                            pca.Mest.sm.noise.obj))
}
# save(list = c("pca.est"),
#      file = "RData/20210819_detect_boente-4.RData")

round(colMeans(sens), 3)
round(colMeans(spec), 3)



sens <- matrix(NA, num_sim, 6)
spec <- matrix(NA, num_sim, 6)
colnames(sens) <- c("Yao","Boente",
                    "M-est","M-est-noise",
                    "M-est(smooth)","M-est(smooth)-noise")
colnames(spec) <- c("Yao","Boente",
                    "M-est","M-est-noise",
                    "M-est(smooth)","M-est(smooth)-noise")
for (sim in 1:num_sim) {
  x.2 <- pca.est[[sim]]$x.2
  x <- list2matrix(x.2)
  k <- pca.est[[sim]]$pca.obj[[1]]$K
  
  ### Outlier detection
  for (i in 1:6) {
    # Outlier map (Score dist VS Orthogonal dist), See Hubert(2005)
    SD <- score_dist(pca.est[[sim]]$pca.obj[[i]])   # score distance
    OD <- orthogonal_dist(pca.est[[sim]]$pca.obj[[i]], x)   # orthogonal distance
    # mcd_fit <- covMcd(OD^(2/3))
    mcd_fit <- MASS::cov.mcd(matrix(OD^(2/3)))
    cut_y <- (mcd_fit$center + sqrt(as.numeric(mcd_fit$cov))*qnorm(0.975))^(3/2)
    cut_x <- sqrt(qchisq(0.975, k))
    
    y_hat <- rep(0, n)
    y_hat[which(SD > cut_x & OD > cut_y)] <- 1

    # # Boxplot for Score distance    
    # SD <- score_dist(pca.est[[sim]]$pca.obj[[i]])   # score distance
    # cut_off <- boxplot(SD, plot = F)$stats[5, 1]
    # 
    # y_hat <- rep(0, n)
    # y_hat[which(SD > cut_off)] <- 1
    
    sens[sim, i] <- caret::sensitivity(factor(y_hat, levels = c(1, 0)), 
                                             factor(x.2$y, levels = c(1, 0)))
    spec[sim, i] <- caret::specificity(factor(y_hat, levels = c(1, 0)), 
                                             factor(x.2$y, levels = c(1, 0)))
  }
}
round(colMeans(sens), 3)
round(colMeans(spec), 3)



sim <- 11
x.2 <- pca.est[[sim]]$x.2
x <- list2matrix(x.2)
k <- pca.est[[sim]]$pca.obj[[1]]$K
mname <- c("Yao","Boente",
           "M-est","M-est-noise",
           "M-est(smooth)","M-est(smooth)-noise")
par(mfrow = c(2, 3))
for (i in 1:6) {
  SD <- score_dist(pca.est[[sim]]$pca.obj[[i]])   # score distance
  OD <- orthogonal_dist(pca.est[[sim]]$pca.obj[[i]], x)   # orthogonal distance
  # mcd_fit <- covMcd(OD^(2/3))
  mcd_fit <- cov.mcd(matrix(OD^(2/3)))
  cut_y <- (mcd_fit$center + sqrt(as.numeric(mcd_fit$cov))*qnorm(0.975))^(3/2)
  cut_x <- sqrt(qchisq(0.975, k))
  
  plot(SD, OD,
       col = 3*x.2$y + 1,
       xlab = "Score distance",
       ylab = "Orthogonal distance",
       main = mname[i],
       cex.lab = 1.5,
       cex.main = 1.5)
  grid()
  abline(v = cut_x, col = 2)
  abline(h = cut_y, col = 2)
}


cov.mcd(matrix(OD))
covMcd(OD)

plot(sqrt(SD*OD), col = x.2$y+1)
boxplot(SD*OD)
adjbox(SD*OD)





### Completion한 다음에 outlier detection
# Mah, Hu 제외한 경우에 sensitivity = 0
library(rainbow)
sens_list <- list()
spec_list <- list()
for (sim in 1:num_sim) {
  x.2 <- pca.est[[sim]]$x.2
  x <- list2matrix(x.2)
  k <- pca.est[[sim]]$pca.obj[[1]]$K
  y <- x.2$y   # true outlier index
  
  # repeat with each PCA methods
  for (m in 1:6) {
    print(paste(sim, "-", m))
    # reconstruction
    recon_mat <- predict(pca.est[[sim]]$pca.obj[[m]], K = NULL)
    
    # make matrix of complete curve
    X_comp <- x
    rownames(X_comp) <- 1:n
    colnames(X_comp) <- gr
    for (i in 1:n) {
      X_i <- x[i, ]
      X_i_hat <- recon_mat[i, ]
      
      idx_miss <- which(is.na(X_i))
      if (length(idx_miss) > 0) {
        X_comp[i, idx_miss] <- X_i_hat[idx_miss]
      }
    }
    
    # outlier detection using complete curves
    fds.obj <- fds(x = gr,
                   y = t(X_comp), 
                   xname = "Time", yname = "Value")
    fout <- list()
    fout[[1]] <- foutliers(data = fds.obj, method = "robMah")$outliers
    fout[[2]] <- foutliers(data = fds.obj, method = "lrt")$outliers
    fout[[3]] <- foutliers(data = fds.obj, method = "depth.trim")$outliers
    fout[[4]] <- foutliers(data = fds.obj, method = "depth.pond")$outliers
    fout[[5]] <- as.numeric( foutliers(data = fds.obj, method = "HUoutliers")$outliers )
    
    # summary outlier detection
    sens <- numeric(5)
    spec <- numeric(5)
    for (j in 1:5) {
      y_hat <- rep(0, n)
      y_hat[fout[[j]]] <- 1
      
      sens[j] <- caret::sensitivity(factor(y_hat, levels = c(1, 0)), 
                                    factor(y, levels = c(1, 0)))
      spec[j] <- caret::specificity(factor(y_hat, levels = c(1, 0)), 
                                    factor(y, levels = c(1, 0)))
    }

    # save result
    if (sim == 1) {
      names(sens) <- c("Mah","LRT","DTR","DWE","HU")
      names(spec) <- c("Mah","LRT","DTR","DWE","HU")
      sens_list[[m]] <- sens
      spec_list[[m]] <- spec
    } else {
      sens_list[[m]] <- rbind(sens_list[[m]],
                              sens)
      spec_list[[m]] <- rbind(spec_list[[m]],
                              spec)
    }
  }
}

sapply(sens_list, colMeans)[c(1, 5), ]
sapply(spec_list, colMeans)[c(1, 5), ]




### Outlier detection
library(rainbow)
sens_list <- list()
spec_list <- list()
for (sim in 1:num_sim) {
  x.2 <- pca.est[[sim]]$x.2
  x <- list2matrix(x.2)
  k <- pca.est[[sim]]$pca.obj[[1]]$K
  y <- x.2$y   # true outlier index
  
  # repeat with each PCA methods
  for (m in 1:6) {
    print(paste(sim, "-", m))
    # reconstruction
    recon_mat <- predict(pca.est[[sim]]$pca.obj[[m]], K = NULL)
    
    # make matrix of complete curve
    X_comp <- x
    rownames(X_comp) <- 1:n
    colnames(X_comp) <- gr
    for (i in 1:n) {
      X_i <- x[i, ]
      X_i_hat <- recon_mat[i, ]
      
      idx_miss <- which(is.na(X_i))
      if (length(idx_miss) > 0) {
        X_comp[i, idx_miss] <- X_i_hat[idx_miss]
      }
    }
    
    # outlier detection using complete curves
    fds.obj <- fds(x = gr,
                   y = t(X_comp), 
                   xname = "Time", yname = "Value")
    fout <- list()
    fout[[1]] <- foutliers(data = fds.obj, method = "robMah")$outliers
    # fout[[2]] <- foutliers(data = fds.obj, method = "lrt")$outliers
    # fout[[3]] <- foutliers(data = fds.obj, method = "depth.trim")$outliers
    # fout[[4]] <- foutliers(data = fds.obj, method = "depth.pond")$outliers
    fout[[2]] <- as.numeric( foutliers(data = fds.obj, method = "HUoutliers")$outliers )
    
    
    # Outlier map (Score dist VS Orthogonal dist), See Hubert(2005)
    SD <- score_dist(pca.est[[sim]]$pca.obj[[m]])   # score distance
    OD <- orthogonal_dist(pca.est[[sim]]$pca.obj[[m]], x)   # orthogonal distance
    # mcd_fit <- covMcd(OD^(2/3))
    mcd_fit <- MASS::cov.mcd(matrix(OD^(2/3)))
    cut_y <- (mcd_fit$center + sqrt(as.numeric(mcd_fit$cov))*qnorm(0.975))^(3/2)
    cut_x <- sqrt(qchisq(0.975, k))
    fout[[3]] <- which(SD > cut_x & OD > cut_y)
    
    
    # summary outlier detection
    sens <- numeric(length(fout))
    spec <- numeric(length(fout))
    for (j in 1:length(fout)) {
      y_hat <- rep(0, n)
      y_hat[fout[[j]]] <- 1
      
      sens[j] <- caret::sensitivity(factor(y_hat, levels = c(1, 0)), 
                                    factor(y, levels = c(1, 0)))
      spec[j] <- caret::specificity(factor(y_hat, levels = c(1, 0)), 
                                    factor(y, levels = c(1, 0)))
    }
    
    # save result
    if (sim == 1) {
      names(sens) <- c("Mah","HU","SOmap")
      names(spec) <- c("Mah","HU","SOmap")
      sens_list[[m]] <- sens
      spec_list[[m]] <- spec
    } else {
      sens_list[[m]] <- rbind(sens_list[[m]],
                              sens)
      spec_list[[m]] <- rbind(spec_list[[m]],
                              spec)
    }
  }
}
sapply(sens_list, colMeans)
sapply(spec_list, colMeans)


