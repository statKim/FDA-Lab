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
source("R/sim_delaigle.R")
source("R/sim_Lin_Wang(2020).R")
source("R/sim_kraus.R")
source("R/sim_boente.R")
source("Kraus(2015)/pred.missfd.R")
source("Kraus(2015)/simul.missfd.R")
# source("robust_Kraus.R")
source("Boente_cov.R")
source("sig2_yao_rob.R")
source("cov_gk.R")
source("cov_pm.R")
source("function_outlier.R")



#####################################
### Simulation Parameters
#####################################
# num_sim <- 30   # number of simulations
# out_prop <- 0.2   # proportion of outliers
# model <- 4   # type of outliers
# data_type <- "partial"   # type of functional data
# kernel <- "epanechnikov"   # kernel function for local smoothing
# # kernel <- "gauss"   # kernel function for local smoothing
# bw_boente <- 0.1   # bandwidth for Boente(2020) - Error occurs for small bw
# n_cores <- 12   # number of threads for parallel computing
# pve <- 0.95   # Not used if K is given
# # fixed number of PCs (If NULL, it is selected by PVE)
# if (model == 3) {
#   K <- 5
# } else {
#   K <- 2
# }

### simulation type
sim_type <- "delaigle"
# sim_type <- "kraus"
# sim_type <- "boente"

### Overall setting
num_sim <- 30   # number of simulations
out_prop <- 0.2   # proportion of outliers
data_type <- "partial"   # type of functional data
n_cores <- 12   # number of threads for parallel computing
kernel <- "epanechnikov"   # kernel function for local smoothing
pve <- 0.95   # Not used if K is given



# simulation result
pca.est <- list()
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
  if (sim_type == "delaigle") {
    out_type <- 2   # type of outliers
    sig <- 0.1   # true noise variance
    x.2 <- sim_delaigle(n = n, 
                        model = 2,
                        type = data_type,
                        out.prop = out_prop, 
                        out.type = out_type,
                        noise = sig)
    K <- 4   # True number of PCs
  } else if (sim_type == "kraus") {
    out_type <- 2   # type of outliers
    sig <- 0.01   # true noise variance
    x.2 <- sim_kraus(n = n, 
                     type = data_type,
                     out.prop = out_prop, 
                     out.type = out_type,
                     noise = sig)
    K <- 5   # True number of PCs
  } else if (sim_type == "boente") {
    model <- 4   # outlier type for Boente(2020) setting
    x.2 <- sim_boente(n = n, 
                      model = model,
                      type = data_type,
                      eps1 = out_prop,
                      outlier = TRUE)
    # fixed number of PCs (If NULL, it is selected by PVE)
    if (model == 3) {
      K <- 5
    } else {
      K <- 2
    }
  }
  
  x <- list2matrix(x.2)
  # matplot(t(x), type = "l")
  
  #############################
  ### Covariance estimation
  #############################
  skip_sim <- FALSE   # if skip_sim == TRUE, pass this seed
  work.grid <- seq(0, 1, length.out = n.grid)
  
  ### Product moment(PM) correlation estimate
  start_time <- Sys.time()
  registerDoRNG(seed)
  tryCatch({
    noise.var.pm <- noise_var_pm(x)
    
    # # Not smoothed  
    # cov.obj <- cov_pm(x, noise.var = noise.var.pm)
    # mu.pm <- cov.obj$mean
    # cov.pm <- cov.obj$cov
    
    # Smoothed
    cov.obj <- cov_pm(x, smooth = TRUE, noise.var = noise.var.pm)
    mu.pm.sm <- cov.obj$mean
    cov.pm.sm <- cov.obj$cov
  }, error = function(e) { 
    print("PM cov error")
    print(e)
    skip_sim <<- TRUE
  })
  if (skip_sim == TRUE) {
    next
  }
  end_time <- Sys.time()
  print(paste0("PM : ", 
               round(difftime(end_time, start_time, units = "secs"), 3),
               " secs"))
  
  
  ### Gnanadesikan-Kettenring(GK) estimate
  start_time <- Sys.time()
  registerDoRNG(seed)
  tryCatch({
    noise.var.gk <- sigma2.rob.yao.gk(x)   # Yao(2005) like noise variance estimator
    
    # # Not smoothed GK
    # cov.obj <- cov_gk(x,
    #                   noise.var = noise.var.gk)
    # mu.gk <- cov.obj$mean
    # cov.gk <- cov.obj$cov
    
    # Smoothed GK
    cov.obj <- cov_gk(x, 
                      smooth = T,
                      noise.var = noise.var.gk)
    mu.gk.sm <- cov.obj$mean
    cov.gk.sm <- cov.obj$cov
  }, error = function(e) { 
    print("GK cov error")
    print(e)
    skip_sim <<- TRUE
  })
  if (skip_sim == TRUE) {
    next
  }
  end_time <- Sys.time()
  print(paste0("GK : ", 
               round(difftime(end_time, start_time, units = "secs"), 3),
               " secs"))
  
  
  ### M-estimator
  start_time <- Sys.time()
  registerDoRNG(seed)
  tryCatch({
    # noise_var <- sigma2.rob(x.2$Lt, x.2$Ly)   # robust noise variance estimator
    noise.var.Mest <- sigma2.rob.yao(x)   # Yao(2005) like noise variance estimator
    
    # # Not smoothed M-est
    # mu.Mest <- mean_Mest(x)
    # cov.Mest <- cov_Mest(x, noise.var = noise.var.Mest)
    
    # smoothed M-est
    mu.Mest.sm <- mean_Mest(x, smooth = TRUE)
    cov.Mest.sm <- cov_Mest(x, smooth = T,
                            noise.var = noise.var.Mest)
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
  start_time <- Sys.time()
  registerDoRNG(seed)
  kern <- ifelse(kernel == "epanechnikov", "epan", kernel)
  optns <- list(methodXi = "CE", dataType = "Sparse", kernel = kern, verbose = FALSE,
                # userBwMu = bw, userBwCov = bw)
                kFoldMuCov = 5, methodBwMu = "CV", methodBwCov = "CV", useBinnedCov = FALSE)
  tryCatch({
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
  
  
  
  ### Principal component analysis
  # Yao
  pca.yao.obj <- funPCA(x.2$Lt, x.2$Ly, 
                        mu.yao, cov.yao, sig2 = cov.yao.obj$sigma2, 
                        work.grid, PVE = pve, K = K)
  # M-est
  pca.Mest.sm.obj <- funPCA(x.2$Lt, x.2$Ly,
                            mu.Mest.sm, cov.Mest.sm, sig2 = noise.var.Mest,
                            work.grid, PVE = pve, K = K)
  # GK
  pca.gk.sm.obj <- funPCA(x.2$Lt, x.2$Ly,
                          mu.gk.sm, cov.gk.sm, sig2 = noise.var.gk,
                          work.grid, PVE = pve, K = K)
  # PM
  pca.pm.sm.obj <- funPCA(x.2$Lt, x.2$Ly,
                          mu.pm.sm, cov.pm.sm, sig2 = noise.var.pm,
                          work.grid, PVE = pve, K = K)
  
  
  num.sim <- num.sim + 1 
  print(paste0("Total # of simulations: ", num.sim))
  
  pca.est[[num.sim]] <- list(seed = seed,
                             x.2 = x.2,
                             work.grid = work.grid,
                             pca.obj = list(pca.yao.obj,
                                            pca.Mest.sm.obj,
                                            pca.gk.sm.obj,
                                            pca.pm.sm.obj))
}
# save(list = c("pca.est"),
#      file = "RData/20210911_detect_boente-4.RData")
# save(list = c("pca.est"),
#      file = "RData/20210911_detect_delaigle-2.RData")


### Outlier detection - PC score based method
## 1. 1st PC - adjbox
## 2. 1st PC - boxplot
## 3. SD - adjbox
## 4. SD - boxplot
## 5. robMah
## 6. Outlier map
mat <- matrix(NA, num_sim, 7)
colnames(mat) <- c("PC1-adjbox","PC1-box",
                   "SD-adjbox","SD-box",
                   "robMah-comp",
                   "robMah","outmap")
sens <- rep(list(mat), 4)
spec <- rep(list(mat), 4)
names(sens) <- c("Yao","Mest","GK","PM")
names(spec) <- c("Yao","Mest","GK","PM")
for (sim in 1:num_sim) {
  x.2 <- pca.est[[sim]]$x.2
  x <- list2matrix(x.2)
  
  ### Outlier detection
  for (i in 1:4) {
    out <- list()
    
    pc1 <- pca.est[[sim]]$pca.obj[[i]]$pc.score[, 1]
    y_hat <- rep(0, n)
    y_hat[which(SD %in% adjbox(pc1, plot = F)$out)] <- 1
    out[[1]] <- y_hat
    
    y_hat <- rep(0, n)
    y_hat[which(SD %in% boxplot(pc1, plot = F)$out)] <- 1
    out[[2]] <- y_hat
    
    
    # Boxplot for Score distance
    SD <- score_dist(pca.est[[sim]]$pca.obj[[i]])   # score distance
    y_hat <- rep(0, n)
    y_hat[which(SD %in% adjbox(SD, plot = F)$out)] <- 1
    out[[3]] <- y_hat
    
    y_hat <- rep(0, n)
    y_hat[which(SD %in% boxplot(SD, plot = F)$out)] <- 1
    out[[4]] <- y_hat
    

    # foutlier - robMah : completion based
    recon_mat <- predict(pca.est[[sim]]$pca.obj[[i]], K = NULL)   # reconstruction
    X_comp <- x
    rownames(X_comp) <- 1:n
    colnames(X_comp) <- gr
    for (row in 1:n) {
      X_i <- x[row, ]
      X_i_hat <- recon_mat[row, ]
      idx_miss <- which(is.na(X_i))
      if (length(idx_miss) > 0) {
        X_comp[row, idx_miss] <- X_i_hat[idx_miss]
      }
    }
    # outlier detection using complete curves
    fds.obj <- rainbow::fds(x = gr,
                            y = t(X_comp),
                            xname = "Time", yname = "Value")
    y_hat <- rep(0, n)
    y_hat[rainbow::foutliers(data = fds.obj, method = "robMah")$outliers] <- 1
    out[[5]] <- y_hat
    
    # Proposed method
    out[[6]] <- fun_outlier(pca.est[[sim]]$pca.obj[[i]], method = "robMah")
    out[[7]] <- fun_outlier(pca.est[[sim]]$pca.obj[[i]], method = "outmap")
    
    
    # Sensitivity and Specificity
    sens[[i]][sim, ] <- sapply(out, function(y_hat){
      caret::sensitivity(factor(y_hat, levels = c(1, 0)), 
                         factor(x.2$y, levels = c(1, 0)))
    })
    spec[[i]][sim, ] <- sapply(out, function(y_hat){
      caret::specificity(factor(y_hat, levels = c(1, 0)), 
                         factor(x.2$y, levels = c(1, 0)))
    })
  }
}

mapply(
  function(x, y) {
    rbind("sensitivity" = x,
          "specififity" = y) %>% 
      round(3) %>% 
      list()
  },
  lapply(sens, colMeans),
  lapply(spec, colMeans)
)







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
  OD <- orthogonal_dist(pca.est[[sim]]$pca.obj[[i]])   # orthogonal distance
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



################################################################
### Complete curve 사용하는 방법
################################################################
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
  for (m in 1:4) {
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

sapply(sens_list, colMeans)
sapply(spec_list, colMeans)
# > sapply(sens_list, colMeans)
# [,1]       [,2]       [,3]       [,4]
# Mah 0.06166667 0.06333333 0.07166667 0.07666667
# LRT 0.15500000 0.16500000 0.16833333 0.16666667
# DTR 0.00000000 0.00000000 0.00000000 0.00000000
# DWE 0.00000000 0.00000000 0.00000000 0.00000000
# HU  0.70333333 0.71000000 0.71500000 0.71666667
# > sapply(spec_list, colMeans)
# [,1]      [,2]      [,3]      [,4]
# Mah 0.9808333 0.9854167 0.9850000 0.9850000
# LRT 0.9991667 0.9991667 0.9991667 0.9987500
# DTR 1.0000000 1.0000000 1.0000000 1.0000000
# DWE 1.0000000 1.0000000 1.0000000 1.0000000
# HU  0.8054167 0.7991667 0.7912500 0.7966667



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


