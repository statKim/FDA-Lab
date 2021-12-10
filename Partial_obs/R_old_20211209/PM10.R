################################################
### Real data - PM10
################################################
library(GA)   # persp plot
library(fdapace)
# library(mcfda)   # 7
# library(synfd)   # 7
# library(doParallel)   # parallel computing
# library(doRNG)   # set.seed for foreach
# library(MASS)   # huber, rlm
library(tidyverse)
# library(latex2exp)
# library(xtable)
library(robfpca)
source("R/sim_boente.R")
source("Kraus(2015)/pred.missfd.R")
source("Kraus(2015)/simul.missfd.R")
source("robust_Kraus.R")
source("Boente_cov.R")
source("sig2_yao_rob.R")
source("cov_gk.R")
source("cov_pm.R")
source("R/sim_delaigle.R")
source("R/sim_Lin_Wang(2020).R")


#########################################
### Data Load
#########################################
load("/Users/hyunsung/Google 드라이브/Lab/KHS/partiall_obs/real_data/PM10/pm_data_korea.RData")
head(new_data)
dim(new_data)

var <- 3

loc <- unique(new_data[, 1])
loc <- as.numeric(loc)

### 3,4,5월 자료만 뽑고, missing이 많은 location만 뽑음
A <- vector()
for (loc_num in 1:length(loc)) {
  data <- new_data %>% 
    filter(new_data[, 1] == loc[loc_num])
  if (length(which(data[, 2] == '2017030101')) != 0) {
    time_index <- c(which(data[, 2] == '2017030101'):which(data[, 2] == '2017053124'))
    
    data <- data[time_index, ]
    A[loc_num] <- length(which(is.na(data[, var]) == T)) / nrow(data) * 100
  }
}

### missing 15% 이상인 location
rbind(which(A > 15), A[which(A > 15)])

### missing 15% 이상인 location 4개에 대해서 full_data 만듬
### 각 지역은 list에 들어있고, 각 list에는 day x hour matrix 들어있음
full_data <- list()
for (i in 1:4) {
  loc_num <- c(139, 240, 168, 228)[i]
  data <- new_data %>% 
    filter(new_data[, 1] == loc[loc_num])
  
  time_index <- c(which(data[, 2] == '2017030101'):which(data[, 2] == '2017053124'))
  
  full_data[[i]] <- matrix(data[time_index, var],
                           ncol = 24,
                           byrow = T)
}
dim(full_data[[1]])   # 92 x 24

### trajectories
par(mfrow = c(2, 2))
for (i in 1:4) {
  matplot(t(full_data[[i]]), 
          type = "l", 
          ylim = c(0, 300),
          col = gray(0.8), 
          lwd = 2,
          lty = 1,
          main = paste('location ID:', i),
          ylab = 'PM10',
          xlab = 'Hour')
  abline(h = 100, lty = 2)
}

# number of missing
sapply(full_data, function(x){ sum(is.na(x)) })







#############################
### Data of ith region
#############################
for (region in 1:4) {   # 4 regions
  # region <- 4
  x <- full_data[[region]]
  all_na_ind <- which( apply(x, 1, function(row){ sum(is.na(row)) == ncol(x) }) )
  x <- x[-all_na_ind, ]   # remove not observed day
  dim(x)
  x.2 <- matrix2list(x)
  
  n <- nrow(x)
  p <- ncol(x)
  work.grid <- seq(0, 1, length.out = p)
  
  #############################
  ### Covariance estimation
  #############################
  seed <- 100
  data_type <- "partial"   # type of functional data
  kernel <- "epanechnikov"   # kernel function for local smoothing
  # kernel <- "gauss"   # kernel function for local smoothing
  bw_boente <- 0.1   # bandwidth for Boente(2020) - Error occurs for small bw
  n_cores <- 12   # number of threads for parallel computing
  pve <- 0.95   # Not used if K is given
  K <- NULL   # fixed number of PCs (If NULL, it is selected by PVE)
  
  
  ### Yao, Müller, and Wang (2005)
  start_time <- Sys.time()
  set.seed(seed)
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
  end_time <- Sys.time()
  print(paste0("Yao et al. : ", 
               round(difftime(end_time, start_time, units = "secs"), 3),
               " secs"))
  
  
  ### Boente et al. (2020)
  start_time <- Sys.time()
  set.seed(seed)
  cov.boente.obj <- cov_boente(x.2, bw.mu = bw_boente, bw.cov = bw_boente)
  mu.boente <- cov.boente.obj$mu
  cov.boente <- cov.boente.obj$cov
  # noise var from source code of sparseFPCA package
  noise_boente <- eigen(cov.boente)$values[1] / (1e3 - 1)
  end_time <- Sys.time()
  print(paste0("Boente (2020) : ", 
               round(difftime(end_time, start_time, units = "secs"), 3),
               " secs"))
  
  
  
  ### M-estimator
  start_time <- Sys.time()
  set.seed(seed)
  # noise_var_Mest <- sigma2.rob(x.2$Lt, x.2$Ly)   # robust noise variance estimator
  noise.var.Mest <- sigma2.rob.yao(x)   # Yao(2005) like noise variance estimator
  print(noise.var.Mest)
  
  # Not smoothed
  mu.Mest <- mean_Mest(x)
  cov.Mest.noise <- cov_Mest(x, 
                             smooth = F,
                             noise.var = noise.var.Mest)
  
  # Smoothed
  # noise.var.Mest <- 1
  mu.Mest.sm <- mean_Mest(x, smooth = TRUE)
  cov.Mest.sm.noise <- cov_Mest(x, 
                                smooth = T,
                                noise.var = noise.var.Mest)
  end_time <- Sys.time()
  print(paste0("M-est : ", 
               round(difftime(end_time, start_time, units = "secs"), 3),
               " secs"))
  
  
  ### GK
  start_time <- Sys.time()
  set.seed(seed)
  noise.var.gk <- sigma2.rob.yao.gk(x)   # Yao(2005) like noise variance estimator
  print(noise.var.gk)
  # noise.var.gk <- 1
  
  # Not smoothed GK
  cov.obj <- cov_gk(x, 
                    noise.var = noise.var.gk)
  mu.gk <- cov.obj$mean
  cov.gk.noise <- cov.obj$cov
  
  # Smoothed GK
  cov.obj <- cov_gk(x, 
                    smooth = T,
                    noise.var = noise.var.gk)
  mu.gk.sm <- cov.obj$mean
  cov.gk.sm.noise <- cov.obj$cov
  end_time <- Sys.time()
  print(paste0("GK : ", 
               round(difftime(end_time, start_time, units = "secs"), 3),
               " secs"))
  
  
  ### PM
  start_time <- Sys.time()
  set.seed(seed)
  noise.var.pm <- noise_var_pm(x)
  print(noise.var.pm)
  
  # Not smoothed  
  cov.obj <- cov_pm(x, noise.var = noise.var.pm)
  mu.pm <- cov.obj$mean
  cov.pm.noise <- cov.obj$cov
  
  # Smoothed
  cov.obj <- cov_pm(x, smooth = TRUE, noise.var = noise.var.pm)
  mu.pm.sm <- cov.obj$mean
  cov.pm.sm.noise <- cov.obj$cov
  end_time <- Sys.time()
  print(paste0("PM : ", 
               round(difftime(end_time, start_time, units = "secs"), 3),
               " secs"))
  
  
  
  ### Principal component analysis
  pca.obj <- list()
  # Yao
  pca.obj[[1]] <- funPCA(x.2$Lt, x.2$Ly, 
                         mu.yao, cov.yao, sig2 = cov.yao.obj$sigma2, 
                         work.grid, PVE = pve, K = K)
  # Boente
  pca.obj[[2]] <- funPCA(x.2$Lt, x.2$Ly, 
                         mu.boente, cov.boente, sig2 = noise_boente, 
                         work.grid, PVE = pve, K = K)
  # # Mest
  # pca.Mest.obj <- funPCA(x.2$Lt, x.2$Ly,
  #                        mu.Mest, cov.Mest, sig2 = 0,
  #                        work.grid, PVE = pve, K = K)
  # pca.Mest.sm.obj <- funPCA(x.2$Lt, x.2$Ly,
  #                           mu.Mest.sm, cov.Mest.sm, sig2 = 0,
  #                           work.grid, PVE = pve, K = K)
  # Mest-noise
  pca.obj[[3]] <- funPCA(x.2$Lt, x.2$Ly,
                         mu.Mest, cov.Mest.noise, sig2 = noise.var.Mest,
                         work.grid, PVE = pve, K = K)
  pca.obj[[4]] <- funPCA(x.2$Lt, x.2$Ly,
                         mu.Mest.sm, cov.Mest.sm.noise, sig2 = noise.var.Mest,
                         work.grid, PVE = pve, K = K)
  # GK
  pca.obj[[5]] <- funPCA(x.2$Lt, x.2$Ly,
                         mu.gk, cov.gk.noise, sig2 = noise.var.gk,
                         work.grid, PVE = pve, K = K)
  pca.obj[[6]] <- funPCA(x.2$Lt, x.2$Ly,
                         mu.gk.sm, cov.gk.sm.noise, sig2 = noise.var.gk,
                         work.grid, PVE = pve, K = K)
  # PM
  pca.obj[[7]] <- funPCA(x.2$Lt, x.2$Ly,
                         mu.pm, cov.pm.noise, sig2 = noise.var.pm,
                         work.grid, PVE = pve, K = K)
  pca.obj[[8]] <- funPCA(x.2$Lt, x.2$Ly,
                         mu.pm.sm, cov.pm.sm.noise, sig2 = noise.var.pm,
                         work.grid, PVE = pve, K = K)
  
  # ### Eigen function
  # min_K <- sapply(pca.obj, function(obj){ obj$K }) %>% 
  #   min()
  # p_list <- list()
  # p_list2 <- list()
  # # par(mfcol = c(2, min_K))
  # for (i in 1:min_K) {
  #   # ith eigenfunction
  #   df <- sapply(pca.obj, function(obj){ obj$eig.fun[, i] })
  #   # align eigen sign with Mest-sm
  #   df <- apply(df, 2, function(col){ check_eigen_sign(col, df[, 4]) })
  #   
  #   df <- data.frame(df)
  #   colnames(df) <- c("Yao","Boente","Mest","Mest-sm",
  #                     "GK","GK-sm","PM","PM-sm")
  #   df$time <- 1:24
  #   p_list[[i]] <- df %>% 
  #     gather(key = "method",
  #            value = "val",
  #            -time) %>% 
  #     mutate(method = factor(method,
  #                            levels = c("Yao","Boente","Mest","Mest-sm",
  #                                       "GK","GK-sm","PM","PM-sm"))) %>% 
  #     ggplot(aes(x = time,
  #                y = val,
  #                color = method)) +
  #     geom_line(size = 0.7) +
  #     theme_bw() +
  #     labs(x = "Hour", y = "", title = paste0("Eigenfunction ", i)) +
  #     theme(legend.position = "bottom",
  #           legend.title = element_blank(),
  #           plot.title = element_text(hjust = 0.5),
  #           aspect.ratio = 1)
  #   
  #   p_list2[[i]] <- df[, c(1,2,4,6,8,9)] %>% 
  #     gather(key = "method",
  #            value = "val",
  #            -time) %>% 
  #     mutate(method = factor(method,
  #                            levels = c("Yao","Boente","Mest","Mest-sm",
  #                                       "GK","GK-sm","PM","PM-sm"))) %>% 
  #     ggplot(aes(x = time,
  #                y = val,
  #                color = method)) +
  #     geom_line(size = 0.7) +
  #     theme_bw() +
  #     labs(x = "Hour", y = "", title = paste0("Eigenfunction ", i)) +
  #     theme(legend.position = "bottom",
  #           legend.title = element_blank(),
  #           plot.title = element_text(hjust = 0.5),
  #           aspect.ratio = 1)
  #  
  #   
  #   # matplot(1:24, 
  #   #         df, 
  #   #         type = "l",
  #   #         col = 1:8,
  #   #         lty = 1,
  #   #         lwd = rep(2, 8),
  #   #         xlab = "Hour", 
  #   #         ylab = "", 
  #   #         main = paste0("Eigenfunction ", i))
  #   # grid()
  #   # if (i == 1) {
  #   #   legend("topleft",
  #   #          c("Yao","Boente","Mest","Mest-sm",
  #   #            "GK","GK-sm","PM","PM-sm"),
  #   #          col = 1:8,
  #   #          lty = 1,
  #   #          lwd = rep(2, 8))
  #   # }
  #   # 
  #   # matplot(1:24, 
  #   #         df[, c(1,2,4,6,8)], 
  #   #         type = "l",
  #   #         col = c(1,2,4,6,8),
  #   #         lty = 1,
  #   #         lwd = rep(2, 5),
  #   #         xlab = "Hour", 
  #   #         ylab = "", 
  #   #         main = paste0("Eigenfunction ", i))
  #   # grid()
  #   # if (i == 1) {
  #   #   legend("topleft",
  #   #          c("Yao","Boente","Mest-sm","GK-sm","PM-sm"),
  #   #          col = c(1,2,4,6,8),
  #   #          lty = 1,
  #   #          lwd = rep(2, 5))
  #   # }
  # }
  # 
  # library(ggpubr)
  # p1 <- ggarrange(plotlist = p_list, 
  #                 nrow = 1, ncol = 3,
  #                 common.legend = TRUE, legend = "bottom")
  # p2 <- ggarrange(plotlist = p_list2, 
  #                 nrow = 1, ncol = 3,
  #                 common.legend = TRUE, legend = "bottom")
  # ggarrange(p1, p2,
  #           nrow = 2)
  # 
  # 
  # 
  # 
  # gr <- work.grid
  # par(mfrow = c(2, 4))
  # GA::persp3D(gr, gr, cov.yao,
  #             theta = -70, phi = 30, expand = 1,
  #             xlab = "", ylab = "", zlab = "",
  #             main = "Yao")
  # GA::persp3D(gr, gr, cov.boente,
  #             theta = -70, phi = 30, expand = 1,
  #             xlab = "", ylab = "", zlab = "",
  #             main = "Boente")
  # GA::persp3D(gr, gr, cov.Mest.noise,
  #             theta = -70, phi = 30, expand = 1,
  #             xlab = "", ylab = "", zlab = "",
  #             main = "Mest")
  # GA::persp3D(gr, gr, cov.Mest.sm.noise,
  #             theta = -70, phi = 30, expand = 1,
  #             xlab = "", ylab = "", zlab = "",
  #             main = "Mest-sm")
  # GA::persp3D(gr, gr, cov.gk.noise,
  #             theta = -70, phi = 30, expand = 1,
  #             xlab = "", ylab = "", zlab = "",
  #             main = "GK")
  # GA::persp3D(gr, gr, cov.gk.sm.noise,
  #             theta = -70, phi = 30, expand = 1,
  #             xlab = "", ylab = "", zlab = "",
  #             main = "GK-sm")
  # GA::persp3D(gr, gr, cov.pm.noise,
  #             theta = -70, phi = 30, expand = 1,
  #             xlab = "", ylab = "", zlab = "",
  #             main = "PM")
  # GA::persp3D(gr, gr, cov.pm.sm.noise,
  #             theta = -70, phi = 30, expand = 1,
  #             xlab = "", ylab = "", zlab = "",
  #             main = "PM-sm")
  # 
  
  
  
  
  ################################################
  ### Outlier Detection
  ################################################
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
  
  
  # outlier <- list()
  # mname <- c("Yao","Boente","Mest","Mest-sm",
  #            "GK","GK-sm","PM","PM-sm")
  # # par(mfrow = c(2, 4))
  # p_list <- list()
  # for (i in 1:8) {
  #   k <- pca.obj[[i]]$K   # number of PCs
  #   SD <- score_dist(pca.obj[[i]])   # score distance
  #   OD <- orthogonal_dist(pca.obj[[i]], x)   # orthogonal distance
  #   # mcd_fit <- covMcd(OD^(2/3))
  #   mcd_fit <- cov.mcd(matrix(OD^(2/3)))
  #   cut_y <- (mcd_fit$center + sqrt(as.numeric(mcd_fit$cov))*qnorm(0.975))^(3/2)
  #   cut_x <- sqrt(qchisq(0.975, k))
  #   
  #   pt_col <- ifelse(SD >= cut_x | OD >= cut_y, 4, 1)
  #   pt_pch <- ifelse(SD >= cut_x | OD >= cut_y, 19, 1)
  #   
  #   outlier[[i]] <- which(SD >= cut_x & OD >= cut_y)
  #   
  #   
  #   p_list[[i]] <- qplot(SD, OD, 
  #         color = factor(pt_col)) +
  #     geom_point(size = 0.5) +
  #     scale_color_manual(values = c("grey","blue")) +
  #     geom_hline(yintercept = cut_y, color = "red") +
  #     geom_vline(xintercept = cut_x, color = "red") +
  #     theme_bw() +
  #     labs(x = "Score distance", y = "Orthogonal distance", title = mname[i]) +
  #     theme(legend.position = "none",
  #           plot.title = element_text(hjust = 0.5),
  #           aspect.ratio = 1) %>% 
  #     ggplotGrob()
  #     
  #   # plot(SD, OD,
  #   #      col = pt_col,
  #   #      pch = pt_pch,
  #   #      cex = 1.5,
  #   #      xlab = "Score distance",
  #   #      ylab = "Orthogonal distance",
  #   #      main = mname[i],
  #   #      cex.lab = 1.5,
  #   #      cex.main = 1.5)
  #   # grid()
  #   # abline(v = cut_x, col = 2)
  #   # abline(h = cut_y, col = 2)
  # }
  # 
  # library(ggpubr)
  # ggarrange(plotlist = p_list, 
  #           nrow = 2, ncol = 4)
  # 
  # 
  # # par(mfrow = c(2, 4))
  # p_list <- list()
  # for (i in 1:8) {
  #   outlier_col <- rep(8, nrow(x))
  #   outlier_col[ outlier[[i]] ] <- 4
  #   
  #   outlier_width <- rep(1, nrow(x))
  #   outlier_width[ outlier[[i]] ] <- 3
  #   
  #   
  #   p <- as.data.frame(cbind(t(x),
  #                                      time = 1:24)) %>% 
  #     gather(key = "ind",
  #            value = "val",
  #            -time) %>% 
  #     ggplot(aes(x = time,
  #                y = val,
  #                group = ind,
  #                color = rep(factor(outlier_col, levels = c(8, 4)), each = 24))) +
  #     geom_line(size = 0.7) +
  #     scale_color_manual(values = c("grey","blue")) +
  #     theme_bw() +
  #     labs(x = "Hour", 
  #          y = "PM10", 
  #          title = paste0(mname[i],
  #                         " (", sum(outlier_col == 4), ")")) +
  #     theme(legend.position = "none",
  #           plot.title = element_text(hjust = 0.5),
  #           aspect.ratio = 1)
  #   
  #   p_list[[i]] <- ggplotGrob(p)
  #   
  #   # matplot(t(x), 
  #   #         type = "l", 
  #   #         ylim = c(0, 300),
  #   #         col = outlier_col, 
  #   #         lwd = outlier_width,
  #   #         lty = 1,
  #   #         main = mname[i],
  #   #         ylab = 'PM10',
  #   #         xlab = 'Hour')
  #   # abline(h = 100, lty = 2)
  #   # grid()
  # }
  # ggarrange(plotlist = p_list, 
  #           nrow = 2, ncol = 4)
  # 
  
  
  
  
  ### Outlier detection
  library(rainbow)
  for (m in 1:8) {
    # reconstruction
    recon_mat <- predict(pca.obj[[m]], K = NULL)
    
    # make matrix of complete curve
    X_comp <- x
    rownames(X_comp) <- 1:n
    colnames(X_comp) <- work.grid
    for (i in 1:n) {
      X_i <- x[i, ]
      X_i_hat <- recon_mat[i, ]
      
      idx_miss <- which(is.na(X_i))
      if (length(idx_miss) > 0) {
        X_comp[i, idx_miss] <- X_i_hat[idx_miss]
      }
    }
    k <- pca.obj[[m]]$K   # number of PCs
    
    # outlier detection using complete curves
    fds.obj <- fds(x = work.grid,
                   y = t(X_comp), 
                   xname = "Time", yname = "Value")
    fout <- list()
    fout[[1]] <- foutliers(data = fds.obj, 
                           method = "robMah")$outliers
    fout[[2]] <- foutliers(data = fds.obj, 
                           method = "lrt")$outliers
    # fout[[3]] <- foutliers(data = fds.obj, method = "depth.trim")$outliers
    # fout[[4]] <- foutliers(data = fds.obj, method = "depth.pond")$outliers
    fout[[3]] <- as.numeric( foutliers(data = fds.obj, 
                                       method = "HUoutliers",
                                       order = k)$outliers )
    
    # PCA distance based outlier detection
    
    SD <- score_dist(pca.obj[[m]])   # score distance
    OD <- orthogonal_dist(pca.obj[[m]], x)   # orthogonal distance
    # mcd_fit <- covMcd(OD^(2/3))
    mcd_fit <- cov.mcd(matrix(OD^(2/3)))
    cut_y <- (mcd_fit$center + sqrt(as.numeric(mcd_fit$cov))*qnorm(0.975))^(3/2)
    cut_x <- sqrt(qchisq(0.975, k))
    fout[[4]] <- which(SD >= cut_x & OD >= cut_y)
    
    outlier[[m]] <- fout
  }
  
  
  # ### c("robMah","LRT","HU","PCA_dist") for each row
  # p_list2 <- list()
  # # par(mfcol = c(4, 8))
  # mname <- c("Yao","Boente","Mest","Mest-sm",
  #            "GK","GK-sm","PM","PM-sm")
  # # mname <- c("robMah","LRT","Depth_trim","Depth_pond","HU","PCA_dist")
  # # mname <- c("robMah","LRT","HU","PCA_dist")
  # for (m in 1:8) {
  #   fout <- outlier[[m]]
  #   
  #   p_list <- list()
  #   for (i in 1:4) {
  #     outlier_col <- rep(8, n)   # grey
  #     outlier_width <- rep(1, n)
  #     
  #     outlier_col[ fout[[i]] ] <- 4
  #     # outlier_width[ fout[[i]] ] <- 2
  #     
  #     p <- as.data.frame(cbind(t(x),
  #                                        time = 1:24)) %>% 
  #       gather(key = "ind",
  #              value = "val",
  #              -time) %>% 
  #       ggplot(aes(x = time,
  #                  y = val,
  #                  group = ind,
  #                  color = rep(factor(outlier_col, levels = c(8, 4)), each = 24))) +
  #       geom_line(size = 0.7) +
  #       geom_hline(yintercept = 100, color = "black", linetype = "dashed") +
  #       scale_color_manual(values = c("grey","blue")) +
  #       theme_bw() +
  #       labs(x = "Hour", 
  #            y = "PM10", 
  #            title = paste0(mname[m],
  #                           " (", sum(outlier_col != 8), ")")) +
  #       theme(legend.position = "none",
  #             plot.title = element_text(hjust = 0.5),
  #             aspect.ratio = 1)
  #     
  #     p_list[[i]] <- ggplotGrob(p)
  #     # matplot(t(x), 
  #     #         type = "l", 
  #     #         ylim = c(0, 300),
  #     #         col = outlier_col, 
  #     #         lwd = outlier_width,
  #     #         lty = 1,
  #     #         main = paste0(mname[m],
  #     #                       " (",
  #     #                       sum(outlier_col != 8),
  #     #                       ")"),
  #     #         ylab = 'PM10',
  #     #         xlab = 'Hour')
  #     # abline(h = 100, lty = 2)
  #     # grid()
  #   }
  #   
  #   p_list2[[m]] <- ggarrange(plotlist = p_list, 
  #                             nrow = 4, ncol = 1) %>% 
  #     ggplotGrob()
  # }
  # 
  # ggarrange(plotlist = p_list2, 
  #           nrow = 1, ncol = 8)
  
  ### save results
  save(list = c("full_data","pca.obj","outlier"),
       file = paste0("RData/20210906_pm10_", region, ".RData"))
}





