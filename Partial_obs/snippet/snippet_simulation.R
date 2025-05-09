################################################
### Simulation using Delaigle(2020) setting
### - Partially observed case
### - 5-fold CV is performed for hyperparameters
################################################
library(mvtnorm)
library(fdapace)   # 1, 2
library(mcfda)   # 7
library(synfd)   # 7
library(doParallel)   # parallel computing
library(doRNG)   # set.seed for foreach
library(MASS)   # huber, rlm
library(tidyverse)
# library(latex2exp)
# library(xtable)
library(robfpca)
source("R/sim_delaigle.R")
source("R/sim_delaigle2.R")
source("R/sim_lin.R")
source("Boente_cov.R")



#####################################
### Simulation Parameters
#####################################
num_sim <- 30   # number of simulations
out_prop <- 0.2   # proportion of outliers
out_type <- 2   # type of outliers
data_type <- "snippet"   # type of functional data
# kernel <- "epanechnikov"   # kernel function for local smoothing
kernel <- "gauss"   # kernel function for local smoothing
bw_boente <- 0.2   # bandwidth for Boente(2020) - Error occurs for small bw
n_cores <- 12   # number of threads for parallel computing
pve <- 0.95   # Not used if K is given

sim_type <- "delaigle"
# sim_type <- "lin"

### t dist case
dist <- "t"
out_prop <- 0

### Outlier 2 case
dist <- "normal"
out_prop <- 0

#####################################
### Simulation
#####################################
pca.est <- list()   # list containing PCA objects
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
    # x.2 <- sim_delaigle(n = n,
    #                     model = 2,
    #                     type = data_type,
    #                     out.prop = 0.2,
    #                     out.type = out_type,
    #                     noise = 0.1)
    x.2 <- sim_delaigle2(n = n,
                         model = 2,
                         type = data_type,
                         dist = dist,
                         out.prop = out_prop,
                         out.type = out_type,
                         noise = 0.1)
    K <- 4   # fixed number of PCs (If NULL, it is selected by PVE)
  } else if (sim_type == "lin") {
    x.2 <- sim_lin(n = n,
                   model = 1,
                   out.prop = out_prop,
                   out.type = 1,
                   noise = 0.1)
    K <- 5   # fixed number of PCs (If NULL, it is selected by PVE)
  }

  x <- list2matrix(x.2)
  # matplot(t(x), type = "l")
  # plot(x.2$Lt[[1]], x.2$Ly[[1]], type = "l", xlim = c(0, 1), ylim = range(unlist(x.2$Ly)))
  # for (i in 2:n) {
  #   lines(x.2$Lt[[i]], x.2$Ly[[i]])
  # }
  
  
  #############################
  ### Covariance estimation
  #############################
  skip_sim <- FALSE   # if skip_sim == TRUE, pass this seed
  # work.grid <- seq(0, 1, length.out = n.grid)
  t.range <- range(unlist(x.2$Lt))
  work.grid <- seq(t.range[1], t.range[2], length.out = n.grid)
  
  
  ### Lin & Wang (2020)
  start_time <- Sys.time()
  registerDoRNG(seed)
  tryCatch({
    # 5-fold CV (It took very long time when we use CV option in mcfda package.)
    cv.mu.lin.obj <- meanfunc.rob(x.2$Lt, x.2$Ly, method = "L2", kernel = kernel,
                                  cv_bw_loss = "L2", ncores = n_cores,
                                  cv = TRUE)
    cv.var.lin.obj <- varfunc.rob(x.2$Lt, x.2$Ly, mu = cv.mu.lin.obj, kernel = kernel,
                                  method = "L2",  cv_bw_loss = "L2", ncores = n_cores,
                                  cv = TRUE)
    # estimate mean, variance, covariance
    mu.lin.obj <- meanfunc(x.2$Lt, x.2$Ly, method = "PACE", kernel = kernel,
                           bw = cv.mu.lin.obj$bw)   # It occurs error or very slow.
    var.lin.obj <- varfunc(x.2$Lt, x.2$Ly, method = "PACE", kernel = kernel,
                           mu = mu.lin.obj, bw = cv.var.lin.obj$obj$bw)
    cov.lin.obj <- covfunc(x.2$Lt, x.2$Ly, method = "SP",
                           mu = mu.lin.obj, sig2x = var.lin.obj)
    mu.lin <- predict(mu.lin.obj, work.grid)
    cov.lin <- predict(cov.lin.obj, work.grid)
    
    # # too large noise variance estimate...
    # if (0 %in% diag(cov.lin)) {
    #   stop("0 variance occur!")
    # }
    
    # # estimate mean by local polynomial method
    # mu.lin.obj <- meanfunc(x$Lt, x$Ly, method = "PACE", 
    #                        kernel = kernel, bw = bw)
    # cov.lin.obj <- covfunc(x$Lt, x$Ly, mu = mu.lin.obj, method = "SP")
  }, error = function(e) { 
    print("Lin cov error")
    print(e)
    skip_sim <<- TRUE
  })
  if (skip_sim == TRUE) {
    next
  }
  end_time <- Sys.time()
  print(paste0("Lin & Wang : ", 
               round(difftime(end_time, start_time, units = "secs"), 3),
               " secs"))
  
  
  ### Robust + Lin
  start_time <- Sys.time()
  registerDoRNG(seed)
  tryCatch({
    ## Huber loss
    mu.huber.obj <- meanfunc.rob(x.2$Lt, x.2$Ly, 
                                 method = "huber", kernel = kernel, 
                                 cv = TRUE, ncores = n_cores)
    var.huber.obj <- varfunc.rob(x.2$Lt, x.2$Ly,
                                 method = "huber", kernel = kernel,
                                 mu = mu.huber.obj, 
                                 cv = TRUE, ncores = n_cores)
    cov.huber.obj <- covfunc.rob(x.2$Lt, x.2$Ly, 
                                 method = "huber", 
                                 kernel = kernel, 
                                 mu = mu.huber.obj,
                                 sig2x = var.huber.obj,
                                 cv = TRUE, 
                                 ncores = n_cores)
    print(c(mu.huber.obj$bw,
            var.huber.obj$obj$bw,
            cov.huber.obj$theta))
    mu.huber <- predict(mu.huber.obj, work.grid)
    cov.huber <- predict(cov.huber.obj, work.grid)
    
    ## Bisquare loss
    mu.bisquare.obj <- meanfunc.rob(x.2$Lt, x.2$Ly, 
                                    method = "bisquare", kernel = kernel, 
                                    cv = TRUE, ncores = n_cores)
    var.bisquare.obj <- varfunc.rob(x.2$Lt, x.2$Ly,
                                    method = "bisquare", kernel = kernel,
                                    mu = mu.bisquare.obj, 
                                    cv = TRUE, ncores = n_cores)
    cov.bisquare.obj <- covfunc.rob(x.2$Lt, x.2$Ly, 
                                    method = "bisquare", 
                                    kernel = kernel, 
                                    mu = mu.bisquare.obj,
                                    sig2x = var.bisquare.obj,
                                    cv = TRUE, 
                                    ncores = n_cores)
    print(c(mu.bisquare.obj$bw,
            var.bisquare.obj$obj$bw,
            cov.bisquare.obj$theta))
    mu.bisquare <- predict(mu.bisquare.obj, work.grid)
    cov.bisquare <- predict(cov.bisquare.obj, work.grid)
  }, error = function(e) { 
    print("Huber cov error")
    print(e)
    skip_sim <<- TRUE
  })
  if (skip_sim == TRUE) {
    next
  }
  end_time <- Sys.time()
  print(paste0("Proposed : ", 
               round(difftime(end_time, start_time, units = "secs"), 3),
               " secs"))
  
  
  
  ### Boente et al. (2020)
  start_time <- Sys.time()
  registerDoRNG(seed)
  tryCatch({
    # bw_boente <- 0.03
    # cov.boente.obj <- cov_boente(x.2, bw.mu = bw_boente, bw.cov = bw_boente)
    cov.boente.obj <- cov_boente(x.2, cv = TRUE, kern = kernel, seed = seed)   # 5-fold CV
    mu.boente <- cov.boente.obj$mu
    cov.boente <- cov.boente.obj$cov
    # noise var from source code of sparseFPCA package
    noise_boente <- cov.boente.obj$noise_var
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
  
  
  
  ### Yao, Müller, and Wang (2005)
  start_time <- Sys.time()
  registerDoRNG(seed)
  kern <- ifelse(kernel == "epanechnikov", "epan", kernel)
  optns <- list(methodXi = "CE", dataType = "Sparse", kernel = kern, verbose = FALSE,
                # userBwMu = 0.3, userBwCov = 0.4)
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
  
  
  # gr <- work.grid
  # if (sim_type == "delaigle") {
  #   cov.true <- get_cov_fragm(gr)
  # } else if (sim_type == "lin") {
  #   cov.true <- get_cov_lin(gr)
  # }
  # par(mfrow = c(2, 3))
  # GA::persp3D(gr, gr, cov.true,
  #             main = "True",
  #             theta = -70, phi = 30, expand = 1)
  # GA::persp3D(gr, gr, cov.yao,
  #             theta = -70, phi = 30, expand = 1)
  # GA::persp3D(gr, gr, cov.lin,
  #             theta = -70, phi = 30, expand = 1)
  # GA::persp3D(gr, gr, cov.boente,
  #             theta = -70, phi = 30, expand = 1)
  # GA::persp3D(gr, gr, cov.huber,
  #             theta = -70, phi = 30, expand = 1)
  # GA::persp3D(gr, gr, cov.bisquare,
  #             theta = -70, phi = 30, expand = 1)
  # par(mfrow = c(1, 1))
  
  
  # if some covariances is a not finite value
  if (!is.finite(sum(cov.yao)) | !is.finite(sum(cov.lin)) | 
      !is.finite(sum(cov.boente)) |
      !is.finite(sum(cov.huber)) | !is.finite(sum(cov.bisquare))) {
    cat("Estimated covariances do not have finite values. \n")
    next
  }
  
  # if all covariances are 0
  if ((sum(cov.yao) == 0) | (sum(cov.lin) == 0) | 
      (sum(cov.boente) == 0) |
      (sum(cov.huber) == 0) | (sum(cov.bisquare) == 0)) {
    cat("Estimated covariance have all 0 values. \n")
    next
  }
  
  
  ### Principal component analysis
  # Yao
  pca.yao.obj <- funPCA(x.2$Lt, x.2$Ly, 
                        mu.yao, cov.yao, sig2 = cov.yao.obj$sigma2, 
                        work.grid, PVE = pve, K = K)
  # Lin
  pca.lin.obj <- funPCA(x.2$Lt, x.2$Ly, 
                        mu.lin, cov.lin, sig2 = cov.lin.obj$sig2e, 
                        work.grid, PVE = pve, K = K)
  # Boente
  pca.boente.obj <- funPCA(x.2$Lt, x.2$Ly, 
                           mu.boente, cov.boente, sig2 = noise_boente, 
                           work.grid, PVE = pve, K = K)
  # Huber
  pca.huber.obj <- funPCA(x.2$Lt, x.2$Ly, 
                          mu.huber, cov.huber, sig2 = cov.huber.obj$sig2e, 
                          work.grid, PVE = pve, K = K)
  # Bisquare
  pca.bisquare.obj <- funPCA(x.2$Lt, x.2$Ly, 
                             mu.bisquare, cov.bisquare, sig2 = cov.bisquare.obj$sig2e, 
                             work.grid, PVE = pve, K = K)
  
  ### Save PCA object
  num.sim <- num.sim + 1 
  print(paste0("Total # of simulations: ", num.sim))
  pca.est[[num.sim]] <- list(seed = seed,
                             x.2 = x.2,
                             work.grid = work.grid,
                             pca.obj = list(pca.yao.obj = pca.yao.obj,
                                            pca.lin.obj = pca.lin.obj,
                                            pca.boente.obj = pca.boente.obj,
                                            pca.huber.obj = pca.huber.obj,
                                            pca.bisquare.obj = pca.bisquare.obj))
}
if (sim_type == "delaigle") {
  save(list = c("pca.est"),
       file = paste0("RData/pca-snippet-", dist, "-", out_prop, ".RData"))
} else if (sim_type == "lin") {
  save(list = c("pca.est"),
       file = paste0("RData/pca-snippet-lin-", out_prop, ".RData"))
}





###############################
### Result summary
###############################
mse_eigen <- matrix(NA, num_sim, 5)
colnames(mse_eigen) <- c("Yao","Lin","Boente","Huber","Bisquare")
mse_score <- mse_eigen
# mse_score1 <- mse_eigen
mse_var <- mse_eigen
mse_cov <- mse_eigen
mse_intra <- mse_eigen
mse_extra <- mse_eigen
m_pve <- mse_eigen

for (num.sim in 1:num_sim) {
  ### Get generated data
  x.2 <- pca.est[[num.sim]]$x.2
  x <- list2matrix(x.2)
  work.grid <- pca.est[[num.sim]]$work.grid
  
  ### Get funPCA.obj
  pca.yao.obj <- pca.est[[num.sim]]$pca.obj$pca.yao.obj
  pca.lin.obj <- pca.est[[num.sim]]$pca.obj$pca.lin.obj
  pca.boente.obj <- pca.est[[num.sim]]$pca.obj$pca.boente.obj
  pca.huber.obj <- pca.est[[num.sim]]$pca.obj$pca.huber.obj
  pca.bisquare.obj <- pca.est[[num.sim]]$pca.obj$pca.bisquare.obj
  
  ### Combine PCA objects
  pca.obj <- list(
    pca.yao.obj,
    pca.lin.obj,
    pca.boente.obj,
    pca.huber.obj,
    pca.bisquare.obj
  )
  
  ### Average of PVE
  m_pve[num.sim, ] <- sapply(pca.obj, function(method){ method$PVE })
  
  ### Eigen function - Compute for fixed K
  if (sim_type == "delaigle") {
    cov.true <- get_cov_fragm(work.grid)
    eig.true <- get_delaigle_eigen(work.grid, model = 2)
  } else if (sim_type == "lin") {
    cov.true <- get_cov_lin(work.grid)
    eig.true <- get_eigen(cov.true, work.grid)$phi[, 1:K]
  }
  # calculate MSE
  mse_eigen[num.sim, ] <- sapply(pca.obj, function(method){
    mean((check_eigen_sign(method$eig.fun, eig.true) - eig.true)^2)
  })
  
  ### PC score
  if (sim_type == "delaigle") {
    pc.true <- x.2$xi
    if (dist == "t") {
      out.ind <- apply(pc.true, 2, function(pc){ 
        which(pc %in% boxplot(pc, plot = F)$out) 
      }) %>% 
        unlist() %>% 
        unique()
      ind <- which(1:n %in% out.ind)
    } else {
      ind <- which(x.2$y == 0)
    }
    mse_score[num.sim, ] <- sapply(pca.obj, function(method){
      mean((method$pc.score - pc.true)[ind, ]^2)
    })
    # mse_score1[num.sim, ] <- sapply(pca.obj, function(method){
    #   mean((method$pc.score[1, ] - pc.true[1, ])^2)
    # })
  } else if (sim_type == "lin") {
    mse_score[num.sim, ] <- 0
  }

  
  ### Variance and Covariance
  mse_var[num.sim, ] <- sapply(pca.obj, function(method){
    mean((diag(method$cov) - diag(cov.true))^2)
  })
  mse_cov[num.sim, ] <- sapply(pca.obj, function(method){
    mean((method$cov - cov.true)^2)
  })
  
  ### Intra and Extrapolated region - Delaigle(2020)
  design.idx <- get_design_index(x.2$Lt, work.grid)
  mse_intra[num.sim, ] <- sapply(pca.obj, function(method){
    get_ise(cov_inter(cov.true, design.idx), 
            cov_inter(method$cov, design.idx), 
            work.grid)
  })
  mse_extra[num.sim, ] <- sapply(pca.obj, function(method){
    get_ise(cov_extra(cov.true, design.idx), 
            cov_extra(method$cov, design.idx), 
            work.grid)
  })
  
  
  # par(mfrow = c(2, 2))
  # for (i in 1:4) {
  #   matplot(work.grid,
  #           cbind(
  #             eig.true[, i],
  #             check_eigen_sign(pca.boente.obj$eig.fun[, i], eig.true[, i]),
  #             check_eigen_sign(pca.huber.obj$eig.fun[, i], eig.true[, i]),
  #             check_eigen_sign(pca.bisquare.obj$eig.fun[, i], eig.true[, i])
  #           ),
  #           type = "l")
  # }
}


data.frame(Method = c("Yao","Lin","Boente",
                      "Huber","Bisquare")) %>% 
  left_join(data.frame(
    Method = colnames(mse_eigen),
    PVE = format(round(colMeans(m_pve), 2), 2)
  ), by = "Method") %>%
  left_join(data.frame(
    Method = colnames(mse_eigen),
    Eigen = paste0(
      format(round(sqrt( colMeans(mse_eigen) ), 2), 2),
      " (",
      format(round(sqrt( apply(mse_eigen, 2, sd) ), 2), 2),
      ")"
    )
  ), by = "Method") %>% 
  left_join(data.frame(
    Method = colnames(mse_score),
    Score = paste0(
      format(round(sqrt( colMeans(mse_score) ), 2), 2),
      " (",
      format(round(sqrt( apply(mse_score, 2, sd) ), 2), 2),
      ")"
    )
  ), by = "Method") %>% 
  left_join(data.frame(
    Method = colnames(mse_var),
    Var = paste0(
      format(round(sqrt( colMeans(mse_var) ), 2), 2),
      " (",
      format(round(sqrt( apply(mse_var, 2, sd) ), 2), 2),
      ")"
    )
  ), by = "Method") %>% 
  left_join(data.frame(
    Method = colnames(mse_cov),
    Cov = paste0(
      format(round(sqrt( colMeans(mse_cov) ), 2), 2),
      " (",
      format(round(sqrt( apply(mse_cov, 2, sd) ), 2), 2),
      ")"
    )
  ), by = "Method") %>% 
  left_join(data.frame(
    Method = colnames(mse_intra),
    Intra = paste0(
      format(round(sqrt( colMeans(mse_intra) ), 2), 2),
      " (",
      format(round(sqrt( apply(mse_intra, 2, sd) ), 2), 2),
      ")"
    )
  ), by = "Method") %>% 
  left_join(data.frame(
    Method = colnames(mse_extra),
    Extra = paste0(
      format(round(sqrt( colMeans(mse_extra) ), 2), 2),
      " (",
      format(round(sqrt( apply(mse_extra, 2, sd) ), 2), 2),
      ")"
    )
  ), by = "Method") %>% 
  print()





### Eigen function trajectories
par(mfrow = c(2, ceiling(K/2)))
for (k in 1:K) {
  matplot(work.grid,
          cbind(
            eig.true[, k],
            check_eigen_sign(pca.yao.obj$eig.fun, eig.true)[, k],
            check_eigen_sign(pca.lin.obj$eig.fun, eig.true)[, k],
            check_eigen_sign(pca.boente.obj$eig.fun, eig.true)[, k],
            check_eigen_sign(pca.huber.obj$eig.fun, eig.true)[, k],
            check_eigen_sign(pca.bisquare.obj$eig.fun, eig.true)[, k]
          ),
          type = "l",
          col = 1:6,
          lty = 1:6,
          main = paste("Eigenfunction", k),
          xlab = "", ylab = "",
          lwd = rep(2, 6))
  if (k == 1) {
    legend("topleft",
           c("True","Yao","Lin","Boente","Huber","Bisquare"),
           col = 1:6,
           lty = 1:6,
           lwd = rep(2, 6))
  }
}




source("function_outlier.R")
sim <- 30
x.2 <- pca.est[[sim]]$x.2
x <- list2matrix(x.2)
k <- pca.est[[sim]]$pca.obj[[1]]$K
mname <- c("Yao","Lin","Boente","Huber","Bisquare")
par(mfrow = c(2, 3))
for (i in 1:5) {
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
       # xlim = c(0, 1), ylim = c(0, 60),
       cex.lab = 1.5,
       cex.main = 1.5)
  grid()
  abline(v = cut_x, col = 2)
  abline(h = cut_y, col = 2)
}



pca.yao.obj$sig2
pca.lin.obj$sig2
pca.boente.obj$sig2
pca.huber.obj$sig2
pca.bisquare.obj$sig2



domain <- range(unlist(Lt))   # range of timepoints
min_bw <- 2*max(unlist(lapply(Lt, diff)))   # minimun candidate of bw
bw_cand <- seq(min_bw, diff(domain)/3, length.out = 10)
