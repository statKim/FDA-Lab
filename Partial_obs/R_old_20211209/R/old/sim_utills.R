##########################################
### Simulation summary functions
##########################################
source("R/sim_Delaigle(2020).R")
source("R/sim_Lin_Wang(2020).R")


### Summary the ISE for simulations
summary_ise <- function(data.list, cov.est, method = "var") {
  
  if (method %in% c("var","cov")) {   ## validation for Lin & Wang(2020)
    if (method == "var") {   # variance
      ise.cov <- sapply(cov.est, function(x) {
        num_method <- length(x$cov) - 1
        ise <- numeric(num_method)
        for (j in 1:num_method) {
          ise[j] <- get_ise(diag(x$cov$true), diag(x$cov[[j+1]]), x$work.grid)
        }
        return(ise)
        # c(get_ise(diag(x$cov$true), diag(x$cov$yao), x$work.grid),
        #   get_ise(diag(x$cov$true), diag(x$cov$lin), x$work.grid),
        #   get_ise(diag(x$cov$true), diag(x$cov$huber), x$work.grid))
      })
    } else if (method == "cov") {   # covariance
      ise.cov <- sapply(cov.est, function(x) {
        num_method <- length(x$cov) - 1
        ise <- numeric(num_method)
        for (j in 1:num_method) {
          ise[j] <- get_ise(x$cov$true, x$cov[[j+1]], x$work.grid)
        }
        return(ise)
        # c(get_ise(x$cov$true, x$cov$yao, x$work.grid),
        #   get_ise(x$cov$true, x$cov$lin, x$work.grid),
        #   get_ise(x$cov$true, x$cov$huber, x$work.grid))
      })
    } 
  } else if (method %in% c("intra","extra")) {   ## validation for Delaigle(2020)
    ise.cov <- mapply(function(x, y) {
      cov.list <- x$cov
      ind <- get_design_index(y$x$Lt)
      num_method <- length(cov.list) - 1
      ise <- numeric(num_method)
      
      if (method == "intra") {   # intrapolated part
        for (j in 1:num_method) {
          ise[j] <- get_ise(cov_inter(cov.list$true, ind), 
                            cov_inter(cov.list[[j+1]], ind), 
                            x$work.grid)
        }
        # cov.true <- cov_inter(cov.list$true, ind)
        # cov.yao <- cov_inter(cov.list$yao, ind)
        # cov.lin <- cov_inter(cov.list$lin, ind)
        # cov.huber <- cov_inter(cov.list$huber, ind)
      } else if (method == "extra") {   # extrapolated part
        for (j in 1:num_method) {
          ise[j] <- get_ise(cov_extra(cov.list$true, ind), 
                            cov_extra(cov.list[[j+1]], ind), 
                            x$work.grid)
        }
        # cov.true <- cov_extra(cov.list$true, ind)
        # cov.yao <- cov_extra(cov.list$yao, ind)
        # cov.lin <- cov_extra(cov.list$lin, ind)
        # cov.huber <- cov_extra(cov.list$huber, ind)
      }
      
      return(ise)
      # c(get_ise(cov.true, cov.yao, x$work.grid),
      #   get_ise(cov.true, cov.lin, x$work.grid),
      #   get_ise(cov.true, cov.huber, x$work.grid))
    },
    cov.est, data.list)
  }
  
  return(ise.cov)
}



### Get PCA results for each simulation datasets
sim_eigen_result <- function(cov.est, num.sim, seed = 1000) {
  # list of packages
  packages <- c("fdapace","mcfda","synfd")
  
  registerDoRNG(seed)
  pca.est <- foreach(sim = 1:num.sim, .packages = packages, 
                     .export = c("get_eigen","check_eigen_sign")) %dopar% {
    # estimated covariances from Simulation 3
    work.grid <- cov.est[[sim]]$work.grid
    cov.true <- cov.est[[sim]]$cov$true
    cov.yao <- cov.est[[sim]]$cov$yao
    cov.lin <- cov.est[[sim]]$cov$lin
    cov.huber <- cov.est[[sim]]$cov$huber
    cov.wrm <- cov.est[[sim]]$cov$wrm
    
    # eigen analysis
    eig.true <- get_eigen(cov = cov.true, grid = work.grid)
    eig.yao <- get_eigen(cov = cov.yao, grid = work.grid)
    eig.lin <- get_eigen(cov = cov.lin, grid = work.grid)
    eig.huber <- get_eigen(cov = cov.huber, grid = work.grid)
    eig.wrm <- get_eigen(cov = cov.wrm, grid = work.grid)
    
    # change eigen direction(sign) for first K eigenvectors
    K <- min(ncol(eig.true$phi),
             ncol(eig.yao$phi),
             ncol(eig.lin$phi),
             ncol(eig.huber$phi),
             ncol(eig.wrm$phi))
    eig.yao$phi[, 1:K] <- check_eigen_sign(eig.yao$phi[, 1:K], eig.true$phi[, 1:K])
    eig.lin$phi[, 1:K] <- check_eigen_sign(eig.lin$phi[, 1:K], eig.true$phi[, 1:K])
    eig.huber$phi[, 1:K] <- check_eigen_sign(eig.huber$phi[, 1:K], eig.true$phi[, 1:K])
    eig.wrm$phi[, 1:K] <- check_eigen_sign(eig.wrm$phi[, 1:K], eig.true$phi[, 1:K])
    
    # output list
    out <- list(work.grid = work.grid,
                true = eig.true,
                yao = eig.yao,
                lin = eig.lin,
                huber = eig.huber,
                wrm = eig.wrm)
    
    return(out)
  }
  return(pca.est)
}



##########################################
### Visualization functions for simulations
##########################################
require(tidyverse)
## Visualize the variance functions
ggplot_var <- function(cov.est, sim, main = "Outlier 1") {
  fig <- list()   # figure list
  
  # remove list contating "null"  
  ind <- which(!sapply(cov.est, is.null))
  # num.sim <- length(ind)   # number of simulations
  # data.list <- data.list[ind]
  cov.est <- cov.est[ind]

  # time points for simulated data  
  work.grid <- cov.est[[sim]]$work.grid
  cov.list <- cov.est[[sim]]$cov
  
  # dataframe for ggplot2
  df <- data.frame(
    x = rep(work.grid, 5),
    y = c(diag(cov.list$true),
          diag(cov.list$yao),
          diag(cov.list$lin),
          diag(cov.list$huber),
          diag(cov.list$wrm)),
    method = factor(rep(c("True","Yao","Lin","Huber","WRM"),
                        each = length(work.grid)),
                    levels = c("True","Yao","Lin","Huber","WRM"))
    # x = rep(work.grid, 4),
    # y = c(diag(cov.list$true),
    #       diag(cov.list$yao),
    #       diag(cov.list$lin),
    #       diag(cov.list$huber)),
    # method = factor(rep(c("True","Yao(2005)","Lin(2020)","Lin + Huber loss"),
    #                     each = length(work.grid)),
    #                 levels = c("True","Yao(2005)","Lin(2020)","Lin + Huber loss"))
  )
  
  # figure of all methods
  fig[[1]] <- ggplot(df, aes(x, y, group = method, color = method, linetype = method)) +
    geom_line(size = 1) +
    labs(x = TeX("$t$"), y = TeX("$\\sigma_X^2(t)$"), title = main) +
    # geom_hline(yintercept = 0, size = 0.8) +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.position = "bottom")
  
  fig[[2]] <- fig[[1]] + 
    ylim(0, 5) +
    labs(title = "")
  # # figure of true and a proposed method
  # fig[[2]] <- df %>% 
  #   # filter(method %in% c("True","Huber")) %>% 
  #   filter(method %in% c("True","Huber","WRM")) %>% 
  #   ggplot(aes(x, y, group = method, color = method)) +
  #   geom_line(size = 1) +
  #   labs(x = TeX("$t$"), y = TeX("$\\sigma_X^2(t)$"), title = "") +
  #   geom_hline(yintercept = 0, size = 0.8) +
  #   theme_bw() +
  #   theme(legend.title = element_blank(),
  #         legend.position = "bottom")
  
  return(fig)
}


## Visualize the covaraince surfaces
plot_cov_surf <- function(cov.est, sim, title = TRUE, lab = "Outlier 1") {
  # remove list contating "null"  
  ind <- which(!sapply(cov.est, is.null))
  # num.sim <- length(ind)   # number of simulations
  # data.list <- data.list[ind]
  cov.est <- cov.est[ind]
  
  # time points for simulated data  
  work.grid <- cov.est[[sim]]$work.grid
  cov.list <- cov.est[[sim]]$cov
  
  if (isTRUE(title)) {
    # main <- c("True","Yao et al. (2005)","Lin & Wang (2020)","Lin + Huber loss")
    main <- c("True","Yao et al. (2005)","Lin & Wang (2020)","Huber Loss","WRM")
  } else {
    main <- c("","","","")
  }
  
  ### Covariance surface
  # par(mfrow = c(1, 3),
  #     mar = c(0, 2, 7, 2))
  persp3D(work.grid, work.grid, cov.list$true,
          theta = -70, phi = 30, expand = 1,
          main = main[1],
          xlab = "s", ylab = "t", zlab = "C(s,t)")
  mtext(lab, side = 2)
  persp3D(work.grid, work.grid, cov.list$yao, 
          theta = -70, phi = 30, expand = 1,
          main = main[2], 
          xlab = "s", ylab = "t", zlab = "C(s,t)")
  # mtext(lab, side = 2)
  persp3D(work.grid, work.grid, cov.list$lin, 
          theta = -70, phi = 30, expand = 1,
          main = main[3],
          xlab = "s", ylab = "t", zlab = "C(s,t)")
  persp3D(work.grid, work.grid, cov.list$huber, 
          theta = -70, phi = 30, expand = 1,
          main = main[4],
          xlab = "s", ylab = "t", zlab = "C(s,t)")
  persp3D(work.grid, work.grid, cov.list$wrm, 
          theta = -70, phi = 30, expand = 1,
          main = main[5],
          xlab = "s", ylab = "t", zlab = "C(s,t)")
}

## Visualize the first 3 eigenfunctions
ggplot_eig <- function(cov.est, sim, main = "Outlier 1") {
  fig <- list()   # figure list
  
  # remove list contating "null"  
  ind <- which(!sapply(cov.est, is.null))
  cov.est <- cov.est[ind]
  
  # estimated covariances from Simulation 3
  work.grid <- cov.est[[sim]]$work.grid
  cov.true <- cov.est[[sim]]$cov$true
  cov.yao <- cov.est[[sim]]$cov$yao
  cov.lin <- cov.est[[sim]]$cov$lin
  cov.huber <- cov.est[[sim]]$cov$huber  
  cov.wrm <- cov.est[[sim]]$cov$wrm  
  
  # eigen analysis
  eig.true <- get_eigen(cov = cov.true, grid = work.grid)
  eig.yao <- get_eigen(cov = cov.yao, grid = work.grid)
  eig.lin <- get_eigen(cov = cov.lin, grid = work.grid)
  eig.huber <- get_eigen(cov = cov.huber, grid = work.grid)
  eig.wrm <- get_eigen(cov = cov.wrm, grid = work.grid)
  
  # change eigen direction(sign) for first K eigenvectors
  K <- 3
  eig.yao$phi[, 1:K] <- check_eigen_sign(eig.yao$phi[, 1:K], eig.true$phi[, 1:K])
  eig.lin$phi[, 1:K] <- check_eigen_sign(eig.lin$phi[, 1:K], eig.true$phi[, 1:K])
  eig.huber$phi[, 1:K] <- check_eigen_sign(eig.huber$phi[, 1:K], eig.true$phi[, 1:K])
  eig.wrm$phi[, 1:K] <- check_eigen_sign(eig.wrm$phi[, 1:K], eig.true$phi[, 1:K])
  
  # fitst 3 eigenfunctions
  for (i in 1:K) {
    if (i == 1) {
      # p_title <- paste("Outlier", k)
      p_title <- main
    } else {
      p_title <- ""
    }
    
    fig.data <- data.frame(
      work.grid = rep(work.grid, 5),
      phi = c(eig.true$phi[, i],
              eig.yao$phi[, i],
              eig.lin$phi[, i],
              eig.huber$phi[, i],
              eig.wrm$phi[, i]),
      method = factor(rep(c("True","Yao","Lin","Huber","WRM"),
                          each = length(work.grid)),
                      levels = c("True","Yao","Lin","Huber","WRM"))
      # work.grid = rep(work.grid, 4),
      # phi = c(eig.true$phi[, i],
      #         eig.yao$phi[, i],
      #         eig.lin$phi[, i],
      #         eig.huber$phi[, i]),
      # method = rep(c("True","Yao (2005)","Lin (2020)","Lin + Huber"), 
      #              each = length(work.grid))
    )
    fig[[i]] <- ggplot(data = fig.data, 
                             mapping = aes(work.grid, phi, color = method, linetype = method)) +
      geom_line(size = 1) +
      labs(x = TeX("$t$"), y = TeX(paste0("$\\phi_", i, "(t)$")), title = p_title) +
      # scale_color_discrete(breaks = c("True","Yao","Lin","Huber","WRM")) +
      theme_bw() +
      theme(legend.title = element_blank(),
            legend.position = "bottom")
  }
  
  return(fig)
}





# ##########################################
# ### Functions for simulations
# ##########################################
# cov_est_sim <- function(data.list, num.sim = NULL, seed = 1000, kernel = "epanechnikov") {
#   if (is.null(num.sim)) {
#     num.sim <- length(data.list)   # number of simulations
#   }
#   
#   # list of manual functions and packages
#   ftns <- fun2char()
#   packages <- c("fdapace","mcfda","synfd")
#   
#   model.cov <- 2   # covariance function setting of the paper (1, 2)
#   
#   registerDoRNG(seed)
#   # cov.est <- foreach(sim = 1:num.sim, .packages = packages, .export = ftns) %do% {
#   cov.est <- list()
#   for (sim in 1:num.sim) {
#     print(paste0(sim, "th simulation:"))
#     start_time <- Sys.time()
#     
#     # Get simulation data
#     x <- data.list[[sim]]$x
#     gr <- data.list[[sim]]$gr
#     
#     ### Covariance estimation
#     work.grid <- seq(min(gr), max(gr), length.out = 51)
#     cov.true <- get_cov_fragm(gr, model = model.cov)   # true covariance
#     cov.true <- ConvertSupport(fromGrid = gr, toGrid = work.grid, Cov = cov.true)   # transform to the observed grid
#     
#     ## 1. Yao, Müller, and Wang (2005)
#     ## 2. Liu and Müller (2009) - fitted.FPCA()
#     x.2 <- list(Ly = x$y,
#                 Lt = x$t)
#     kern <- ifelse(kernel == "epanechnikov", "epan", kernel)
#     optns <- list(methodXi = 'CE', dataType = 'Sparse', kernel = kern, verbose = FALSE)
#     mu.yao.obj <- GetMeanCurve(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)   # get bandwidth of mean estimation
#     cov.yao.obj <- GetCovSurface(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)
#     cov.yao <- cov.yao.obj$cov
#     if (length(work.grid) != 51) {
#       cov.yao <- ConvertSupport(fromGrid = cov.yao.obj$workGrid, toGrid = work.grid,
#                                 Cov = cov.yao.obj$cov)   # transform to the observed grid
#     }
#     
#     
#     ## 7. Lin & Wang (2020)
#     # estimate mean by local polynomial method
#     mu.lin.obj <- meanfunc(x.2$Lt, x.2$Ly, method = "PACE", 
#                            kernel = kernel, bw = mu.yao.obj$optns$userBwMu)
#     cov.lin.obj <- tryCatch({
#       covfunc(x.2$Lt, x.2$Ly, mu = mu.lin.obj, method = "SP")
#     },  error = function(e) { 
#       print("Lin cov error")
#       print(e)
#       return(NA) 
#     })
#     if (is.na(cov.lin.obj)) {
#       # return(NULL)
#       next
#     }
#     cov.lin <- predict(cov.lin.obj, work.grid)
#     
#     
#     ## Huber loss
#     # For computation times, we specified bw_mu.
#     mu.huber.obj <- tryCatch({
#       # meanfunc(x.2$Lt, x.2$Ly, method = "HUBER", kernel = kernel, bw = mu.yao.obj$optns$userBwMu)
#       meanfunc.rob(x.2$Lt, x.2$Ly, method = "Huber", kernel = kernel, bw = mu.yao.obj$optns$userBwMu)
#     }, error = function(e) { 
#       print("Huber mean error")
#       print(e) 
#       return(NA) 
#     })
#     if (is.na(mu.huber.obj)) {
#       # return(NULL)
#       next
#     }
#     # bandwidth are selected from 5-fold CV (almost 3 minutes)
#     cov.huber.obj <-  tryCatch({
#       # covfunc(x.2$Lt, x.2$Ly, mu = mu.huber.obj, kernel = kernel, method = "HUBER")   # CV
#       covfunc.rob(x.2$Lt, x.2$Ly, mu = mu.huber.obj, kernel = kernel, method = "Huber")
#       # covfunc(x.2$Lt, x.2$Ly, mu = mu.huber.obj, kernel = kernel, method = "HUBER",
#       #         bw = cov.lin.obj$sig2x$obj$bw)   # CV X
#     }, error = function(e) { 
#       print("Huber cov error")
#       print(e) 
#       return(NA) 
#     })
#     if (is.na(cov.huber.obj)) {
#       # return(NULL)
#       next
#     }
#     cov.huber <- predict(cov.huber.obj, work.grid)
#     
#     end_time <- Sys.time()
#     print(end_time - start_time)
#     
#     
#     # if some covariances is a not finite value
#     if (!is.finite(sum(cov.yao)) | !is.finite(sum(cov.lin)) | !is.finite(sum(cov.huber))) {
#       # return(NULL)
#       next
#     }
#     # if all covariances are 0
#     if ((sum(cov.yao) == 0) | (sum(cov.lin) == 0) | (sum(cov.huber) == 0)) {
#       # return(NULL) 
#       next
#     }
#     
#     
#     # output list
#     out <- list(work.grid = work.grid,
#                 mu.obj = list(yao = mu.yao.obj,
#                               lin = mu.lin.obj,
#                               huber = mu.huber.obj),
#                 cov.obj = list(yao = cov.yao.obj,
#                                lin = cov.lin.obj,
#                                huber = cov.huber.obj),
#                 cov = list(true = cov.true,
#                            yao = cov.yao,
#                            lin = cov.lin,
#                            huber = cov.huber))
#     
#     # return(out)
#     cov.est[[sim]] <- out
#   }
#   
#   return(cov.est)
# }

