load_sources <- function() {
  library(Rcpp)
  path <- "../../mcfda/src/"
  flist <- list.files(path)
  flist <- c("cov.cpp","mean.cpp","Rhessian.cpp" )
  for (fname in flist) {
    print(fname)
    sourceCpp(paste0(path, fname))
  }
  
  path <- "../../mcfda/R/"
  flist <- list.files(path)
  for (fname in flist) {
    print(fname)
    source(paste0(path, fname))
  }
}

##########################################
### Functions for simulations
##########################################
cov_est_sim <- function(data.list, num.sim = NULL, seed = 1000, kernel = "epanechnikov") {
  if (is.null(num.sim)) {
    num.sim <- length(data.list)   # number of simulations
  }
  
  # list of manual functions and packages
  ftns <- fun2char()
  packages <- c("fdapace","mcfda","synfd")
  
  model.cov <- 2   # covariance function setting of the paper (1, 2)
  
  registerDoRNG(seed)
  # cov.est <- foreach(sim = 1:num.sim, .packages = packages, .export = ftns) %do% {
  cov.est <- list()
  for (sim in 1:num.sim) {
    print(paste0(sim, "th simulation:"))
    start_time <- Sys.time()
    
    # Get simulation data
    x <- data.list[[sim]]$x
    gr <- data.list[[sim]]$gr
    
    ### Covariance estimation
    work.grid <- seq(min(gr), max(gr), length.out = 51)
    cov.true <- get_cov_fragm(gr, model = model.cov)   # true covariance
    cov.true <- ConvertSupport(fromGrid = gr, toGrid = work.grid, Cov = cov.true)   # transform to the observed grid
    
    ## 1. Yao, Müller, and Wang (2005)
    ## 2. Liu and Müller (2009) - fitted.FPCA()
    x.2 <- list(Ly = x$y,
                Lt = x$t)
    kern <- ifelse(kernel == "epanechnikov", "epan", kernel)
    optns <- list(methodXi = 'CE', dataType = 'Sparse', kernel = kern, verbose = FALSE)
    mu.yao.obj <- GetMeanCurve(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)   # get bandwidth of mean estimation
    cov.yao.obj <- GetCovSurface(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)
    cov.yao <- cov.yao.obj$cov
    if (length(work.grid) != 51) {
      cov.yao <- ConvertSupport(fromGrid = cov.yao.obj$workGrid, toGrid = work.grid,
                                Cov = cov.yao.obj$cov)   # transform to the observed grid
    }
    
    
    ## 7. Lin & Wang (2020)
    # estimate mean by local polynomial method
    mu.lin.obj <- meanfunc(x.2$Lt, x.2$Ly, method = "PACE", 
                           kernel = kernel, bw = mu.yao.obj$optns$userBwMu)
    cov.lin.obj <- tryCatch({
      covfunc(x.2$Lt, x.2$Ly, mu = mu.lin.obj, method = "SP")
    },  error = function(e) { 
      print("Lin cov error")
      print(e)
      return(NA) 
    })
    if (is.na(cov.lin.obj)) {
      # return(NULL)
      next
    }
    cov.lin <- predict(cov.lin.obj, work.grid)
    
    
    ## Huber loss
    # For computation times, we specified bw_mu.
    mu.huber.obj <- tryCatch({
      # meanfunc(x.2$Lt, x.2$Ly, method = "HUBER", kernel = kernel, bw = mu.yao.obj$optns$userBwMu)
      meanfunc.rob(x.2$Lt, x.2$Ly, method = "Huber", kernel = kernel, bw = mu.yao.obj$optns$userBwMu)
    }, error = function(e) { 
      print("Huber mean error")
      print(e) 
      return(NA) 
    })
    if (is.na(mu.huber.obj)) {
      # return(NULL)
      next
    }
    # bandwidth are selected from 5-fold CV (almost 3 minutes)
    cov.huber.obj <-  tryCatch({
      # covfunc(x.2$Lt, x.2$Ly, mu = mu.huber.obj, kernel = kernel, method = "HUBER")   # CV
      covfunc.rob(x.2$Lt, x.2$Ly, mu = mu.huber.obj, kernel = kernel, method = "Huber")
      # covfunc(x.2$Lt, x.2$Ly, mu = mu.huber.obj, kernel = kernel, method = "HUBER",
      #         bw = cov.lin.obj$sig2x$obj$bw)   # CV X
    }, error = function(e) { 
      print("Huber cov error")
      print(e) 
      return(NA) 
    })
    if (is.na(cov.huber.obj)) {
      # return(NULL)
      next
    }
    cov.huber <- predict(cov.huber.obj, work.grid)
    
    end_time <- Sys.time()
    print(end_time - start_time)
    
    
    # if some covariances is a not finite value
    if (!is.finite(sum(cov.yao)) | !is.finite(sum(cov.lin)) | !is.finite(sum(cov.huber))) {
      # return(NULL)
      next
    }
    # if all covariances are 0
    if ((sum(cov.yao) == 0) | (sum(cov.lin) == 0) | (sum(cov.huber) == 0)) {
      # return(NULL) 
      next
    }
    
    
    # output list
    out <- list(work.grid = work.grid,
                mu.obj = list(yao = mu.yao.obj,
                              lin = mu.lin.obj,
                              huber = mu.huber.obj),
                cov.obj = list(yao = cov.yao.obj,
                               lin = cov.lin.obj,
                               huber = cov.huber.obj),
                cov = list(true = cov.true,
                           yao = cov.yao,
                           lin = cov.lin,
                           huber = cov.huber))
    
    # return(out)
    cov.est[[sim]] <- out
  }
  
  return(cov.est)
}







##########################################
### Simulation summary functions
##########################################

### Summary the ISE for simulations
summary_ise <- function(data.list, cov.est, method = "var") {
  
  if (method %in% c("var","cov")) {   ## validation for Lin & Wang(2020)
    if (method == "var") {   # variance
      ise.cov <- sapply(cov.est, function(x) {
        c(get_ise(diag(x$cov$true), diag(x$cov$yao), x$work.grid),
          get_ise(diag(x$cov$true), diag(x$cov$lin), x$work.grid),
          get_ise(diag(x$cov$true), diag(x$cov$huber), x$work.grid))
      })
    } else if (method == "cov") {   # covariance
      ise.cov <- sapply(cov.est, function(x) {
        c(get_ise(x$cov$true, x$cov$yao, x$work.grid),
          get_ise(x$cov$true, x$cov$lin, x$work.grid),
          get_ise(x$cov$true, x$cov$huber, x$work.grid))
      })
    } 
  } else if (method %in% c("intra","extra")) {   ## validation for Delaigle(2020)
    ise.cov <- mapply(function(x, y) {
      cov.list <- x$cov
      ind <- get_design_index(y$x$t)
      
      if (method == "intra") {   # intrapolated part
        cov.true <- cov_inter(cov.list$true, ind)
        cov.yao <- cov_inter(cov.list$yao, ind)
        cov.lin <- cov_inter(cov.list$lin, ind)
        cov.huber <- cov_inter(cov.list$huber, ind)
      } else if (method == "extra") {   # extrapolated part
        cov.true <- cov_extra(cov.list$true, ind)
        cov.yao <- cov_extra(cov.list$yao, ind)
        cov.lin <- cov_extra(cov.list$lin, ind)
        cov.huber <- cov_extra(cov.list$huber, ind)
      }
      
      c(get_ise(cov.true, cov.yao, x$work.grid),
        get_ise(cov.true, cov.lin, x$work.grid),
        get_ise(cov.true, cov.huber, x$work.grid))
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
    
    # eigen analysis
    eig.true <- get_eigen(cov = cov.true, grid = work.grid)
    eig.yao <- get_eigen(cov = cov.yao, grid = work.grid)
    eig.lin <- get_eigen(cov = cov.lin, grid = work.grid)
    eig.huber <- get_eigen(cov = cov.huber, grid = work.grid)
    
    # change eigen direction(sign) for first K eigenvectors
    K <- min(ncol(eig.true$phi),
             ncol(eig.yao$phi),
             ncol(eig.lin$phi),
             ncol(eig.huber$phi))
    eig.yao$phi[, 1:K] <- check_eigen_sign(eig.yao$phi[, 1:K], eig.true$phi[, 1:K])
    eig.lin$phi[, 1:K] <- check_eigen_sign(eig.lin$phi[, 1:K], eig.true$phi[, 1:K])
    eig.huber$phi[, 1:K] <- check_eigen_sign(eig.huber$phi[, 1:K], eig.true$phi[, 1:K])
    
    # output list
    out <- list(work.grid = work.grid,
                true = eig.true,
                yao = eig.yao,
                lin = eig.lin,
                huber = eig.huber)
    
    return(out)
  }
  return(pca.est)
}


