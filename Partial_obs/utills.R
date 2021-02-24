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


