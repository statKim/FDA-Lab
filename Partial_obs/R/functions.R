### Load pacakges and source codes
require(fdapace)
require(tidyverse)
require(tidyverse)
source("R/utills.R")
load_sources()

########################################################
### Robust mean estimation for functional snippets
########################################################

### Local polynomial kernel smoothing with huber loss (mean estimator)
# Lt : a list of vectors or a vector containing time points for all curves
# Ly : a list of vectors or a vector containing observations for all curves
# newt : a vector containing time points to estimate
# method : "Huber"
# bw : bandwidth
# kernel : a kernel function for kernel smoothing ("epanechnikov", "gauss" are supported.)
# deg : a numeric scalar of polynomial degrees for the kernel smoother
# k2 : If method == "Huber", it uses for \rho function in Huber loss.
# cv : If cv == TRUE, it performs K-fold cross-validation for optimal bandwidth(bw).
# cv_optns : a option for K-fold cross-validation if cv == TRUE.

# loss : a loss function for kernel smoothing("L2" is squared loss, "Huber" is huber loss.)
#   For loss = "Huber", it uses `rlm()` in `MASS` and fits the robust regression with Huber loss. 
#   So additional parameters of `rlm()` can be applied. (k2, maxit, ...)
meanfunc.rob <- function(Lt, 
                         Ly, 
                         newt = NULL, 
                         method = c("L2","Huber","WRM","Bisquare"), 
                         bw = NULL,
                         kernel = "epanechnikov",
                         deg = 1, 
                         k2 = 1.345,
                         ncores = 1,
                         cv_delta_loss = "L1",
                         cv_bw_loss = "HUBER",
                         cv_K = 5,
                         # cv = FALSE, 
                         # cv_optns = list(K = 5,
                         #                 ncores = 1,
                         #                 Loss = "Huber"), 
                         ...) {
  method <- toupper(method)
  if (!(method %in% c("L2","HUBER","WRM","BISQUARE"))) {
    stop(paste0(method, " is not provided. Check method parameter."))
  }
  
  R <- NULL
  
  others <- list(...)
  domain <- c(0,0)
  domain[1] <- min(unlist(Lt))
  domain[2] <- max(unlist(Lt))
  domain <- c(domain[1]-0.01*(domain[2]-domain[1]),
              domain[2]+0.01*(domain[2]-domain[1]))
  
  if (!is.list(Lt)) {
    stop("Lt and Ly must be list type.")
  }
  
  # kernel <- tolower(get.optional.param('kernel',others,'epanechnikov'))
  # bw <- get.optional.param('bw',others,NULL)
  
  cv <- FALSE
  # 5-fold CV for delta in Huber function
  if ((method %in% c("HUBER","BISQUARE")) && (is.null(k2) | ncores > 1)) {
    print(paste0(cv_K, "-fold CV is performed for delta in Huber function."))
    delta_cv_obj <- delta.local_kern_smooth(Lt = Lt, 
                                            Ly = Ly, 
                                            method = method,
                                            kernel = kernel, 
                                            deg = deg,
                                            # k2 = k2,
                                            bw = bw,
                                            cv_loss = cv_delta_loss, 
                                            K = cv_K, 
                                            ncores = ncores,
                                            ...)
    k2 <- delta_cv_obj$selected_delta
    cv <- TRUE
  }
  
  # 5-fold CV for bandwidth
  if (is.null(bw) | ncores > 1) {
    print(paste0(cv_K, "-fold CV is performed for bandwidth."))
    bw_cv_obj <- bw.local_kern_smooth(Lt = Lt, 
                                      Ly = Ly, 
                                      method = method,
                                      kernel = kernel, 
                                      deg = deg,
                                      # k2 = k2,
                                      cv_loss = cv_bw_loss, 
                                      K = cv_K, 
                                      ncores = ncores,
                                      ...)
    bw <- bw_cv_obj$selected_bw
    cv <- TRUE
  }
  
  n <- length(Lt)
  t <- unlist(Lt)
  y <- unlist(Ly)
  
  ord <- sort(t, index.return = T)$ix
  t <- t[ord]
  y <- y[ord]
  
  R <- list(bw = bw,
            t = t,
            y = y,
            n = n,
            domain = domain,
            method = method,
            kernel = kernel,
            deg = deg,
            k2 = k2,
            yend = c(NULL,NULL))
  class(R) <- 'meanfunc.rob'
  
  L0 <- domain[2]-domain[1]
  yend <- predict(R, c(domain[1]+L0/100, domain[2]-L0/100))
  R$yend <- yend
  
  if (!is.null(newt)) {
    R$fitted <- predict(R, newt)
  }
  
  if (cv == TRUE) {
    # R$cv_optns <- cv_optns
    R$cv_optns <- list(K = cv_K,
                       ncores = ncores,
                       delta_loss = cv_delta_loss,
                       bw_loss = cv_bw_loss)
  }
  
  return(R)
}


### Predict mean at new time points
predict.meanfunc.rob <- function(meanfunc.obj, newt) {
  pred <- function(newt) {   # newt must be a vector
    idxl <- newt < meanfunc.obj$domain[1]
    idxu <- newt > meanfunc.obj$domain[2]
    idx <- (!idxl) & (!idxu)
    
    newt0 <- newt[idx]
    ord <- sort(newt0, index.return = T)$ix
    
    tmp <- rep(Inf, length(newt0))
    
    # print(meanfunc.obj$method)
    
    tmp[ord] <- local_kern_smooth(Lt = meanfunc.obj$t, 
                                  Ly = meanfunc.obj$y, 
                                  newt = newt0[ord],
                                  method = meanfunc.obj$method,
                                  bw = meanfunc.obj$bw, 
                                  kernel = meanfunc.obj$kernel, 
                                  # loss = "Huber", 
                                  k2 = meanfunc.obj$k2)
    
    yhat <- rep(0, length(newt))
    yhat[idx] <- tmp
    yhat[idxl] <- meanfunc.obj$yend[1]
    yhat[idxu] <- meanfunc.obj$yend[2]
    
    return(yhat)
  }
  
  if (is.list(newt)) {   # Need to test this code for list type.
    mi <- lapply(newt, length)
    newt <- unlist(newt)
    fitted <- pred(newt)
    
    cm <- c(0, cumsum(mi))
    
    R <- sapply(1:length(mi), function(i) {
      res <- list()
      res[[1]] <- fitted[(cm[i]+1):cm[i+1]]
      res
    })
    
    return(R)
  } else if (is.vector(newt)) {
    # obtain unique time point => merge to newt
    df_newt <- data.frame(t = newt)
    newt_unique <- unique(newt)
    pred_unique <- pred(newt_unique)
    df_newt_unique <- data.frame(t = newt_unique,
                                 pred = pred_unique)
    pred_newt <- df_newt %>% 
      left_join(df_newt_unique, by = "t") %>% 
      dplyr::select(pred) %>% 
      unlist()
    
    return(as.numeric(pred_newt))
    # return(pred(newt))
  }
  else stop('newt must be a vector or a list of vectors of real numbers')
}


#' #' estimate the window width of snippets
#' #' @export
#' estimate.delta <- function(Lt) {
#'   if(is.list(Lt)) {
#'     tmp <- lapply(Lt, function(v) { max(v) - min(v) })
#'     return( max(unlist(tmp)) )
#'   } else if (is.vector(Lt)) {
#'     return( max(Lt)-min(Lt) )
#'   } else {
#'     stop('unsupported type of Lt')
#'   }
#' }



#########################################################
### Robust variance estimation for functional snippets
#########################################################
varfunc.rob <- function(Lt,
                        Ly,
                        newt = NULL,
                        sig2 = NULL,
                        method = c("L2","Huber","WRM","Bisquare"),
                        mu = NULL,
                        # weig=NULL,
                        ...) {
  
  method <- toupper(method)
  if (!(method %in% c("L2","HUBER","WRM","BISQUARE"))) {
    stop(paste0(method, " is not provided. Check method parameter."))
  }
  
  n <- length(Lt)
  
  if (!(is.list(Lt) && is.list(Ly))) {
    stop("Lt and Ly must be list type.")
  }
  
  if (is.null(sig2)) {
    sig2 <- sigma2.rob(Lt, Ly)
  }
  
  if (is.null(mu)) {
    mu <- meanfunc.rob(Lt, Ly, method = method)
  }
  
  # calculate sum of squares
  gr <- sort(unique(unlist(Lt)))
  mu_hat <- predict(mu, gr)
  ss <- lapply(1:n, function(i) {
    ind <- match(Lt[[i]], gr)
    if (length(ind) == length(Lt[[i]])) {
      return( (Ly[[i]] - mu_hat[ind])^2 )
    } else {
      mui <- predict(mu, Lt[[i]])
      return( (Ly[[i]] - mui)^2 )
    }
  })
  
  R <- list(obj = meanfunc.rob(Lt, ss, method = method, ...),
            sig2 = sig2)
  class(R) <- 'varfunc.rob'
  
  if (is.null(newt)) {
    newt <- Lt
  }
  
  # R$fitted <- predict(R, newt)
  
  return(R)
}

### Predict variance at new time points
predict.varfunc.rob <- function(R, newt) {
  
  if (is.list(newt)) {
    newt <- unlist(newt)
  }
  res <- predict(R$obj, newt)
  res <- res - R$sig2
  res[res < 0] <- 0
  
  return(res)
}



#########################################################
### Robust covariance estimation for functional snippets
#########################################################
covfunc.rob <- function(Lt, 
                        Ly, 
                        newt = NULL, 
                        mu = NULL, 
                        weig = NULL, 
                        method = c("Huber","WRM","Bisquare"), ...) {
  
  method <- toupper(method)
  if (!(method %in% c("HUBER","WRM","BISQUARE"))) {
    stop(paste0(method, " is not provided. Check method parameter."))
  }
  
  return(covfunc.rob.huber(Lt, Ly, mu = mu, newt = newt, method = method, ...))
  # if (method == 'Huber') {
  #   return(cov.huber(Lt, Ly, mu = mu, newt = newt, ...))
  # } else {
  #   stop("method is not supported.")
  # }
}

covfunc.rob.huber <- function(Lt, 
                              Ly, 
                              newt = NULL,
                              method = "HUBER",
                              domain = NULL,
                              weig = NULL,
                              corf = NULL, # correlation function(theta,x,y)
                              mu = NULL,
                              sig2e = NULL,
                              sig2x = NULL,
                              pfunc = NULL,
                              theta0 = NULL,
                              lb = NULL,
                              ub = NULL,
                              D = NULL,
                              kernel = "epanechnikov", ...) {
  
  if (is.null(corf)) {
    corf <- function(x, y, theta) {
      matern(x, y, nu = theta)
    }
    D <- 1
  } else {
    if (is.null(theta0) && is.null(D)) {
      stop('The dimension D must be specified')
    }
  }
  
  if (is.null(mu)) {
    mu <- meanfunc.rob(Lt, Ly, kernel = kernel, method = method)   # Huber option
  }
  mu.hat <- predict(mu, unlist(Lt))   #  TOO SLOW => Need to improve speed
  
  print("Finish mean estimation!")
  
  if (is.null(sig2e)) {
    sig2e <- sigma2.rob(Lt, Ly)
  }
  
  if(is.null(sig2x)) {
    # sig2x <- varfunc(Lt,Ly,mu=mu,sig2=sig2e)
    sig2x <- varfunc.rob(Lt, Ly, mu = mu, sig2 = sig2e, kernel = kernel, 
                         method = method, ...)   # Huber option
  }
  var.hat <- predict(sig2x, unlist(Lt))   #  TOO SLOW => Need to improve speed
  
  print("Finish variance estimation!")
  
  if (is.null(domain)) {
    t.vec <- unlist(Lt)
    domain <- c(min(t.vec), max(t.vec))
  }
  
  th.est <- estimate.theta(Lt,
                           Ly,
                           D = D,
                           var.hat = var.hat,
                           mu.hat = mu.hat,
                           method = 'LS',
                           rho = corf,
                           weig = weig,
                           pfunc = pfunc,
                           theta.lb = lb,
                           theta.ub = ub,
                           theta0 = theta0,
                           domain = domain)$LS
  
  rslt <- list(sig2e = sig2e,
               theta = th.est,
               mu.hat = mu.hat,
               domain = domain,
               mu = mu,
               sig2x = sig2x,
               rho = function(x, y) { corf(x, y, th.est) },
               method = method)
  class(rslt) <- 'covfunc.rob'
  
  if (!is.null(newt)) {
    rslt$fitted <- predict(rslt, newt)
  }
  print("Finish covariance estimation!")
  
  return(rslt)
}


### Predict covariance at new time points
predict.covfunc.rob <- function(covobj, newt) {
  pred <- function(newt) {  # newt is a vector
    stopifnot(is.vector(newt))
    
    corr.est <- covobj$rho(newt, newt)
    var.hat <- predict(covobj$sig2x, newt)
    sig.hat <- sqrt(var.hat)
    cov.fitted <- corr.est * (sig.hat %*% t(sig.hat))
    return(cov.fitted)
  }
  
  if (is.list(newt)) {
    return(lapply(newt, pred))
  } else if (is.vector(newt)) {
    return(pred(newt))
  } else {
    stop('newt must be a vector or a list of vectors of real numbers')
  }
}




#########################################################
### Robust noise variance estimation for functional snippets
#########################################################
sigma2.rob <- function(t, y, h = NULL) {
  if (is.list(y)) {   # irregular data
    n <- length(t)
    
    t.min <- min(unlist(t))
    t.max <- max(unlist(t))
    if (is.null(h)) {
      h <- select.sig2.bw(t,y)
    } else if (2*h >= t.max-t.min) {
      stop('h is too large')
    }
    
    A0A1 <- sapply(1:n, function(i){
      tobs <- t[[i]]
      y <- y[[i]]
      m <- length(tobs)
      A0 <- 0
      A1 <- 0
      B <- 0
      
      for(j in 1:m) {
        for(k in 1:m) {
          if( (k!=j) && abs(tobs[j]-tobs[k])<h ) {
            A0 <- A0 + (y[j]^2 + y[k]^2)/2
            A1 <- A1 + y[j]*y[k]
            B <- B + 1
          }
        }
      }
      return(c(A0/(m*(m-1)), 
               A1/(m*(m-1)),
               B/(m*(m-1))))
    })
    
    # mean((A0A1[1, ] - A0A1[2, ])/A0A1[3, ])
    # mean( mean(A0A1[1, ] - A0A1[2, ]) / mean(A0A1[3, ]) )
    sig2 <- median(A0A1[1, ] - A0A1[2, ]) / median(A0A1[3, ])
    
  } else {
    stop('unsupported data type')
  } 
  
  return(sig2)
}
