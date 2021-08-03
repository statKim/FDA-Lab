##############################################################
### Robust noise variance estimation for functional snippets
##############################################################

#' Noise variance estimator
#'
#' @param t a list of vectors containing time points for each curve
#' @param y a list of vectors containing observations for each curve
#' @param h a bandwidth
#'
#' @importFrom stats median quantile
#'
#' @export
#' @useDynLib robfpca
sigma2.rob.t <- function(t, y, h = NULL) {
  if (is.list(y)) {   # irregular data
    n <- length(t)
    
    t.min <- min(unlist(t))
    t.max <- max(unlist(t))
    if (is.null(h)) {
      h <- select.sig2.rob.bw(t, y)
    } else if (2*h >= t.max-t.min) {
      stop('h is too large')
    }
    
    gr <- sort(unique(unlist(t)))
    sig2 <- numeric(length(gr))
    
    for (t_ind in 1:length(gr)) {
      A0A1 <- sapply(1:n, function(i){
        tobs <- t[[i]]
        y <- y[[i]]
        m <- length(tobs)
        A0 <- 0
        A1 <- 0
        B <- 0
        
        for (j in 1:m) {
          for (k in 1:m) {
            if ( (k!=j) && (abs(tobs[j]-gr[t_ind]) < h) && (abs(tobs[k]-gr[t_ind]) < h) ) {
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
      
      # ------------------------- test version ---------------------------
      # mcfda
      # mean((A0A1[1, ] - A0A1[2, ])/A0A1[3, ])
      # mean( mean(A0A1[1, ] - A0A1[2, ]) / mean(A0A1[3, ]) )
      
      # wrong version of median
      # sig2 <- median(A0A1[1, ] - A0A1[2, ]) / median(A0A1[3, ])
      
      # median version
      # sig2 <- (median(A0A1[1, ]) - median(A0A1[2, ])) / median(A0A1[3, ])
      
      # 25%-trimmed mean (too large estimate)
      # sig2 <- (mean(A0A1[1, ], trim = 0.25) - mean(A0A1[2, ], trim = 0.25)) / mean(A0A1[3, ], trim = 0.25)
      
      # wrong 25%-trimmed mean
      # sig2 <- mean(A0A1[1, ] - A0A1[2, ], trim = 0.25) / mean(A0A1[3, ], trim = 0.25)
      
      # 25% trimmed estimator (interquartile mean; Q1 ~ Q3)
      # cutoff_trimmed <- quantile(A0A1[1, ], c(0.25, 0.75))
      # ind_trimmed <- which(A0A1[1, ] > cutoff_trimmed[1] & A0A1[1, ] < cutoff_trimmed[2])
      # # sig2 <- (mean(A0A1[1, ind_trimmed]) - mean(A0A1[2, ind_trimmed])) / mean(A0A1[3, ind_trimmed])
      # mean(A0A1[1, ind_trimmed] - A0A1[2, ind_trimmed]) / mean(A0A1[3, ind_trimmed])
      
      # trimmed mean where A_1 <= 0.75 quantile
      cutoff_trimmed <- quantile(A0A1[1, ], 0.75)
      ind_trimmed <- which(A0A1[1, ] < cutoff_trimmed)
      # sig2 <- (mean(A0A1[1, ind_trimmed]) - mean(A0A1[2, ind_trimmed])) / mean(A0A1[3, ind_trimmed])
      sig2[t_ind] <- mean(A0A1[1, ind_trimmed] - A0A1[2, ind_trimmed]) / mean(A0A1[3, ind_trimmed])
      # ------------------------------------------------------------------
    }
    
    
    
  } else {
    stop('unsupported data type')
  }
  
  return(sig2)
}


### select optimal bandwidth
select.sig2.rob.bw <- function(Lt, Ly, ss = NULL, mu = NULL) {
  n <- length(Lt)
  
  t.min <- min(unlist(Lt))
  t.max <- max(unlist(Lt))
  
  delta <- max(sapply(Lt, function(ts){ max(ts)-min(ts) }))
  m <- mean(sapply(Lt, length))
  M <- n * m^2
  #h <- delta * (M^(-1/5)) / 7
  
  # calculate sum of squares - Not recommended (Already it was calculated once!)
  if (is.null(ss)) {
    # # mu.hat <- predict(meanfunc(Lt,Ly),Ly)
    # mu.hat <- predict(meanfunc.rob(Lt, Ly), Lt)
    # ss <- lapply(1:length(Lt),function(i){
    #   rr <- Ly[[i]] - mu.hat[[i]]
    #   rr^2
    # })
    
    gr <- sort(unique(unlist(Lt)))
    if (is.null(mu)) {
      mu <- meanfunc.rob(Lt, Ly, method = "huber", kernel = "epanechnikov",
                         bw = (t.max-t.min)/5, delta = 1.345)
    }
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
  }
  
  vn <- sqrt(mean(unlist(ss)))
  
  h <- 0.29 * delta * vn * (M^(-1/5))
  
  max.it <- 1000
  it <- 0
  while (it < max.it) {  # h0 two small
    it <- it + 1
    # print(paste(it, ":"))
    
    cnt <- sapply(1:n, function(i) {
      tobs <- Lt[[i]]
      # y <- Ly[[i]]
      m <- length(tobs)
      v1 <- 0
      
      if (m < 2) {
        return(0)
      }
      for (j in 1:m) {
        for (k in 1:m) {
          if ( (k!=j) && (abs(tobs[j]-tobs[k]) < h) ) {
            v1 <- v1 + 1
          }
        }
      }
      
      return(v1)
    })
    
    cnt <- sum(cnt)
    if (cnt >= min(50, 0.1*n*m*(m-1))) {
      break
    } else {
      h <- h*1.01
    }
  }
  
  return(h)
}



sigma2.rob(x.2$Lt, x.2$Ly)

sigma2.rob.t(x.2$Lt, x.2$Ly)
