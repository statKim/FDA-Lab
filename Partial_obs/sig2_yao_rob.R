

### Yao version
sigma2.rob.yao <- function(x, gr = NULL, cov = NULL) {
  m <- 51
  h <- 0.02
  # 이부분 cov_Mest로 수정하기
  cov_hat <- cov_Mest(x, 
                      smooth = FALSE, 
                      make.pos.semidef = F)
  
  if (is.null(gr)) {
    gr <- seq(0, 1, length.out = m)
  }
  
  # 1D smoothing
  var_y <- diag(cov_hat)
  var_y <- smooth.spline(gr, var_y)$y
  
  df <- data.frame(
    v = as.numeric(cov_hat),
    s = rep(gr, m),
    t = rep(gr, each = m)
  )
  
  # 2D smoothing
  var_x <- rep(NA, m)
  for (i in 1:m) {
    idx <- which((abs(df$s - gr[i]) <= h + .Machine$double.eps) & 
                   (abs(df$t - gr[i]) <= h + .Machine$double.eps) &
                   (df$s != df$t))
    var_x[i] <- mean(df$v[idx])
  }
  diag(cov_hat) <- var_x
  cov.sm.obj <- refund::fbps(cov_hat, list(x = gr,
                                           z = gr))
  rob.var <- cov.sm.obj$Yhat
  
  int_inf <- min(gr) + (max(gr) - min(gr))/4
  int_sup <- max(gr) - (max(gr) - min(gr))/4
  idx <- which(gr > int_inf & gr < int_sup)
  noise_var <- 2*trapzRcpp(gr[idx], var_y[idx] - diag(rob.var)[idx])
  
  if (noise_var < 0) {
    noise_var <- 1e-6
    warning("noise variance is estimated negative")
  }
  
  return(noise_var)
}



### Winsorization
sigma2.rob3 <- function(Lt, Ly, h = NULL) {
  if (is.list(Ly)) {   # irregular data
    n <- length(Lt)
    
    t.min <- min(unlist(Lt))
    t.max <- max(unlist(Lt))
    if (is.null(h)) {
      h <- select.sig2.rob.bw(Lt, Ly)
    } else if (2*h >= t.max-t.min) {
      stop('h is too large')
    }
    
    A0A1 <- sapply(1:n, function(i){
      t <- Lt[[i]]
      y <- Ly[[i]]
      m <- length(t)
      A0 <- 0
      A1 <- 0
      B <- 0
      
      for (j in 1:m) {
        for (k in 1:m) {
          if ( (k!=j) && (abs(t[j]-t[k]) <= h + .Machine$double.eps) ) {
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
    
    A <- data.frame(t(A0A1))
    colnames(A) <- c("A0","A1","B")
    A$diff <- A$A0 - A$A1
    A$y <- x.2$y
    # 
    # plot(density(A$diff))
    # cut_off <- max( A$diff[which(A$diff < min(boxplot(A$diff, plot = F)$out))] )
    # abline(v = cut_off, lty = 2, lwd = 2)
    # table(A$y,
    #       (A$diff > cut_off))
    
    
    # # trimmed mean where A_1 <= 0.75 quantile
    # cutoff_trimmed <- quantile(A0A1[1, ], 0.75)
    # ind_trimmed <- which(A0A1[1, ] < cutoff_trimmed)
    # # sig2 <- (mean(A0A1[1, ind_trimmed]) - mean(A0A1[2, ind_trimmed])) / mean(A0A1[3, ind_trimmed])
    # sig2 <- mean(A0A1[1, ind_trimmed] - A0A1[2, ind_trimmed]) / mean(A0A1[3, ind_trimmed])
    # # ------------------------------------------------------------------
    
    
    A_diff <- A0A1[1, ] - A0A1[2, ]
    
    # winsorization
    cut_off <- quantile(A_diff, 0.75)
    # cut_off <- max( A_diff[which(A_diff < min(boxplot(A_diff, plot = F)$out))] )
    # cut_off <- max( A_diff[which(A_diff < max(adjbox(A_diff, plot = F)$stats))] )
    ind_trimmed <- which(A_diff > cut_off)
    A_diff[ind_trimmed] <- cut_off
    # A_diff[ind_trimmed] <- median(A_diff)
    sig2 <- mean(A_diff) / mean(A0A1[3, ])
    
    
    # # weighted average
    # density_A <- density(A_diff)
    # w <- sapply(A_diff, function(a){
    #   density_A$y[which.min(abs(a - density_A$x))]
    # })
    # w <- w^2
    # numerator <- sum(w*A_diff) / sum(w)
    # denominator <- sum(w*A0A1[3, ]) / sum(w)
    # sig2 <- numerator / denominator
    # # sig2 <- numerator / mean(A0A1[3, ])
    
    # # density cut-off
    # density_A <- density(A_diff)
    # cut_off <- quantile(density_A$y, 0.99)
    # cut_range <- range( density_A$x[which(density_A$y > cut_off)] )
    # ind_trimmed <- which(A_diff >= cut_range[1] & A_diff <= cut_range[2])
    # sig2 <- mean(A_diff[ind_trimmed]) / mean(A0A1[3, ind_trimmed])
  } else {
    stop('unsupported data type')
  }
  
  return(sig2)
}


### Trimming
sigma2.rob2 <- function(Lt, Ly, h = NULL, alpha = 0.75) {
  if (is.list(Ly)) {   # irregular data
    n <- length(Lt)
    
    t.min <- min(unlist(Lt))
    t.max <- max(unlist(Lt))
    if (is.null(h)) {
      h <- select.sig2.rob.bw(Lt, Ly)
    } else if (2*h >= t.max-t.min) {
      stop('h is too large')
    }
    
    A0A1 <- sapply(1:n, function(i){
      t <- Lt[[i]]
      y <- Ly[[i]]
      m <- length(t)
      A0 <- 0
      A1 <- 0
      B <- 0
      
      for (j in 1:m) {
        for (k in 1:m) {
          if ( (k!=j) && (abs(t[j]-t[k]) <= h + .Machine$double.eps) ) {
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

    A <- data.frame(t(A0A1))
    colnames(A) <- c("A0","A1","B")
    A$diff <- A$A0 - A$A1
    A$y <- x.2$y
    # 
    # plot(density(A$diff))
    # cut_off <- max( A$diff[which(A$diff < min(boxplot(A$diff, plot = F)$out))] )
    # abline(v = cut_off, lty = 2, lwd = 2)
    # table(A$y,
    #       (A$diff > cut_off))
    

    # # trimmed mean where A_1 <= 0.75 quantile
    # cutoff_trimmed <- quantile(A0A1[1, ], quant)
    # ind_trimmed <- which(A0A1[1, ] < cutoff_trimmed)
    # # sig2 <- (mean(A0A1[1, ind_trimmed]) - mean(A0A1[2, ind_trimmed])) / mean(A0A1[3, ind_trimmed])
    # sig2 <- mean(A0A1[1, ind_trimmed] - A0A1[2, ind_trimmed]) / mean(A0A1[3, ind_trimmed])
    # # ------------------------------------------------------------------
    
    # diff(quantile(A_diff, c(0.25, 0.75)))*1.5
    
    
    A_diff <- A0A1[1, ] - A0A1[2, ]
    cut_off <- quantile(A_diff, alpha)
    # cut_off <- max( A_diff[which(A_diff < min(boxplot(A_diff, plot = F)$out))] )
    # cut_off <- max( A_diff[which(A_diff < max(adjbox(A_diff, plot = F)$stats))] )
    ind_trimmed <- which(A_diff <= cut_off)
    sig2 <- mean(A_diff[ind_trimmed]) / mean(A0A1[3, ind_trimmed])
  } else {
    stop('unsupported data type')
  }
  
  return(sig2)
}




### select optimal bandwidth
select.sig2.rob.bw <- function(Lt, Ly, c = 0.29) {
  n <- length(Lt)
  
  t.min <- min(unlist(Lt))
  t.max <- max(unlist(Lt))
  
  delta <- max(sapply(Lt, function(ts){ max(ts)-min(ts) }))
  m <- mean(sapply(Lt, length))
  M <- n * m^2
  #h <- delta * (M^(-1/5)) / 7
  
  # # calculate sum of squares - Not recommended (Already it was calculated once!)
  # if (is.null(ss)) {
  #   # # mu.hat <- predict(meanfunc(Lt,Ly),Ly)
  #   # mu.hat <- predict(meanfunc.rob(Lt, Ly), Lt)
  #   # ss <- lapply(1:length(Lt),function(i){
  #   #   rr <- Ly[[i]] - mu.hat[[i]]
  #   #   rr^2
  #   # })
  #   
  #   gr <- sort(unique(unlist(Lt)))
  #   if (is.null(mu)) {
  #     mu <- meanfunc.rob(Lt, Ly, method = "huber", kernel = "epanechnikov",
  #                        bw = (t.max-t.min)/5, delta = 1.345)
  #   }
  #   mu_hat <- predict(mu, gr)
  #   ss <- lapply(1:n, function(i) {
  #     ind <- match(Lt[[i]], gr)
  #     if (length(ind) == length(Lt[[i]])) {
  #       return( (Ly[[i]] - mu_hat[ind])^2 )
  #     } else {
  #       mui <- predict(mu, Lt[[i]])
  #       return( (Ly[[i]] - mui)^2 )
  #     }
  #   })
  # }
  # 
  # vn <- sqrt(mean(unlist(ss)))
  
  x <- list2matrix(list(Lt = Lt,
                        Ly = Ly))
  cov_est <- cov_Mest(x, smooth = T)
  # vn <- sqrt( mean(diag(cov_est)^2) )
  vn <- sqrt( fdapace::trapzRcpp(gr, diag(cov_est)^2) )
  h <- c * delta * vn * (M^(-1/5))
  
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
