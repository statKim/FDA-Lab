### Generate functional data with outliers
# model = 1~2 avaliable (In the paper, 1~4 models are shown)
# type
# out.prop
# out.type : 1~3 are available
sim_delaigle2 <- function(n = 100, 
                          model = 2, 
                          type = c("partial","snippet","snippet-long","sparse","dense"),
                          dist = "normal",
                          out.prop = 0, 
                          out.type = 1,
                          noise = 0) {
  
  gr <- seq(0, 1, length.out = 51)   # equispaced points
  
  # generate dense curves
  m <- length(gr)   # legnth of observed grids
  lambda <- 0.5^(0:3)
  phi <- get_delaigle_eigen(gr, model = model)
  
  if (dist == "normal") {
    xi <- mvtnorm::rmvnorm(n, sigma = diag(lambda))
  } else if (dist == "t") {
    xi <- mvtnorm::rmvt(n, sigma = diag(lambda), df = 3)
  }
  y <- xi %*% t(phi)
  # if (dist == "normal") {
  #   xi <- mvtnorm::rmvnorm(n, sigma = diag(1, length(lambda)))
  # } else if (dist == "t") {
  #   xi <- mvtnorm::rmvt(n, sigma = diag(1, length(lambda)), df = 3)
  # }
  # y <- xi %*% diag(sqrt(lambda)) %*% t(phi)
  cov_sim <- get_cov_fragm(gr, model = model)
  
  # random noise
  if (noise > 0) {
    y <- y + matrix(rnorm(n*m, 0, sqrt(noise)), n, m)
  }
  
  x <- list()
  x$Ly <- lapply(1:n, function(i) { y[i, ] })
  x$Lt <- lapply(1:n, function(i) { gr })
  x$y <- rep(0, n)   # indicator of outlier
  x$xi <- xi
  x.full <- t(sapply(x$Ly, cbind))   # matrix containing the fully observed data
  
  # Check type option
  if (type == "dense") {   # Nothing do
    x$x.full <- x.full
  } else if (type == "partial") {   # generate partially obseved data
    # generate observation periods (Kraus(2015) setting)
    # curve 1 will be missing on (.4,.7), other curves on random subsets
    x.obs <- rbind((gr <= .4) | (gr >= .7), 
                   simul.obs(n = n-1, grid = gr)) # TRUE if observed
    # remove missing periods 
    x.partial <- x.full
    x.partial[!x.obs] <- NA
    
    x <- list(Ly = apply(x.partial, 1, function(y){ y[!is.na(y)] }),
              Lt = apply(x.obs, 1, function(y){ gr[y] }),
              y = x$y,
              xi = x$xi,
              x.full = x.full)
  } else if (type %in% c("snippet","snippet-long")) {   # generate functional snippets
    # Delaigle(2020) setting
    if (type == "snippet") {
      len.frag = c(0.1, 0.3)   # length of domain for each curves
    } else {
      len.frag = c(0.4, 0.6)   # length of domain for each curves
    }
    
    a_l <- len.frag[1]
    b_l <- len.frag[2]
    
    Ly <- list()
    Lt <- list()
    for (n_i in 1:n) {
      l_i <- runif(1, a_l, b_l)
      M_i <- runif(1, a_l/2, 1-a_l/2)
      
      A_i <- max(0, M_i-l_i/2)
      B_i <- min(1, M_i+l_i/2)
      
      t <- gr[which(gr >= A_i & gr <= B_i)]   # observed grids
      m <- length(t)   # legnth of observed grids
      idx <- which(gr %in% t)
      
      Ly[[n_i]] <- x$Ly[[n_i]][idx]
      Lt[[n_i]] <- t
      
      # # Exact Lin & Wang(2020) setting
      # cov_sim <- matrix(NA, m, m)
      # for (i in 1:m) {
      #   for (j in 1:m) {
      #     # If upper triangular than compute, else substitute transposed value
      #     if (i <= j) {
      #       cov_sim[i, j] <- get_K(t[i], t[j], model = model)
      #     } else {
      #       cov_sim[i, j] <- cov_sim[j, i]
      #     }
      #   }
      # }
      # 
      # x$Ly[[n_i]] <- as.numeric(rmvnorm(1, rep(0, m), cov_sim))
      # x$Lt[[n_i]] <- t
    }
    x <- list(Ly = Ly,
              Lt = Lt,
              y = x$y,
              xi = x$xi,
              x.full = x.full)
  } else if (type == "sparse") {
    
  } else {
    stop(paste(type, "is not an appropriate argument of type"))
  }
  
  # no outliers
  if (out.prop == 0) {
    return(x)
  }
  
  # generate outlier curves
  n.outlier <- ceiling(n*out.prop)   # number of outliers
  if (out.type %in% 1:3) {
    x.outlier <- list(Ly = x$Ly[(n-n.outlier+1):n],
                      Lt = x$Lt[(n-n.outlier+1):n])
    x.outlier <- make_outlier(x.outlier, out.type = out.type)
    x$Ly[(n-n.outlier+1):n] <- x.outlier$Ly
    x$y[(n-n.outlier+1):n] <- 1   # outlier indicator
    # x$Lt[(n-n.outlier+1):n] <- x.outlier$Lt
  } else if (out.type == 4) {
    # x.outlier <- mvtnorm::rmvt(n = n.outlier, 
    #                            delta = rep(0, m), 
    #                            sigma = cov_sim, 
    #                            df = 1)
    # idx <- lapply(x$Lt[(n-n.outlier+1):n], function(t) { 
    #   which(t %in% gr) 
    # })
    # x$Ly[(n-n.outlier+1):n] <- lapply(1:n.outlier, function(i) { 
    #   x.outlier[i, idx[[i]]] 
    # })
    x.outlier <- rep(gr, each = n) + mvtnorm::rmvt(n = n, 
                                                   sigma = cov_sim, 
                                                   df = 1)
    idx <- lapply(x$Lt, function(t) { 
      which(t %in% gr) 
    })
    x$Ly <- lapply(1:n, function(i) { 
      x.outlier[i, idx[[i]]] 
    })
    x$x.full <- x.outlier
    x$y[(n-n.outlier+1):n] <- 1   # outlier indicator
  } else if (out.type == 5) {
    eps2 <- 0.3
    x.outlier <- lapply(x$Ly[(n-n.outlier+1):n], function(y) {
      m <- length(y)
      W_s <- rbinom(m, 1, eps2)
      z_s <- rnorm(m, 30, 0.1)
      Y <- W_s*z_s
      
      X_i <- y + Y
      return(X_i)
    })
    x$Ly[(n-n.outlier+1):n] <- x.outlier
    x$y[(n-n.outlier+1):n] <- 1   # outlier indicator
  } else if (out.type == 6) {
    # l <- 1/15
    l <- 1/10
    M <- 30
    for (i in (n-n.outlier+1):n) {
      Ly <- x$Ly[[i]]
      Lt <- x$Lt[[i]]
      
      D_i <- sample(c(1, -1), 1)
      
      out_ind <- integer(0)   # just make the variable
      while (length(out_ind) == 0) {
        T_i <- runif(1, 0, 1-l)
        out_ind <- which(Lt > T_i & Lt < T_i+l)
      }
      
      Ly[out_ind] <- Ly[out_ind] + D_i*M
      x$Ly[[i]] <- Ly
      x$y[i] <- 1   # outlier indicator
    }
  } else if (out.type == 7) {
    # l <- 1/15
    l <- 1/10
    M <- 30
    for (i in (n-n.outlier+1):n) {
      Ly <- x$Ly[[i]]
      m <- length(Ly)
      
      eps2 <- 0.3
      # V <- rbinom(1, 1, eps1)
      V <- 1
      
      W_s <- rbinom(m, 1, eps2)
      z_s <- rnorm(m, 30, 0.1)
      pm <- sample(c(-1, 1), 1)
      Y <- W_s*z_s*pm
      
      Ly <- Ly + V*Y
      x$Ly[[i]] <- Ly
      x$y[i] <- 1   # outlier indicator
    }
  } else {
    stop(paste(out.type, "is not correct value of argument out.type! Just integer value between 1~3."))
  }
  
  return(x)
}


# # set.seed(100)
# x <- sim_delaigle2(n = 100, 
#                    model = 2, 
#                    type = "snippet1",
#                    out.prop = 0, 
#                    out.type = 1,
#                    noise = 0)
# x <- list2matrix(x)
# matplot(t(x), type = "l")


