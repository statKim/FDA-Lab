############################################
### Simulation data generation functions
### Boente et al. (2015) simulation
############################################

### Generate functional data with outliers
# model = 1 avaliable
# type
sim_boente <- function(n = 100, 
                       model = 1, 
                       outlier = TRUE,
                       eps1 = 0.2,
                       type = c("partial","snippet","sparse","dense")) {
  gr <- seq(0, 1, length.out = 51)   # equispaced points
  m <- length(gr)   # legnth of observed grids
  
  # generate dense curves
  if (model == 1) {
    obj <- boente_model_1(n, m, gr, outlier, eps1)
  } else if (model == 2) {
    obj <- boente_model_2(n, m, gr, outlier, eps1)
  } else if (model == 3) {
    obj <- boente_model_3(n, m, gr, outlier, eps1)
  } else {
    stop(paste("model =", model, "is not supported."))
  }

  Ly <- lapply(1:n, function(i) { obj$X[i, ] })
  Lt <- lapply(1:n, function(i) { gr })
  y <- obj$y
  phi <- obj$phi
  x.full <- obj$X   # matrix containing the fully observed data
  
  # Check type option
  if (type == "dense") {   # Nothing do
    x <- list(Ly = Ly,
              Lt = Lt,
              y = y,
              phi = phi,
              x.full = x.full)
  } else if (type == "partial") {   # generate partially obseved data
    # generate observation periods (Kraus(2015) setting)
    # curve 1 will be missing on (.4,.7), other curves on random subsets
    x.obs <- rbind((gr <= .4) | (gr >= .7), 
                   simul.obs(n = n-1, grid = gr)) # TRUE if observed
    # remove missing periods 
    x <- x.full
    x[!x.obs] <- NA
    
    x <- list(Ly = apply(x, 1, function(y){ y[!is.na(y)] }),
              Lt = apply(x.obs, 1, function(y){ gr[y] }),
              y = y,
              phi = phi,
              x.full = x.full)
  } else if (type == "snippet") {   # generate functional snippets
    # Lin & Wang(2020) setting
    len.frag = c(0.1, 0.3)   # length of domain for each curves
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
              y = y,
              phi = phi,
              x.full = x.full)
  } else if (type == "sparse") {
    
  } else {
    stop(paste(type, "is not an appropriate argument of type"))
  }
  
  return(x)
}



boente_model_1 <- function(n, m, gr, outlier, eps1 = 0.2) {
  X <- matrix(0, n, m)
  y <- rep(0, m)   # outlier indicator (if outlier, then 1, else 0)
  for (i in 1:n) {
    xi_1 <- rnorm(1, 0, 5/2)
    xi_2 <- rnorm(1, 0, 1/2)
    mu <- sapply(gr, function(t) {
      5 + 10*sin(4*pi*t)*exp(-2*t) + 5*sin(pi*t/3) + 2*cos(pi*t/2)
    })
    phi_1 <- sapply(gr, function(t) {
      sqrt(2)*cos(2*pi*t)
    })
    phi_2 <- sapply(gr, function(t) {
      sqrt(2)*sin(2*pi*t)
    })
    z <- rnorm(1, 0, 1)
    
    X_i <- 10 + mu + xi_1*phi_1 + xi_2*phi_2 + z
    
    # Contaminated trajectories
    if (outlier == TRUE) {
      # eps1 <- 0.2
      eps2 <- 0.3
      V <- rbinom(1, 1, eps1)
      y[i] <- V   # outlier indicator
      
      W_s <- rbinom(m, 1, eps2)
      z_s <- rnorm(m, 30, 0.1)
      Y <- W_s*z_s
      
      X_i <- X_i + V*Y
    }
    
    X[i, ] <- X_i
  }
  
  return(list(X = X,
              y = y,
              phi = cbind(phi_1, 
                          phi_2)))
}



boente_model_2 <- function(n, m, gr, outlier, eps1 = 0.2) {
  X <- matrix(0, n, m)
  y <- rep(0, m)   # outlier indicator (if outlier, then 1, else 0)
  for (i in 1:n) {
    xi_1 <- rnorm(1, 0, 5/2)
    xi_2 <- rnorm(1, 0, 1/2)
    mu <- sapply(gr, function(t) {
      5 + 10*sin(4*pi*t)*exp(-2*t) + 5*sin(pi*t/3) + 2*cos(pi*t/2)
    })
    phi_1 <- sapply(gr, function(t) {
      sqrt(2)*cos(2*pi*t)
    })
    phi_2 <- sapply(gr, function(t) {
      sqrt(2)*sin(2*pi*t)
    })
    z <- rnorm(1, 0, 1)
    
    X_i <- 150 - 2*mu + xi_1*phi_1 + xi_2*phi_2 + z
    
    # Contaminated trajectories
    if (outlier == TRUE) {
      # eps1 <- 0.2
      eps2 <- 0.9
      V <- rbinom(1, 1, eps1)
      y[i] <- V   # outlier indicator
      
      out_ind <- which(gr < 0.4)
      mu_c <- -5 -2*mu[out_ind]
      W_s <- rbinom(length(out_ind), 1, eps2)
      z_s <- mvtnorm::rmvnorm(1, mean = mu_c, sigma = diag(0.01, length(mu_c)))
      Y <- W_s * as.numeric(z_s)
      
      X_i[out_ind] <- X_i[out_ind] + V*Y
    }
    
    X[i, ] <- X_i
  }
  
  return(list(X = X,
              y = y,
              phi = cbind(phi_1, 
                          phi_2)))
}


boente_model_3 <- function(n, m, gr, outlier, eps1 = 0.2) {
  y <- rep(0, m)   # outlier indicator (if outlier, then 1, else 0)
  
  # covariance kernel
  Gamma_x <- 10 * matrix(apply(expand.grid(gr, gr), 1, min),
                         m, m)
  phi <- sapply(1:m, function(j) {
    sqrt(2) * sin((2*j-1)*(pi/2)*gr)
  })
  # lambda <- 10*(2/(d*(2*j-1)*pi))^2
  
  X <- mvtnorm::rmvnorm(n, sigma = Gamma_x)
  
  # Contaminated trajectories
  if (outlier == TRUE) {
    for (i in 1:n) {
      # eps1 <- 0.2
      l <- 1/15
      M <- 30
      V <- rbinom(1, 1, eps1)
      D_i <- sample(c(1, -1), 1)
      T_i <- runif(1, 0, 1-l)
      y[i] <- V   # outlier indicator
      
      out_ind <- which(gr > T_i & gr < T_i+l)
      X[i, out_ind] <- X[i, out_ind] + V*D_i*M
    }
  }
  
  return(list(X = X,
              y = y,
              phi = phi))
}

