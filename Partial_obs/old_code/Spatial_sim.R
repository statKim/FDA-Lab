sim_corr <- function(n = 403,
                     type = c("partial","snippet","dense"),
                     out.prop = 0.2,
                     out.type = 1,
                     dist = "normal",
                     noise = 0,
                     dist.mat = dist.mat) {
  
  gr <- seq(0, 1, length.out = 51)   # equispaced points
  t <- gr
  n.loc <- nrow(dist.mat)   # number of locations
  
  # PC functions
  # phi <- cbind(
  #   sqrt(2)*cos(2*pi*t),
  #   sqrt(2)*sin(2*pi*t),
  #   sqrt(2)*cos(4*pi*t)
  # )
  phi <- get_corr_eigen(gr)
  
  r.par <- 200
  cor.sim <- matrix(fields::Matern(as.vector(dist.mat),
                                   range = r.par,
                                   smoothness = 1),
                    nrow = n.loc)
  
  # FPC scores
  # independent in traditional FPCA. Spatially correlated in spatial FPCA
  if (dist == 'normal') {
    xi <- rbind(
      # mvtnorm::rmvnorm(n = 1, mean = rep(0, n.loc), sigma = 9*cor.sim),
      mvtnorm::rmvnorm(n = 1, mean = rep(0, n.loc), sigma = 4*cor.sim),
      mvtnorm::rmvnorm(n = 1, mean = rep(0, n.loc), sigma = 1*cor.sim)
    )
  } else if (dist == 'tdist') {
    out.prop <- 0   # for heavy-tailed distrubution, we do not set outlier index
    xi <- rbind(
      # LaplacesDemon::rmvt(n = 1,
      #                     mu = rep(0, n.loc),
      #                     S = 9*cor.sim,
      #                     df = 3),
      LaplacesDemon::rmvt(n = 1,
                          mu = rep(0, n.loc),
                          S = 4*cor.sim,
                          df = 3),
      LaplacesDemon::rmvt(n = 1,
                          mu = rep(0, n.loc),
                          S = 1*cor.sim,
                          df = 3)
    )
  }
  
  # 403 observations over p grid points
  y <- t(xi) %*% t(phi)
  
  # random noise
  if (noise > 0) {
    y <- y + matrix(stats::rnorm(n*m, 0, sqrt(noise)), n, m)
  }
  
  x <- list()
  x$Ly <- lapply(1:n, function(i) { y[i, ] })
  x$Lt <- lapply(1:n, function(i) { gr })
  x$out.ind <- rep(0, n)   # indicator of outlier
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
              out.ind = x$out.ind,
              x.full = x.full)
  } else if (type == "snippet") {   # generate functional snippets
    # Lin & Wang(2020) setting
    len.frag = c(0.1, 0.3)   # length of domain for each curves
    a_l <- len.frag[1]
    b_l <- len.frag[2]
    
    Ly <- list()
    Lt <- list()
    for (n_i in 1:n) {
      l_i <- stats::runif(1, a_l, b_l)
      M_i <- stats::runif(1, a_l/2, 1-a_l/2)
      
      A_i <- max(0, M_i-l_i/2)
      B_i <- min(1, M_i+l_i/2)
      
      t <- gr[which(gr >= A_i & gr <= B_i)]   # observed grids
      m <- length(t)   # legnth of observed grids
      idx <- which(gr %in% t)
      
      Ly[[n_i]] <- x$Ly[[n_i]][idx]
      Lt[[n_i]] <- t
    }
    x <- list(Ly = Ly,
              Lt = Lt,
              out.ind = x$out.ind,
              x.full = x.full)
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
    x$out.ind[(n-n.outlier+1):n] <- 1   # outlier indicator
    # x$Lt[(n-n.outlier+1):n] <- x.outlier$Lt
  } else {
    stop(paste(out.type, "is not correct value of argument out.type! Just integer value between 1~3."))
  }
  
  return(x)
}




### Generate outlying curves
make_outlier <- function(x, out.type = 1) {
  n <- length(x$Lt)   # number of outlying curves
  d <- 0.3
  sigma.exp <- 1
  for (k in 1:n) {
    t <- x$Lt[[k]]
    m <- length(t)   # length of time points
    tmp.mat <- matrix(NA, m, m)
    for (j in 1:m){
      tmp.mat[j, ] <- abs(t - t[j])
    }
    Sigma <- exp(-tmp.mat/d) * sigma.exp^2

    mu <- rep(0, m)
    I <- matrix(0, m, m)
    diag(I) <- rep(1, m)
    Sig_norm <- matrix(0, m, m)
    diag(Sig_norm) <- rep(100, m)

    if (out.type == 1) {
      err.out <- LaplacesDemon::rmvt(1, mu, I, df = 3) * LaplacesDemon::rmvn(1, rep(2, m), Sig_norm)   # t with df=3
    } else if (out.type == 2) {
      err.out <- LaplacesDemon::rmvc(1, mu, I)   # cauchy
    } else if (out.type == 3) {
      err.out <- LaplacesDemon::rmvc(1, mu, Sigma)   # cauchy
    }

    # x_i <- rmvn(1, mu, Sigma) * 2 + err.out
    x_i <- err.out
    x$Ly[[k]] <- as.numeric(x_i)
  }

  return(x)
}



#' Get True eigenfunction for simulation of the spatially correlated case.
#'
#' True eigenfunctions are obtained for given grids from the spatially correlated case.
#'
#' @param grid a vector containing the observed timepoints.
#'
#' @return a p x 4 matrix containing the eigenfunctions where p is the length of grids and 4 is the number of true components.
#'
#' @examples
#' gr <- seq(0, 1, length.out = 51)
#' eig.true <- get_corr_eigen(gr)
#' matplot(eig.true, type = "l")
#'
#' @export
get_corr_eigen <- function(grid) {
  t <- grid
  eig_ftn <- cbind(
    1,
    sin(2*pi*t)
    # sqrt(2)*cos(2*pi*t),
    # sqrt(2)*sin(2*pi*t),
    # sqrt(2)*cos(4*pi*t)
  )

  for (j in 1:ncol(eig_ftn)){
    xx <- eig_ftn[, j]
    eig_ftn[, j] <- xx/sqrt(fdapace::trapzRcpp(t, xx^2))
  }
  
  # eig_ftn <- get_delaigle_eigen(grid)[, -1]
  # eig_ftn <- get_delaigle_eigen(grid, model = 1)
  
  return(eig_ftn)
}
