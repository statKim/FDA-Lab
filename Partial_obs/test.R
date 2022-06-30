library(tidyverse)
### Data Load
load("/Users/hyunsung/GoogleDrive/Lab/KHS/partiall_obs/real_data/PM10/pm_data_korea.RData")
head(new_data)
dim(new_data)

data <- new_data[, c(1, 5, 6)] %>% 
  drop_na() %>% 
  distinct()
dim(data)
data
# mat <- dist(data[, -1]) %>% 
#   as.matrix()
# dim(mat)
# mat

library(geosphere)
dist.mat <- distm(data[, -1])
dim(dist.mat)
dist.mat[1:5, 1:5]



############################################
### Simulation data generation functions
### from spatially correlated data
############################################

#' Generate partially observed functional data with spatially correlated case
#'
#' Partially observed functional data is generated with 51 regular grids from the spatially correlated model
#'
#' @param n a number of curves
#' @param type the type of generated data. "partial" means the option for partially observed data, "snippet" is short fragmented data, and "dense" means the fully observed curves.
#' @param out.prop a proportion of outlying curves of total n curves. Only used for dist = "normal".
#' @param out.type a outlier type, 1~3 are supported. Only used for dist = "normal".
#' @param dist a distribution which the data is generated. "normal"(Normal distribution) and "tdist"(t-distribution) are supported. If dist = "tdist", the option of \code{out.prop} and \code{out.type} are ignored.
#' @param noise a numeric value which is added random gaussian noises. Default is 0(No random noise).
#' @param dist.mat a 403 spatial locations and their distances matrix
#' @param d a parameter for missingness when \code{type} is "partial" (See Kraus(2015))
#' @param f a parameter for missingness when \code{type} is "partial" (See Kraus(2015))
#'
#' @return a list contatining as follows:
#' \item{Ly}{a list of n vectors containing the observed values for each individual.}
#' \item{Lt}{a list of n vectors containing the observation time points for each individual corresponding to \code{Ly}}
#' \item{out.ind}{a vector containing outlier index. 0 is non-outlier and 1 is the outlier.}
#' \item{x.full}{a n x 51 dense matrix with n observations per 51 timepoints before making partially observed.}
#'
#' @examples
#' set.seed(100)
#' n <- 100
#' x.list <- sim_corr(n = n,
#'                    type = "partial",
#'                    out.prop = 0.2,
#'                    out.type = 1,
#'                    dist = "normal",
#'                    dist.mat = dist.mat[1:n, 1:n])
#' x <- list2matrix(x.list)
#' matplot(t(x), type = "l")
#'
#' @export
sim_corr <- function(n = 403,
                     type = c("partial","snippet","dense"),
                     out.prop = 0.2,
                     out.type = 1,
                     dist = "normal",
                     noise = 0,
                     dist.mat = dist.mat,
                     d = 1.4,
                     f = 0.2,
                     r.par = 200) {
  
  gr <- seq(0, 1, length.out = 51)   # equispaced points
  t <- gr
  n.loc <- nrow(dist.mat)   # number of locations
  
  # PC functions
  phi <- cbind(
    sqrt(2)*cos(2*pi*t),
    sqrt(2)*cos(4*pi*t),
    sqrt(2)*cos(6*pi*t)
    # sqrt(2)*sin(2*pi*t),
    # sqrt(2)*cos(4*pi*t)
  )
  
  # r.par <- 200
  cor.sim <- matrix(fields::Matern(as.vector(dist.mat),
                                   range = r.par,
                                   smoothness = 1),
                    nrow = n.loc)
  
  # FPC scores
  # independent in traditional FPCA. Spatially correlated in spatial FPCA
  if (dist == 'normal') {
    xi <- rbind(
      mvtnorm::rmvnorm(n = 1, mean = rep(0, n.loc), sigma = 9*cor.sim),
      mvtnorm::rmvnorm(n = 1, mean = rep(0, n.loc), sigma = 4*cor.sim),
      mvtnorm::rmvnorm(n = 1, mean = rep(0, n.loc), sigma = 1*cor.sim)
    )
  } else if (dist == 'tdist') {
    out.prop <- 0   # for heavy-tailed distrubution, we do not set outlier index
    xi <- rbind(
      LaplacesDemon::rmvt(n = 1,
                          mu = rep(0, n.loc),
                          S = 9*cor.sim,
                          df = 3),
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
                   simul.obs(n = n-1, grid = gr,
                             d = d, f = f)) # TRUE if observed
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
    sqrt(2)*cos(2*pi*t),
    sqrt(2)*cos(4*pi*t),
    sqrt(2)*cos(6*pi*t)
    # sqrt(2)*sin(2*pi*t),
    # sqrt(2)*cos(4*pi*t)
  )
  
  for (j in 1:3){
    xx <- eig_ftn[, j]
    eig_ftn[, j] <- xx/sqrt(fdapace::trapzRcpp(t, xx^2))
  }
  
  return(eig_ftn)
}



# ###################################
# ### Functions from Kraus (2015)
# ###################################
# ### Functions for partially observed case
# norm.randcoef = function(n) stats::rnorm(n,0,1)
# unif.randcoef = function(n) stats::runif(n,-1,1)*sqrt(3)
# t5.randcoef = function(n) stats::rt(n,5)/sqrt(5/3)
# simul.obs <- function(n = 100, grid = seq(0, 1, len = 200), d = 1.4, f = .2) {
#   out = matrix(TRUE,n,length(grid))
#   for (i in 1:n) {
#     cen = d*sqrt(stats::runif(1))
#     e = f*stats::runif(1)
#     out[i,(cen-e<grid)&(grid<cen+e)] = FALSE # missing interval = (cen-u,cen+u)
#   }
#   out
# }
