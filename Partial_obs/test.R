############################################
### Simulation data generation functions
### Kraus (2015) simulation with outliers
### - Assume true components is 5
############################################


#' Get True eigenfunction for simulation in Kraus(2015).
#'
#' True eigenfunctions are obtained for given grids from the setting of Kraus(2015).
#'
#' @param grid a vector containing the observed timepoints.
#'
#' @return a p x k matrix containing the eigenfunctions where p is the length of grids and k is the number of true components.
#'
#' @examples
#' gr <- seq(0, 1, length.out = 51)
#' eig.true <- get_kraus_eigen(gr)
#' matplot(eig.true, type = "l")
#'
#' @references
#' \cite{Kraus, D. (2015). Components and completion of partially observed functional data. Journal of the Royal Statistical Society: Series B: Statistical Methodology, 777-801.}
#'
#' @export
#'
### Get True eigenfunction for simulation in Kraus(2015) setting
get_kraus_eigen <- function(grid) {
  t <- grid
  eig_ftn <- matrix(0, length(t), 3)
  j <- 1
  eig_ftn[, 1] <-cos(2*pi*j*t)
  eig_ftn[, 2] <- sin(2*pi*j*t)
  j <- 2
  eig_ftn[, 3] <- cos(2*pi*j*t)
  
  for (j in 1:3){
    xx <- eig_ftn[, j]
    eig_ftn[, j] <- xx/sqrt(fdapace::trapzRcpp(t, xx^2))
  }
  
  return(eig_ftn)
}


#' Generate partially observed functional data from Kraus. (2015).
#'
#' Partially observed functional data is generated with 51 regular grids by using the setting of Kraus(2015).
#'
#' @param n a number of curves
#' @param type the type of generated data. "partial" means the option for partially observed data, "snippet" is short fragmented data, and "dense" means the fully observed curves.
#' @param num.comp the number of components when the data are generated. See Kraus(2015).
#' @param out.prop a proportion of outlying curves of total n curves. Only used for dist = "normal".
#' @param out.type a outlier type, 1~3 are supported. Only used for dist = "normal".
#' @param dist a distribution which the data is generated. "normal"(Normal distribution) and "tdist"(t-distribution) are supported. If dist = "tdist", the option of \code{out.prop} and \code{out.type} are ignored.
#' @param noise a numeric value which is added random gaussian noises. Default is 0(No random noise).
#'
#' @return a list contatining as follows:
#' \item{Ly}{a list of n vectors containing the observed values for each individual.}
#' \item{Lt}{a list of n vectors containing the observation time points for each individual corresponding to \code{Ly}}
#' \item{out.ind}{a vector containing outlier index. 0 is non-outlier and 1 is the outlier.}
#' \item{x.full}{a n x 51 dense matrix with n observations per 51 timepoints before making partially observed.}
#'
#' @examples
#' set.seed(100)
#' x.list <- sim_kraus(n = 100,
#'                     type = "partial",
#'                     num.comp = 5,
#'                     out.prop = 0.2,
#'                     out.type = 1,
#'                     dist = "normal")
#' x <- list2matrix(x.list)
#' matplot(t(x), type = "l")
#'
#' @references
#' \cite{Kraus, D. (2015). Components and completion of partially observed functional data. Journal of the Royal Statistical Society: Series B: Statistical Methodology, 777-801.}
#'
#' @export
sim_kraus <- function(n = 100,
                      type = c("partial","snippet","dense"),
                      num.comp = 100,
                      out.prop = 0.2,
                      out.type = 1,
                      dist = "normal",
                      noise = 0,
                      f = 0.2) {
  gr <- seq(0, 1, length.out = 51)
  
  # generate dense curves
  m <- length(gr)   # legnth of observed grids
  
  if (dist == 'tdist') {
    out.prop <- 0   # for heavy-tailed distrubution, we do not set outlier index
    x.full <- simul.fd(n = n,
                       grid = gr,
                       lambda.cos = 3^(-(2*(1:num.comp)-1)),
                       lambda.sin = 3^(-(2*(1:num.comp))) ,
                       randcoef = t3.randcoef)
  } else if (dist == 'normal') {
    x.full <- simul.fd(n = n,
                       grid = gr,
                       lambda.cos = 3^(-(2*(1:num.comp)-1)),
                       lambda.sin = 3^(-(2*(1:num.comp))) ,
                       randcoef = norm.randcoef)
  }
  
  # random noise
  if (noise > 0) {
    x.full <- x.full + matrix(stats::rnorm(n*m, 0, sqrt(noise)), n, m)
  }
  
  x <- list()
  x$Ly <- lapply(1:n, function(i) { x.full[i, ] })
  x$Lt <- lapply(1:n, function(i) { gr })
  x$out.ind <- rep(0, n)   # indicator of outlier
  
  # Check type option
  if (type == "dense") {   # Nothing do
    x$x.full <- x.full
  } else if (type == "partial") {   # generate partially obseved data
    # generate observation periods (Kraus(2015) setting)
    # curve 1 will be missing on (.4,.7), other curves on random subsets
    x.obs <- rbind((gr <= .4) | (gr >= .7),
                   simul.obs(n = n-1, grid = gr,
                             d = 1.4, f = f)) # TRUE if observed
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
    x.outlier <- make_outlier2(x.outlier,
                               out.type = out.type,
                               gr = gr,
                               num.comp = num.comp)
    x$Ly[(n-n.outlier+1):n] <- x.outlier$Ly
    x$out.ind[(n-n.outlier+1):n] <- 1   # outlier indicator
    # x$Lt[(n-n.outlier+1):n] <- x.outlier$Lt
  } else {
    stop(paste(out.type, "is not correct value of argument out.type! Just integer value between 1~3."))
  }
  
  return(x)
}


### Generate outlying curves
make_outlier2 <- function(x, out.type = 1, gr, num.comp) {
  n <- length(x$Lt)   # number of outlying curves
  
  for (k in 1:n) {
    err.out <- simul.fd(n = 1,
                        grid = gr,
                        lambda.cos = 3^(-(2*(1:num.comp)-1)),
                        lambda.sin = 3^(-(2*(1:num.comp))) ,
                        randcoef = out.norm.randcoef)
    obs_index <- which(gr %in% x$Lt[[k]])
    x_i <- err.out[obs_index]
    x$Ly[[k]] <- as.numeric(x_i)
  }
  # print(range(x$Ly))
  
  return(x)
}





###################################
### Functions from Kraus (2015)
###################################

## functions for generating random functional data and missing periods
simul.fd = function(n = 200, grid = seq(0,1,len=200), lambda.cos = 3^(-(2*(1:300)-1)), lambda.sin = 3^(-(2*(1:300))), randcoef = norm.randcoef)
{
  x = matrix(0,n,length(grid))
  R = matrix(0,length(grid),length(grid))
  for (j in 1:length(lambda.cos)) {
    f = sqrt(lambda.cos[j])*sqrt(2)*cos(2*pi*j*grid)
    x = x + randcoef(n)%*%t(f)
    R = R + f%*%t(f)
  }
  for (j in 1:length(lambda.sin)) {
    f = sqrt(lambda.sin[j])*sqrt(2)*sin(2*pi*j*grid)
    x = x + randcoef(n)%*%t(f)
    R = R + f%*%t(f)
  }
  attr(x,"R") = R
  x
}

norm.randcoef = function(n) stats::rnorm(n,0,1)
unif.randcoef = function(n) stats::runif(n,-1,1)*sqrt(3)
t5.randcoef = function(n) stats::rt(n,5)/sqrt(5/3)
t3.randcoef = function(n) stats::rt(n,3)
out.norm.randcoef = function(n) stats::rnorm(n,12,10)


simul.obs = function(n = 100, grid = seq(0,1,len=200),d=1.4,f=.2)
{
  out = matrix(TRUE,n,length(grid))
  for (i in 1:n) {
    cen = d*sqrt(stats::runif(1))
    e = f*stats::runif(1)
    out[i,(cen-e<grid)&(grid<cen+e)] = FALSE # missing interval = (cen-u,cen+u)
  }
  out
}




############################################
### Simulation data generation functions
### Delaigle et al. (2020) simulation
############################################

#' Generate partially observed functional data from Delaigle et al.(2021).
#'
#' Partially observed functional data is generated with 51 regular grids by using the setting of Delaigle et al.(2021) with model 2.
#'
#' @param n a number of curves
#' @param type the type of generated data. "partial" means the option for partially observed data, "snippet" is short fragmented data, and "dense" means the fully observed curves.
#' @param out.prop a proportion of outlying curves of total n curves. Only used for dist = "normal".
#' @param out.type a outlier type, 1~3 are supported. Only used for dist = "normal".
#' @param dist a distribution which the data is generated. "normal"(Normal distribution) and "tdist"(t-distribution) are supported. If dist = "tdist", the option of \code{out.prop} and \code{out.type} are ignored.
#' @param noise a numeric value which is added random gaussian noises. Default is 0(No random noise).
#'
#' @return a list contatining as follows:
#' \item{Ly}{a list of n vectors containing the observed values for each individual.}
#' \item{Lt}{a list of n vectors containing the observation time points for each individual corresponding to \code{Ly}}
#' \item{out.ind}{a vector containing outlier index. 0 is non-outlier and 1 is the outlier.}
#' \item{x.full}{a n x 51 dense matrix with n observations per 51 timepoints before making partially observed.}
#'
#' @examples
#' set.seed(100)
#' x.list <- sim_delaigle(n = 100,
#'                        type = "partial",
#'                        out.prop = 0.2,
#'                        out.type = 1,
#'                        dist = "normal")
#' x <- list2matrix(x.list)
#' matplot(t(x), type = "l")
#'
#' @references
#' \cite{Delaigle, A., Hall, P., Huang, W., & Kneip, A. (2021). Estimating the covariance of fragmented and other related types of functional data. Journal of the American Statistical Association, 116(535), 1383-1401.}
#'
#' @export
sim_delaigle <- function(n = 100,
                         type = c("partial","snippet","dense"),
                         out.prop = 0.2,
                         out.type = 1,
                         dist = "normal",
                         noise = 0,
                         f = 0.2) {
  
  gr <- seq(0, 1, length.out = 51)   # equispaced points
  
  # generate dense curves
  m <- length(gr)   # legnth of observed grids
  cov_sim <- get_delaigle_cov(gr, model = 2)
  
  if (dist == 'normal') {
    y <- mvtnorm::rmvnorm(n, rep(0, m), cov_sim)
  } else if (dist == 'tdist') {
    out.prop <- 0   # for heavy-tailed distrubution, we do not set outlier index
    y <- LaplacesDemon::rmvt(n = n,
                             mu = rep(0, m),
                             S = lqmm::make.positive.definite(cov_sim),
                             df = 3)
  } else if (dist == 'laplace') {
    out.prop <- 0   # for heavy-tailed distrubution, we do not set outlier index
    y <- LaplacesDemon::rmvl(n = n,
                             mu = rep(0, m),
                             Sigma = lqmm::make.positive.definite(cov_sim))
  }
  
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
                             d = 1.4, f = f)) # TRUE if observed
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



### Get covariance function for simulation
# model = 1~2 avaliable (In the paper, 1~4 models are shown)
get_K <- function(s, t, model = 2) {
  if (model == 1) {
    K <- sqrt(5)*(6*t^2-6*t+1) * sqrt(5)*(6*s^2-6*s+1) +
      0.8 * sqrt(2)*log(t+0.5)*ifelse(t <= 0.5, 1, 0) * sqrt(2)*log(s+0.5)*ifelse(s <= 0.5, 1, 0) +
      0.3 * sqrt(22)*(252*t^5-630*t^4+560*t^3-210*t^2+30*t-1)*ifelse(t > 0.5, 1, 0) *
      sqrt(22)*(252*s^5-630*s^4+560*s^3-210*s^2+30*s-1)*ifelse(s > 0.5, 1, 0)
  } else if (model == 2) {
    K <- 1 + 0.5*(2*t-1)*(2*s-1)*3 + 0.5^2*(6*t^2-6*t+1)*(6*s^2-6*s+1)*5 +
      0.5^3*(20*t^3-30*t^2+12*t-1)*(20*s^3-30*s^2+12*s-1)*7
    
  }else if (model == 3) {
    K <-  2*sin(2*1*pi*t)*sin(2*1*pi*s)+ (2^{-2})*2*cos(2*1*pi*t)*cos(2*1*pi*s)
    + (3^{-2})* 2*sin(2*2*pi*t)*sin(2*2*pi*s)+ 2*cos(2*2*pi*t)*cos(2*2*pi*s)
    +(4^{-2})* 2*sin(2*2*pi*t)*sin(2*2*pi*s)+ 2*cos(2*2*pi*t)*cos(2*2*pi*s)
    +(5^{-2})* 2*sin(2*3*pi*t)*sin(2*3*pi*s)+ 2*cos(2*3*pi*t)*cos(2*3*pi*s)
    +(6^{-2})* 2*sin(2*3*pi*t)*sin(2*3*pi*s)+ 2*cos(2*3*pi*t)*cos(2*3*pi*s)
  }
  
  return(K)
}



#' Get True eigenfunction for simulation in Delaigle et al.(2021).
#'
#' True eigenfunctions are obtained for given grids from the setting of Delaigle et al.(2021) with model 2.
#'
#' @param grid a vector containing the observed timepoints.
#' @param model The model option for Delaigle et al.(2021).
#'
#' @return a p x k matrix containing the eigenfunctions where p is the length of grids and k is the number of true components.
#'
#' @examples
#' gr <- seq(0, 1, length.out = 51)
#' eig.true <- get_delaigle_eigen(gr, model = 2)
#' matplot(eig.true, type = "l")
#'
#' @references
#' \cite{Delaigle, A., Hall, P., Huang, W., & Kneip, A. (2021). Estimating the covariance of fragmented and other related types of functional data. Journal of the American Statistical Association, 116(535), 1383-1401.}
#'
#' @export
get_delaigle_eigen <- function(grid, model = 2) {
  t <- grid
  if (model == 1) {
    eig_ftn <- matrix(0, length(t), 3)
    eig_ftn[, 1] <-sqrt(5)*(6*t^2-6*t+1)
    eig_ftn[, 2] <-  sqrt(2)*log(t+0.5)*ifelse(t <= 0.5, 1, 0)
    eig_ftn[, 3] <-  sqrt(22)*(252*t^5-630*t^4+560*t^3-210*t^2+30*t-1)*ifelse(t > 0.5, 1, 0)
  } else if (model == 2) {
    eig_ftn <- matrix(0, length(t), 4)
    eig_ftn[, 1] <- rep(1, length(t))
    eig_ftn[, 2] <- (2*t-1)*sqrt(3)
    eig_ftn[, 3] <- (6*t^2-6*t+1)*sqrt(5)
    eig_ftn[, 4] <- (20*t^3-30*t^2+12*t-1)*sqrt(7)
  }else if (model == 3) {
    eig_ftn <- matrix(0, length(t), 6)
    k=1;eig_ftn[, 1] <- sqrt(2)*sin(2*k*pi*t)
    eig_ftn[, 2] <- sqrt(2)*cos(2*k*pi*t)
    k=2;eig_ftn[, 3] <- sqrt(2)*sin(2*k*pi*t)
    eig_ftn[, 4] <- sqrt(2)*cos(2*k*pi*t)
    k=3;eig_ftn[, 5] <- sqrt(2)*sin(2*k*pi*t)
    eig_ftn[, 6] <- sqrt(2)*cos(2*k*pi*t)
  }
  return(eig_ftn)
}


#' Get true covariance function from Delaigle et al.(2021).
#'
#' True covariance function is obtained for given grids from the setting of Delaigle et al.(2021) with model 2.
#'
#' @param grid a vector containing the observed timepoints.
#' @param model The model option for Delaigle et al.(2021).
#'
#' @return a covariance matrix with the dimension, length of grids.
#'
#' @examples
#' gr <- seq(0, 1, length.out = 51)
#' cov.true <- get_delaigle_cov(gr, model = 2)
#' library(GA)
#' persp3D(gr, gr, cov.true,
#'         theta = -70, phi = 30, expand = 1)
#'
#' @references
#' \cite{Delaigle, A., Hall, P., Huang, W., & Kneip, A. (2021). Estimating the covariance of fragmented and other related types of functional data. Journal of the American Statistical Association, 116(535), 1383-1401.}
#'
#' @export
get_delaigle_cov <- function(grid, model = 2) {
  m <- length(grid)   # legnth of observed grids
  cov_sim <- matrix(NA, m, m)
  for (i in 1:m) {
    for (j in 1:m) {
      # If upper triangular than compute, else substitute transposed value
      if (i <= j) {
        cov_sim[i, j] <- get_K(grid[i], grid[j], model = model)
      } else {
        cov_sim[i, j] <- cov_sim[j, i]
      }
    }
  }
  
  return(cov_sim)
}


###################################
### Functions from Kraus (2015)
###################################
### Functions for partially observed case
norm.randcoef = function(n) stats::rnorm(n,0,1)
unif.randcoef = function(n) stats::runif(n,-1,1)*sqrt(3)
t5.randcoef = function(n) stats::rt(n,5)/sqrt(5/3)
simul.obs <- function(n = 100, grid = seq(0, 1, len = 200), d = 1.4, f = .2) {
  out = matrix(TRUE,n,length(grid))
  for (i in 1:n) {
    cen = d*sqrt(stats::runif(1))
    e = f*stats::runif(1)
    out[i,(cen-e<grid)&(grid<cen+e)] = FALSE # missing interval = (cen-u,cen+u)
  }
  out
}


# ### get design grids
# get_ind_inter <- function(data.list) {
#   gr <- data.list$gr
#   tt <- lapply(data.list$x$t, function(t) {
#     val <- cbind(rep(t, length(t)),
#                  rep(t, each = length(t)))
#     ind <- cbind(rep(which(gr %in% val[, 1]), length(t)),
#                  rep(which(gr %in% val[, 2]), each = length(t)))
#     return(ind)
#   })
#   tt <- do.call("rbind", tt)   # rbind the argument(matrix type) in list
#   tt <- unique(tt)
#
#   return(tt)
# }

### extrapolation parts of covariance
cov_extra <- function(cov, ind) {
  cov[ind] <- 0
  
  return(cov)
}

### intrapolation parts of covariance
cov_inter <- function(cov, ind) {
  cov_ext <- cov_extra(cov, ind)
  cov <- cov - cov_ext
  
  return(cov)
}


