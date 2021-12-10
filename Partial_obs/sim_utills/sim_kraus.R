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
                      noise = 0) {
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
