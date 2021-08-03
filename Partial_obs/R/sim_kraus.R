############################################
### Simulation data generation functions
### Kraus (2015) simulation with outliers
### - Assume true components is 5
############################################

sim_kraus <- function(n = 100,
                      type = c("partial","snippet","sparse","dense"),
                      num.comp = 100,
                      out.prop = 0.2, out.type = 1) {
  gr <- seq(0, 1, length.out = 51)
  
  # generate dense curves
  x.full <- simul.fd(n = n, grid = gr,
                     lambda.cos = 3^(-(2*(1:num.comp)-1)), 
                     lambda.sin = 3^(-(2*(1:num.comp))))
  x <- list()
  x$Ly <- lapply(1:n, function(i) { x.full[i, ] })
  x$Lt <- lapply(1:n, function(i) { gr })
  
  # Check type option
  if (type == "dense") {   # Nothing do
    x$x.full <- x.full
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
    # x$Lt[(n-n.outlier+1):n] <- x.outlier$Lt
  } else {
    stop(paste(out.type, "is not correct value of argument out.type! Just integer value between 1~3."))
  }
  
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

norm.randcoef = function(n) rnorm(n,0,1)
unif.randcoef = function(n) runif(n,-1,1)*sqrt(3)
t5.randcoef = function(n) rt(n,5)/sqrt(5/3)

simul.obs = function(n = 100, grid = seq(0,1,len=200),d=1.4,f=.2)
{
  out = matrix(TRUE,n,length(grid))
  for (i in 1:n) {
    cen = d*sqrt(runif(1))
    e = f*runif(1)
    out[i,(cen-e<grid)&(grid<cen+e)] = FALSE # missing interval = (cen-u,cen+u)
  }
  out
}
