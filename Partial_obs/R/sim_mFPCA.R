
# simul.obs.2 <- function(n = 100, grid = seq(0, 1, len = 200), d = 1.4, f = .2) {
#   out <- matrix(TRUE, n, length(grid))
#   for (i in 1:n) {
#     cen <- d*runif(1, 0, 10)
#     e <- f*runif(1, 0, 10)
#     out[i, (cen-e < grid)&(grid < cen+e)] <- FALSE # missing interval = (cen-u,cen+u)
#   }
#   out
# }



sim_mFPCA <- function(n, type = 1, out.prop = 0.2, out.type = 1, m = 51) {
  gr <- seq(0, 1, length.out = m)
  k <- 10
  
  x.full <- matrix(0, n, m)
  t <- gr
  for (i in 1:n) {
    psi <- (1/sqrt(5))*sin(10*pi*k*t/10)
    
    if (type == 1) {
      xi <- rnorm(m, 0, sqrt(45.25*(2*(k-1))^(-2)))
      x.full[i, ] <- 10*t + sin(10*t*3) + as.numeric(xi %*% psi)
    } else if (type == 2) {
      xi <- rnorm(m, 0, sqrt(45.25*(2*k)^(-2)))
      x.full[i, ] <- 10*t + cos(10*t*3) + as.numeric(xi %*% psi)   
    }
  }
  
  
  # generate observation periods
  # curve 1 will be missing on (.4,.7), other curves on random subsets
  x.obs <- rbind((gr <= .4) | (gr >= .7), 
                 simul.obs(n = n-1, grid = gr)) # TRUE if observed
  
  # gr <- seq(0, 10, length.out = m)
  # k <- 10
  # 
  # x.full <- matrix(0, n, m)
  # t <- gr
  # for (i in 1:n) {
  #   psi <- (1/sqrt(5))*sin(pi*k*t/10)
  #   
  #   if (type == 1) {
  #     xi <- rnorm(m, 0, sqrt(45.25*(2*(k-1))^(-2)))
  #     x.full[i, ] <- t + sin(t*3) + as.numeric(xi %*% psi)
  #   } else if (type == 2) {
  #     xi <- rnorm(m, 0, sqrt(45.25*(2*k)^(-2)))
  #     x.full[i, ] <- t + cos(t*3) + as.numeric(xi %*% psi)   
  #   }
  # }
  # 
  # 
  # # generate observation periods
  # # curve 1 will be missing on (.4,.7), other curves on random subsets
  # x.obs <- rbind((gr <= 4) | (gr >= 7), 
  #                simul.obs.2(n = n-1, grid = gr)) # TRUE if observed
  # remove missing periods 
  x <- x.full
  x[!x.obs] <- NA
    
  x <- list(Ly = apply(x, 1, function(y){ y[!is.na(y)] }),
            Lt = apply(x.obs, 1, function(y){ gr[y] }),
            x.full = x.full)
  
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
    x$Lt[(n-n.outlier+1):n] <- x.outlier$Lt
  } else {
    stop(paste(out.type, "is not correct value of argument out.type! Just integer value between 1~3."))
  }
  
  return(x)
}


# ## Example
# set.seed(1000)
# m <- 101
# x.2 <- sim_mFPCA(n = 100, type = 1, out.prop = 0, out.type = 1, m = m)
# df <- data.frame(
#   id = factor(unlist(sapply(1:length(x.2$Lt),
#                             function(id) {
#                               rep(id, length(x.2$Lt[[id]]))
#                             })
#   )),
#   y = unlist(x.2$Ly),
#   t = unlist(x.2$Lt)
# )
# x <- list2matrix(x.2)
# gr <- seq(0, 1, length.out = m)
# matplot(gr, t(x), type = "l")


