
### Local covariance
local_cov_Mest <- function(x, h = 0.02) {
  # n <- nrow(x)
  
  Lt <- x$Lt
  Ly <- x$Ly
  
  t <- unlist(Lt)
  y <- unlist(Ly)
  
  gr <- seq(0, 1, length.out = 51)
  p <- length(gr)
  
  cov_mat <- matrix(NA, p, p)
  for (i in 1:p) {
    i_neighbor <- which(t < gr[i] + h & t > gr[i] - h)
    mu_i <- mean(y[i_neighbor])
    for (j in 1:p) {
      j_neighbor <- which(t < gr[j] + h & t > gr[j] - h)
      mu_j <- mean(y[j_neighbor])
      
      ind <- intersect(i_neighbor, j_neighbor)
      # tmp <- (y[i_neighbor] - mu_i) * (y[j_neighbor] - mu_j)
      tmp <- (y[ind] - mu_i) * (y[ind] - mu_j)
      tmp <- as.numeric(tmp)
      
      cov_mat[i, j] <- mean(tmp)
    }
  }
  
  return(cov_mat)  
}


x.3 <- list(Lt = x.2$Lt[1:80],
            Ly = x.2$Ly[1:80])
system.time({
  cc <- local_cov_Mest(x.3, h = 0.02)  
})
dim(cc)


gr <- work.grid
GA::persp3D(gr, gr, cc,
            theta = -70, phi = 30, expand = 1)

hist(unlist(x.2$Ly), probability = T, ylim = c(0, 0.4), breaks = 100, xlim = c(-10, 10))
lines(density(unlist(x.2$Ly)), col = 2, lwd = 3)
lines(seq(-6, 6, length.out = 101),
      dnorm(seq(-6, 6, length.out = 101)), col = 3, lwd = 3)







get_raw_cov <- function(x) {
  gr <- seq(0, 1, length.out = 51)
  n <- nrow(x)
  p <- ncol(x)
  
  mu <- mean_Mest(x)
  
  for (i in 1:p) {
    ind <- which(!is.na(x[i, ]))
    gr_sub <- gr[ind]
    x_sub <- x[i, ind]
    
    exp_grid <- expand.grid(gr_sub, gr_sub)
    exp_x <- expand.grid(x_sub, x_sub)
    
    # eliminate diagonal parts
    idx <- which(exp_grid[, 1] == exp_grid[, 2])
    cov_i <- exp_x[-idx, 1] * exp_x[-idx, 2]
    
    if (i == 1) {
      raw_cov <- cbind(cov_i,
                       exp_grid[-idx, ])
    } else {
      raw_cov <- rbind(raw_cov,
                       cbind(cov_i,
                             exp_grid[-idx, ]))
    }
  }
  
  return(raw_cov)  
}


get_smooth_cov <- function(x, bw) {
  cov_raw <- data.frame(get_raw_cov(x))
  
  fit <- loess(cov_i ~ Var1 + Var2, 
               data = cov_raw,
               span = bw, 
               family = "symmetric",
               degree = 1)
  
  gr <- seq(0, 1, length.out = 51)
  exp_grid <- expand.grid(gr, gr)
  cov_sm <- matrix(predict(fit, exp_grid), 51, 51)
  
  return(cov_sm)
}


cov_sm <- get_smooth_cov(x, 0.02)
GA::persp3D(gr, gr, cov_sm,
            theta = -70, phi = 30, expand = 1)


### Eigen function trajectories
eig.true <- get_delaigle_eigen(work.grid, model = 2)
# eig.true <- get_eigen(cov.true, gr)$phi[, 1:4]
eig.huber <- get_eigen(cov.huber, gr)$phi[, 1:4]
eig.boente <- get_eigen(cov.boente, gr)$phi[, 1:4]
eig.Mest <- get_eigen(cov.Mest, gr)$phi[, 1:4]
eig.Mest.sm <- get_eigen(cov.Mest.sm, gr)$phi[, 1:4]
eig.cc.sm <- get_eigen(cov_sm, gr)$phi[, 1:4]
par(mfrow = c(2, 2))
for (k in 1:4) {
  matplot(work.grid,
          cbind(
            eig.true[, k],
            check_eigen_sign(eig.huber, eig.true)[, k],
            check_eigen_sign(eig.boente, eig.true)[, k],
            check_eigen_sign(eig.Mest, eig.true)[, k],
            check_eigen_sign(eig.Mest.sm, eig.true)[, k],
            check_eigen_sign(eig.cc.sm, eig.true)[, k]
          ),
          type = "l",
          col = 1:6,
          lty = 1:6,
          main = paste("Eigenfunction", k),
          xlab = "", ylab = "",
          lwd = rep(2, 4))
  if (k == 1) {
    legend("bottomleft",
           c("True","Huber","Boente","M-est","M-est(smooth)","cc_sm"),
           col = 1:7,
           lty = 1:7,
           lwd = rep(2, 7))
  }
}
