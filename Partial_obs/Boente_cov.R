library(sparseFPCA)

### mean and covariance in Boente (2020)
cov_boente <- function(x, bw.mu, bw.cov) {
  X <- list(x = x$Ly,
            pp = x$Lt)
  gr <- sort(unique(unlist(x$Lt)))
  
  mh <- mu.hat3.lin(X=X, h=bw.mu)
  ma <- matrixx(X, mh)
  
  # Compute the estimated cov function
  cov.fun2 <- cov.fun.hat2.ls(X=X, h=bw.cov, mh=mh, ma=ma, ncov=length(gr), trace=FALSE)
  # smooth it
  yy <- as.vector(cov.fun2$G)
  xx <- cov.fun2$grid
  tmp <- fitted(mgcv::gam(yy ~ s(xx[,1], xx[,2]), family='gaussian'))
  cov.fun2$G <- matrix(tmp, length(unique(xx[,1])), length(unique(xx[,1])))
  cov.fun2$G <- ( cov.fun2$G + t(cov.fun2$G) ) / 2
  
  # obtain mean function
  df <- data.frame(t = unlist(X$pp),
                   mu = unlist(mh))
  df <- unique(df)
  idx <- sort(df$t, index.return = T)$ix
  mu <- df$mu[idx]
  
  return(list(mu = mu,
              cov = cov.fun2$G))
}
