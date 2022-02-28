############################################################
### Principal expectile component analysis
### Tran et al. (2019)
### - PrincipalExpectile algorithm
############################################################

### Principal expectile component analysis using PrincipalExpectile algorithm
### - Y : nxp data matrix (dense functional data)
### - tau : tau percent quantile between 0 and 1
### - k : the number of PECs
### - max.iter : the number of maximum iterations for each PEC
peca <- function(Y, tau = 0.5, k = 2, max.iter = 30) {
  n <- nrow(Y)   # number of curves
  p <- ncol(Y)   # number of timepoints
  
  ### Each PEC is obtained by sequentially
  pec <- matrix(0, nrow = p, ncol = k)   # principal expectile component
  iter.pec <- rep(0, k)   # number of iteration for each component
  resid <- Y
  for (i in 1:k) {
    # Obtain residuals
    if (i > 1) {
      resid <- resid - (resid %*% v + mu_tau) %*% t(v)
    }
    
    # Initialize weights
    w <- ifelse(runif(n) > 0.5, tau, 1-tau)
    w_before <- rep(0, n)
    
    ### Do PrincipalExpectile algorithm (See, Algorithm 1 in Tran et al.(2019))
    iter <- 0   # initial iteration
    while (sum(w - w_before) != 0) {
      # Obtain index set
      tau_pos <- which(w == tau)
      tau_neg <- which(w != tau)
      
      # Compute e_tau
      num <- tau * colSums(resid[tau_pos, ]) + (1-tau) * colSums(resid[tau_neg, ])
      den <- tau * length(tau_pos) + (1-tau) * length(tau_neg)
      e_tau <- num / den
      
      # Compute C_tau
      Y_center <- resid - matrix(rep(e_tau, n), ncol = p, byrow = T)
      C_tau <- (tau/n) * (t(Y_center[tau_pos, ]) %*% Y_center[tau_pos, ] ) + ((1-tau)/n) * (t(Y_center[tau_neg, ]) %*% Y_center[tau_neg, ] )
      
      # Set v into largest eivenvector of C_tau %*% C_tau^T
      eig <- eigen(C_tau %*% t(C_tau), symmetric = T)
      v <- matrix(eig$vectors[, 1], ncol = 1)
      
      # Set mu_tau into tau-expectile of v^T %*% Y_i
      # t(v) %*% Y[1, ]
      v_Y <- as.numeric(resid %*% v)
      mu_tau <- expectile(v_Y, probs = tau)
      
      # Update w_i
      w_before <- w
      w <- ifelse(v_Y > mu_tau, tau, 1-tau)
      
      # Update iteration number
      iter <- iter + 1
      
      # Too many iterations, breaks
      if (iter > max.iter) {
        break
      }
    }
    
    if (i == 1) {
      tau_expect <- e_tau   # tau-expectile
    }
    pec[, i] <- v
    iter.pec[i] <- iter
  }
  
  # Output
  obj <- list(
    center = tau_expect,  # tau-expectile (corresponding to mean function)
    pec = pec,            # first k PECs
    num.iter = iter.pec   # the iteration number of the algorithm
  )
  
  return(obj)
}


### Get expectile
### - x : vector
### - probs : quantile between 0 and 1
### - max.iter : the maximum iteration numbers
expectile <- function(x, probs = 0.5, max.iter = 30) {
  m <- mean(x)
  expect <- 0 * probs
  tol <- max(abs(x)) * 1e-6
  
  for (i in 1:length(probs)) {
    p <- probs[i]
    
    if (p == 0) {
      expect[i] <- min(x, na.rm = TRUE)
    } else if (p == 1) {
      expect[i] <- max(x, na.rm = TRUE)
    } else {
      for (iter in 1:max.iter) {
        w <- ifelse(x < m, 1 - p, p)
        m_new <- sum(w * x) / sum(w)
        m_diff <- max(abs(m_new - m))
        m <- m_new
        if (m_diff < tol) {
          break
        }
      }
      expect[i] <- m
    }
  }
  
  names(expect) <- probs
  
  return(expect)
}




############################################
### Example
############################################

### Generate example data
set.seed(100)
x <- robfpca::sim_delaigle(n = 100,
                           type = "dense",
                           out.prop = 0,
                           dist = "tdist")
Y <- robfpca::list2matrix(x)   # transform list to matrix
matplot(t(Y), type = "l")

### Check the method is same with conventional PCA
peca.obj <- peca(Y, tau = 0.5, k = 3)    # PEC
pca.obj <- eigen(cov(Y))$vectors[, 1:3]  # conventional PCA 

matplot(peca.obj$pec, type = "l", col = 1)
matlines(pca.obj, col = 2)


### 80% expectile components
peca.obj2 <- peca(Y, tau = 0.9, k = 3)

matplot(peca.obj$pec, type = "l", col = 1)   # PCA
matlines(peca.obj2$pec, col = 2)   # 80% PEC


### Mean and expectile function
plot(colMeans(Y), type = "l", ylim = c(-1, 1))   # meanfunction
lines(peca.obj$center, col = 2)   # 50% expectile function (same with mean)
lines(peca.obj2$center, col = 3)  # 80% expectile function




