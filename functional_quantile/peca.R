############################################
### Principal expectile component analysis
### Tran et al. (2019)
### - PrincipalExpectile algorithm
############################################

### Example
# devtools::install_github("obleeker/quant.pca")
library(quant.pca)

# generate data
set.seed(100)
n = 100
X = data.frame(cbind(rnorm(n),rnorm(n),rnorm(n)))

# running PQC on the 0.9-quantile 
pec.fit <- pqcomp(data = Y, 
                  projDim = 2,
                  tau = 0.9, 
                  lambda = 0.1, 
                  muEst = T)
matplot(pec.fit$loadings, type = "l")


peca.obj <- peca(Y, tau = 0.9)

plot(pec.fit$loadings[, 1], type = "l")
lines(peca.obj$pec, col = 2)


set.seed(100)
# n <- 100
# p <- 51
# Y <- mvtnorm::rmvnorm(n = n, mean = rep(0, p))
x <- robfpca::sim_delaigle(n = 100,
                           type = "dense",
                           out.prop = 0,
                           dist = "normal")
Y <- robfpca::list2matrix(x)   # transform list to matrix
matplot(t(Y), type = "l")

peca.obj <- peca(Y, tau = 0.7)

matplot(t(Y), type = "l", col = "gray")
lines(peca.obj$pec, type = "l", lwd = 2)


### Principal expectile component analysis using PrincipalExpectile algorithm
### 확인이 필요
peca <- function(Y, tau = 0.5, max.iter = 30) {
  n <- nrow(Y)
  p <- ncol(Y)
  
  # tau <- 0.7   # tau percent quantile
  # max.iter <- 30   # the number of maximum iterations
  
  iter <- 0   # initial iteration
  # w <- matrix(ifelse(runif(n*p) > tau, tau, 1-tau), ncol = p)   # initial weights (다시 체크하기!!)
  w <- ifelse(runif(n) > tau, tau, 1-tau)
  w_update <- rep(0, n)
  
  # Iteration
  while (sum(w - w_update) != 0) {
    tau_pos <- which(w == tau)
    tau_neg <- which(w != tau)
    
    # Compute e_tau
    num <- tau * colSums(Y[tau_pos, ]) + (1-tau) * colSums(Y[tau_neg, ])
    den <- tau * length(tau_pos) + (1-tau) * length(tau_neg)
    e_tau <- num / den
    
    # Compute C_tau
    Y_center <- Y - matrix(rep(e_tau, n), ncol = p, byrow = T)
    C_tau <- (tau/n) * (t(Y_center[tau_pos, ]) %*% Y_center[tau_pos, ] ) + ((1-tau)/n) * (t(Y_center[tau_neg, ]) %*% Y_center[tau_neg, ] )
    
    # Set v into largest eivenvector of C_tau %*% C_tau^T
    eig <- eigen(C_tau %*% t(C_tau), symmetric = T)
    v <- matrix(eig$vectors[, 1], ncol = 1)
    
    # Set mu_tau into tau-expectile of v^T %*% Y_i
    # t(v) %*% Y[1, ]
    v_Y <- as.numeric(Y %*% v)
    mu_tau <- expectile(v_Y, probs = tau)
    
    # Update w_i
    w_update <- ifelse(as.numeric(Y %*% v) > mu_tau, tau, 1-tau)
    
    iter <- iter + 1
    
    if (iter > max.iter) {
      break
    }
  }
  
  obj <- list(
    pec = v,   # 1st principal expectile component of Y
    num.iter = iter
  )
  
  return(obj)
}



### From "expectreg" 수정해야됨!!
expectile <- function(x, probs = 0.5, dec = 4) {
  if (!is.vector(x)) 
    stop("observations are needed in vector form.")
  if (min(probs) < 0 || max(probs) > 1) 
    stop("only asymmetries between 0 and 1 allowed.")
  e = mean(x)
  ee = 0 * probs
  g = max(abs(x)) * 1e-06
  for (k in 1:length(probs)) {
    p = probs[k]
    if (p == 0) 
      ee[k] = min(x, na.rm = TRUE)
    else if (p == 1) 
      ee[k] = max(x, na.rm = TRUE)
    else {
      for (it in 1:20) {
        w = ifelse(x < e, 1 - p, p)
        enew = sum(w * x)/sum(w)
        de = max(abs(enew - e))
        e = enew
        if (de < g) 
          break
      }
      ee[k] = e
    }
  }
  names(ee) = probs
  ee = round(ee, dec)
  return(ee)
}
  
  
