###############################################################
### k-Centres Riemannian Functional Clustering (kCRFC)
### - RFPCA (Dai and Muller, 2018) + kCFC (Chiou and Li, 2007)
###############################################################

### k-Centres Riemannian Functional Clustering (kCRFC)
### (i.e. kCFC for Riemannian functional data)
#### Algorithm
#### 1. Initial clustering of functional principal component scores
####  1-1. FPCA using overall data
####  1-2. k-means clustering using FPC scores
#### 2. Iterative updating via reclassification
####  2-1. FPCA for each cluster
####       For a cluster in ith curve, FPCA is performed without ith observation.
####  2-2. Predict ith curve for each cluster
####  2-3. Assign new cluster which minimizes the L2 norm between ith curve and prediction for each cluster
#### 3. Repeat 2 until no more curves are classified.
kCRFC <- function(y, 
                  t, 
                  k = 3,
                  kSeed = 123, 
                  maxIter = 125, 
                  fast = TRUE,
                  optnsSW = list(mfdName = "Sphere",
                                 FVEthreshold = 0.90,
                                 userBwMu = "GCV", 
                                 userBwCov = "GCV"),
                  optnsCS = list(mfdName = "Sphere",
                                 FVEthreshold = 0.70, 
                                 userBwMu = 'GCV', 
                                 userBwCov = 'GCV')) { 
  
  if( (k <2) || (floor(length(y)*0.5) < k) ){
    warning("The value of 'k' is outside [2, 0.5*N]; reseting to 3.")
  } 
  if(maxIter <1){
    stop("Please allow at least 1 iteration.")
  }
  
  ## First RFPCA and threshold by FVE
  fpcaObjY <- RFPCA(Ly = y, 
                    Lt = t, 
                    optns = optnsSW)
  # fpcaObjY <- RFPCA.FVE(RFPCA.obj = fpcaObjY, 
  #                       Lt = t, Ly = y, 
  #                       FVEthreshold = optnsSW$FVEthreshold)
  N <- length(y)
  if( fpcaObjY$optns$dataType == 'Sparse' ){
    stop(paste0("The data has to be 'Dense' for kCFC to be relevant; the current dataType is : '", fpcaObjY$optns$dataType,"'!") )
  }
  
  ## Initial clustering and cluster-associated FPCAs
  if(!is.null(kSeed)){
    set.seed(kSeed)
  }
  initialClustering <- kmeans(fpcaObjY$xi, centers = k, algorithm = "MacQueen", 
                              iter.max = maxIter, nstart = 50)
  clustConf0 <- as.factor(initialClustering$cluster)
  indClustIds <- lapply(levels(clustConf0), function(u) which(clustConf0 == u) )
  if( any( min( sapply( indClustIds, length)) <= c(3)) ){
    stop(paste0("kCFC stopped during the initial k-means step. The smallest cluster has three (or less) curves. " ,
                "Consider using a smaller number of clusters (k) or a different random seed (kSeed)."))
  }
  listOfFPCAobjs <- lapply(indClustIds, function(u){
    RFPCA(Ly = y[u], 
          Lt = t[u],
          optns = optnsCS) 
    # obj <- RFPCA(Ly = y[u], 
    #              Lt = t[u],
    #              optns = optnsCS) 
    # RFPCA.FVE(RFPCA.obj = obj,
    #           Lt = t[u], Ly = y[u], 
    #           FVEthreshold = optnsCS$FVEthreshold)
  })
  
  ## Iterative clustering
  convInfo <- "None"
  clustConf <- list() 
  
  for(j in seq_len(maxIter)){ 
    
    # Get new costs and relevant cluster configuration
    if (isTRUE(fast)) {
      iseCosts <- sapply(listOfFPCAobjs, function(u){ GetISEfromRFPCA_fast(u, y, t) })
    } else {
      iseCosts <- sapply(1:k, function(u){ GetISEfromRFPCA(u, listOfFPCAobjs[[u]], y, t,
                                                           indClustIds, optnsCS) })
    }
    clustConf[[j]] <- as.factor(apply(iseCosts, 1, which.min))
    
    # Check that clustering progressed reasonably 
    #ie. Still having k clster AND the minimum cluster size is reasonable 
    if( (length(unique(clustConf[[j]])) < k) || any( min(summary(clustConf[[j]])) <= c(0.01 * N,3)) ){
      convInfo <- ifelse( length(unique(clustConf[[j]])) < k , "LostCluster", "TinyCluster")
      break;
    }
    # Check if algorithm converged
    if( (j>1) && any(sapply(clustConf[1:(j-1)], function(u) all(u == clustConf[[j]]))) ){
      convInfo <- "WeMadeIt!"
      break;
    } 
    
    indClustIds       <- lapply(levels(clustConf[[j]]), function(u) which(clustConf[[j]] == u) )
    listOfFPCAobjs    <- lapply(indClustIds, function(u){ RFPCA(Ly = y[u], Lt = t[u], optns = optnsCS) })
    curvesThatChanged <- ifelse(j > 1, sum(!( as.numeric(clustConf[[j]])  == as.numeric(clustConf[[j-1]] ))),
                                sum(!( as.numeric(clustConf[[j]])  == as.numeric(clustConf0))))
  } 
  
  if(convInfo == 'None'){
    warning(paste0( 'FkC did not converge after maxIter = ', maxIter, ' iterations. ', curvesThatChanged, ' curve(s) are undecided.'))
  }
  if(convInfo == 'TinyCluster'){
    warning(paste0("kCFC did not fully converge. It stopped because the smallest cluster has ",
                   "less than 1% of the samples' curves. Consider using a smaller number of clusters."))
  } 
  if(convInfo == 'LostCluster'){
    warning(paste0("kCFC did not fully converge. It stopped because it 'lost a cluster'. Consider using a smaller number of clusters."))
  }
  
  kCFCobj <-  list(cluster = clustConf[[j]], fpcaList = listOfFPCAobjs, iterToConv = j-1, prevConf = clustConf, clustConf0 = clustConf0)
  class(kCFCobj) <- 'kCFCobj'
  return( kCFCobj )
}




### Obtain ISE using geodesic distance
GetISEfromRFPCA <- function(cluster, fpcaObj, y, t,
                            indClustIds, optnsCS) {
  obs.grid <- fpcaObj$obsGrid
  mfd <- fpcaObj$mfd
  n <- length(t)
  d <- nrow(y[[1]])
  p <- ncol(y[[1]])
  
  pred <- array(0, dim = c(n, d, p))
  # Reconstruction for same cluster
  idx <- indClustIds[[cluster]]
  for (i in idx) {
    ind <- setdiff(idx, i)   # remove ith curve on cluster
    fpcaObj_ind <- RFPCA(Ly = y[ind], Lt = t[ind], optns = optnsCS)
    pred[i, , ] <- predict(object = fpcaObj_ind,
                           newLt = t[i],
                           newLy = y[i],
                           # K = k,
                           xiMethod = "IN",
                           type = "traj")
  }
  
  # Reconstruction using K components
  idx_other_cluster <- unlist(indClustIds[-cluster])
  pred[idx_other_cluster, , ] <- predict(object = fpcaObj,
                                         newLt = t[idx_other_cluster],
                                         newLy = y[idx_other_cluster],
                                         # K = k,
                                         xiMethod = "IN",
                                         type = "traj")
  
  ise <- sapply(1:n, function(i){
    d_0 <- distance(mfd = mfd,
                    X = y[[i]],
                    Y = pred[i, , ])
    trapzRcpp(X = obs.grid, 
              Y = d_0^2)
  })
  
  return(ise)
}

GetISEfromRFPCA_fast <- function(fpcaObj, y, t) {
  obs.grid <- fpcaObj$obsGrid
  mfd <- fpcaObj$mfd
  n <- length(t)
  
  # Reconstruction using K components
  pred <- predict(object = fpcaObj,
                  newLt = t,
                  newLy = y,
                  # K = k,
                  xiMethod = "IN",
                  type = "traj")
  
  ise <- sapply(1:n, function(i){
    d_0 <- distance(mfd = mfd,
                    X = y[[i]],
                    Y = pred[i, , ])
    trapzRcpp(X = obs.grid, 
              Y = d_0^2)
  })
  
  return(ise)
}


### Get fraction of varianc explained (FVE)
### - It does not necessary. (RFPCA에서 옵션 주면 eigenanalysis에서 해줌)
RFPCA.FVE <- function(RFPCA.obj, 
                      Lt, Ly, 
                      FVEthreshold = 0.95) {
  mfd <- RFPCA.obj$mfd
  # work.grid <- RFPCA.obj$workGrid
  obs.grid <- RFPCA.obj$obsGrid
  K <- RFPCA.obj$K
  
  # Null residual variance, U_0
  U_0 <- sapply(Ly, function(y) {
    d_0 <- distance(mfd = mfd,
                    X = y,
                    Y = RFPCA.obj$muObs)
    return( trapzRcpp(obs.grid, d_0^2) )
    # d_0 <- distance(mfd = mfd,
    #                 X = y,
    #                 Y = RFPCA.obj$muWork)
    # return( trapzRcpp(work.grid, d_0^2) )
  })
  U_0 <- mean(U_0)
  
  FVE <- rep(0, K)
  for (k in 1:K) {
    # Reconstruction using K components
    pred <- predict(object = RFPCA.obj,
                    newLt = Lt,
                    newLy = Ly,
                    K = k,
                    xiMethod = "IN",
                    type = "traj")
    
    # Residual variance using K components, U_K
    U_k <- sapply(1:n, function(i) {
      d_k <- distance(mfd = mfd,
                      X = Ly[[i]],
                      Y = pred[i, , ])
      return( trapzRcpp(obs.grid, d_k^2) )
    })
    U_k <- mean(U_k)
    
    FVE[k] <- (U_0 - U_k) / U_0
  }
  
  K <- min( which(FVE > FVEthreshold) )
  RFPCA.obj$phi <- RFPCA.obj$phi[, , 1:K]
  RFPCA.obj$lam <- RFPCA.obj$lam[1:K]
  RFPCA.obj$xi <- RFPCA.obj$xi[, 1:K]
  RFPCA.obj$FVE <- FVE[1:K]
  RFPCA.obj$FVEthreshold <- FVEthreshold
  RFPCA.obj$K <- K
  
  return(RFPCA.obj)
}


### Utility functions
array2list <- function(X, t){
  n <- dim(X)[1]
  Ly <- list()
  Lt <- list()
  for (i in 1:n) {
    Ly[[i]] <- X[i, , ]
    Lt[[i]] <- t
  }

  return(list(Ly = Ly,
              Lt = Lt))
}


list2array <- function(Ly) {
  n <- length(Ly)   # number of curves
  m <- nrow(Ly[[1]])   # number of axis of manifolds
  p <- ncol(Ly[[1]])   # number of timepoints
  Y <- array(0, dim = c(n, m, p))
  
  for (i in 1:n) {
    Y[i, , ] <- Ly[[i]]
  }
  
  return(Y)
}


# list2matrix <- function(Ly, Lt){
#   n = length(Ly)
#   obsGrid = sort(unique(unlist(Lt)))
#   ymat = matrix( rep(NA, n * length(obsGrid)), nrow = n, byrow = TRUE)
#   
#   for (i in 1:n){
#     ymat[i, is.element(obsGrid, Lt[[i]])] = Ly[[i]]   
#   }
#   return(ymat)
# }


### Generate simulated data
# n : the number of curves
# type : the type of spaces. (Default is "Sphere")
# k : the number of clusters
simul_data <- function(n = 100, type = "Sphere", k = 2) {
  
  # n <- 100  # number of curves
  m <- 51   # number of different time points
  K <- 20   # number of components
  # k <- 2    # number of clusters
  
  # Generate for each cluster
  Lt <- list()
  Ly <- list()
  mu_list <- list()
  cluster <- rep(1:k, each = n/k)
  for (i in 1:k) {
    lambda <- (i*0.07)^(seq_len(K) / 2)
    # D <- 3
    basisType <- 'legendre01'
    sigma2 <- 0
    muList <- list(
      function(x) x * 2,
      function(x) sin(x * 1 * pi * i) * pi / 2 * 0.6,
      function(x) rep(0, length(x))
      # function(x) x * 2, 
      # function(x) sin(x * 1 * pi) * pi / 2 * 0.6,
      # function(x) rep(0, length(x))
    )
    pts <- seq(0, 1, length.out = m)
    mfd <- structure(1, class = type)
    mu <- Makemu(mfd, muList, c(0, 0, 1), pts)
    
    # Generate samples
    samp <- MakeMfdProcess(mfd = mfd, 
                           n = n/k, 
                           mu = mu, 
                           pts = pts, 
                           K = K, 
                           lambda = lambda, 
                           basisType = basisType, 
                           sigma2 = sigma2)
    # sparsity <- m
    # # spSamp <- SparsifyM(samp$X, samp$T, sparsity)
    spSamp <- array2list(samp$X, samp$T)
    Ly <- c(Ly, spSamp$Ly)
    Lt <- c(Lt, spSamp$Lt)
    mu_list <- c(mu_list, list(mu))
  }
  
  data <- list(
    Lt = Lt,
    Ly = Ly,
    mu = mu_list
  )
  
  # theta_1 <- c(0.4, 0.3)
  # theta_2 <- c(0.2, 0.1)
  # phi <- list(
  #   cbind(sqrt(2)*sin(pi*gr),
  #         sqrt(2)*cos(pi*gr)),
  #   cbind(sqrt(2)*sin(2*pi*gr),
  #         sqrt(2)*cos(2*pi*gr))
  # )
  
  
  return(data)
}




################################################
### Simulation setting from Dai (2022)
### - Only different with mu_1 = (0, 0, 1)
################################################

# 지금 만드는 function 완성하고 체크해보기

function(mfd) {
    mfd <- structure(1, class = 'Sphere')
    
    
    n_1 <- 50
    n_2 <- 50
    m <- 51   # number of timepoints
    gr <- seq(0, 1, length.out = m)   # timepoints
    
    K_mu <- 5   # 1 or 3 or 5 (number of mean components)
    K_X <- 50   # number of components
    
    theta <- 3^{-(1:K_X)}   # variance
    xi <- sapply(theta, function(th){ rnorm(n_1, sd = th) })   # PC score
    
    psi <- cbind(
        rep(1, m),
        sapply(2:K_X, function(k){
            if (k %% 2 == 0) { 
                sqrt(2) * sin(2*k*pi*gr)
            } else {
                sqrt(2) * cos(2*(k-1)*pi*gr)
            }
        })
    )
    
    
    # mu_1 <- c(0, 0, 1)
    mu_1 <- sqrt(rbeta(3, 2, 1))
    mu_1 <- mu_1 / sqrt(sum(mu_1^2))
    
    inner_prod_hilbert <- function(f, g, gr) { fdapace::trapzRcpp(gr, f*g) }
    rho <- function(f, g, gr){ acos( inner_prod_hilbert(f, g, gr) ) }
    R <- function(q, p) {
        rho_q <- rho(psi[, 1], q, gr)
        u <- (psi[, 1] - rho_q*q) / sqrt(1 - rho_q^2)
        up <- inner_prod_hilbert(u, p, gr)
        qp <- inner_prod_hilbert(q, p, gr)
        return( p + sin(rho_q)*(up*q - qp*u) + (cos(rho_q) - 1)*(qp*q + up*u) )
    }
    phi_1 <- sapply(1:K_X, function(k){ R(mu_1, psi[, k]) })
    
    
    rieExp(mfd, mu_1, xi %*% t(phi_1))
    
    
    
    ### Group 2
    v <- K_mu^{-1/2} * rowSums(phi_1[, 1:K_mu])
    delta <- runif(1, -0.4, 0.4)
    # mu_2 <- RFPCA::rieExp(mfd, mu_1, delta * v)
    
    mu_2 <- sapply(mu_1, function(mu_1i){ 
        RFPCA::rieExp(mfd, mu_1i, matrix(delta * v, nrow = 1))
    })
    RFPCA::rieExp(mfd, 0, matrix(delta * v, nrow = 1))
    
    rieExp(createM("Euclidean"), mu_1, matrix(rep(delta * v, 3), nrow = 3, byrow = T))
    
    V <- matrix(rep(delta * v, 3), nrow = 3, byrow = T)
    RFPCA::rieExp(mfd, mu_1, V)
    
    phi_1 <- d
}



    
