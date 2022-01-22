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
                  optnsSW = list(mfdName = "Sphere",
                                 FVEthreshold = 0.90,
                                 userBwMu = "GCV", 
                                 userBwCov = "GCV"),
                  # kernel = "epan", 
                  # maxK=K, 
                  # error = FALSE),
                  # optnsSW = list(methodMuCovEst = 'smooth', 
                  #                FVEthreshold = 0.90, 
                  #                methodBwCov = 'GCV', 
                  #                methodBwMu = 'GCV'), 
                  # optnsCS = list(methodMuCovEst = 'smooth', 
                  #                FVEthreshold = 0.70, 
                  #                methodBwCov = 'GCV', 
                  #                methodBwMu = 'GCV')
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
  initialClustering <- kmeans(fpcaObjY$xi, centers = k, algorithm = "MacQueen", iter.max = maxIter)
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
    iseCosts <- sapply(listOfFPCAobjs, function(u){ GetISEfromRFPCA(u, y, t) })
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
GetISEfromRFPCA <- function(fpcaObj, y, t) {
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



