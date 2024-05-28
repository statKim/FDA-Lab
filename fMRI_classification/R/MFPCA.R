# Univariate Functional Principal Component Analysis
# Inputs:
#    X: A n-m matrix (n curves observed from m timepoints)
#    grid: A vector containing m timepoints
#    K: Number of FPCs
#    FVE: Fraction of variance explained (Default is 0.95)
#    centered: If TRUE, centering of curves is not conducted. (Default is FALSE)
uFPCA <- function(X, grid = NULL, K = NULL, FVE = 0.95, centered = FALSE) {

    n <- nrow(X)   # number of curves
    m <- ncol(X)   # number of timepoints

    # Observed grid points
    if (is.null(grid)) {
        grid <- seq(0, 1, length.out = m)
    }

    # Centering the curves
    if (!isTRUE(centered)) {
        mu <- colMeans(X)
        X <- X - do.call(rbind, rep(list(mu), n))
    } else {
        mu <- rep(0, length.out = m)
    }

    # SVD
    svd.obj <- svd(t(X) / sqrt(n))
    positive_ind <- which(svd.obj$d > 0)   # indices of positive eigenvalues
    lambda <- svd.obj$d[positive_ind]^2   # eigenvalues
    phi <- svd.obj$u[, positive_ind]   # eigenvectors

    # # Eigen analysis
    # cov.mat <- crossprod(X, X) / n
    # eig.obj <- eigen(cov.mat)
    # positive_ind <- which(eig.obj$values > 0)   # indices of positive eigenvalues
    # lambda <- eig.obj$values[positive_ind]^2   # eigenvalues
    # phi <- eig.obj$vectors[, positive_ind]   # eigenvectors


    # Normalizing eigenvalues and eigenfunctions
    grid_size <- grid[2] - grid[1]
    work.grid <- seq(min(grid), max(grid), length.out = m)
    lambda <- lambda * grid_size
    phi <- apply(phi, 2, function(x) {
        x <- x / sqrt(grid_size)
        # x <- x / sqrt(sum(x^2))   # normalize
        x <- x / sqrt(fdapace::trapzRcpp(work.grid, x^2))
        if ( 0 <= sum(x * work.grid) )
            return(x)
        else
            return(-x)
    })


    # Select the number of FPCs
    cum_FVE <- cumsum(lambda) / sum(lambda)   # cumulative FVE
    if (is.null(K)) {
        K <- which(cum_FVE >= FVE)[1]
    }
    if (K > length(lambda)) {
        stop("K is too large.")
    }


    # FPC scores
    lambda <- lambda[1:K]
    phi <- matrix(phi[, 1:K], ncol = K)
    xi <- get_fpc_scores(X, phi, work.grid)   # C++ function
    # # xi <- matrix(X %*% phi, ncol = K)
    # xi <- sapply(1:K, function(k) {
    #     apply(t(X) * phi[, k], 2, function(x_i_phi_k) {
    #         fdapace::trapzRcpp(work.grid, x_i_phi_k)
    #     })
    # })
    # # xi <- sapply(1:n, function(i){
    # #     apply(phi * X[i, ], 2, function(x_i_phi_k){
    # #         fdapace::trapzRcpp(work.grid, x_i_phi_k)
    # #     })
    # # })
    # # xi <- t(xi)


    res <- list(eig.val = lambda,
                eig.ftn = phi,
                fpc.score = xi,
                K = K,
                FVE = cum_FVE[K],
                work.grid = work.grid,
                mean.ftn = mu)
    class(res) <- "uFPCA"

    return(res)
}

# Predict the FPC scores
predict.uFPCA <- function(obj, newdata, K = NULL, type = "fpc.score") {
  if (is.null(K)) {
    K <- obj$K
  }
  if (K > obj$K) {
    stop("Only K is lower than `uFPCA` object is required!")
  }
  
  # FPC scores
  lambda <- obj$eig.val[1:K]
  phi <- obj$eig.ftn[, 1:K]
  xi <- get_fpc_scores(newdata, phi, work.grid)   # C++ function
  # xi <- sapply(1:K, function(k) {
  #   apply(t(newdata) * phi[, k], 2, function(x_i_phi_k) {
  #     fdapace::trapzRcpp(obj$work.grid, x_i_phi_k)
  #   })
  # })
  
  if (type == "fpc.score") {
    res <- xi
  } else if (type == "raw") {
    res <- obj$mean.ftn + t(xi %*% t(phi))
    res <- t(res)
  }
  
  return(res)
}


# Multivariate Functional Principal Component Analysis
#   See "3.2. Estimation of Multivariate FPCA" from
#       Clara Happ & Sonja Greven (2018) "Multivariate Functional Principal Component Analysis for Data Observed on Different (Dimensional) Domains", JASA.
# Inputs:
#    X: A n-m-d array (d-variate functional data; each functional data consists of n curves observed from m timepoints)
#    grid: A vector containing m timepoints
#    K: Number of FPCs
#    FVE: Fraction of variance explained (Default is 0.95)
#    centered: If TRUE, centering of curves is not conducted. (Default is FALSE)
mFPCA <- function(X, grid = NULL, K = NULL, FVE = 0.95, centered = FALSE) {

    n <- dim(X)[1]   # number of curves
    m <- dim(X)[2]   # number of timepoints
    d <- dim(X)[3]   # number of variables

    # Observed grid points
    if (is.null(grid)) {
        grid <- seq(0, 1, length.out = m)
    }

    # Univariate FPCA
    ufpca.obj <- lapply(1:d, function(i){
        uFPCA(X[, , i], grid = grid, K = NULL, FVE = 1, centered = centered)
        # uFPCA(X[, , i], grid = grid, K = K, FVE = FVE, centered = centered)
    })
    mu <- t( sapply(ufpca.obj, function(obj){ obj$mean.ftn }) )   # d-m
    phi_u <- lapply(ufpca.obj, function(obj){ obj$eig.ftn })   # list( m-K_j )
    xi_u <- do.call(cbind,
                    lapply(ufpca.obj, function(obj){ obj$fpc.score }))   # n-(K_1 + ... + K_d)
    M_list <- sapply(phi_u, ncol)   # number of FPCs for each variable
    # max_K <- sum(M_list)   # maximal dimension of FPCs = K_1 + ... + K_d

    # Define K-K matrix Z
    Z <- crossprod(xi_u, xi_u) / n   # max_K-max_K matrix
    # Z <- cov(xi_u, xi_u)
    eig.obj <- eigen(Z, symmetric = TRUE)
    positive_ind <- which(eig.obj$values > 0)   # indices of positive eigenvalues
    lambda <- eig.obj$values[positive_ind]   # Eigenvalues of Z
    C <- eig.obj$vectors[, positive_ind]   # Eigenvectors of Z
    cum_FVE <- cumsum(lambda) / sum(lambda)   # FVE
    if (is.null(K)) {
        K <- which(cum_FVE > FVE)[1]
    }
    FVE <- cum_FVE[K]

    # Truncation
    lambda <- lambda[1:K]
    C <- C[, 1:K]


    # # Normalizing eigenvalues and eigenfunctions
    # grid_size <- grid[2] - grid[1]
    # work.grid <- seq(min(grid), max(grid), length.out = nrow(C))
    # lambda <- lambda * grid_size
    # C <- apply(C, 2, function(x) {
    #     x <- x / sqrt(grid_size)
    #     # x <- x / sqrt(sum(x^2))   # normalize
    #     x <- x / sqrt(fdapace::trapzRcpp(work.grid, x^2))
    #     if ( 0 <= sum(x * work.grid) )
    #         return(x)
    #     else
    #         return(-x)
    # })

    # d-m-K Eigenfunctions and n-K mFPC scores
    psi <- array(0, dim = c(d, m, K))
    # mfpc.score <- array(0, dim = c(d, n, K))
    for (j in 1:d) {
        M_j <- M_list[j]   # number of PCs of jth variable
        cum_idx <- sum(M_list[0:(j-1)])
        idx <- (cum_idx + 1):(cum_idx + M_j)   # index of jth variables for C
        psi[j, , ] <- phi_u[[j]] %*% C[idx, ]
        # mfpc.score[j, , ] <- xi_u[, idx] %*% C[idx, ]
    }
    # mfpc.score <- apply(mfpc.score, c(2, 3), sum)   # n-K matrix
    mfpc.score <- matrix(xi_u %*% C, ncol = K)   # mFPC scores



    # # Normalization (make the variance = lambda)
    # normalizing_const <- 1/sqrt(diag(crossprod(mfpc.score, mfpc.score)))
    # mfpc.score <- mfpc.score %*% diag(sqrt(lambda) * normalizing_const)
    # for (j in 1:d) {
    #     psi[j, , ] <- psi[j, , ] %*% diag(1/(sqrt(lambda) * normalizing_const))
    # }
    # # cov(mfpc.score)
    # # crossprod(mfpc.score, mfpc.score)
    # # crossprod(mfpc.score, mfpc.score) %>% diag
    # # lambda
    # # psi2 <- aperm(psi, c(2,1,3))
    # # dim(psi2)
    # # psi2 <- matrix(psi3, 3*10, 21)
    # # crossprod(psi2[, j], psi2[, j])


    res <- list(eig.val = lambda,
                eig.ftn = psi,
                fpc.score = mfpc.score,
                K = K,
                FVE = FVE,
                mean.ftn = mu)
    class(res) <- "mFPCA"

    return(res)
}




# ##################################################
# ### Comparison between "MFPCA" and my function
# ##################################################
# 
# library(robfpca)
# library(MFPCA)
# library(tidyverse)
# set.seed(1)
# ### simulate data (one-dimensional domains)
# sim <- simMultiFunData(type = "split", argvals = list(seq(0,1,0.01),
#                                                       seq(0,1,0.01),
#                                                       seq(0,1,0.01)),
#                        M = 5, eFunType = "Poly", eValType = "linear", N = 100)
# X <- array(0, c(100, 101, 3))
# X[, , 1] <- sim$simData[[1]]@X
# X[, , 2] <- sim$simData[[2]]@X
# X[, , 3] <- sim$simData[[3]]@X
# 
# 
# # MFPCA based on univariate FPCA
# uFPCA.obj <- MFPCA(sim$simData, M = 5, uniExpansions = list(list(type = "uFPCA"),
#                                                             list(type = "uFPCA"),
#                                                             list(type = "uFPCA")))
# summary(uFPCA.obj)
# sim$simData[[1]]@X
# 
# uFPCA.obj$functions[[1]]@X
# 
# 
# # 내가 짠 MFPCA
# mFPCA.obj <- mFPCA(X, K = 5)
# psi <- mFPCA.obj$eig.ftn
# dim(psi)
# 
# psi2 <- aperm(psi, c(2,1,3))
# dim(psi2)
# psi2 <- matrix(psi2, 3*101, 5)
# crossprod(psi2, psi2)
# 
# fdapace::trapzRcpp(seq(0, 1, length.out = 101), psi[1, , 1]^2)
# fdapace::trapzRcpp(seq(0, 1, length.out = 101), uFPCA.obj$functions[[1]]@X[1, ]^2)
# sum(psi[1, , 1]^2)
# sum(uFPCA.obj$functions[[1]]@X[1, ]^2)
# 
# 
# par(mfrow = c(3, 2))
# matplot(check_eigen_sign(psi[1, , ],
#                          t(uFPCA.obj$functions[[1]]@X)), type = "l")
# matplot(t(uFPCA.obj$functions[[1]]@X), type = "l")
# matplot(check_eigen_sign(psi[2, , ],
#                          t(uFPCA.obj$functions[[2]]@X)), type = "l")
# matplot(t(uFPCA.obj$functions[[2]]@X), type = "l")
# matplot(check_eigen_sign(psi[3, , ],
#                          t(uFPCA.obj$functions[[3]]@X)), type = "l")
# matplot(t(uFPCA.obj$functions[[3]]@X), type = "l")
# 
# 
# uFPCA.obj$scores %>% head()
# mFPCA.obj$fpc.score %>% head()
# 
# colMeans(uFPCA.obj$scores) %>% round(3)
# colMeans(mFPCA.obj$fpc.score) %>% round(3)
# 
# cov(uFPCA.obj$scores) %>% round(3)
# cov(mFPCA.obj$fpc.score) %>% round(3)
# 
# 
# sum( ( X[, , 1] - mFPCA.obj$fpc.score %*% t(mFPCA.obj$eig.ftn[1, , ]) )^2 )
# sum( ( X[, , 1] - uFPCA.obj$scores %*% uFPCA.obj$functions[[1]]@X )^2 )
# 
# d <- 1
# i <- 15
# plot(X[i, , d], type = "l")
# lines((mFPCA.obj$fpc.score %*% t(mFPCA.obj$eig.ftn[d, , ]))[i, ] + mFPCA.obj$mean.ftn[d, ], col = 2)
# lines((uFPCA.obj$scores %*% uFPCA.obj$functions[[d]]@X)[i, ], col = 3) + uFPCA.obj$meanFunction[[d]]@X[1, ]
# 
# par(mfrow = c(3, 3))
# for (d in 1:3) {
#     for (i in c(15, 16, 17)) {
#         plot(X[i, , d], type = "l")
#         lines((mFPCA.obj$fpc.score %*% t(mFPCA.obj$eig.ftn[d, , ]))[i, ] + mFPCA.obj$mean.ftn[d, ], col = 2)
#         lines((uFPCA.obj$scores %*% uFPCA.obj$functions[[d]]@X)[i, ] + uFPCA.obj$meanFunction[[d]]@X[1, ], col = 3)
#     }
# }
# 
# 
# 
# 
# apply(uFPCA.obj$scores, 2, var)
# uFPCA.obj$values %>% round(3)
# mFPCA.obj$eig.val %>% round(3)
# 
# 
# dim(uFPCA.obj$functions[[1]]@X)
# 
# 
# 
# 
# ##################################################
# ### Comparison between "fdapace" and my function
# ##################################################
# 
# library(fdapace)
# X2 <- X[, , 1]
# dim(X2)
# 
# # fdapace
# X2_list <- matrix2list(X2)
# Lt <- X2_list$Lt
# Ly <- X2_list$Ly
# fdapace.obj <- FPCA(Ly, Lt, optns = list(maxK = 3, FVEthreshold = 1))
# fdapace.obj$selectK
# 
# # 내가 짠 univariate FPCA
# uFPCA.obj <- uFPCA(X2, K = 3)
# psi <- uFPCA.obj$eig.ftn
# dim(psi)
# 
# 
# fdapace::trapzRcpp(seq(0, 1, length.out = 101), psi[, 1]^2)
# fdapace::trapzRcpp(seq(0, 1, length.out = 101), fdapace.obj$phi[, 1]^2)
# sum(psi[, 1]^2)
# sum(fdapace.obj$phi[, 1]^2)
# 
# 
# uFPCA.obj$eig.val %>% round(3)
# fdapace.obj$lambda %>% round(3)
# 
# par(mfrow = c(1, 2))
# matplot(check_eigen_sign(psi, fdapace.obj$phi), type = "l")
# matplot(fdapace.obj$phi, type = "l")
# 
# 
# 
# uFPCA.obj$fpc.score %>% head() %>% round(3)
# fdapace.obj$xiEst %>% head() %>% round(3)
# 
# par(mfrow = c(1, 2))
# plot(uFPCA.obj$fpc.score)
# plot(fdapace.obj$xiEst)
# 
# uFPCA.obj$eig.val
# apply(uFPCA.obj$fpc.score, 2, var)
# apply(fdapace.obj$xiEst, 2, var)
# 
# 
