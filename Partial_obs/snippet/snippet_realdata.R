######################################
### Real data - CD4 counts
######################################
library(tidyverse)
library(catdata)
data(aids)
head(aids)

length(unique(aids$person))   # 369 individuals
aids$cd4 <- log(aids$cd4)   # log-trasnformation

# Remove individuals containing only 1 observation
id <- aids %>% 
  group_by(person) %>% 
  summarise(n = n()) %>% 
  filter(n > 1) %>% 
  select(person) %>% 
  unlist() %>% 
  as.numeric()
aids <- aids[which(aids$person %in% id), ]
length(unique(aids$person))   # 364 individuals

# Normalize time points into (0, 1)
aids$time <- (aids$time - min(aids$time)) / (max(aids$time) - min(aids$time))

# aids <- aids[aids$time >= 0, ]
Ly <- split(aids$cd4, aids$person)
Lt <- split(aids$time, aids$person)

# Design plot and scatter
par(mfrow = c(1, 2))
CreateDesignPlot(Lt)
plot(1,
     type = "n",
     xlim = range(unlist(Lt)), ylim = range(unlist(Ly)),
     xlab = "", ylab = "")
for (i in 1:length(id)) {
  lines(Lt[[i]], Ly[[i]], lty = i, col = i, type = "o")
}

# # Only consider postiive time grids
# id <- aids %>% 
#   group_by(person) %>% 
#   summarize(n = n()) %>% 
#   filter(n > 1) %>% 
#   select(person) %>% 
#   unlist() %>% 
#   as.numeric()
# Lt <- list()
# Ly <- list()
# for (i in 1:length(id)) {
#   idx <- which(aids$person == id[i])
#   Ly[[i]] <- aids$cd4[idx]
#   Lt[[i]] <- aids$time[idx]
# }

# Individual trajectories
par(mfrow = c(2, 3))
for (j in 1:6) {
  plot(1,
       type = "n",
       xlim = range(unlist(Lt)), ylim = range(unlist(Ly)),
       xlab = "", ylab = "")
  for (i in ((j-1)*10 + 1):(j*10)) {
    lines(Lt[[i]], Ly[[i]], lty = i, col = i, type = "o")
  }
}

# Combine dataset
x.2 <- list(Lt = Lt,
            Ly = Ly)



#############################
### Covariance estimation
#############################
library(fdapace)
library(mcfda)
library(tidyverse)
library(doParallel)   # parallel computing
library(doRNG)   # set.seed for foreach
library(robfpca)
source("Boente_cov.R")

data_type <- "snippet"   # type of functional data
# kernel <- "epanechnikov"   # kernel function for local smoothing
kernel <- "gauss"   # kernel function for local smoothing
# bw_boente <- 0.2   # bandwidth for Boente(2020) - Error occurs for small bw
n_cores <- 12   # number of threads for parallel computing
pve <- 0.95   # Not used if K is given

seed <- 100
t.range <- range(unlist(x.2$Lt))
work.grid <- seq(t.range[1], t.range[2], length.out = 51)


### Lin & Wang (2020)
start_time <- Sys.time()
registerDoRNG(seed)
tryCatch({
  # # 5-fold CV (It took very long time when we use CV option in mcfda package.)
  # cv.mu.lin.obj <- meanfunc.rob(x.2$Lt, x.2$Ly, method = "L2", kernel = kernel,
  #                               cv_bw_loss = "L2", ncores = n_cores,
  #                               cv = TRUE)
  # cv.var.lin.obj <- varfunc.rob(x.2$Lt, x.2$Ly, mu = cv.mu.lin.obj, kernel = kernel,
  #                               method = "L2",  cv_bw_loss = "L2", ncores = n_cores,
  #                               cv = TRUE)
  # estimate mean, variance, covariance
  mu.lin.obj <- meanfunc(x.2$Lt, x.2$Ly, method = "PACE", kernel = kernel,
                         bw = 0.2)
                         # bw = cv.mu.lin.obj$bw)   # It occurs error or very slow.
  var.lin.obj <- varfunc(x.2$Lt, x.2$Ly, method = "PACE", kernel = kernel,
                         mu = mu.lin.obj, 
                         bw = 0.2)
                         # bw = cv.var.lin.obj$obj$bw)
  cov.lin.obj <- covfunc(x.2$Lt, x.2$Ly, method = "SP",
                         mu = mu.lin.obj, sig2x = var.lin.obj)
  mu.lin <- predict(mu.lin.obj, work.grid)
  cov.lin <- predict(cov.lin.obj, work.grid)
  
  # # too large noise variance estimate...
  # if (0 %in% diag(cov.lin)) {
  #   stop("0 variance occur!")
  # }
  
  # # estimate mean by local polynomial method
  # mu.lin.obj <- meanfunc(x$Lt, x$Ly, method = "PACE", 
  #                        kernel = kernel, bw = bw)
  # cov.lin.obj <- covfunc(x$Lt, x$Ly, mu = mu.lin.obj, method = "SP")
}, error = function(e) { 
  print("Lin cov error")
  print(e)
})
end_time <- Sys.time()
print(paste0("Lin & Wang : ", 
             round(difftime(end_time, start_time, units = "secs"), 3),
             " secs"))


### Robust + Lin
start_time <- Sys.time()
registerDoRNG(seed)
tryCatch({
  ## Huber loss
  mu.huber.obj <- meanfunc.rob(x.2$Lt, x.2$Ly, 
                               method = "huber", kernel = kernel, 
                               cv = F)
                               # cv = TRUE, ncores = n_cores)
  
  Lt <- x.2$Lt
  Ly <- x.2$Ly
  mu <- mu.huber.obj
  n <- length(Lt)
  # calculate sum of squares
  gr <- sort(unique(unlist(Lt)))
  mu_hat <- predict(mu, gr)
  ss <- lapply(1:n, function(i) {
    ind <- match(Lt[[i]], gr)
    if (length(ind) == length(Lt[[i]])) {
      return( (Ly[[i]] - mu_hat[ind])^2 )
    } else {
      mui <- predict(mu, Lt[[i]])
      return( (Ly[[i]] - mui)^2 )
    }
  })
  
  # obtain noise variance
  h.sig2 <- select.sig2.rob.bw(Lt, Ly, ss)
  sig2 <- sigma2.rob(Lt, Ly, h = h.sig2)
  
  
  ### bw
  t.min <- min(unlist(Lt))
  t.max <- max(unlist(Lt))
  
  delta <- max(sapply(Lt, function(ts){ max(ts)-min(ts) }))
  m <- mean(sapply(Lt, length))
  M <- n * m^2
  #h <- delta * (M^(-1/5)) / 7
  
  vn <- sqrt(mean(unlist(ss)))
  
  ss <- unlist(ss)
  cut_off <- quantile(ss, probs = 0.75)
  t <- unlist(Lt)
  idx <- order(t)
  vn <- sqrt(trapzRcpp(t[idx][ss[idx] < cut_off], ss[idx][ss[idx] < cut_off]))
  zeta <- meanfunc.rob(Lt = Lt,
                       Ly = ss,
                       newt = work.grid,
                       method = "huber",
                       kernel = kernel,
                       bw = 0.2)
  sqrt(mean(zeta$fitted))
  # ss <- unlist(ss)
  # cut_off <- quantile(ss, probs = 0.75)
  # vn <- sqrt(mean(ss[ss < cut_off]))
  
  h <- 0.29 * delta * vn * (M^(-1/5))
  
  max.it <- 1000
  it <- 0
  while (it < max.it) {  # h0 two small
    it <- it + 1
    # print(paste(it, ":"))
    
    cnt <- sapply(1:n, function(i) {
      tobs <- Lt[[i]]
      # y <- Ly[[i]]
      m <- length(tobs)
      v1 <- 0
      
      if (m < 2) {
        return(0)
      }
      for (j in 1:m) {
        for (k in 1:m) {
          if ( (k!=j) && (abs(tobs[j]-tobs[k]) < h) ) {
            v1 <- v1 + 1
          }
        }
      }
      
      return(v1)
    })
    
    cnt <- sum(cnt)
    if (cnt >= min(50, 0.1*n*m*(m-1))) {
      break
    } else {
      h <- h*1.01
    }
  }
  
  h
  
  
  
  
  
  var.huber.obj <- varfunc.rob(x.2$Lt, x.2$Ly,
                               method = "huber", kernel = kernel,
                               mu = mu.huber.obj, 
                               cv = F)
                               # cv = TRUE, ncores = n_cores)
  cov.huber.obj <- covfunc.rob(x.2$Lt, x.2$Ly, 
                               method = "huber", 
                               kernel = kernel, 
                               mu = mu.huber.obj,
                               sig2x = var.huber.obj,
                               # cv = TRUE, 
                               ncores = n_cores)
  print(c(mu.huber.obj$bw,
          var.huber.obj$obj$bw,
          cov.huber.obj$theta))
  mu.huber <- predict(mu.huber.obj, work.grid)
  cov.huber <- predict(cov.huber.obj, work.grid)
  
  ## Bisquare loss
  mu.bisquare.obj <- meanfunc.rob(x.2$Lt, x.2$Ly, 
                                  method = "bisquare", kernel = kernel, 
                                  cv = TRUE, ncores = n_cores)
  var.bisquare.obj <- varfunc.rob(x.2$Lt, x.2$Ly,
                                  method = "bisquare", kernel = kernel,
                                  mu = mu.bisquare.obj, 
                                  cv = TRUE, ncores = n_cores)
  cov.bisquare.obj <- covfunc.rob(x.2$Lt, x.2$Ly, 
                                  method = "bisquare", 
                                  kernel = kernel, 
                                  mu = mu.bisquare.obj,
                                  sig2x = var.bisquare.obj,
                                  cv = TRUE, 
                                  ncores = n_cores)
  print(c(mu.bisquare.obj$bw,
          var.bisquare.obj$obj$bw,
          cov.bisquare.obj$theta))
  mu.bisquare <- predict(mu.bisquare.obj, work.grid)
  cov.bisquare <- predict(cov.bisquare.obj, work.grid)
}, error = function(e) { 
  print("Huber cov error")
  print(e)
})
end_time <- Sys.time()
print(paste0("Proposed : ", 
             round(difftime(end_time, start_time, units = "secs"), 3),
             " secs"))



### Boente et al. (2020)
start_time <- Sys.time()
registerDoRNG(seed)
tryCatch({
  bw_boente <- 0.2
  cov.boente.obj <- cov_boente(x.2, bw.mu = bw_boente, bw.cov = bw_boente)
  # cov.boente.obj <- cov_boente(x.2, cv = TRUE, kern = kernel, seed = seed)   # 5-fold CV
  mu.boente <- cov.boente.obj$mu
  cov.boente <- cov.boente.obj$cov
  # noise var from source code of sparseFPCA package
  noise_boente <- cov.boente.obj$noise_var
}, error = function(e) {
  print("Boente (2020) cov error")
  print(e)
})
end_time <- Sys.time()
print(paste0("Boente (2020) : ", 
             round(difftime(end_time, start_time, units = "secs"), 3),
             " secs"))



### Yao, Müller, and Wang (2005)
start_time <- Sys.time()
registerDoRNG(seed)
kern <- ifelse(kernel == "epanechnikov", "epan", kernel)
optns <- list(methodXi = "CE", dataType = "Sparse", kernel = kern, verbose = FALSE,
              userBwMu = 0.2, userBwCov = 0.2)
              # kFoldMuCov = 5, methodBwMu = "CV", methodBwCov = "CV", useBinnedCov = FALSE)
tryCatch({
  mu.yao.obj <- GetMeanCurve(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)
  cov.yao.obj <- GetCovSurface(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)
  mu.yao <- mu.yao.obj$mu
  cov.yao <- cov.yao.obj$cov
  if (length(work.grid) != 51) {
    mu.yao <- ConvertSupport(fromGrid = cov.yao.obj$workGrid, toGrid = work.grid,
                             mu = mu.yao.obj$mu)
    cov.yao <- ConvertSupport(fromGrid = cov.yao.obj$workGrid, toGrid = work.grid,
                              Cov = cov.yao.obj$cov)
  }
}, error = function(e) { 
  print("Yao cov error")
  print(e)
})
end_time <- Sys.time()
print(paste0("Yao et al. : ", 
             round(difftime(end_time, start_time, units = "secs"), 3),
             " secs"))

# Covariance surfaces
gr <- work.grid
par(mfrow = c(2, 3))
GA::persp3D(gr, gr, cov.yao,
            theta = -70, phi = 30, expand = 1)
GA::persp3D(gr, gr, cov.lin,
            theta = -70, phi = 30, expand = 1)
GA::persp3D(gr, gr, cov.boente,
            theta = -70, phi = 30, expand = 1)
GA::persp3D(gr, gr, cov.huber,
            theta = -70, phi = 30, expand = 1)
GA::persp3D(gr, gr, cov.bisquare,
            theta = -70, phi = 30, expand = 1)
par(mfrow = c(1, 1))


#######################################
### Principal component analysis
#######################################
pve <- NULL
K <- 2
# Yao
pca.yao.obj <- funPCA(x.2$Lt, x.2$Ly, 
                      mu.yao, cov.yao, sig2 = cov.yao.obj$sigma2, 
                      work.grid, PVE = pve, K = K)
# Lin
pca.lin.obj <- funPCA(x.2$Lt, x.2$Ly, 
                      mu.lin, cov.lin, sig2 = cov.lin.obj$sig2e, 
                      work.grid, PVE = pve, K = K)
# Boente
pca.boente.obj <- funPCA(x.2$Lt, x.2$Ly, 
                         mu.boente, cov.boente, sig2 = noise_boente, 
                         work.grid, PVE = pve, K = K)
# Huber
pca.huber.obj <- funPCA(x.2$Lt, x.2$Ly, 
                        mu.huber, cov.huber, sig2 = cov.huber.obj$sig2e, 
                        work.grid, PVE = pve, K = K)
# Bisquare
pca.bisquare.obj <- funPCA(x.2$Lt, x.2$Ly, 
                           mu.bisquare, cov.bisquare, sig2 = cov.bisquare.obj$sig2e, 
                           work.grid, PVE = pve, K = K)

### Save PCA object
pca.est <- list(seed = seed,
                x.2 = x.2,
                work.grid = work.grid,
                pca.obj = list(pca.yao.obj = pca.yao.obj,
                               pca.lin.obj = pca.lin.obj,
                               pca.boente.obj = pca.boente.obj,
                               pca.huber.obj = pca.huber.obj,
                               pca.bisquare.obj = pca.bisquare.obj))

# Eigenfunction
par(mfrow = c(1, 2))
for (k in 1:K) {
  matplot(work.grid, 
          cbind(pca.boente.obj$eig.fun[, k],
                check_eigen_sign(pca.yao.obj$eig.fun[, k],
                                 pca.boente.obj$eig.fun[, k]),
                check_eigen_sign(pca.lin.obj$eig.fun[, k],
                                 pca.boente.obj$eig.fun[, k])),
          type = "l",
          xlab = "", ylab = "")
  if (k == 1) {
    legend("topright", 
           lty = 1:3, 
           col = 1:3,
           legend = c("Boente","Yao","Lin"))
  }
}
































######################################
### Other datasets
######################################

data <- read.table("/Users/hyunsung/Desktop/fat.txt")
colnames(data) <- c("ID","age","ageatmen","timetomen","bfat")
head(data)

Ly <- split(data$bfat, data$ID)
Lt <- split(data$timetomen, data$ID)
CreateDesignPlot(Lt)




data <- read.table("/Users/hyunsung/Desktop/Homework/Longitudinal_Data_Analysis/HW2/fev1.txt", header = T)
head(data)

Ly <- split(data$logFev1, data$ID)
Lt <- split(data$Age, data$ID)

idx <- which(sapply(Lt, length) <= 7 & sapply(Lt, length) > 1)
Lt <- Lt[idx]
Ly <- Ly[idx]
CreateDesignPlot(Lt)

plot(1,
     type = "n",
     xlim = range(unlist(Lt)), ylim = range(unlist(Ly)),
     xlab = "", ylab = "")
for (i in 1:length(Ly)) {
  lines(Lt[[i]], Ly[[i]], lty = i, col = i, type = "o")
}

par(mfrow = c(2, 3))
for (j in 1:6) {
  plot(1,
       type = "n",
       xlim = range(unlist(Lt)), ylim = range(unlist(Ly)),
       xlab = "", ylab = "")
  for (i in ((j-1)*10 + 1):(j*10)) {
    lines(Lt[[i]], Ly[[i]], lty = i, col = i, type = "o")
  }
}