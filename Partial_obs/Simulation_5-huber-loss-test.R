############################################
### Simulation 5
### - Huber loss for estimation on simulation 3
### - Obtain estimations for mean and variances 
###   using huber loss
############################################

# library(dplyr)
library(GA)   # persp plot
library(mvtnorm)
library(fdapace)   # 1, 2
library(mcfda)   # 7
library(synfd)   # 7
library(doParallel)   # parallel computing
library(doRNG)   # set.seed for foreach
library(MASS)   # huber, rlm
source("functions.R")
# source("Kraus(2015)/pred.missfd.R")   # 3
# source("Kraus(2015)/simul.missfd.R")  # 3


load("RData/sim3-1_20210204.RData")
model.cov <- 2   # covariance function setting of the paper (1, 2)
sim <- 1

# Get simulation data
x <- data.list.outlier[[sim]]$x
gr <- data.list.outlier[[sim]]$gr

### Covariance estimation
work.grid <- seq(min(gr), max(gr), length.out = 51)
cov.true <- get_cov_fragm(gr, model = model.cov)   # true covariance
cov.true <- ConvertSupport(fromGrid = gr, toGrid = work.grid, Cov = cov.true)   # transform to the observed grid

## 1. Yao, Müller, and Wang (2005)
## 2. Liu and Müller (2009) - fitted.FPCA()
x.2 <- list(Ly = x$y,
            Lt = x$t)
optns <- list(methodXi = 'CE', dataType = 'Sparse', kernel = 'gauss', verbose = FALSE)
mu.yao.obj <- GetMeanCurve(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)   # get bandwidth of mean estimation
# with(mu.yao.obj, lines(workGrid, mu, lwd = 2))   # draw mean curve
cov.yao.obj <- GetCovSurface(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)
# system.time({ 
#   cov.yao.obj <- GetCovSurface(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)   # 31 sec
# })
cov.yao <- cov.yao.obj$cov
if (length(work.grid) != 51) {
  cov.yao <- ConvertSupport(fromGrid = cov.yao.obj$workGrid, toGrid = work.grid,
                            Cov = cov.yao.obj$cov)   # transform to the observed grid
}


## 7. Lin & Wang (2020)
# estimate mean by local polynomial method
mu.lin.obj <- meanfunc(x.2$Lt, x.2$Ly, method = "PACE", 
                       kernel = "gauss", bw = mu.yao.obj$optns$userBwMu)
cov.lin.obj <- covfunc(x.2$Lt, x.2$Ly, mu = mu.lin.obj, method = "SP")
# system.time({ 
#   cov.lin.obj <- covfunc(x.2$Lt, x.2$Ly, mu = mu.lin.obj, method="SP")
# })
cov.lin <- predict(cov.lin.obj, work.grid)


## Lin % Wang (2020) with Huber loss
# devtools::install_github("statKim/mcfda")
start_time <- Sys.time()
mu.huber.obj <- meanfunc(x.2$Lt, x.2$Ly, method = "HUBER", kernel = "gauss", bw = mu.yao.obj$optns$userBwMu)
cov.huber.obj <- covfunc(x.2$Lt, x.2$Ly, mu = mu.huber.obj, kernel = "gauss", method = "HUBER")
cov.huber <- predict(cov.huber.obj, work.grid)
end_time <- Sys.time()
end_time - start_time


### Variance trajectories
par(mfrow = c(1, 1))
matplot(work.grid, 
        cbind(diag(cov.true),
              diag(cov.yao),
              diag(cov.lin),
              diag(cov.huber)),
        type = "l", lwd = 2,
        xlab = "", ylab = "")
abline(h = 0)
legend("topright",
       c("True","Yao et al. (2005)","Lin & Wang (2020)","Huber"),
       col = 1:4,
       lty = 1:4)

### Covariance surface
par(mfrow = c(2, 2))
persp3D(work.grid, work.grid, cov.true, 
        theta = -40, phi = 30, expand = 1,
        xlab = "s", ylab = "t", zlab = "C(s,t)", main = "True")
persp3D(work.grid, work.grid, cov.yao, 
        theta = -40, phi = 30, expand = 1,
        xlab = "s", ylab = "t", zlab = "C(s,t)", main = "Yao et al. (2005)")
persp3D(work.grid, work.grid, cov.lin, 
        theta = -40, phi = 30, expand = 1,
        xlab = "s", ylab = "t", zlab = "C(s,t)", main = "Lin & Wang (2020)")
persp3D(work.grid, work.grid, cov.huber, 
        theta = -40, phi = 30, expand = 1,
        xlab = "s", ylab = "t", zlab = "C(s,t)", main = "Huber")




############################################
### Huber loss
### - LSE for huber regression
############################################
Ly <- unlist(x.2$Ly)
Lt <- unlist(x.2$Lt)

### mean
mu.loc <- predict(mu.lin.obj, work.grid)
mu.huber.pt <- sapply(gr, function(t) {
  ind <- which(Lt == t)
  return( huber(Ly[ind], k = 1.345)$mu )
})

# local kernel smoothing
bw <- mu.yao.obj$optns$userBwMu   # bandwidth
kernel <- "gauss"
bw_mu_huber <- cv.local_kern_smooth(Lt = Lt, Ly = Ly, newt = NULL, 
                                    kernel = kernel, loss = "Huber", K = 5, parallel = TRUE)
mu.huber <- local_kern_smooth(Lt = Lt, Ly = Ly, newt = gr,
                              bw = bw_mu_huber$selected_bw, kernel = kernel, loss = "Huber", k2 = 1.345)
matplot(work.grid, 
        cbind(mu.loc,
              mu.huber.pt,
              mu.huber),
        type = "l",
        xlab = "", ylab = "")
legend("topleft",
       c("Yao","huber()","rlm()"),
       col = 1:3,
       lty = 1:3)


### Variance
ind <- sapply(Lt, function(t) { match(t, gr) })
length(ind)
length(Ly)
sig2e <- sigma2(x.2$Lt, x.2$Ly)   # the estimate for sigma^2

ss <- (Ly - mu.huber[ind])^2
bw_cov_huber <- cv.local_kern_smooth(Lt = Lt, Ly = ss, newt = NULL, 
                                     kernel = kernel, loss = "Huber", K = 5, parallel = TRUE)
var.huber <- local_kern_smooth(Lt = Lt, Ly = ss, newt = gr, 
                               bw = bw_cov_huber$selected_bw, kernel = kernel, loss = "Huber", k2 = 1.345)
# ss <- (Ly - mu.loc[ind])^2
# var.huber <- local_kern_smooth(Lt = Lt, Ly = ss,
#                                newt = gr, bw = bw, kernel = kernel, loss = "L2")
# var.huber <- ifelse(var.huber - sig2e < 0, 0, var.huber - sig2e)   # 이렇게 하면 다 0 나옴...ㅠㅠ

var.huber.pt <- sapply(gr, function(t) {
  ind <- which(Lt == t)
  sigma2 <- huber(Ly[ind], k = 1.345)$mu
  sigma2 <- ifelse(sigma2 < 0, 0, sigma2)   # set 0 for negative variances
  
  return(sigma2)
})

matplot(work.grid, 
        cbind(diag(cov.true),
              diag(cov.yao),
              diag(cov.lin),
              var.huber,
              var.huber.pt),
        type = "l", lwd = 2,
        xlab = "", ylab = "")
abline(h = 0)
legend("topright",
       c("True","Yao","Lin","rlm()","huber()"),
       col = 1:5,
       lty = 1:5)


library(tidyverse)
library(latex2exp)
df <- data.frame(
  x = rep(gr, 5),
  y = c(diag(cov.true),
        diag(cov.yao),
        diag(cov.lin),
        var.huber,
        var.huber.pt),
  method = factor(rep(c("True","Yao(2005)","Lin(2020)","Lin+rlm()","Lin+huber()"),
                      each = length(gr)),
                  levels = c("True","Yao(2005)","Lin(2020)","Lin+rlm()","Lin+huber()"))
)
ggplot(df, aes(x, y, group = method, color = method)) +
  geom_line(size = 1) +
  labs(x = TeX("$t$"), y = TeX("$\\sigma_X^2(t)$")) +
  geom_hline(yintercept = 0, size = 0.8) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "bottom")

df %>% 
  filter(method %in% c("True","Lin+rlm()","Lin+huber()")) %>% 
  ggplot(aes(x, y, group = method, color = method)) +
  geom_line(size = 1) +
  labs(x = TeX("$t$"), y = TeX("$\\sigma_X^2(t)$")) +
  geom_hline(yintercept = 0, size = 0.8) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "bottom")

