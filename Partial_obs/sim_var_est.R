library(GA)   # persp plot
library(mvtnorm)
library(fdapace)   # 1, 2
library(mcfda)   # 7
library(synfd)   # 7
library(doParallel)   # parallel computing
library(doRNG)   # set.seed for foreach
library(MASS)   # huber, rlm
library(latex2exp)
library(tidyverse)
library(robfilter)
source("R/functions.R")


### IRLS test for Huber regression
# test 1
IRLS(Y = stackloss$stack.loss, 
     X = stackloss[, 1:3],
     method = "huber")
rlm(stack.loss ~ .-1, stackloss, method = "M")

# test 2
require(foreign)
require(MASS)
cdata <- read.dta("https://stats.idre.ucla.edu/stat/data/crime.dta")
rlm(crime ~ poverty + single, data = cdata, scale.est = "Huber")
IRLS(Y = cdata$crime, 
     X = cbind(rep(1, nrow(cdata)),
               cdata[, c("poverty","single")]))


### IRLS test for Tukey's bisquare regression(MM-estimation)
# test 1
IRLS(Y = stackloss$stack.loss, 
     X = stackloss[, 1:3],
     method = "bisquare", k = 4.685)
rlm(stack.loss ~ .-1, stackloss, method = "MM")



#####################################
### Load example data
#####################################
sim <- 1
data <- data.list.outlier[[sim]]
test_var <- function(data, kernel = "epanechnikov", bw = 0.1, k2 = 1.345) {
  model.cov <- 2   # covariance function setting of the paper (1, 2)
  
  # Get simulation data
  x <- data$x
  gr <- data$gr
  
  if ( !identical(range(unlist(x$t)), c(0, 1)) ) {
    warning("It can be more time consumming!")
  }
  
  ### Covariance estimation
  work.grid <- seq(min(gr), max(gr), length.out = 51)
  cov.true <- get_cov_fragm(gr, model = model.cov)   # true covariance
  cov.true <- ConvertSupport(fromGrid = gr, toGrid = work.grid, Cov = cov.true)   # transform to the observed grid
  
  ## 1. Yao, Müller, and Wang (2005)
  ## 2. Liu and Müller (2009) - fitted.FPCA()
  x.2 <- list(Ly = x$y,
              Lt = x$t)
  cat("PACE")
  system.time({
    if (kernel == "epanechnikov") {
      kern <- "epan"
    } else {
      kern <- kernel
    }
    optns <- list(methodXi = 'CE', dataType = 'Sparse', kernel = kern, verbose = FALSE,
                  userBwMu = bw, userBwCov = bw)
    mu.yao.obj <- GetMeanCurve(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)   # get bandwidth of mean estimation
    cov.yao.obj <- GetCovSurface(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)
    cov.yao <- cov.yao.obj$cov
    if (length(work.grid) != 51) {
      cov.yao <- ConvertSupport(fromGrid = cov.yao.obj$workGrid, toGrid = work.grid,
                                Cov = cov.yao.obj$cov)   # transform to the observed grid
    }
  })
  
  
  ## 7. Lin & Wang (2020)
  cat("Lin & Wang")
  system.time({
    # estimate mean by local polynomial method
    mu.lin.obj <- meanfunc(x.2$Lt, x.2$Ly, method = "PACE", 
                           kernel = kernel, bw = bw)
    cov.lin.obj <- covfunc(x.2$Lt, x.2$Ly, mu = mu.lin.obj, method = "SP")
    cov.lin <- predict(cov.lin.obj, work.grid)
  })
  
  
  ## Huber loss
  cat("Huber loss")
  system.time({
    # k2 <- 1 / max(abs(unlist(x.2$Ly)))
    mu.huber.obj <- meanfunc.rob(x.2$Lt, x.2$Ly, method = "Huber", kernel = kernel, bw = bw, k2 = k2)
    var.huber.obj <- varfunc.rob(x.2$Lt, x.2$Ly, mu = mu.huber.obj,  
                                 method = "Huber", kernel = kernel, bw = bw, k2 = k2)
    var.huber <- predict(var.huber.obj, work.grid)
  })
  
  
  ## WRM
  cat("WRM")
  mu.wrm.obj <- meanfunc.rob(x.2$Lt, x.2$Ly, method = "WRM", kernel = kernel, bw = bw)
  var.wrm.obj <- varfunc.rob(x.2$Lt, x.2$Ly, mu = mu.wrm.obj,  
                             method = "WRM", kernel = kernel, bw = bw)
  var.wrm <- predict(var.wrm.obj, work.grid)
  
  end_time <- Sys.time()
  print(end_time - start_time)
  
  ## Visualization
  df <- data.frame(
    x = rep(work.grid, 5),
    y = c(diag(cov.true),
          diag(cov.yao),
          diag(cov.lin),
          var.huber,
          # var.huber + var.huber.obj$sig2 - cov.lin.obj$sig2e,
          var.wrm),
    method = factor(rep(c("True","Yao(2005)","Lin(2020)","Huber","WRM"),
                        each = length(work.grid)),
                    levels = c("True","Yao(2005)","Lin(2020)","Huber","WRM"))
  )
  df$y[df$y < 0] <- 0
  fig <- ggplot(df, aes(x, y, group = method, color = method)) +
    geom_line(size = 1) +
    labs(x = TeX("$t$"), y = TeX("$\\sigma_X^2(t)$")) +
    geom_hline(yintercept = 0, size = 0.8) +
    ylim(0, 5) +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.position = "bottom")
  
  return(fig)
}


load("RData/sim3-2_20210204.RData")
sim <- 4
system.time({
  test_var(data.list.outlier[[sim]], kernel = "gauss", bw = 0.1, k2 = 1.345)  
})
system.time({
  test_var(data.list.outlier[[sim]], kernel = "gauss", bw = 0.1, k2 = 0.001)  
})



#####################################
### Load example data
#####################################
load("RData/sim3-2_20210204.RData")
# sim <- 20
sim <- 1
model.cov <- 2   # covariance function setting of the paper (1, 2)
kernel <- "gauss"
bw <- 0.1

# Get simulation data
x <- data.list.outlier[[sim]]$x
gr <- data.list.outlier[[sim]]$gr
range(unlist(x$t))

### Covariance estimation
work.grid <- seq(min(gr), max(gr), length.out = 51)
cov.true <- get_cov_fragm(gr, model = model.cov)   # true covariance
cov.true <- ConvertSupport(fromGrid = gr, toGrid = work.grid, Cov = cov.true)   # transform to the observed grid

x.2 <- list(Ly = x$y,
            Lt = x$t)

## 1. Yao, Müller, and Wang (2005)
## 2. Liu and Müller (2009) - fitted.FPCA()
system.time({
  if (kernel == "epanechnikov") {
    kern <- "epan"
  } else {
    kern <- kernel
  }
  optns <- list(methodXi = 'CE', dataType = 'Sparse', kernel = kern, verbose = FALSE,
                kFoldMuCov = 5, methodBwMu = "CV", methodBwCov = "CV", useBinnedCov = FALSE)
                # userBwMu = bw, userBwCov = bw)
  mu.yao.obj <- GetMeanCurve(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)   # get bandwidth of mean estimation
  cov.yao.obj <- GetCovSurface(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)
  cov.yao <- cov.yao.obj$cov
  if (length(work.grid) != 51) {
    cov.yao <- ConvertSupport(fromGrid = cov.yao.obj$workGrid, toGrid = work.grid,
                              Cov = cov.yao.obj$cov)   # transform to the observed grid
  }
})


## 7. Lin & Wang (2020)
system.time({
  # 5-fold CV (It took very long time when we use CV option in mcfda package.)
  cv.mu.lin.obj <- meanfunc.rob(x.2$Lt, x.2$Ly, method = "L2", kernel = kernel,
                                cv_bw_loss = "L2", ncores = 9,
                                bw = NULL)
  cv.var.lin.obj <- varfunc.rob(x.2$Lt, x.2$Ly, mu = cv.mu.lin.obj, kernel = kernel,
                                method = "L2",  cv_bw_loss = "L2", ncores = 9,
                                bw = NULL)
  # estimate mean, variance, covariance
  mu.lin.obj <- meanfunc(x.2$Lt, x.2$Ly, method = "PACE", kernel = kernel,
                         bw = cv.mu.lin.obj$bw)   # It occurs error or very slow.
  var.lin.obj <- varfunc(x.2$Lt, x.2$Ly, method = "PACE", kernel = kernel,
                         mu = mu.lin.obj, bw = cv.var.lin.obj$obj$bw)
  cov.lin.obj <- covfunc(x.2$Lt, x.2$Ly, method = "SP",
                         mu = mu.lin.obj, sig2x = var.lin.obj)
  cov.lin <- predict(cov.lin.obj, work.grid)
})


# trajectories of generated data
df <- data.frame(id = factor(unlist(sapply(1:length(x.2$Lt), 
                                           function(id) { 
                                             rep(id, length(x.2$Lt[[id]])) 
                                           }) 
                                    )),
                 t = unlist(x.2$Lt),
                 y = unlist(x.2$Ly))
ggplot(df, aes(t, y, color = id)) +
  geom_line(size = 1) +
  theme_bw() +
  # ylim(-10, 10) +
  theme(legend.position = "none")



#####################################
### Find optimal delta for Huber loss
### - k <- 1 / max(abs(y))
#####################################
# k_cand <- c(1 / max(abs(unlist(x.2$Ly))), 1e-4, 0.001, 0.01, 0.1, 0.5, 0.7, 1, 1.345)
k_cand <- union(1.345,
                10^seq(-3, 2, length.out = 10) * (1 - 0)/3)
var_est <- matrix(0, length(gr), length(k_cand)+1)
var_est[, 1] <- diag(cov.true)
system.time({
  for (i in 1:length(k_cand)) {
    print(k_cand[i])
    
    bw_fixed <- 0.1
    kernel <- "gauss"
    mu.huber.obj <- meanfunc.rob(x.2$Lt, x.2$Ly, method = "huber", kernel = kernel, 
                                 bw = bw_fixed, k2 = k_cand[i])
    var.huber.obj <- varfunc.rob(x.2$Lt, x.2$Ly, mu = mu.huber.obj,  
                                 method = "huber", kernel = kernel, bw = bw_fixed, k2 = k_cand[i])
    var_est[, i+1] <- predict(var.huber.obj, gr)
  }
}) 
df_delta <- data.frame(t = rep(gr, length(k_cand)+3),
                       y = c(as.numeric(var_est),
                             diag(cov.yao),
                             diag(cov.lin)),
                       k_huber = rep(factor(c("True", round(k_cand, 4), "Yao", "Lin"), 
                                            levels = c("True", round(k_cand, 4), "Yao", "Lin")), 
                                     each = length(gr)))
p1 <- ggplot(df_delta, 
             aes(t, y, color = k_huber, linetype = k_huber)) +
  geom_line(size = 1) +
  scale_linetype_manual(breaks = c("True", round(k_cand, 4), "Yao", "Lin"),
                        values = c("solid", rep("dashed", 11), rep("solid", 2))) +
  theme_bw() +
  theme(legend.position = c(0.8, 0.7),
        legend.title = element_blank())
p2 <- p1 + 
  ylim(0, 5) +
  theme(legend.position = "none")
gridExtra::grid.arrange(p1, p2, 
                        nrow = 1)

# Compare CV results and true optimal delta
var.huber.obj <- varfunc.rob(x.2$Lt, x.2$Ly, mu = mu.huber.obj,  
                             method = "huber", kernel = kernel, bw = bw_fixed, k2 = NULL)
ise_var <- df_delta %>% 
  group_by(k_huber) %>% 
  summarise(ise = get_ise(y, diag(cov.true), gr))

var.huber.obj$obj$k2
k_cand[which.min(ise_var$ise[-1])]



# # par(mfrow = c(1, 2))
# # matplot(var_est, type = "l")
# # legend("topright",
# #        c("True", "Naive", round(k_cand, 5)),
# #        col = 1:ncol(var_est),
# #        lty = 1:ncol(var_est))
# # matplot(var_est, type = "l",
# #         ylim = c(0, 5))
# # # matplot(var_est[, 1:4], type= "l")
# # lines(diag(cov.yao), col = 2)
# # lines(diag(cov.lin), col = 3)
# 
# ise_var <- df_delta %>% 
#   group_by(k_huber) %>% 
#   summarise(ise = get_ise(y, diag(cov.true), gr))
# # ise_var[-1, ] %>% 
# #   filter(ise == min(ise)) %>% 
# #   select(k_huber)
# k_cand[which.min(ise_var$ise[-1])]
#   
# # ise_var <- apply(var_est[, -1], 2, function(cov){ get_ise(cov, diag(cov.true), gr) })
# # k_cand[which.min(ise_var)]


#####################################
### bandwidth test for Huber loss
#####################################
bw_cand <- 10^seq(-2, 0, length.out = 10) * (1 - 0)/3
var_est <- matrix(0, length(gr), length(bw_cand)+1)
var_est[, 1] <- diag(cov.true)
system.time({
  for (i in 1:length(bw_cand)) {
    print(bw_cand[i])
    
    # k <- 1 / max(abs(unlist(x.2$Ly)))
    k <- k_cand[which.min(ise_var$ise[-1])]   # optimal delta for Huber loss from above procedure
    kernel <- "gauss"
    tryCatch({
      mu.huber.obj <- meanfunc.rob(x.2$Lt, x.2$Ly, method = "Huber", kernel = kernel, bw = bw_cand[i], k2 = k)
      var.huber.obj <- varfunc.rob(x.2$Lt, x.2$Ly, mu = mu.huber.obj,  
                                   method = "Huber", kernel = kernel, bw = bw_cand[i], k2 = k)
      var_est[, i+1] <- predict(var.huber.obj, gr)
    }, error = function(e){
      print(e)
    })
  }
}) 
df_bw <- data.frame(t = rep(gr, length(bw_cand)+1),
                    y = as.numeric(var_est),
                    bw = rep(factor(c("True",round(bw_cand, 3)),
                                    levels = c("True",round(bw_cand, 3))), 
                             each = length(gr)))
ggplot(df_bw, aes(t, y, color = bw)) +
  geom_line(size = 1) +
  scale_linetype_manual(breaks = c("True", round(bw_cand, 3)),
                        values = c("solid", rep("dashed", 10))) +
  theme_bw() +
  theme(legend.position = c(0.5, 0.75),
        legend.title = element_blank())

ise_var <- apply(var_est[, -1], 2, function(cov){ get_ise(cov, diag(cov.true), gr) })
bw_cand[which.min(ise_var)]


# Compare CV results and true optimal bandwidth
var.huber.obj <- varfunc.rob(x.2$Lt, x.2$Ly, mu = mu.huber.obj,  
                             method = "huber", kernel = kernel, bw = NULL, k2 = k)
ise_var <- apply(var_est[, -1], 2, function(cov){ get_ise(cov, diag(cov.true), gr) })

var.huber.obj$obj$bw
bw_cand[which.min(ise_var)]

cbind(bw_cand, ise_var)


#####################################
### bandwidth test for WRM
#####################################
bw_cand <- 10^seq(-2, 0, length.out = 10) * (1 - 0)/3
var_est <- matrix(0, length(gr), length(bw_cand)+1)
var_est[, 1] <- diag(cov.true)
system.time({
  for (i in 1:length(bw_cand)) {
    print(bw_cand[i])
    
    kernel <- "gauss"
    mu.wrm.obj <- meanfunc.rob(x.2$Lt, x.2$Ly, method = "WRM", kernel = kernel, bw = bw_cand[i])
    var.wrm.obj <- varfunc.rob(x.2$Lt, x.2$Ly, mu = mu.wrm.obj,  
                                 method = "WRM", kernel = kernel, bw = bw_cand[i])
    var_est[, i+1] <- predict(var.wrm.obj, gr)
  }
}) 
df_bw_wrm <- data.frame(t = rep(gr, length(bw_cand)+1),
                        y = as.numeric(var_est),
                        bw = rep(factor(c("True", round(bw_cand, 3)),
                                        levels = c("True", round(bw_cand, 3))),
                                 each = length(gr)))
ggplot(df_bw_wrm, aes(t, y, color = bw, linetype = bw)) +
  geom_line(size = 1) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank())

# save(list = c("k_cand","df_delta","bw_cand","df_bw","df_bw_wrm"),
#      file = "RData/20210310_fig_1st.RData")
# save(list = c("k_cand","df_delta","bw_cand","df_bw","df_bw_wrm"),
#      file = "RData/20210310_fig_20th.RData")


#####################################
### CV test - bandwidth
#####################################
# 아 왜 오래 걸려 ㅡㅡ 첫번째 부분에서 에러뜨긴 하는듯... k2 값이 작아서 오래 걸리는건가??
# C++로 짠 이후에 매우 빠르게 개선됨!!
source("R/functions.R")
system.time({
  bw_cv_obj <- bw.local_kern_smooth(Lt = x.2$Lt,
                                    Ly = x.2$Ly, 
                                    method = "HUBER",
                                    kernel = "gauss", 
                                    k2 = 0.001,
                                    cv_loss = "L1",
                                    K = 5, 
                                    ncores = 1)
})
bw_cv_obj

# WRM
system.time({
  bw_cv_obj <- bw.local_kern_smooth(Lt = x.2$Lt,
                                    Ly = x.2$Ly, 
                                    method = "WRM",
                                    kernel = "gauss", 
                                    cv_loss = "L1",
                                    K = 5,
                                    ncores = 9)
})
bw_cv_obj
