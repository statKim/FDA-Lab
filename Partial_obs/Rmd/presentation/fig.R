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
source("Kraus(2015)/pred.missfd.R")
source("Kraus(2015)/simul.missfd.R")


### sample trajectories of partially obseved functional data
set.seed(10)
n.grid <- 51
x.2 <- sim.kraus(n = 100, out.prop = 0, out.type = 4, len.grid = n.grid)
df <- data.frame(
  id = factor(unlist(sapply(1:length(x.2$Lt), 
                            function(id) { 
                              rep(id, length(x.2$Lt[[id]])) 
                            }) 
  )),
  y = unlist(x.2$Ly),
  t = unlist(x.2$Lt)
)
ggplot(df, aes(t, y, color = id)) +
  geom_line() +
  theme_bw() +
  # ylim(-10, 10) +
  theme(legend.position = "none")

x <- df %>% 
  spread(key = "t", value = "y")
x <- x[, -1] %>% 
  as.matrix
x[1, ]
dim(x)
matplot(t(x), type = "l")

num.NA <- apply(x, 1, function(x){sum(is.na(x))})
ind <- which(num.NA > 0)[c(1,3,5,7,9)]
gr <- seq(0, 1, length.out = n.grid)

postscript("partial_obs.eps",
           width=6, height=6, onefile=FALSE, horizontal=FALSE)
matplot(gr, t(x[ind, ]), 
        type = "l", xlab = "", ylab = "", lwd = 2)
dev.off()


### Dense functional data
postscript("funcional_obs.eps",
           width=6, height=6, onefile=FALSE, horizontal=FALSE)
matplot(gr, t(x.2$x.full[ind, ]), 
        type = "l", xlab = "", ylab = "", lwd = 2)
dev.off()



### functional snippets
set.seed(10)
x <- fun.fragm(n = 30, model = 2, out.prop = 0, out.type = 1)
df <- data.frame(
  id = factor(unlist(sapply(1:length(x$Lt), 
                            function(id) { 
                              rep(id, length(x$Lt[[id]])) 
                            }) 
  )),
  y = unlist(x$Ly),
  t = unlist(x$Lt)
)
ggplot(df, aes(t, y, color = id)) +
  geom_line(size = 1) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(legend.position = "none")
ggsave("fun_snippets.eps", width = 6, height = 6)

# deign plot
obs_ind <- get_design_index(x$Lt)
df_ind <- data.frame(
  x = gr[obs_ind[, 1]],
  y = gr[obs_ind[, 2]]
)
ggplot(df_ind, aes(x, y)) +
  geom_point() +
  labs(x = "", y = "") +
  theme_bw()
ggsave("fun_snippets_design.eps", width = 6, height = 6)




###########################
### Covariance - RMISE
###########################
library(xtable)
load("RData/sim_20210323.RData")
cname <- c("Yao (2005)","Lin (2020)","Huber","WRM")
res.mat <- matrix(NA, 8, 5)
colnames(res.mat) <- c("", cname)
row.names(res.mat) <- rep(c("Outlier X",
                            paste0("Outlier ", 1:3)),
                          each = 2)
res.mat[, 1] <- rep(c("Var", "Cov"), 4)
res.mat2 <- res.mat
num.sim <- 100   # number of simulations

for (i in 1:4) {
  data.list <- sim.obj[[i]]$data.list
  cov.est <- sim.obj[[i]]$cov.est
  
  ### variance
  ise.var <- summary_ise(data.list, cov.est, method = "var")
  ### covariance
  ise.cov <- summary_ise(data.list, cov.est, method = "cov")
  
  res.mat[2*(i-1)+1, 2:5] <- paste0(round(sqrt(rowMeans(ise.var)), 2),
                                    " (",
                                    round(apply(ise.var, 1, sd), 2),
                                    ")")
  
  res.mat[2*(i-1)+2, 2:5] <- paste0(round(sqrt(rowMeans(ise.cov)), 2),
                                    " (",
                                    round(apply(ise.cov, 1, sd), 2),
                                    ")")
  
  
  ### Intrapolation parts (D_0)
  ise.intra <- summary_ise(data.list, cov.est, method = "intra")
  ### Extrapolation parts (S_0 \ D_0)
  ise.extra <- summary_ise(data.list, cov.est, method = "extra")
  
  res.mat2[2*(i-1)+1, 2:5] <- paste0(round(sqrt(rowMeans(ise.intra)), 2),
                                     " (",
                                     round(apply(ise.intra, 1, sd), 2),
                                     ")")
  
  res.mat2[2*(i-1)+2, 2:5] <- paste0(round(sqrt(rowMeans(ise.extra)), 2),
                                     " (",
                                     round(apply(ise.extra, 1, sd), 2),
                                     ")")
}
xtable(res.mat)
xtable(res.mat2)


###########################
### PCA - RMISE
###########################
# Parallel computing setting
ncores <- detectCores() - 3
cl <- makeCluster(ncores)
registerDoParallel(cl)

res <- list()

### Eigen analysis
for (i in 1:4) {
  data.list <- sim.obj[[i]]$data.list
  cov.est <- sim.obj[[i]]$cov.est
  num.sim <- length(cov.est)   # number of simulations
  
  res[[i]] <- sim_eigen_result(cov.est, num.sim, seed = 1000)  
}
stopCluster(cl)


### Calculate ISE for 1~3 eigenfunctions
cname <- c("Yao (2005)","Lin (2020)","Huber","WRM")
res.mat <- matrix(NA, 8, 4)
colnames(res.mat) <- cname
row.names(res.mat) <- rep(c("Outlier X",
                            paste0("Outlier ", 1:3)),
                          each = 2)
K <- 3   # consider first 3 eigenfunctions

for (i in 1:4) {
  pca.est <- res[[i]]
  num.sim <- length(pca.est)
  ise <- matrix(NA, num.sim, 4)
  pve <- matrix(NA, num.sim, 4)
  
  for (sim in 1:num.sim) {
    work.grid <- pca.est[[sim]]$work.grid
    eig.true <- pca.est[[sim]]$true
    eig.yao <- pca.est[[sim]]$yao
    eig.lin <- pca.est[[sim]]$lin
    eig.huber <- pca.est[[sim]]$huber
    eig.wrm <- pca.est[[sim]]$wrm
    
    
    # calculate ISE for k eigenfunctions
    ise_eig <- matrix(NA, K, 4)
    for (k in 1:K) {
      ise_eig[k, ] <- c(
        get_ise(eig.true$phi[, k], eig.yao$phi[, k], work.grid),
        get_ise(eig.true$phi[, k], eig.lin$phi[, k], work.grid),
        get_ise(eig.true$phi[, k], eig.huber$phi[, k], work.grid),
        get_ise(eig.true$phi[, k], eig.wrm$phi[, k], work.grid)
      )
    }
    
    ise[sim, ] <- colSums(ise_eig)
    
    # PVE for first 3 PCs
    pve[sim, ] <- c(
      # eig.true$PVE[3],
      eig.yao$PVE[3],
      eig.lin$PVE[3],
      eig.huber$PVE[3],
      eig.wrm$PVE[3]
    )
  }
  
  res.mat[2*(i-1)+1, ] <- paste0(round(sqrt(colMeans(ise)), 2),
                                 " (",
                                 round(apply(ise, 2, sd), 2),
                                 ")")
  
  res.mat[2*(i-1)+2, ] <- round(colMeans(pve)*100, 2)
  
  # # number of PCs w.r.t. PVE >= 0.99
  # selected_K <- c(
  #   which(eig.true$PVE >= 0.99)[1],
  #   which(eig.yao$PVE >= 0.99)[1],
  #   which(eig.lin$PVE >= 0.99)[1],
  #   which(eig.huber$PVE >= 0.99)[1],
  #   which(eig.wrm$PVE >= 0.99)[1]
  # )
  
  
}
xtable(res.mat)





###########################
### Figures of simulations
###########################
rname <- c("Outlier X",
           paste0("Outlier ", 1:3))
set.seed(123)
ind <- sample(1:100, 4)

### Sample trajectories
fig_traj <- list()   # variance trajectories
for (i in 1:4) {
  x <- sim.obj[[i]]$data.list[[ ind[i] ]]
  x.2 <- x$x
  gr <- x$gr
  
  df <- data.frame(
    id = factor(unlist(sapply(1:length(x.2$Lt), 
                              function(id) { 
                                rep(id, length(x.2$Lt[[id]])) 
                              }) 
    )),
    y = unlist(x.2$Ly),
    t = unlist(x.2$Lt)
  )
  
  # fig_traj[[i]] <- ggplot(df, aes(t, y, color = id)) +
  #   geom_line() +
  #   theme_bw() +
  #   labs(title = paste(rname[i], "-", ind[i], "th")) +
  #   # ylim(-10, 10) +
  #   theme(legend.position = "none")
  fig_traj[[2*(i-1)+1]] <- ggplot(df, aes(t, y, color = id)) +
    geom_line() +
    theme_bw() +
    labs(y = "", title = paste(rname[i], "-", ind[i], "th")) +
    # ylim(-10, 10) +
    theme(legend.position = "none")
  fig_traj[[2*(i-1)+2]] <- ggplot(df, aes(t, y, color = id)) +
    geom_line() +
    theme_bw() +
    labs(y = "", title = "") +
    ylim(-20, 20) +
    theme(legend.position = "none")
}

p <- gridExtra::grid.arrange(grobs = fig_traj, 
                             nrow = 4)
ggsave("samp_traj.eps", p, width = 12, height = 12)


### Estimated variance and eigenfunctions
fig_var <- list()   # variance trajectories
fig_pca <- list()   # pca results
for (i in 1:4) {
  data.list <- sim.obj[[i]]$data.list
  cov.est <- sim.obj[[i]]$cov.est
  
  fig_var <- c(fig_var, 
               ggplot_var(cov.est, ind[i], main = paste(rname[i], "-", ind[i], "th"))) 
  fig_pca <- c(fig_pca,
               ggplot_eig(cov.est, ind[i], main = paste(rname[i], "-", ind[i], "th"))) 
}

p <- gridExtra::grid.arrange(grobs = fig_var, 
                             nrow = 4)
ggsave("var_traj.eps", p, width = 12, height = 12)
p <- gridExtra::grid.arrange(grobs = fig_pca, 
                             nrow = 4)
ggsave("eig_traj.eps", p, width = 12, height = 12)


### Covariance surfaces
# postscript("cov_surf.eps",
#            width = 8, height = 8, onefile = FALSE, horizontal = FALSE)
par(mfrow = c(4, 5),
    mar = c(2, 2, 2, 2))
for (i in 1:4) {
  data.list <- sim.obj[[i]]$data.list
  cov.est <- sim.obj[[i]]$cov.est
  title <- FALSE
  
  if (i == 1) {
    title <- TRUE
  }
  
  plot_cov_surf(cov.est, ind[i], title = title, lab = paste(rname[i], "-", ind[i], "th"))
}
# dev.off()
