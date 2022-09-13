#####################################################################
### Medfly data in Chiou and Muller (2014)
### - Refer paper : Chiou and Muller (2014), 
###                 Linear manifold modelling of multivariate functional data
#####################################################################
library(tidyverse)
library(RFPCA)    # RFPCA and MFPCA
library(mclust)   # clustering measure
library(Funclustering)   # funclust (Currently, it is not supported by cran.)
library(funHDDC)   # funHDDC
library(gmfd)   # gmfd
source("functions.R")

### Load data
data <- read.table("~/GoogleDrive/Lab/KHS/manifold_clust/real_data/fly_log_130521.txt", header = T)
head(data)

id <- unique(data$id)
length(id)   # 62 individuals

data %>% 
    group_by(id) %>% 
    summarise(m = max(day)) %>% 
    as.data.frame()

# Only use times less than 37 weeks
data <- data %>%
    filter(day <= 37)


Lt <- list()
Ly <- list()
for (i in 1:length(id)) {
    ind <- which(data$id == id[i])
    Lt[[i]] <- data$day[ind]
    Ly[[i]] <- as.matrix(data[ind, 3:6])
}
Lt
sapply(Ly, class)

data[, -(1:2)] %>% 
    rowSums()


ggplot(data,
       aes(x = day, y = resting, group = id, color = factor(id))) +
    geom_line(size = 0.3) +
    theme_bw() +
    theme(legend.position = "none")

par(mfrow = c(4, 4))
Y <- matrix(0, 62, 51)
for (i in 1:length(id)){
    plot(Lt[[i]], Ly[[i]][, 1], type = "l", main = i, ylim = c(-2, 1)) 
    lines(ksmooth(Lt[[i]], Ly[[i]][, 1],
                  kernel = "normal",
                  n.points = 51, 
                  bandwidth = 5),
          col = 2)
    Y[i, ] <- ksmooth(Lt[[i]], Ly[[i]][, 1],
                      kernel = "normal",
                      n.points = 51, 
                      bandwidth = max(diff(Lt[[i]])))$y
}
plot(colMeans(Y), type = "l", main = "Mean function")




### using raw count data
data2 <- exp(data[, -(1:2)]) - 0.5
data2[which(data2 < 1e-6, arr.ind = T)] <- 0

rowSums(data2)
rowSums(data[, -(1:2)])
data2 <- cbind(data[, 1:2], data2)

Lt <- list()
Ly <- list()
for (i in 1:length(id)) {
    ind <- which(data2$id == id[i])
    Lt[[i]] <- data2$day[ind]
    Ly[[i]] <- as.matrix(data2[ind, 3:6])
}
Lt
sapply(Ly, class)

data2[, -(1:2)] %>% 
    rowSums()

# # Filter the max of timepoints == 40
# id <- id[sapply(Lt, max) == 40]
# Ly <- Ly[sapply(Lt, max) == 40]
# Lt <- Lt[sapply(Lt, max) == 40]

### Pre-smoothing for regular grids using local linear smoother
n <- length(id)
num_grid <- 101   # number of timepoints
Ly <- lapply(1:n, function(i) {
    y <- Ly[[i]]
    t <- Lt[[i]]
    # bw <- max(diff(t))/2   # very small bandwidth
    bw <- 5

    # kernel smoothing with 51 regular grids
    y <- apply(y, 2, function(col) {
        stats::ksmooth(x = t,
                       y = col,
                       kernel = "normal",
                       bandwidth = bw,
                       n.points = num_grid)$y
    })
    
    # make spherical data
    apply(y, 1, function(row){ sqrt(row / sum(row)) })
})
Lt <- rep(list(seq(1, 37, length.out = num_grid)), n)

apply(Ly[[1]], 2, function(x){ sum(x^2) })   # check that it is the unit sphere


### Example of data
y_name <- c("flying","feeding","walking","resting")
fig_list <- list()
set.seed(100)
for (i in sample(1:62, 3)) {
    df <- data.frame(
        time = Lt[[i]],
        t(Ly[[i]])
    ) %>% 
        gather(var_name, val, -time)
    # head(df)
    p <- ggplot(
        data = df, 
        aes(
            x = time,
            y = val,
            group = var_name,
            color = var_name
        )
    ) +
        geom_line() +
        labs(x = "Days", y = "") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5),
              legend.title = element_blank(),
              legend.text = element_text(size = 12),
              legend.position = "bottom")
    # p
    fig_list <- c(fig_list,
                  list(p))
}
# fig_list2 <- fig_list[c(1,3,5,7,2,4,6,8)]
ggpubr::ggarrange(plotlist = fig_list, 
                  nrow = 1, ncol = 3,
                  common.legend = TRUE, legend = "bottom")
ggsave("./figure/medfly.eps", width = 9, height = 3, dpi = 600)


# ### Example of data - base R version
# y_name <- c("flying","feeding","walking","resting")
# postscript("./figure/medfly.eps", horizontal = F, onefile = F, paper = "special",
#            width = 10, height = 7)
# par(mfrow = c(3, 4))
# # for (i in c(12, 3)) {
# set.seed(100)
# for (i in sample(1:62, 3)) {
#     for (j in 1:4) {
#         plot(Lt[[i]], Ly[[i]][j, ], 
#              type = "l", lwd = 3, col = j,
#              cex.lab = 1.5, cex.axis = 1.5,
#              xlab = "Days", ylab = y_name[j],
#              xlim = c(0, 37), ylim = c(0, 1))
#         grid()
#     }
# }
# dev.off()



### Riemannian FPCA and multivariate FPCA
rfpca.obj <- RFPCA(Lt = Lt,
                   Ly = Ly,
                   optns = list(mfdName = "Sphere",
                                FVEthreshold = 1,
                                userBwMu = "GCV", 
                                userBwCov = "GCV"))
mfpca.obj <- RFPCA(Lt = Lt,
                   Ly = Ly,
                   optns = list(mfdName = "Euclidean",
                                FVEthreshold = 1,
                                userBwMu = "GCV", 
                                userBwCov = "GCV"))
rfpca.obj$K
mfpca.obj$K
cumsum(rfpca.obj$lam[1:5]) / sum(rfpca.obj$lam)
cumsum(mfpca.obj$lam[1:5]) / sum(mfpca.obj$lam)

par(mfrow = c(1, 2))
plot(rfpca.obj$xi[, 1:2], main = "RFPCA")
plot(mfpca.obj$xi[, 1:2], main = "MFPCA")


### Estimated trajectories
par(mfrow = c(2, 2))
for (i in c(12, 3)) {
    # RFPCA
    matplot(t(Ly[[i]]), type = "l", lty = 1, lwd = 2, 
            main = paste0("RFPCA - ", i, "th medfly"), xlab = "t")
    pred <- predict(object = rfpca.obj,
                    newLt = Lt[i],
                    newLy = Ly[i],
                    K = 10,
                    xiMethod = "IN",
                    type = "traj")[1, , ]
    matlines(t(pred), lty = 2, lwd = 2)
    grid()
    if (i == 12) {
        legend("right", 
               c("flying","feeding","walking","resting"), 
               lty = 1, col = 1:4)
    }
    
    # MFPCA
    matplot(t(Ly[[i]]), type = "l", lty = 1, lwd = 2, 
            main = paste0("MFPCA - ", i, "th medfly"), xlab = "t")
    pred <- predict(object = mfpca.obj,
                    newLt = Lt[i],
                    newLy = Ly[i],
                    K = 10,
                    xiMethod = "IN",
                    type = "traj")[1, , ]
    matlines(t(pred), lty = 2, lwd = 2)
    grid()
}



######################################################
### Clustering
######################################################

# devtools::install_github('CrossD/RFPCA')
# devtools::install_url("https://cran.r-project.org/src/contrib/Archive/Funclustering/Funclustering_1.0.2.tar.gz")
library(RFPCA)    # RFPCA and MFPCA
library(mclust)   # clustering measure
library(Funclustering)   # funclust (Currently, it is not supported by cran.)
library(funHDDC)   # funHDDC
library(gmfd)   # gmfd
source("functions.R")

### Model parameters
seed <- 1000
k <- 2    # number of clusters (the number of airlines)
num.pc.method <- "FVE"   # using FVE thresholds
# num.pc.method <- 2     # fixed number
if (num.pc.method == "FVE") {
    FVEthresholdSW <- 0.90
    FVEthresholdCS <- 0.70
    maxK <- Inf
} else if (as.integer(num.pc.method)) {
    FVEthresholdSW <- 1
    FVEthresholdCS <- 1
    maxK <- num.pc.method
}

### kCFC with Riemannian metric
t1 <- Sys.time()
fit.kCFC.Riemann <- kCRFC(y = Ly, 
                          t = Lt, 
                          k = k,
                          kSeed = seed, 
                          maxIter = 125, 
                          optnsSW = list(mfdName = "Sphere",
                                         FVEthreshold = FVEthresholdSW,
                                         maxK = maxK,
                                         # error = T,
                                         userBwMu = "GCV", 
                                         userBwCov = "GCV"),
                          optnsCS = list(mfdName = "Sphere",
                                         FVEthreshold = FVEthresholdCS,
                                         maxK = maxK,
                                         # error = T,
                                         userBwMu = 'GCV', 
                                         userBwCov = 'GCV'))
clust.kCFC.Riemann <- fit.kCFC.Riemann$cluster   # clustering index
clust.kmeans.Riemann <- fit.kCFC.Riemann$clustConf0   # initial cluster
# fit.kCFC.Riemann$clustConf0   # initial clustering index from k-means
t2 <- Sys.time()
print(paste("kCFC (R):", round(t2 - t1, 2)))   # 14.31 secs


### kCFC with Euclidean metric (multivariate FPCA)
t1 <- Sys.time()
fit.kCFC.L2 <- kCRFC(y = Ly, 
                     t = Lt, 
                     k = k,
                     kSeed = seed, 
                     maxIter = 125, 
                     optnsSW = list(mfdName = "Euclidean",
                                    FVEthreshold = FVEthresholdSW,
                                    maxK = maxK,
                                    # error = T,
                                    userBwMu = "GCV", 
                                    userBwCov = "GCV"),
                     optnsCS = list(mfdName = "Euclidean",
                                    FVEthreshold = FVEthresholdCS,
                                    maxK = maxK,
                                    # error = T,
                                    userBwMu = 'GCV', 
                                    userBwCov = 'GCV'))
clust.kCFC.L2 <- fit.kCFC.L2$cluster   # clustering index
clust.kmeans.L2 <- fit.kCFC.L2$clustConf0   # initial cluster
t2 <- Sys.time()
print(paste("kCFC (M):", round(t2 - t1, 2)))   # 7.11 secs


### funclust - set.seed does not working!!
t1 <- Sys.time()
set.seed(seed)
CWtime <- Lt[[1]]
CWfd <- lapply(1:3, function(mdim){
    data <- sapply(Ly, function(y){ y[mdim, ] })
    fda::smooth.basisPar(CWtime, data, lambda = 1e-2)$fd   # B-spline basis
})
# set.seed(seed)
fit.funclust <- funclust(CWfd, K = k, increaseDimension = T)
clust.funclust <- fit.funclust$cls
t2 <- Sys.time()
print(paste("funclust:", round(t2 - t1, 2)))   # 2.86 mins


### funHDDC
t1 <- Sys.time()
set.seed(seed)
fit.funHDDC <- funHDDC(CWfd, 
                       K = k,
                       model = "AkjBQkDk",
                       init = "kmeans",
                       threshold = 0.2)
clust.funHDDC <- fit.funHDDC$class
t2 <- Sys.time()
print(paste("funHDDC:", round(t2 - t1, 2)))   # 0.76 secs


### gmfd
t1 <- Sys.time()
set.seed(seed)
FD <- funData(Lt[[1]], list(
    t( sapply(Ly, function(y){ y[1, ] }) ),
    t( sapply(Ly, function(y){ y[2, ] }) ),
    t( sapply(Ly, function(y){ y[3, ] }) )
))
fit.gmfd <- gmfd_kmeans(FD, n.cl = k, metric = "mahalanobis", p = 10^5)
graphics.off()   # remove plot panel
clust.gmfd <- fit.gmfd$cluster
t2 <- Sys.time()
print(paste("gmfd:", round(t2 - t1, 2)))   # 1.53 mins


table(clust.kCFC.Riemann, clust.kCFC.L2)
table(clust.kCFC.Riemann, clust.kmeans.L2)
table(clust.kCFC.Riemann, clust.kmeans.Riemann)
table(clust.kCFC.Riemann, clust.funclust)
table(clust.kCFC.Riemann, clust.funHDDC)
table(clust.kCFC.Riemann, clust.gmfd)


### Mean functions for each methods
y_name <- c("Cluster 1", "Cluster 2")
m_name <- c("kCRFC","funclust","funHDDC","gmfd")
# m_name <- c("kCRFC","kCFC","funclust","funHDDC","gmfd")
clust_list <- list(clust.kCFC.Riemann,
                   # clust.kCFC.L2,
                   clust.funclust,
                   clust.funHDDC,
                   clust.gmfd)
fig_list <- list()
for (i in 1:length(clust_list)) {
    clust <- clust_list[[i]]
    # # match cluster
    # if (length( mclust::classError(clust.kCFC.Riemann, clust)$misclassified ) > n/3) {
    #     clust <- ifelse(clust == 1, 2, 1)
    # }
    
    mean_ftn <- lapply(1:2, function(cl) {
        ind <- which(clust == cl)
        Reduce("+", Ly[ind]) / length(ind)
    })
    time_points <- Lt[[1]]
    
    for (k in 1:2) {
        df <- data.frame(
            time = time_points,
            t(mean_ftn[[k]])
        ) %>% 
            gather(var_name, val, -time)
        # head(df)
        p <- ggplot(
            data = df, 
            aes(
                x = time,
                y = val,
                group = var_name,
                color = var_name
            )
        ) +
            geom_line() +
            labs(x = "Days", y = "", title = paste0(m_name[i], " - ", y_name[k])) +
            theme_bw() +
            theme(plot.title = element_text(hjust = 0.5),
                  legend.title = element_blank(),
                  legend.text = element_text(size = 12),
                  legend.position = "bottom")
        # p
        fig_list <- c(fig_list,
                      list(p))
    }
}
# fig_list %>% 
#     length
fig_list2 <- fig_list[c(1,3,5,7,2,4,6,8)]
ggpubr::ggarrange(plotlist = fig_list2, 
                  nrow = 2, ncol = 4,
                  common.legend = TRUE, legend = "bottom")
ggsave("./figure/medfly_mean.eps", width = 10, height = 5, dpi = 600)



##########################
### Test plot codes
##########################
### Plot for each variable
par(mfrow = c(2, 4))
y_name <- c("flying","feeding","walking","resting")
# RFPCA
for (j in 1:4) {
    plot("n",
         xlab = "Days", ylab = y_name[j], main = paste0("kCFC(R) - ", y_name[j]),
         xlim = c(0, 37), ylim = c(0, 1))
    for (i in 1:n) {
        # l_col <- ifelse(clust.kCFC.Riemann[i] == 1, "gray", "black")
        # lines(Lt[[i]], Ly[[i]][j, ], col = l_col)
        lines(Lt[[i]], Ly[[i]][j, ], col = clust.kCFC.Riemann[i])
    }
}
# MFPCA
for (j in 1:4) {
    plot("n",
         xlab = "Days", ylab = y_name[j], main = paste0("kCFC(M) - ", y_name[j]),
         xlim = c(0, 37), ylim = c(0, 1))
    # match cluster
    if (length( mclust::classError(clust.kCFC.Riemann, clust.kCFC.L2)$misclassified ) < n/3) {
        clust <- clust.kCFC.L2
    } else {
        clust <- ifelse(clust.kCFC.L2 == 1, 2, 1)
    }
    for (i in 1:n) {
        lines(Lt[[i]], Ly[[i]][j, ], col = clust[i])
    }
}

### Curves of different cluster result
par(mfrow = c(4, 4))
y_name <- c("flying","feeding","walking","resting")
m_name <- c("kCFC","funclust","funHDDC","gmfd")
clust_list <- list(clust.kCFC.L2,
                   clust.funclust,
                   clust.funHDDC,
                   clust.gmfd)
## kCRFC vs kCFC
for (i in 1:length(clust_list)) {
    clust <- clust_list[[i]]
    # match cluster
    if (length( mclust::classError(clust.kCFC.Riemann, clust)$misclassified ) > n/3) {
        clust <- ifelse(clust == 1, 2, 1)
    }
    ind <- which(clust.kCFC.Riemann != clust)   # different cluster result
    for (j in 1:4) {
        plot("n",
             xlab = "Days", ylab = y_name[j], main = paste0("kCRFC vs ", m_name[i], " - ", y_name[j]),
             xlim = c(0, 37), ylim = c(0, 1))
        for (k in ind) {
            lines(Lt[[k]], Ly[[k]][j, ], col = clust.kCFC.Riemann[k])
        }
    }
}
# # match cluster
# if (length( mclust::classError(clust.kCFC.Riemann, clust.kCFC.L2)$misclassified ) < n/3) {
#     clust <- clust.kCFC.L2
# } else {
#     clust <- ifelse(clust.kCFC.L2 == 1, 2, 1)
# }
# ind <- which(clust.kCFC.Riemann != clust)   # different cluster result
# for (j in 1:4) {
#     plot("n",
#          xlab = "Days", ylab = y_name[j], main = paste0("kCRFC - ", y_name[j]),
#          xlim = c(0, 37), ylim = c(0, 1))
#     for (i in ind) {
#         lines(Lt[[i]], Ly[[i]][j, ], col = clust.kCFC.Riemann[i])
#     }
# }


### Mean functions for each methods - base R version
postscript("./figure/medfly_mean.eps", horizontal = F, onefile = F, paper = "special",
           width = 10, height = 9)
par(mfrow = c(4, 4))
y_name <- c("flying","feeding","walking","resting")
m_name <- c("kCRFC","funclust","funHDDC","gmfd")
# m_name <- c("kCRFC","kCFC","funclust","funHDDC","gmfd")
clust_list <- list(clust.kCFC.Riemann,
                   # clust.kCFC.L2,
                   clust.funclust,
                   clust.funHDDC,
                   clust.gmfd)
for (i in 1:length(clust_list)) {
    clust <- clust_list[[i]]
    # match cluster
    if (length( mclust::classError(clust.kCFC.Riemann, clust)$misclassified ) > n/3) {
        clust <- ifelse(clust == 1, 2, 1)
    }
    
    mean_ftn <- lapply(1:2, function(cl) {
        ind <- which(clust == cl)
        Reduce("+", Ly[ind]) / length(ind)
    })
    
    time_points <- Lt[[1]]
    
    for (j in 1:4) {
        plot("n",
             xlab = "Days", ylab = y_name[j], main = paste0(m_name[i], " - ", y_name[j]),
             cex.lab = 1.5, cex.axis = 1.5,
             xlim = c(0, 37), ylim = c(0, 1))
        for (k in 1:2) {
            lines(time_points, mean_ftn[[k]][j, ], col = k, lwd = 2)
        }
        grid()
    }
}
dev.off()


### Plot for each variable
par(mfrow = c(4, 4))
y_name <- c("flying","feeding","walking","resting")
# RFPCA
for (j in 1:4) {
    plot("n",
         xlab = "Days", ylab = y_name[j], main = paste0("kCRFC - ", y_name[j]),
         xlim = c(0, 37), ylim = c(0, 1))
    for (i in which(clust.kCFC.Riemann == 1)) {
        lines(Lt[[i]], Ly[[i]][j, ], col = clust.kCFC.Riemann[i])
    }
}
for (j in 1:4) {
    plot("n",
         xlab = "Days", ylab = y_name[j], main = paste0("kCRFC - ", y_name[j]),
         xlim = c(0, 37), ylim = c(0, 1))
    for (i in which(clust.kCFC.Riemann == 2)) {
        lines(Lt[[i]], Ly[[i]][j, ], col = clust.kCFC.Riemann[i])
    }
}
# funHDDC
for (j in 1:4) {
    plot("n",
         xlab = "Days", ylab = y_name[j], main = paste0("funHDDC - ", y_name[j]),
         xlim = c(0, 37), ylim = c(0, 1))
    for (i in which(clust.funHDDC == 1)) {
        lines(Lt[[i]], Ly[[i]][j, ], col = clust.funHDDC[i])
    }
}
for (j in 1:4) {
    plot("n",
         xlab = "Days", ylab = y_name[j], main = paste0("funHDDC - ", y_name[j]),
         xlim = c(0, 37), ylim = c(0, 1))
    for (i in which(clust.funHDDC == 2)) {
        lines(Lt[[i]], Ly[[i]][j, ], col = clust.funHDDC[i])
    }
}



### Scatter plot
par(mfrow = c(2, 2))
plot(rfpca.obj$xi[, 1:2], main = "kCFC(R)", col = clust.kCFC.Riemann)
plot(rfpca.obj$xi[, 1:2], main = "kmeans(R)", col = clust.kmeans.Riemann)
plot(mfpca.obj$xi[, 1:2], main = "kCFC(M)", col = clust.kCFC.L2)
plot(mfpca.obj$xi[, 1:2], main = "kmeans(M)", col = clust.kmeans.L2)
