### Boxplot of eigenfunction angle
par(mfrow = c(2, 2))
load("RData/Delaigle-normal-prop0.RData")
boxplot(mse_eigen2[, c(1,2,3,4,7)],
        ylab = "Eigenfunction angle",
        main = "Gaussian")

load("RData/Delaigle-tdist.RData")
boxplot(mse_eigen2[, c(1,2,3,4,7)],
        ylab = "Eigenfunction angle",
        main = "t(3)")

load("RData/Delaigle-normal-prop1.RData")
boxplot(mse_eigen2[, c(1,2,3,4,7)],
        ylab = "Eigenfunction angle",
        main = "10% contaminated")

load("RData/Delaigle-normal-prop2.RData")
boxplot(mse_eigen2[, c(1,2,3,4,7)],
        ylab = "Eigenfunction angle",
        main = "20% contaminated")


subspace(c(1,1), c(-1,-1))
work.grid <- seq(0, 1, length.out = 51)
eig.true <- get_delaigle_eigen(work.grid, model = 2) 
eig.est <- pca.est[[1]]$pca.obj$pca.ogk.sm.MM.obj$eig.fun
subspace(eig.true, eig.est)
subspace(eig.true, check_eigen_sign(eig.est, eig.true))

i <- 1
subspace(eig.true[, i], eig.est[, i])


### 각 eigenfunction마다 계산하면 작게 나오는데, 전체 subspace로 하면 나쁘게 나오네...
load("RData/Delaigle-normal-prop0.RData")
k <- 4
sapply(pca.est[[1]]$pca.obj, function(obj) {
  subspace(eig.true[, k], obj$eig.fun[, k])
})



# black : True
# red : yao
# green : kraus
# blue : proposed
par(mfrow = c(2, 2))
for (i in 1:4) {
  matplot(cbind(eig.true[, i],
                check_eigen_sign(pca.est[[1]]$pca.obj$pca.yao.obj$eig.fun[, i], eig.true[, i]),
                check_eigen_sign(pca.est[[1]]$pca.obj$pca.kraus.obj$eig.fun[, i], eig.true[, i]),
                check_eigen_sign(pca.est[[1]]$pca.obj$pca.ogk.sm.MM.obj$eig.fun[, i], eig.true[, i])),
          type = "l", 
          ylab = "", main = paste("Eigen", i))
}



par(mfrow = c(2, 2))
x.2 <- robfpca::sim_corr(n = n, 
                         type = data_type,  
                         out.prop = 0, 
                         out.type = out_type, 
                         dist = "normal",
                         dist.mat = dist.mat[1:n, 1:n])
x <- list2matrix(x.2)
matplot(t(x), type = "l")
x.2 <- robfpca::sim_corr(n = n, 
                         type = data_type,  
                         out.prop = out_prop, 
                         out.type = out_type, 
                         dist = "tdist",
                         dist.mat = dist.mat[1:n, 1:n])
x <- list2matrix(x.2)
matplot(t(x), type = "l")
x.2 <- sim_corr(n = n, 
                type = data_type,  
                out.prop = 0, 
                out.type = out_type, 
                dist = "normal",
                dist.mat = dist.mat[1:n, 1:n])
x <- list2matrix(x.2)
matplot(t(x), type = "l")
x.2 <- sim_corr(n = n, 
                type = data_type,  
                out.prop = out_prop, 
                out.type = out_type, 
                dist = "tdist",
                dist.mat = dist.mat[1:n, 1:n])
x <- list2matrix(x.2)
matplot(t(x), type = "l")



# matplot(get_corr_eigen(work.grid), type = "l")

work.grid <- seq(0, 1, length.out = 21)

phi <- get_delaigle_eigen(work.grid)

r.par <- 200
cor.sim <- matrix(fields::Matern(as.vector(dist.mat),
                                 range = r.par,
                                 smoothness = 1),
                  nrow = 403)
cor.sim[1:5, 1:5]

cc <- get_delaigle_cov(work.grid)
dim(cc)

n <- 50
q <- 21
cov_cross <- matrix(0, n*q, n*q)
dim(cov_cross)

cov_cross <- kronecker(cor.sim[1:n, 1:n], cc)
dim(cov_cross)
# heatmap(cov_cross)


par(mfrow = c(2, 2))
system.time({
  x <- mvtnorm::rmvnorm(n = 1, mean = rep(0, n*q), sigma = cov_cross)
})
dim(x)
x <- matrix(x, nrow = n, byrow = T)
matplot(t(x), type = "l")

system.time({
  x <- LaplacesDemon::rmvt(n = 1,
                           mu = rep(0, n*q),
                           S = lqmm::make.positive.definite(cov_cross),
                           df = 3)
})
dim(x)
x <- matrix(x, nrow = n, byrow = T)
matplot(t(x), type = "l")






plot(x.2$x.full[ind, ], type = "l", ylim = c(-10, 4))
lines(pred_reconstr[[6]][ind, ], col = 2)
# lines(pred_comp[[5]], col = 2)
lines(x[ind, ], col = 3)




## bw 줄여서 다시 해볼 것
robfpca::cv.cov_ogk(x,  
                    K = 5, 
                    bw_cand = seq(0.01, 0.1, length.out = 10),
                    MM = TRUE,
                    type = 'huber')
cv.cov_ogk(x,  
           K = 5, 
           bw_cand = seq(0.01, 0.1, length.out = 10),
           MM = TRUE,
           type = 'huber')


GA::persp3D(work.grid, work.grid,
            cov(x.2$x.full),
            theta = -70, phi = 30, expand = 1,
            xlab = 't', ylab = 's', zlab = '', main = 'True')
GA::persp3D(work.grid, work.grid,
            cov.ogk.sm,
            theta = -70, phi = 30, expand = 1,
            xlab = 't', ylab = 's', zlab = '', main = 'Proposed FPCA')  

GA::persp3D(work.grid, work.grid,
            cov.tdist,
            theta = -70, phi = 30, expand = 1,
            xlab = 't', ylab = 's', zlab = '', main = 'Proposed FPCA')  



cov.obj <- cov_ogk(x,   
                   type = "huber",
                   MM = TRUE,
                   smooth = T, 
                   bw = 0.3)
mu.ogk.sm <- cov.obj$mean
cov.ogk.sm <- cov.obj$cov
pca.ogk.sm.obj <- funPCA(x.2$Lt, x.2$Ly,
                         mu.ogk.sm, cov.ogk.sm, sig2 = 0,
                         work.grid, PVE = pve, K = K)
pred_reconstr <- list(
  predict(pca.yao.obj, K = K),
  NA,   # Kraus does not do reconstruction
  NA,   # R-Kraus does not do reconstruction
  predict(pca.boente.obj, K = K),
  predict(pca.ogk.obj, K = K),
  predict(pca.ogk.sm.obj, K = K)
)



cv.cov_ogk(x,  
           K = 5, 
           bw_cand = seq(0.3, 0.4, length.out = 10),
           MM = TRUE,
           type = 'huber')


p <- c()
for (i in 1:100) {
  set.seed(i)
  x.2 <- sim_delaigle(n = n,  
                      type = data_type, 
                      out.prop = out_prop, 
                      out.type = out_type, 
                      dist = dist_type,
                      f = 0.0) 
  x <- list2matrix(x.2)
  p[i] <- sum(is.na(x)) / length(x)
}
mean(p)
# d=1.4, f=0.0 => 0.002745098
# d=1.4, f=0.05 => 0.02788235
# d=1.4, f=0.1 => 0.05305686
# d=1.4, f=0.15 => 0.07828431
# d=1.4, f=0.2 => 0.1035333
# d=1.4, f=0.25 => 0.1289373
# d=1.4, f=0.3 => 0.1545
# d=1.4, f=0.35 => 0.1805157
# d=1.4, f=0.4 => 0.2064902
# d=1.4, f=0.5 => 0.2588078

par(mfrow = c(2, 2))
boxplot(p)
hist(p)
matplot(t(x), type = "l")






#############################################
### Calculate distance matrix from real data
### PM10
### - Example of partially observed data
#############################################
library(tidyverse)
### Data Load
load("/Users/hyunsung/GoogleDrive/Lab/KHS/partiall_obs/real_data/PM10/pm_data_korea.RData")
head(new_data)
dim(new_data)

data <- new_data[, c(1, 5, 6)] %>% 
  drop_na() %>% 
  distinct()
dim(data)
data
mat <- dist(data[, -1]) %>% 
  as.matrix()
dim(mat)
mat

library(geosphere)
mat <- distm(data[, -1])
dim(mat)
mat[1:5, 1:5]

hav.dist <- function(long1, lat1, long2, lat2) {
  R <- 6371
  diff.long <- (long2 - long1)
  diff.lat <- (lat2 - lat1)
  a <- sin(diff.lat/2)^2 + cos(lat1) * cos(lat2) * sin(diff.long/2)^2
  b <- 2 * asin(pmin(1, sqrt(a))) 
  d = R * b
  return(d)
}
hav.dist(data[1, 2], data[1, 3], data[2, 2], data[2, 3])

distm(data[1:3, -1], fun = distVincentyEllipsoid)




var <- 3

loc <- unique(new_data[, 1])
loc <- as.numeric(loc)

### 3,4,5월 자료만 뽑고, missing이 많은 location만 뽑음
A <- vector()
for (loc_num in 1:length(loc)) {
  data <- new_data %>% 
    filter(new_data[, 1] == loc[loc_num])
  if (length(which(data[, 2] == '2017030101')) != 0) {
    time_index <- c(which(data[, 2] == '2017030101'):which(data[, 2] == '2017053124'))
    
    data <- data[time_index, ]
    A[loc_num] <- length(which(is.na(data[, var]) == T)) / nrow(data) * 100
  }
}

### missing 15% 이상인 location
rbind(which(A > 15), A[which(A > 15)])

### missing 15% 이상인 location 4개에 대해서 full_data 만듬
### 각 지역은 list에 들어있고, 각 list에는 day x hour matrix 들어있음
full_data <- list()
for (i in 1:4) {
  loc_num <- c(139, 240, 168, 228)[i]
  data <- new_data %>% 
    filter(new_data[, 1] == loc[loc_num])
  
  time_index <- c(which(data[, 2] == '2017030101'):which(data[, 2] == '2017053124'))
  
  full_data[[i]] <- matrix(data[time_index, var],
                           ncol = 24,
                           byrow = T)
}
dim(full_data[[1]])   # 92 x 24

### trajectories
par(mfrow = c(2, 2))
for (i in 1:4) {
  matplot(t(full_data[[i]]), 
          type = "l", 
          ylim = c(0, 300),
          col = gray(0.8), 
          lwd = 2,
          lty = 1,
          main = paste('location ID:', i),
          ylab = 'PM10',
          xlab = 'Hour')
  abline(h = 100, lty = 2)
}

# number of missing
sapply(full_data, function(x){ sum(is.na(x)) })

# different color between fully observed and with missing
col_ind <- 1:nrow(full_data[[3]])
idx <- which(apply(full_data[[3]], 1, function(x){ sum(is.na(x)) }) > 0)
col_ind[-idx] <- 0

df <- as.data.frame(full_data[[3]])
colnames(df) <- 1:24
df <- df %>% 
  mutate(id = 1:nrow(df)) %>% 
  gather("time", "value", -id) %>% 
  mutate(time = as.numeric(time),
         col_ind = ifelse(id %in% idx, 1, 0),
         id = factor(id)) %>% 
  mutate(col_ind = factor(col_ind))
p <- ggplot(df, aes(x = time,
                    y = value,
                    color = id,
                    group = id)) +
  geom_line(size = 1) +
  # geom_line(size = 1,
  #           color = "darkgray") +
  geom_hline(yintercept = 100,
             linetype = "dashed",
             size = 1) +
  theme_bw() +
  labs(x = "Hour", 
       y = "PM10") +
  scale_x_continuous(breaks = c(0,6,12,18,24)) +
  scale_color_manual(values = c("darkgray","blue")) +
  # guides(color = guide_legend(nrow = 1)) +   # legend 1 line
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    text = element_text(size = 30),   # text size
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    panel.border = element_rect(colour = "black", fill = NA, size = 2),   # panel border
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    # panel.grid.major = element_blank(), # get rid of major grid
    # panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    # legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )
p






x.2 <- sim_corr(n = n,
                type = data_type,
                out.prop = 0,
                out.type = out_type,
                dist = "tdist",
                dist.mat = dist.mat[1:n, 1:n])
x <- list2matrix(x.2)
matplot(t(x), type = "l")

dist.mat[1:5, 1:5]


# 원래거
cor.sim <- fields::Matern(as.vector(dist.mat),
                          range = 500,
                          smoothness = 1)
cor.sim <- cor.sim[-which(cor.sim == 1)]
summary(cor.sim)



cor.sim <- matrix(fields::Matern(as.vector(dist.mat),
                                 range = 800,
                                 smoothness = 1),
                  nrow = nrow(dist.mat))
cor.sim[1:5, 1:5]
cor.sim[c(1,100,200,300,310), c(1,100,200,300,310)]

library(RColorBrewer)
cc <- rainbow(nrow(cor.sim))
# heatmap(cor.sim, Colv = NA, Rowv = NA, scale="column")
heatmap(cor.sim, Rowv = NA, Colv = NA, symm = TRUE, 
        # RowSideColors = cc, ColSideColors = cc,
        labRow = FALSE, labCol = FALSE,
        col = colorRampPalette(brewer.pal(8, "Blues"))(25))

# ggplot으로 그리기
library(hrbrthemes)
df <- expand.grid(X = as.character(1:nrow(cor.sim)),
                  Y = as.character(1:nrow(cor.sim)))
df$Z <- as.numeric(cor.sim)
ggplot(df, aes(X, Y, fill = Z)) + 
  geom_tile() +
  labs(x = "", y = "") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())



cor.sim %>% round(1) %>% sum

r.par <- 350
cor.sim <- fields::Matern(as.vector(dist.mat),
                          range = r.par,
                          smoothness = 1)
# cor.sim <- cor.sim[which(cor.sim > 0.01)]
cor.sim <- cor.sim[which(cor.sim < 1)]
summary(cor.sim)
hist(cor.sim, main = paste("zeta =", r.par))


par(mfrow = c(2, 3))
r.par <- 800
loc_ind <- sample(1:ncol(dist.mat), 100)
x.2 <- sim_corr(n = 100,
                type = "partial",
                dist = "tdist",
                dist.mat = dist.mat[loc_ind, loc_ind],
                r.par = r.par,
                nu = 1)
x <- list2matrix(x.2)
matplot(t(x), type = "l")




load("RData/휴지통/Corr-tdist_zeta300.RData")

### Summary results
if (is.null(K)) {
  PVE_K <- K_res
} else {
  PVE_K <- pve_res
}

res <- data.frame(Method = c("Yao","Kraus","R-Kraus","Boente",
                             "OGK(non-smooth)","OGK(smooth)")) %>% 
  # PVE
  left_join(data.frame(
    Method = colnames(PVE_K),
    "PVE" = format(round(colMeans(PVE_K), 3), 3)
  ), by = "Method") %>% 
  # Eigen MISE
  left_join(data.frame(
    Method = colnames(mse_eigen),
    "Eigen MISE" = paste0(
      format(round(colMeans(mse_eigen), 3), 3),
      " (",
      format(round(apply(mse_eigen, 2, sd), 3), 3),
      ")"
    )
  ), by = "Method") %>% 
  # Eigen Angle
  left_join(data.frame(
    Method = colnames(mse_eigen2),
    "Eigen angle" = paste0(
      format(round(colMeans(mse_eigen2), 3), 3),
      " (",
      format(round(apply(mse_eigen2, 2, sd), 3), 3),
      ")"
    )
  ), by = "Method") %>% 
  # Reconstruction MISE
  left_join(data.frame(
    Method = colnames(mse_reconstr),
    "Recon MISE" = paste0(
      format(round(colMeans(mse_reconstr), 3), 3),
      " (",
      format(round(apply(mse_reconstr, 2, sd), 3), 3),
      ")"
    )
  ), by = "Method") %>% 
  # Completion MISE
  left_join(data.frame(
    Method = colnames(mse_completion),
    "Comp MISE" = paste0(
      format(round(colMeans(mse_completion), 3), 3),
      " (",
      format(round(apply(mse_completion, 2, sd), 3), 3),
      ")"
    )
  ), by = "Method")
print(res)






matplot(cbind(
  -sqrt(2)*cos(2*pi*t),
  sqrt(2)*sin(4*pi*t),
  sqrt(2)*cos(6*pi*t)
), type = "l")
