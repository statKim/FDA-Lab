# Simulation: mirror statistic centering해도 될까??
# X_0i ~ sqrt(1-rho) * y_i + sqrt(rho) * z,  y_i ~ iid N(0,1), z ~ N(0,1)
# y_i and z are independent
# X_1i ~ N(3, 1)

library(tidyverse)
source("R/homo_test.R")

n <- 10000
p <- 0.9
n1 <- n*p
n2 <- n*(1-p)

rho <- 0.9
z <- rnorm(1)
y <- rnorm(n1)
x0 <- sqrt(1-rho)*y + sqrt(rho)*z
x1 <- rnorm(n2, 2, 1)
x <- c(x0, x1)
z

par(mfrow = c(2, 2))
hist(x0, probability = T, breaks = 30,
     xlab = "", main = latex2exp::TeX("$X_0$"))
lines(density(x0), col = 2, lwd = 2)
hist(x1, probability = T, breaks = 30,
     xlab = "", main = latex2exp::TeX("$X_1$"))
lines(density(x1), col = 2, lwd = 2)
hist(x, probability = T, breaks = 30,
     xlab = "", main = latex2exp::TeX("$X$"))
lines(density(x), col = 2, lwd = 2)
# lines(density(x0), col = 4, lwd = 2)
# lines(density(x1), col = 5, lwd = 2)
# lines(seq(-4, 4, by = 0.1), dnorm(seq(-4, 4, by = 0.1)), col = 3, lwd = 2)
# legend("topright",
#        c(latex2exp::TeX("$X$"), latex2exp::TeX("$X_0$"), latex2exp::TeX("$X_1$"), "N(0,1)"),
#        lty = c(1, 1, 1, 1), lwd = c(2, 2, 2, 2), col = c(2, 4, 5, 3))


# obj <- mirror_stat_cutoff(x, q = 0.1, est = TRUE)
# fdp <- obj$fdp
# cutoff <- obj$cutoff
# 
# fdp <- sum(x0 > cutoff) / sum(x > cutoff)
# tpr <- sum(x1 > cutoff) / n2
# fdp
# tpr

sort(x, decreasing = T) %>% head()
sort(x0, decreasing = T) %>% head()
sort(x1, decreasing = T) %>% head()

# Centering between mode
dens <- density(x)
center <- dens$x[which.max(dens$y)]

hist(x - center, probability = T, breaks = 30,
     xlab = "", main = latex2exp::TeX("$X_{centered}$"))
lines(density(x - center), col = 2, lwd = 2)

cutoff <- mirror_stat_cutoff(x - center, q = 0.1)
cutoff

sum(x0 - center > cutoff) / sum(x - center > cutoff)
sum(x1 - center > cutoff) / n2


####################################
### Simulation - Gaussian
####################################
B <- 1000
delta_list <- c(1, 1.5, 2, 2.5, 3)
fdp_list <- list()
tpr_list <- list()
fdp_res <- data.frame(
  Raw = rep(NA, length(delta_list)),
  Centered = rep(NA, length(delta_list))
)
rownames(fdp_res) <- delta_list
tpr_res <- fdp_res
for (i in 1:length(delta_list)) {
  delta <- delta_list[i]
  print(delta)
  
  fdp <- matrix(NA, B, 2)
  tpr <- matrix(NA, B, 2)
  # mu <- data.frame(
  #   true = rep(NA, B),
  #   mode = rep(NA, B)
  # )
  for (b in 1:B) {
    if (b %% 100 == 0) {
      cat(paste(b, " "))
      if (b == B) {
        cat("\n")
      }
    }
    set.seed(b)
    
    n <- 10000
    p <- 0.9
    n1 <- n*p
    n2 <- n*(1-p)
    
    rho <- 0.9
    z <- rnorm(1)
    y <- rnorm(n1)
    x0 <- sqrt(1-rho)*y + sqrt(rho)*z
    x1 <- rnorm(n2, sqrt(rho)*z + delta, 1)
    x <- c(x0, x1)
    
    # # True center
    # mu$true[b] <- sqrt(rho)*z
    
    # Not centered FDP, TPR
    cutoff <- mirror_stat_cutoff(x, q = 0.1)
    fdp[b, 1] <- sum(x0 > cutoff) / max(1, sum(x > cutoff))
    tpr[b, 1] <- sum(x1 > cutoff) / n2
    
    # Centering between mode
    center <- find_center(x)
    # mu$mode[b] <- center
    cutoff <- mirror_stat_cutoff(x - center, q = 0.1)
    fdp[b, 2] <- sum(x0 - center > cutoff) / max(1, sum(x - center > cutoff))
    tpr[b, 2] <- sum(x1 - center > cutoff) / n2
  }
  # colMeans(fdp)
  # colMeans(tpr)
  # sum((mu$true - mu$mode)^2)
  
  print(colMeans(fdp))
  print(colMeans(tpr))
  
  fdp_list[[i]] <- fdp
  tpr_list[[i]] <- tpr
  
  fdp_res[i, ] <- colMeans(fdp)
  tpr_res[i, ] <- colMeans(tpr)
}
save(fdp_list, tpr_list, fdp_res, tpr_res, file = "RData/sim_gauss.RData")
# fdp_res
# tpr_res
cbind(fdp_res,
      tpr_res) %>% round(3)
cbind(fdp_res,
      tpr_res) %>% xtable::xtable(digits = 3)


# Visualization of estimated FDR and TPR
fig_list <- list()
fig_list[[1]] <- fdp_res %>% 
  rownames_to_column("x") %>% 
  gather("group","value", -x) %>% 
  mutate(group = ifelse(group == "Raw", "Not Centered", "Centered")) %>% 
  ggplot(aes(x, value, group = group, color = group, shape = group)) +
  geom_point(size = 3) +
  geom_line(size = 1) +
  ylim(c(0, 0.5)) +
  labs(x = latex2exp::TeX("$\\delta$"), y = "FDR") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "bottom"
  )
fig_list[[2]] <- tpr_res %>% 
  rownames_to_column("x") %>% 
  gather("group","value", -x) %>% 
  mutate(group = ifelse(group == "Raw", "Not Centered", "Centered")) %>% 
  ggplot(aes(x, value, group = group, color = group, shape = group)) +
  geom_point(size = 3) +
  geom_line(size = 1) +
  ylim(c(0.25, 1)) +
  labs(x = latex2exp::TeX("$\\delta$"), y = "Power") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "bottom"
  )
ggpubr::ggarrange(plotlist = fig_list, 
                  nrow = 1, ncol = 2,
                  common.legend = TRUE, legend = "bottom")



####################################
### Simulation - Gaussian (Block correlation)
####################################
B <- 1000
delta_list <- c(1, 1.5, 2, 2.5, 3)
fdp_list <- list()
tpr_list <- list()
fdp_res <- data.frame(
  Raw = rep(NA, length(delta_list)),
  Centered = rep(NA, length(delta_list))
)
rownames(fdp_res) <- delta_list
tpr_res <- fdp_res
for (i in 1:length(delta_list)) {
  delta <- delta_list[i]
  print(delta)
  
  fdp <- matrix(NA, B, 2)
  tpr <- matrix(NA, B, 2)
  for (b in 1:B) {
    if (b %% 100 == 0) {
      cat(paste(b, " "))
      if (b == B) {
        cat("\n")
      }
    }
    set.seed(b)
    
    n <- 10000
    p <- 0.9
    n1 <- n*p
    n2 <- n*(1-p)
    
    # x0
    rho_list <- c(0.6, 0.8)
    z1 <- rnorm(1)
    y1 <- rnorm(n1/2)
    x0 <- sqrt(1-rho_list[1])*y1 + sqrt(rho_list[1])*z1
    # z2 <- rnorm(1)
    z2 <- z1
    y2 <- rnorm(n1/2)
    x0 <- c(x0,
            sqrt(1-rho_list[2])*y2 + sqrt(rho_list[2])*z2)
    # x1
    x1 <- rnorm(n2, max(sqrt(rho_list[1])*z1, sqrt(rho_list[2])*z2) + delta, 1)
    x <- c(x0, x1)
    
    
    # par(mfrow = c(2, 2))
    # hist(x0, probability = T, breaks = 30,
    #      xlab = "", main = latex2exp::TeX("$X_0$"))
    # lines(density(x0), col = 2, lwd = 2)
    # hist(x1, probability = T, breaks = 30,
    #      xlab = "", main = latex2exp::TeX("$X_1$"))
    # lines(density(x1), col = 2, lwd = 2)
    # hist(x, probability = T, breaks = 30,
    #      xlab = "", main = latex2exp::TeX("$X$"))
    # lines(density(x), col = 2, lwd = 2)
    # 
    # center <- find_center(x, density_cutoff = 0.6)
    # abline(v = center, col = 4, lwd = 3)
    # sum(x0 - center > cutoff) / max(1, sum(x - center > cutoff))
    
    
    # Not centered FDP, TPR
    cutoff <- mirror_stat_cutoff(x, q = 0.1)
    fdp[b, 1] <- sum(x0 > cutoff) / max(1, sum(x > cutoff))
    tpr[b, 1] <- sum(x1 > cutoff) / n2
    
    # Centering between mode
    center <- find_center(x)
    # mu$mode[b] <- center
    cutoff <- mirror_stat_cutoff(x - center, q = 0.1)
    fdp[b, 2] <- sum(x0 - center > cutoff) / max(1, sum(x - center > cutoff))
    tpr[b, 2] <- sum(x1 - center > cutoff) / n2
  }
  # colMeans(fdp)
  # colMeans(tpr)
  
  print(colMeans(fdp))
  print(colMeans(tpr))
  
  fdp_list[[i]] <- fdp
  tpr_list[[i]] <- tpr
  
  fdp_res[i, ] <- colMeans(fdp)
  tpr_res[i, ] <- colMeans(tpr)
}
save(fdp_list, tpr_list, fdp_res, tpr_res, file = "RData/sim_gauss_block.RData")
# fdp_res
# tpr_res
cbind(fdp_res,
      tpr_res) %>% round(3)
cbind(fdp_res,
      tpr_res) %>% xtable::xtable(digits = 3)



# Visualization of estimated FDR and TPR
fig_list <- list()
fig_list[[1]] <- fdp_res %>% 
  rownames_to_column("x") %>% 
  gather("group","value", -x) %>% 
  mutate(group = ifelse(group == "Raw", "Not Centered", "Centered")) %>% 
  ggplot(aes(x, value, group = group, color = group, shape = group)) +
  geom_point(size = 3) +
  geom_line(size = 1) +
  ylim(c(0, 0.5)) +
  labs(x = latex2exp::TeX("$\\delta$"), y = "FDR") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "bottom"
  )
fig_list[[2]] <- tpr_res %>% 
  rownames_to_column("x") %>% 
  gather("group","value", -x) %>% 
  mutate(group = ifelse(group == "Raw", "Not Centered", "Centered")) %>% 
  ggplot(aes(x, value, group = group, color = group, shape = group)) +
  geom_point(size = 3) +
  geom_line(size = 1) +
  ylim(c(0.25, 1)) +
  labs(x = latex2exp::TeX("$\\delta$"), y = "Power") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "bottom"
  )
ggpubr::ggarrange(plotlist = fig_list, 
                  nrow = 1, ncol = 2,
                  common.legend = TRUE, legend = "bottom")



####################################
### Simulation - t(3) + delta
### delta candidate: 2 ~ 4 (1, 1.5일 때 power가 0.1보다도 낮음...)
####################################
df <- 3   # degrees of freedom of t-distribution
B <- 1000
delta_list <- c(2, 2.5, 3, 3.5, 4)
fdp_list <- list()
tpr_list <- list()
fdp_res <- data.frame(
  Raw = rep(NA, length(delta_list)),
  Centered = rep(NA, length(delta_list))
)
rownames(fdp_res) <- delta_list
tpr_res <- fdp_res
for (i in 1:length(delta_list)) {
  delta <- delta_list[i]
  print(delta)
  
  fdp <- matrix(NA, B, 2)
  tpr <- matrix(NA, B, 2)
  for (b in 1:B) {
    if (b %% 100 == 0) {
      cat(paste(b, " "))
      if (b == B) {
        cat("\n")
      }
    }
    set.seed(b)
    
    n <- 10000
    p <- 0.9
    n1 <- n*p
    n2 <- n*(1-p)
    
    rho <- 0.9
    z <- rt(1, df = df)
    y <- rt(n1, df = df)
    x0 <- sqrt(1-rho)*y + sqrt(rho)*z
    x1 <- rt(n2, df = df, ncp = sqrt(rho)*z + delta)
    x <- c(x0, x1)

    # par(mfrow = c(2, 2))
    # hist(x0, probability = T, breaks = 80,
    #      xlab = "", main = latex2exp::TeX("$X_0$"))
    # lines(density(x0), col = 2, lwd = 2)
    # hist(x1, probability = T, breaks = 30,
    #      xlab = "", main = latex2exp::TeX("$X_1$"))
    # lines(density(x1), col = 2, lwd = 2)
    # hist(x, probability = T, breaks = 80, xlim = c(-30, 30),
    #      xlab = "", main = latex2exp::TeX("$X$"))
    # lines(density(x), col = 2, lwd = 2)
    # 
    # center <- find_center(x, density_cutoff = 0.6)
    # abline(v = center, col = 4, lwd = 3)
    # sum(x0 - center > cutoff) / max(1, sum(x - center > cutoff))
    
    # Not centered FDP, TPR
    cutoff <- mirror_stat_cutoff(x, q = 0.1)
    fdp[b, 1] <- sum(x0 > cutoff) / max(1, sum(x > cutoff))
    tpr[b, 1] <- sum(x1 > cutoff) / n2
    
    # Centering between mode
    center <- find_center(x)
    # mu$mode[b] <- center
    cutoff <- mirror_stat_cutoff(x - center, q = 0.1)
    fdp[b, 2] <- sum(x0 - center > cutoff) / max(1, sum(x - center > cutoff))
    tpr[b, 2] <- sum(x1 - center > cutoff) / n2
  }
  # colMeans(fdp)
  # colMeans(tpr)
  # sum((mu$true - mu$mode)^2)
  
  print(colMeans(fdp))
  print(colMeans(tpr))
  
  fdp_list[[i]] <- fdp
  tpr_list[[i]] <- tpr
  
  fdp_res[i, ] <- colMeans(fdp)
  tpr_res[i, ] <- colMeans(tpr)
}
save(fdp_list, tpr_list, fdp_res, tpr_res, file = paste0("RData/sim_t", df, ".RData"))
# fdp_res
# tpr_res
cbind(fdp_res,
      tpr_res) %>% round(3)
cbind(fdp_res,
      tpr_res) %>% xtable::xtable(digits = 3)


# Visualization of estimated FDR and TPR
fig_list <- list()
fig_list[[1]] <- fdp_res %>% 
  rownames_to_column("x") %>% 
  gather("group","value", -x) %>% 
  mutate(group = ifelse(group == "Raw", "Not Centered", "Centered")) %>% 
  ggplot(aes(x, value, group = group, color = group, shape = group)) +
  geom_point(size = 3) +
  geom_line(size = 1) +
  ylim(c(0, 0.5)) +
  labs(x = latex2exp::TeX("$\\delta$"), y = "FDR") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "bottom"
  )
fig_list[[2]] <- tpr_res %>% 
  rownames_to_column("x") %>% 
  gather("group","value", -x) %>% 
  mutate(group = ifelse(group == "Raw", "Not Centered", "Centered")) %>% 
  ggplot(aes(x, value, group = group, color = group, shape = group)) +
  geom_point(size = 3) +
  geom_line(size = 1) +
  ylim(c(0.25, 1)) +
  labs(x = latex2exp::TeX("$\\delta$"), y = "Power") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "bottom"
  )
ggpubr::ggarrange(plotlist = fig_list, 
                  nrow = 1, ncol = 2,
                  common.legend = TRUE, legend = "bottom")




####################################
### Simulation - Cauchy + delta
### delta candidate: 2 ~ 4 (1, 1.5일 때 power가 0.1보다도 낮음...)
### FDR control 안됨
####################################
n <- 10000
p <- 0.9
n1 <- n*p
n2 <- n*(1-p)

rho <- 0.9
z <- rcauchy(1)
y <- rcauchy(n1)
x0 <- sqrt(1-rho)*y + sqrt(rho)*z
x1 <- rcauchy(n2, sqrt(rho)*z + delta, 1)
x <- c(x0, x1)

par(mfrow = c(2, 2))
# hist(x0, probability = T, breaks = 80,
#      xlab = "", main = latex2exp::TeX("$X_0$"))
# lines(density(x0), col = 2, lwd = 2)
# hist(x1, probability = T, breaks = 30,
#      xlab = "", main = latex2exp::TeX("$X_1$"))
# lines(density(x1), col = 2, lwd = 2)
# hist(x, probability = T, breaks = 80, 
#      xlab = "", main = latex2exp::TeX("$X$"))
# lines(density(x), col = 2, lwd = 2)
hist(x0[which(abs(x0) < quantile(abs(x0), 0.99))], probability = T, breaks = 200, xlim = c(-10, 10),
     xlab = "", main = latex2exp::TeX("$X_0$"))
lines(density(x0[which(abs(x0) < quantile(abs(x0), 0.99))]), col = 2, lwd = 2)
hist(x1[which(abs(x1) < quantile(abs(x1), 0.99))], probability = T, breaks = 200, xlim = c(-10, 10)+delta,
     xlab = "", main = latex2exp::TeX("$X_1$"))
lines(density(x1[which(abs(x1) < quantile(abs(x1), 0.99))]), col = 2, lwd = 2)
hist(x[which(abs(x) < quantile(abs(x), 0.99))], probability = T, breaks = 200, xlim = c(-10, 10+delta),
     xlab = "", main = latex2exp::TeX("$X$"))
lines(density(x[which(abs(x) < quantile(abs(x), 0.99))]), col = 2, lwd = 2)

center <- find_center(x)
abline(v = center, col = 4, lwd = 3)
cutoff <- mirror_stat_cutoff(x - center, q = 0.1)
sum(x0 - center > cutoff) / max(1, sum(x - center > cutoff))
sum(x1 > cutoff) / n2




### Centering test
x <- c(rnorm(1000, -1, 1),
       rnorm(1000, 5, 1),
       rnorm(100, 10, 1))
hist(x)
mean(x)
median(x)

idx_trim <- which(x < quantile(x, 0.65) & x > quantile(x, 0.25))
mean(x[idx_trim])
median(x[idx_trim])
hist(x - mean(x[idx_trim]))


sqrt(rho) * z

p_value[, 4] %>% hist
p_value[, 4] %>% qnorm() %>% hist
p_value[, 4] %>% mean
hist(p_value[, 4] - mean(p_value[, 4]) + 1/2)
hist(qnorm(p_value[, 4] - mean(p_value[, 4]) + 1/2))

hist(scale(p_value[, 4]))




### bi-modal centering test
x <- c(rnorm(400, 2, 1),
       rnorm(500, 7, 1),
       rnorm(100, 14, 2))
x <- c(rnorm(400, 2, 1),
       rnorm(500, 7, 1),
       rnorm(400, 13, 1),
       rnorm(100, 17, 2))



# Centering between mode
# dens <- density(x)
# center1 <- dens$x[which.max(dens$y)]
center1 <- find_center(x)

# Centering between trimmed mean
idx_trim <- which(x < median(x) + trans(0.8) & x > median(x) - trans(0.8))
# idx_trim <- which(x < quantile(x, 0.75) & x > quantile(x, 0.25))
# center <- mean(mean(x[idx_trim]), median(x[idx_trim]))
center2 <- mean(x[idx_trim])

# Centering between trimmed mean
idx_trim <- which(x < median(x) + trans(0.7) & x > median(x) - trans(0.7))
center3 <- mean(x[idx_trim])


hist(x, breaks = 40, probability = T, xlim = c(-3, 18),
     xlab = latex2exp::TeX("$F^{-1}(1-p)$"))   
lines(density(x), col = 2, lwd = 2)
abline(v = 0, col = 1, lwd = 3)
abline(v = center1, col = 3, lwd = 3)
abline(v = center2, col = 4, lwd = 3)
abline(v = center3, col = 5, lwd = 3)
# lines(density(x - center1), col = 3, lwd = 2)
# lines(density(x - center2), col = 4, lwd = 2)
# lines(density(x - center3), col = 5, lwd = 2)  





