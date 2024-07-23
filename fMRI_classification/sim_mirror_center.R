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





### Simulation: mirror statistic centering해도 될까??
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
transf_list <- list(
  qnorm,   # Probit transform (inverse Gaussian)
  function(x){ qt(x, df = 3) },   # Inverse t(3)
  function(x){ qt(x, df = 2) },   # Inverse t(2)
  qcauchy   # Inverse Cauchy  
)
fdp_list <- list()
tpr_list <- list()
fdp_res <- data.frame(
  Raw = rep(NA, length(delta_list)),
  Centered = rep(NA, length(delta_list))
)
rownames(fdp_res) <- delta_list
fdp_res <- rep(list(fdp_res), length(transf_list))
tpr_res <- fdp_res
for (i in 1:length(delta_list)) {
  delta <- delta_list[i]
  print(delta)
  
  fdp <- rep(list(matrix(NA, B, 2)), length(transf_list))
  tpr <- rep(list(matrix(NA, B, 2)), length(transf_list))
  for (b in 1:B) {
    if (b %% 100 == 0) {
      cat(paste(b, " "))
      if (b == B) {
        cat("\n")
      }
    }
    
    # Generate correlated p-values
    set.seed(b)
    n <- 10000
    p <- 0.9
    n1 <- n*p
    n2 <- n*(1-p)
    idx_null <- 1:n1
    
    rho <- 0.9
    z <- rnorm(1)
    y <- rnorm(n1)
    x0 <- sqrt(1-rho)*y + sqrt(rho)*z
    x1 <- rnorm(n2, delta, 1)
    x <- c(x0, x1)
    p_value <- 1 - pnorm(x)
    
    # hist(x)
    # hist(1 - qnorm(p_value))
    # hist(1 - qcauchy(p_value))
    
    
    for (j in 1:length(transf_list)) {
      # Inverse Gaussian transform and find the center
      inv_p_value <- qnorm(1 - p_value)
      center <- find_center(inv_p_value)
      
      # Re-transform to Unif(0,1) and transform using symmetric distribution
      transf <- transf_list[[j]]
      mirror_stat <- transf(pnorm(inv_p_value - center))
      mirror_stat_uncentered <- transf(pnorm(inv_p_value))
      
      # Not centered FDP, TPR
      cutoff <- mirror_stat_cutoff(mirror_stat_uncentered, q = 0.1)
      fdp[[j]][b, 1] <- sum(mirror_stat_uncentered[idx_null] > cutoff) / max(1, sum(mirror_stat_uncentered > cutoff))
      tpr[[j]][b, 1] <- sum(mirror_stat_uncentered[-idx_null] > cutoff) / n2
      
      # Centering between mode
      cutoff <- mirror_stat_cutoff(mirror_stat, q = 0.1)
      fdp[[j]][b, 2] <- sum(mirror_stat[idx_null] > cutoff) / max(1, sum(mirror_stat > cutoff))
      tpr[[j]][b, 2] <- sum(mirror_stat[-idx_null] > cutoff) / n2
      
      
      # # Transform using symmetric distribution
      # transf <- transf_list[[j]]
      # 
      # # Mirror statistics
      # mirror_stat <- transf(1 - p_value)
      # 
      # # Not centered FDP, TPR
      # cutoff <- mirror_stat_cutoff(mirror_stat, q = 0.1)
      # fdp[[j]][b, 1] <- sum(mirror_stat[idx_null] > cutoff) / max(1, sum(mirror_stat > cutoff))
      # tpr[[j]][b, 1] <- sum(mirror_stat[-idx_null] > cutoff) / n2
      # 
      # # Centering between mode
      # center <- find_center(mirror_stat)
      # cutoff <- mirror_stat_cutoff(mirror_stat - center, q = 0.1)
      # fdp[[j]][b, 2] <- sum(mirror_stat[idx_null] - center > cutoff) / max(1, sum(mirror_stat - center > cutoff))
      # tpr[[j]][b, 2] <- sum(mirror_stat[-idx_null] - center > cutoff) / n2
    }
    
    # # Not centered FDP, TPR
    # cutoff <- mirror_stat_cutoff(x, q = 0.1)
    # fdp[b, 1] <- sum(x0 > cutoff) / max(1, sum(x > cutoff))
    # tpr[b, 1] <- sum(x1 > cutoff) / n2
    # 
    # # Centering between mode
    # center <- find_center(x)
    # cutoff <- mirror_stat_cutoff(x - center, q = 0.1)
    # fdp[b, 2] <- sum(x0 - center > cutoff) / max(1, sum(x - center > cutoff))
    # tpr[b, 2] <- sum(x1 - center > cutoff) / n2
  }
  # colMeans(fdp)
  # colMeans(tpr)
  
  print(sapply(fdp, colMeans, na.rm = T))
  print(sapply(tpr, colMeans, na.rm = T))
  
  fdp_list[[i]] <- fdp
  tpr_list[[i]] <- tpr
  
  for (j in 1:length(transf_list)) {
    fdp_res[[j]][i, ] <- colMeans(fdp[[j]])
    tpr_res[[j]][i, ] <- colMeans(tpr[[j]])
  }
  # fdp_res[i, ] <- colMeans(fdp)
  # tpr_res[i, ] <- colMeans(tpr)
}
save(fdp_list, tpr_list, fdp_res, tpr_res, file = "RData/sim_gauss.RData")
res <- mapply(function(x, y){
  cbind(x, y) %>% round(3)
}, fdp_res, tpr_res, SIMPLIFY = F)
res
lapply(res, function(x){ xtable::xtable(x, digits = 3) })

# Average with sd
res_sub <- res
for (i in 1:length(fdp_list)) {
  for (j in 1:length(fdp_list[[i]])) {
    fdp_sd <- apply(fdp_list[[i]][[j]], 2, sd)
    tpr_sd <- apply(tpr_list[[i]][[j]], 2, sd)
    
    res_sub[[j]][i, ] <- c(
      paste0(
        paste(format(res[[j]][i, 1:2], nsmall = 3), 
              format(round(fdp_sd, 3), nsmall = 3), 
              sep = " ("), 
        ")"
      ),
      paste0(
        paste(format(res[[j]][i, 3:4], nsmall = 3), 
              format(round(tpr_sd, 3), nsmall = 3), 
              sep = " ("), 
        ")"
      )
    )
    
  }
}
res_sub
lapply(res_sub, function(x){ xtable::xtable(x, digits = 3) })



# Visualization of estimated FDR and TPR
library(ggpubr)
library(patchwork)
title_list <- c("Inv-N(0,1)", "Inv-t(3)", " Inv-t(2)", "Inv-Cauchy(0,1)")
fig_list <- list()
for (i in 1:length(title_list)) {
  fdp_plot <- fdp_res[[i]] %>% 
    rownames_to_column("x") %>% 
    gather("group","value", -x) %>% 
    mutate(group = ifelse(group == "Raw", "Un-centered", "Centered")) %>% 
    ggplot(aes(x, value, group = group, color = group, shape = group)) +
    geom_point(size = 3) +
    geom_line(linewidth = 0.8) +
    ylim(c(0, 0.5)) +
    labs(x = latex2exp::TeX("$\\delta$"), y = "FDR") +
    theme_light() +
    theme(
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      legend.position = "bottom"
    )
  
  tpr_plot <- tpr_res[[i]] %>% 
    rownames_to_column("x") %>% 
    gather("group","value", -x) %>% 
    mutate(group = ifelse(group == "Raw", "Un-centered", "Centered")) %>% 
    ggplot(aes(x, value, group = group, color = group, shape = group)) +
    geom_point(size = 3) +
    geom_line(linewidth = 0.8) +
    ylim(c(0, 1)) +
    labs(x = latex2exp::TeX("$\\delta$"), y = "Power") +
    theme_light() +
    theme(
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      legend.position = "bottom"
    )
  
  # Combine fdp_plot and tpr_plot into one row
  combined_plot <- (fdp_plot | tpr_plot)
  
  fig_list[[i]] <- combined_plot
}

# Create title plots
title_plots <- lapply(title_list, function(title) {
  ggplot() + 
    labs(title = title) + 
    theme_void() + 
    theme(
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5, margin = margin(t = 20))
    )
})

# Combine title plots and figure plots
combined_plots <- map2(title_plots, fig_list, ~ .x / .y + plot_layout(heights = c(1, 10))) 

# Combine all plots into one
wrap_plots(combined_plots, ncol = 2) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")



####################################
### Simulation - Gaussian (Block correlation)
####################################
B <- 1000
delta_list <- c(1, 1.5, 2, 2.5, 3)
transf_list <- list(
  qnorm,   # Probit transform (inverse Gaussian)
  function(x){ qt(x, df = 3) },   # Inverse t(3)
  function(x){ qt(x, df = 2) },   # Inverse t(2)
  qcauchy   # Inverse Cauchy  
)
fdp_list <- list()
tpr_list <- list()
fdp_res <- data.frame(
  Raw = rep(NA, length(delta_list)),
  Centered = rep(NA, length(delta_list))
)
rownames(fdp_res) <- delta_list
fdp_res <- rep(list(fdp_res), length(transf_list))
tpr_res <- fdp_res
for (i in 1:length(delta_list)) {
  delta <- delta_list[i]
  print(delta)
  
  fdp <- rep(list(matrix(NA, B, 2)), length(transf_list))
  tpr <- rep(list(matrix(NA, B, 2)), length(transf_list))
  for (b in 1:B) {
    if (b %% 100 == 0) {
      cat(paste(b, " "))
      if (b == B) {
        cat("\n")
      }
    }
    
    # Generate correlated p-values
    set.seed(b)
    n <- 10000
    p <- 0.9
    n1 <- n*p
    n2 <- n*(1-p)
    idx_null <- 1:n1
    
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
    x1 <- rnorm(n2, delta, 1)
    x <- c(x0, x1)
    p_value <- 1 - pnorm(x)
    
    # hist(x)
    # hist(1 - qnorm(p_value))
    # hist(1 - qcauchy(p_value))
    
    
    for (j in 1:length(transf_list)) {
      # Inverse Gaussian transform and find the center
      inv_p_value <- qnorm(1 - p_value)
      center <- find_center(inv_p_value)
      
      # Re-transform to Unif(0,1) and transform using symmetric distribution
      transf <- transf_list[[j]]
      mirror_stat <- transf(pnorm(inv_p_value - center))
      mirror_stat_uncentered <- transf(pnorm(inv_p_value))
      
      # Not centered FDP, TPR
      cutoff <- mirror_stat_cutoff(mirror_stat_uncentered, q = 0.1)
      fdp[[j]][b, 1] <- sum(mirror_stat_uncentered[idx_null] > cutoff) / max(1, sum(mirror_stat_uncentered > cutoff))
      tpr[[j]][b, 1] <- sum(mirror_stat_uncentered[-idx_null] > cutoff) / n2
      
      # Centering between mode
      cutoff <- mirror_stat_cutoff(mirror_stat, q = 0.1)
      fdp[[j]][b, 2] <- sum(mirror_stat[idx_null] > cutoff) / max(1, sum(mirror_stat > cutoff))
      tpr[[j]][b, 2] <- sum(mirror_stat[-idx_null] > cutoff) / n2
      
      
      # # Transform using symmetric distribution
      # transf <- transf_list[[j]]
      # 
      # # Mirror statistics
      # mirror_stat <- transf(1 - p_value)
      # 
      # # Not centered FDP, TPR
      # cutoff <- mirror_stat_cutoff(mirror_stat, q = 0.1)
      # fdp[[j]][b, 1] <- sum(mirror_stat[idx_null] > cutoff) / max(1, sum(mirror_stat > cutoff))
      # tpr[[j]][b, 1] <- sum(mirror_stat[-idx_null] > cutoff) / n2
      # 
      # # Centering between mode
      # center <- find_center(mirror_stat)
      # cutoff <- mirror_stat_cutoff(mirror_stat - center, q = 0.1)
      # fdp[[j]][b, 2] <- sum(mirror_stat[idx_null] - center > cutoff) / max(1, sum(mirror_stat - center > cutoff))
      # tpr[[j]][b, 2] <- sum(mirror_stat[-idx_null] - center > cutoff) / n2
    }
    
  }
  
  print(sapply(fdp, colMeans, na.rm = T))
  print(sapply(tpr, colMeans, na.rm = T))
  
  fdp_list[[i]] <- fdp
  tpr_list[[i]] <- tpr
  
  for (j in 1:length(transf_list)) {
    fdp_res[[j]][i, ] <- colMeans(fdp[[j]])
    tpr_res[[j]][i, ] <- colMeans(tpr[[j]])
  }
}
save(fdp_list, tpr_list, fdp_res, tpr_res, file = "RData/sim_gauss_block.RData")
res <- mapply(function(x, y){
  cbind(x, y) %>% round(3)
}, fdp_res, tpr_res, SIMPLIFY = F)
res
lapply(res, function(x){ xtable::xtable(x, digits = 3) })

# Average with sd
res_sub <- res
for (i in 1:length(fdp_list)) {
  for (j in 1:length(fdp_list[[i]])) {
    fdp_sd <- apply(fdp_list[[i]][[j]], 2, sd)
    tpr_sd <- apply(tpr_list[[i]][[j]], 2, sd)
    
    res_sub[[j]][i, ] <- c(
      paste0(
        paste(format(res[[j]][i, 1:2], nsmall = 3), 
              format(round(fdp_sd, 3), nsmall = 3), 
              sep = " ("), 
        ")"
      ),
      paste0(
        paste(format(res[[j]][i, 3:4], nsmall = 3), 
              format(round(tpr_sd, 3), nsmall = 3), 
              sep = " ("), 
        ")"
      )
    )
   
  }
}
res_sub
lapply(res_sub, function(x){ xtable::xtable(x, digits = 3) })


# Visualization of estimated FDR and TPR
library(ggpubr)
library(patchwork)
title_list <- c("Inv-N(0,1)", "Inv-t(3)", " Inv-t(2)", "Inv-Cauchy(0,1)")
fig_list <- list()
for (i in 1:length(title_list)) {
  fdp_plot <- fdp_res[[i]] %>% 
    rownames_to_column("x") %>% 
    gather("group","value", -x) %>% 
    mutate(group = ifelse(group == "Raw", "Un-centered", "Centered")) %>% 
    ggplot(aes(x, value, group = group, color = group, shape = group)) +
    geom_point(size = 3) +
    geom_line(linewidth = 0.8) +
    ylim(c(0, 0.5)) +
    labs(x = latex2exp::TeX("$\\delta$"), y = "FDR") +
    theme_light() +
    theme(
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      legend.position = "bottom"
    )
  
  tpr_plot <- tpr_res[[i]] %>% 
    rownames_to_column("x") %>% 
    gather("group","value", -x) %>% 
    mutate(group = ifelse(group == "Raw", "Un-centered", "Centered")) %>% 
    ggplot(aes(x, value, group = group, color = group, shape = group)) +
    geom_point(size = 3) +
    geom_line(linewidth = 0.8) +
    ylim(c(0, 1)) +
    labs(x = latex2exp::TeX("$\\delta$"), y = "Power") +
    theme_light() +
    theme(
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      legend.position = "bottom"
    )
  
  # Combine fdp_plot and tpr_plot into one row
  combined_plot <- (fdp_plot | tpr_plot)
  
  fig_list[[i]] <- combined_plot
}

# Create title plots
title_plots <- lapply(title_list, function(title) {
  ggplot() + 
    labs(title = title) + 
    theme_void() + 
    theme(
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5, margin = margin(t = 20))
    )
})

# Combine title plots and figure plots
combined_plots <- map2(title_plots, fig_list, ~ .x / .y + plot_layout(heights = c(1, 10))) 

# Combine all plots into one
wrap_plots(combined_plots, ncol = 2) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")





########################################################################
### Simulation - Gaussian (Block correlation - different sample sizes)
########################################################################
B <- 1000
delta_list <- c(1, 1.5, 2, 2.5, 3)
transf_list <- list(
  qnorm,   # Probit transform (inverse Gaussian)
  function(x){ qt(x, df = 3) },   # Inverse t(3)
  function(x){ qt(x, df = 2) },   # Inverse t(2)
  qcauchy   # Inverse Cauchy  
)
fdp_list <- list()
tpr_list <- list()
fdp_res <- data.frame(
  Raw = rep(NA, length(delta_list)),
  Centered = rep(NA, length(delta_list))
)
rownames(fdp_res) <- delta_list
fdp_res <- rep(list(fdp_res), length(transf_list))
tpr_res <- fdp_res
for (i in 1:length(delta_list)) {
  delta <- delta_list[i]
  print(delta)
  
  fdp <- rep(list(matrix(NA, B, 2)), length(transf_list))
  tpr <- rep(list(matrix(NA, B, 2)), length(transf_list))
  for (b in 1:B) {
    if (b %% 100 == 0) {
      cat(paste(b, " "))
      if (b == B) {
        cat("\n")
      }
    }
    
    # Generate correlated p-values
    set.seed(b)
    n <- 10000
    p <- 0.9
    n1 <- n*p
    n2 <- n - n1
    idx_null <- 1:n1
    
    # x0
    rho_list <- c(0.2, 0.8)
    z1 <- rnorm(1)
    n11 <- round(0.6*n1)
    y1 <- rnorm(n11)
    x0 <- sqrt(1-rho_list[1])*y1 + sqrt(rho_list[1])*z1
    # z2 <- rnorm(1)
    z2 <- z1
    y2 <- rnorm(n1 - n11)
    x0 <- c(x0,
            sqrt(1-rho_list[2])*y2 + sqrt(rho_list[2])*z2)
    # x1
    x1 <- rnorm(n2, delta, 1)
    x <- c(x0, x1)
    p_value <- 1 - pnorm(x)
    
    # par(mfrow = c(1, 2))
    # hist(x0, breaks = 80)
    # abline(v = sqrt(rho_list[1])*z1, col = 2, lwd = 2)
    # abline(v = sqrt(rho_list[2])*z2, col = 3, lwd = 2)
    # # abline(v = find_center(x0), col = 2, lwd = 2)
    # # abline(v = find_center(x), col = 3, lwd = 2)
    # hist(x, breaks = 80)
    # # hist(1 - qnorm(p_value))
    # # hist(1 - qcauchy(p_value))
    
    
    for (j in 1:length(transf_list)) {
      # Inverse Gaussian transform and find the center
      inv_p_value <- qnorm(1 - p_value)
      center <- find_center(inv_p_value)
      
      # Re-transform to Unif(0,1) and transform using symmetric distribution
      transf <- transf_list[[j]]
      mirror_stat <- transf(pnorm(inv_p_value - center))
      mirror_stat_uncentered <- transf(pnorm(inv_p_value))
      
      # Not centered FDP, TPR
      cutoff <- mirror_stat_cutoff(mirror_stat_uncentered, q = 0.1)
      fdp[[j]][b, 1] <- sum(mirror_stat_uncentered[idx_null] > cutoff) / max(1, sum(mirror_stat_uncentered > cutoff))
      tpr[[j]][b, 1] <- sum(mirror_stat_uncentered[-idx_null] > cutoff) / n2
      
      # Centering between mode
      cutoff <- mirror_stat_cutoff(mirror_stat, q = 0.1)
      fdp[[j]][b, 2] <- sum(mirror_stat[idx_null] > cutoff) / max(1, sum(mirror_stat > cutoff))
      tpr[[j]][b, 2] <- sum(mirror_stat[-idx_null] > cutoff) / n2
      
      
      # # Transform using symmetric distribution
      # transf <- transf_list[[j]]
      # 
      # # Mirror statistics
      # mirror_stat <- transf(1 - p_value)
      # 
      # # Not centered FDP, TPR
      # cutoff <- mirror_stat_cutoff(mirror_stat, q = 0.1)
      # fdp[[j]][b, 1] <- sum(mirror_stat[idx_null] > cutoff) / max(1, sum(mirror_stat > cutoff))
      # tpr[[j]][b, 1] <- sum(mirror_stat[-idx_null] > cutoff) / n2
      # 
      # # Centering between mode
      # center <- find_center(mirror_stat)
      # cutoff <- mirror_stat_cutoff(mirror_stat - center, q = 0.1)
      # fdp[[j]][b, 2] <- sum(mirror_stat[idx_null] - center > cutoff) / max(1, sum(mirror_stat - center > cutoff))
      # tpr[[j]][b, 2] <- sum(mirror_stat[-idx_null] - center > cutoff) / n2
    }
    
  }
  
  print(sapply(fdp, colMeans, na.rm = T))
  print(sapply(tpr, colMeans, na.rm = T))
  
  fdp_list[[i]] <- fdp
  tpr_list[[i]] <- tpr
  
  for (j in 1:length(transf_list)) {
    fdp_res[[j]][i, ] <- colMeans(fdp[[j]])
    tpr_res[[j]][i, ] <- colMeans(tpr[[j]])
  }
}
save(fdp_list, tpr_list, fdp_res, tpr_res, file = "RData/sim_gauss_block_unbalanced.RData")
res <- mapply(function(x, y){
  cbind(x, y) %>% round(3)
}, fdp_res, tpr_res, SIMPLIFY = F)
res
lapply(res, function(x){ xtable::xtable(x, digits = 3) })

# Average with sd
res_sub <- res
for (i in 1:length(fdp_list)) {
  for (j in 1:length(fdp_list[[i]])) {
    fdp_sd <- apply(fdp_list[[i]][[j]], 2, sd)
    tpr_sd <- apply(tpr_list[[i]][[j]], 2, sd)
    
    res_sub[[j]][i, ] <- c(
      paste0(
        paste(format(res[[j]][i, 1:2], nsmall = 3), 
              format(round(fdp_sd, 3), nsmall = 3), 
              sep = " ("), 
        ")"
      ),
      paste0(
        paste(format(res[[j]][i, 3:4], nsmall = 3), 
              format(round(tpr_sd, 3), nsmall = 3), 
              sep = " ("), 
        ")"
      )
    )
    
  }
}
res_sub
lapply(res_sub, function(x){ xtable::xtable(x, digits = 3) })


# Visualization of estimated FDR and TPR
library(ggpubr)
library(patchwork)
title_list <- c("Inv-N(0,1)", "Inv-t(3)", " Inv-t(2)", "Inv-Cauchy(0,1)")
fig_list <- list()
for (i in 1:length(title_list)) {
  fdp_plot <- fdp_res[[i]] %>% 
    rownames_to_column("x") %>% 
    gather("group","value", -x) %>% 
    mutate(group = ifelse(group == "Raw", "Un-centered", "Centered")) %>% 
    ggplot(aes(x, value, group = group, color = group, shape = group)) +
    geom_point(size = 3) +
    geom_line(linewidth = 0.8) +
    ylim(c(0, 0.5)) +
    labs(x = latex2exp::TeX("$\\delta$"), y = "FDR") +
    theme_light() +
    theme(
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      legend.position = "bottom"
    )
  
  tpr_plot <- tpr_res[[i]] %>% 
    rownames_to_column("x") %>% 
    gather("group","value", -x) %>% 
    mutate(group = ifelse(group == "Raw", "Un-centered", "Centered")) %>% 
    ggplot(aes(x, value, group = group, color = group, shape = group)) +
    geom_point(size = 3) +
    geom_line(linewidth = 0.8) +
    ylim(c(0, 1)) +
    labs(x = latex2exp::TeX("$\\delta$"), y = "Power") +
    theme_light() +
    theme(
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      legend.position = "bottom"
    )
  
  # Combine fdp_plot and tpr_plot into one row
  combined_plot <- (fdp_plot | tpr_plot)
  
  fig_list[[i]] <- combined_plot
}

# Create title plots
title_plots <- lapply(title_list, function(title) {
  ggplot() + 
    labs(title = title) + 
    theme_void() + 
    theme(
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5, margin = margin(t = 20))
    )
})

# Combine title plots and figure plots
combined_plots <- map2(title_plots, fig_list, ~ .x / .y + plot_layout(heights = c(1, 10))) 

# Combine all plots into one
wrap_plots(combined_plots, ncol = 2) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")
