library(tidyverse)
source("R/homo_test.R")

load("RData/homo_test_ddplot.RData")
# p-values
p_value <- data.frame(
  FM = sapply(test_obj, function(x){ x$FM$p_value }),
  RT = sapply(test_obj, function(x){ x$RT$p_value }),
  hM = sapply(test_obj, function(x){ x$hM$p_value }),
  fd2 = sapply(test_obj, function(x){ x$fd2$p_value })
) %>% as.matrix()

load("RData/homo_test_kern.RData")
p_value <- data.frame(
  SET_ID = sapply(test_obj, function(x){ x$SET_ID[[2]] }),
  COV = sapply(test_obj, function(x){ x$COV[[2]] }),
  SET_FPCA = sapply(test_obj, function(x){ x$SET_FPCA[[2]] }),
  SET_SQR = sapply(test_obj, function(x){ x$SET_SQR[[2]] })
) %>% as.matrix()


# Benjamini-Hochberg FDR
q <- 0.1
for (i in 1:4) {
  # which(isTRUE(sort(p_value[, 1]) < q/p * (1:p)))
  print( apply(p_value, 2, function(x){ which(isTRUE(sort(x) < q/p * (1:p))) }) )
}


# Transform using symmetric distribution
trans <- qnorm   # Probit transform (inverse Gaussian)
trans <- qcauchy   # Inverse Cauchy  
trans <- function(x){ qt(x, df = 3) }   # Inverse t(3)
trans <- function(x){ qt(x, df = 2) }   # Inverse t(2)

# Mirror statistics
mirror_stat <- trans(1 - p_value)

# Histogram
par(mfrow = c(2, 2))
for (i in 1:4) {
  x <- mirror_stat[, i]
  # x[which(is.infinite(x))] <- sign(x[which(is.infinite(x))]) * (max(abs(x[which(is.finite(x))])) + 1)
  
  
  # Centering between mode
  center <- find_center(x)
  
  # # Centering between trimmed mean
  # idx_trim <- which(x < median(x) + trans(0.8) & x > median(x) - trans(0.8))
  # # idx_trim <- which(x < quantile(x, 0.75) & x > quantile(x, 0.25))
  # # center <- mean(mean(x[idx_trim]), median(x[idx_trim]))
  # center2 <- mean(x[idx_trim])
  # 
  # # Centering between trimmed mean
  # idx_trim <- which(x < median(x) + trans(0.7) & x > median(x) - trans(0.7))
  # center3 <- mean(x[idx_trim])
  
  hist(x, breaks = 40, probability = T, 
       xlab = latex2exp::TeX("$F^{-1}(1-p)$"),
       main = colnames(mirror_stat)[i])   
  lines(density(x), col = 2, lwd = 2)
  abline(v = 0, col = 1, lwd = 3)
  abline(v = center, col = 4, lwd = 3)
  
  print(colnames(mirror_stat)[i])
  cutoff <-  mirror_stat_cutoff(x, q = 0.1)
  print(paste0("Not centered: ", sum(x > cutoff)))
  cutoff <-  mirror_stat_cutoff(x - center, q = 0.1)
  print(paste0("Centered: ", sum(x - center > cutoff)))
}




####################################
### Variable selection results
####################################
library(tidyverse)
source("R/homo_test.R")

load("RData/homo_test_ddplot.RData")
# p-values
p_value <- data.frame(
  FM = sapply(test_obj, function(x){ x$FM$p_value }),
  RT = sapply(test_obj, function(x){ x$RT$p_value }),
  hM = sapply(test_obj, function(x){ x$hM$p_value }),
  fd2 = sapply(test_obj, function(x){ x$fd2$p_value })
) %>% as.matrix()

load("RData/homo_test_kern.RData")
p_value <- data.frame(
  SET_ID = sapply(test_obj, function(x){ x$SET_ID[[2]] }),
  COV = sapply(test_obj, function(x){ x$COV[[2]] }),
  SET_FPCA = sapply(test_obj, function(x){ x$SET_FPCA[[2]] }),
  SET_SQR = sapply(test_obj, function(x){ x$SET_SQR[[2]] })
) %>% as.matrix()


# Transform using symmetric distribution
transf_list <- list(
  qnorm,   # Probit transform (inverse Gaussian)
  function(x){ qt(x, df = 3) },   # Inverse t(3)
  function(x){ qt(x, df = 2) },   # Inverse t(2)
  qcauchy   # Inverse Cauchy  
)
df <- as.data.frame(matrix(NA, ncol(p_value), length(transf_list)*2))
# rownames(df) <- c("N(0,1)","t(3)","t(2)","Cauchy")
rownames(df) <- colnames(p_value)
colnames(df) <- rep(c("Un-centered", "Centered"), ncol(p_value))
for (i in 1:ncol(p_value)) {
  num_rej <- rep(NA, ncol(df))
  
  for (j in 1:length(transf_list)) {
    # Mirror statistics
    transf <- transf_list[[j]]    # inverse transform function
    mirror_stat_uncentered <- transf(1 - p_value[, i])
    # mirror_stat_uncentered[which(is.infinite(mirror_stat_uncentered))] <- sign(mirror_stat_uncentered[which(is.infinite(mirror_stat_uncentered))]) * (max(abs(mirror_stat_uncentered[which(is.finite(mirror_stat_uncentered))])) + 1)
    
    # Centering between mode
    center <- find_center(mirror_stat_uncentered)
    mirror_stat <- mirror_stat_uncentered - center
    
    # # Inverse Gaussian transform and find the center
    # inv_p_value <- qnorm(1 - p_value[, i])
    # center <- find_center(inv_p_value)
    # 
    # # Re-transform to Unif(0,1) and transform using symmetric distribution
    # transf <- transf_list[[j]]
    # mirror_stat <- transf(pnorm(inv_p_value - center))
    # mirror_stat_uncentered <- transf(pnorm(inv_p_value))
    
    # Not centered FDP, TPR
    cutoff <- mirror_stat_cutoff(mirror_stat_uncentered, q = 0.1)
    num_rej[(j-1)*2 + 1] <- sum(mirror_stat_uncentered > cutoff)
    
    # Centering between mode
    cutoff <- mirror_stat_cutoff(mirror_stat, q = 0.1)
    num_rej[(j-1)*2 + 2] <- sum(mirror_stat > cutoff)
  }
  df[i, ] <- num_rej
}
df
xtable::xtable(df[, c(1,2,4,6,8)])

