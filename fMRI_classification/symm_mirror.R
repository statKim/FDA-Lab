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
  x[which(is.infinite(x))] <- sign(x[which(is.infinite(x))]) * (max(abs(x[which(is.finite(x))])) + 1)
  
  
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


