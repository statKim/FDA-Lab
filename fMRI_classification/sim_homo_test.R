##################################################
### Load fMRI data
##################################################
library(tidyverse)

# Class label of 191 subjects
y <- read.csv("./fMRI_Classification/class_label.csv", header = T)
y <- y$y

# Functional covariates from 82 retions
file_list <- list.files("./fMRI_Classification/AlignedSubject/")
ord <- sapply(strsplit(file_list, ".csv"), function(x){ strsplit(x, "AlignedSub")[[1]][2] }) %>% 
  as.integer() %>% 
  order()
file_list <- file_list[ord]
for (i in 1:length(file_list)) {
  df <- read.csv(paste0("./fMRI_Classification/AlignedSubject/", file_list[i]), header = T)
  df <- df[, -1] %>% as.matrix()
  
  if (i == 1) {
    X <- array(0, c(length(file_list), ncol(df), nrow(df)))
  }
  X[i, , ] <- t(df)
}
dim(X)

n <- dim(X)[1]   # number of curves
m <- dim(X)[2]   # number of timepoints
p <- dim(X)[3]  # number of functional variables
gr <- seq(0, 1, length.out = m)


##################################################
### Homogeneity test for each region
##################################################

source("R/homo_test.R")

# Parallel computing
cl <- makePSOCKcluster(10)
registerDoParallel(cl)

B <- 500   # number of bootstrap
test_obj <- list()
for (i in 1:p) {
  print(paste(i, "th region"))
  
  set.seed(i)
  X1 <- X[y == 0, , i]
  X2 <- X[y == 1, , i]
  
  test_obj[[i]] <- list()
  
  # FM - 0.6 secs for B=10
  print("FM")
  start_time <- Sys.time()
  test_obj[[i]]$FM <- fun_homo_test(X1, X2, method = "ddplot", depth = "FM", B = B)
  end_time <- Sys.time()
  print(end_time - start_time)
  
  # # RPD - 43 secs for B=10
  # print("RPD")
  # start_time <- Sys.time()
  # test_obj[[i]]$RPD <- fun_homo_test(X1, X2, method = "ddplot", depth = "RPD", B = B)
  # end_time <- Sys.time()
  # print(end_time - start_time)
  
  # RT - 2 secs for B=10
  print("RT")
  start_time <- Sys.time()
  test_obj[[i]]$RT <- fun_homo_test(X1, X2, method = "ddplot", depth = "RT", B = B)
  end_time <- Sys.time()
  print(end_time - start_time)
  
  # # ABD - 23 secs for B=10
  # print("ABD")
  # start_time <- Sys.time()
  # test_obj[[i]]$ABD <- fun_homo_test(X1, X2, method = "ddplot", depth = "ABD", B = B)
  # end_time <- Sys.time()
  # print(end_time - start_time)
  
  # hM - 0.5 secs for B=10
  print("hM")
  start_time <- Sys.time()
  test_obj[[i]]$hM <- fun_homo_test(X1, X2, method = "ddplot", depth = "hM", B = B)
  end_time <- Sys.time()
  print(end_time - start_time)
  
  # fd2 - 6.5 secs for B=10
  print("fd2")
  start_time <- Sys.time()
  test_obj[[i]]$fd2 <- fun_homo_test(X1, X2, method = "ddplot", depth = "fd2", B = B)
  end_time <- Sys.time()
  print(end_time - start_time)
  
  # print(test_obj[[i]]$p_value)
  # print(end_time - start_time)
}
stopCluster(cl)
unregister()
# save(test_obj, file = "RData/homo_test_ddplot.RData")


###############################################################
### Variable selection using mirror statistics controlling FDR
###############################################################

# p-values
p_value <- data.frame(
  FM = sapply(test_obj, function(x){ x$FM$p_value }),
  RT = sapply(test_obj, function(x){ x$RT$p_value }),
  hM = sapply(test_obj, function(x){ x$hM$p_value }),
  fd2 = sapply(test_obj, function(x){ x$fd2$p_value })
) %>% as.matrix()


# Transform using symmetric distribution
transf <- "norm"
transf <- "cauchy"

if (transf == "norm") {
  mirror_stat <- qnorm(1 - p_value)   # probit transform  
} else if (transf == "cauchy") {
  mirror_stat <- qcauchy(1 - p_value)   # Inverse Cauchy transform  
}

# Histogram
par(mfrow = c(2, 2))
for (i in 1:4) {
  hist(mirror_stat[, i], breaks = 20, probability = T, 
       xlab = latex2exp::TeX("$F^{-1}(1-p)$"),
       main = colnames(mirror_stat)[i])   
  lines(density(mirror_stat[, i]), col = 2, lwd = 2)
  abline(v = 0, col = 1, lwd = 3)
  
  if (transf == "norm") {
    lines(seq(-3, 3, by = 0.1), dnorm(seq(-3, 3, by = 0.1)), col = 3, lwd = 2)
  } else if (transf == "cauchy") {
    lines(seq(-30, 30, by = 0.1), dnorm(seq(-30, 30, by = 0.1)), col = 3, lwd = 2)
  }
}
legend("topright",
       c("Density", "N(0,1)"),
       lty = c(1, 1), lwd = c(2, 2), col = c(2, 3))


mirror_stat_cutoff(mirror_stat[, 1], q = 0.1)
mirror_stat_cutoff(mirror_stat[, 2], q = 0.1)
mirror_stat_cutoff(mirror_stat[, 3], q = 0.1)
mirror_stat_cutoff(mirror_stat[, 4], q = 0.1)

mirror_stat_cutoff(mirror_stat[, 1], q = 0.2)
mirror_stat_cutoff(mirror_stat[, 2], q = 0.2)
mirror_stat_cutoff(mirror_stat[, 3], q = 0.2)
mirror_stat_cutoff(mirror_stat[, 4], q = 0.2)

# Selected variables
cutoff <- mirror_stat_cutoff(mirror_stat[, 3], q = 0.2)
sum(mirror_stat[, 3] > cutoff)

cutoff <- mirror_stat_cutoff(mirror_stat[, 4], q = 0.1)
sum(mirror_stat[, 4] > cutoff)

cutoff <- mirror_stat_cutoff(mirror_stat[, 4], q = 0.2)
sum(mirror_stat[, 4] > cutoff)


# Histogram with FDR cutoff
par(mfrow = c(1, 2))
i <- 3
hist(mirror_stat[, i], breaks = 50, probability = T, 
     xlab = latex2exp::TeX("$F^{-1} (p)$"),
     main = colnames(mirror_stat)[i])   
lines(density(mirror_stat[, i]), col = 2, lwd = 2)
abline(v = 0, col = 1, lwd = 3)
abline(v = mirror_stat_cutoff(mirror_stat[, i], q = 0.2),
       col = 4, lwd = 3)
i <- 4
hist(mirror_stat[, i], breaks = 50, probability = T, 
     xlab = latex2exp::TeX("$F^{-1} (p)$"),
     main = colnames(mirror_stat)[i])   
lines(density(mirror_stat[, i]), col = 2, lwd = 2)
abline(v = 0, col = 1, lwd = 3)
abline(v = mirror_stat_cutoff(mirror_stat[, i], q = 0.1),
       col = 3, lwd = 3)
abline(v = mirror_stat_cutoff(mirror_stat[, i], q = 0.2),
       col = 4, lwd = 3)
legend("topright",
       c("FDR = 0.1", "FDR = 0.2"),
       lty = c(1, 1), lwd = c(3, 3), col = c(3, 4))




##################################################
### Load EEG data
##################################################
library(tidyverse)

# Load EEG dataset
X1 <- data.table::fread("eeg/alcoholic_data.txt")
X2 <- data.table::fread("eeg/control_data.txt")
# dim(X1)
# dim(X2)

# Combine 2 datasets
X <- rbind(X1, X2)
X <- as.matrix(X)
# dim(X)
# par(mfrow = c(2, 2))
# for (p in c(1,2,20,30)-1) {
#   matplot(t(X[90:120, (p*256+1):(p*256+256)]), type = "l", col = 1)
#   matlines(t(X[1:30, (p*256+1):(p*256+256)]), type = "l", col = 2)  
# }

# Transform to 3D array
n <- 122
m <- 256
p <- 64
X <- array(X, c(n, m, p))
# dim(X)
# par(mfrow = c(2, 2))
# for (p in c(1,2,20,30)) {
#   matplot(t(X[90:120, , p]), type = "l", col = 1)
#   matlines(t(X[1:30, , p]), type = "l", col = 2)
# }

# Class labels
y <- c(rep(0, 77), rep(1, 45))


### Homogeneity test for each region
# Parallel computing
cl <- makePSOCKcluster(10)
registerDoParallel(cl)

B <- 500   # number of bootstrap
test_obj <- list()
for (i in 1:p) {
  print(paste(i, "th region"))
  
  set.seed(i)
  X1 <- X[y == 0, , i]
  X2 <- X[y == 1, , i]
  
  test_obj[[i]] <- list()
  
  # FM - 0.6 secs for B=10
  print("FM")
  start_time <- Sys.time()
  test_obj[[i]]$FM <- fun_homo_test(X1, X2, method = "ddplot", depth = "FM", B = B)
  end_time <- Sys.time()
  print(end_time - start_time)
  
  # # RPD - 43 secs for B=10
  # print("RPD")
  # start_time <- Sys.time()
  # test_obj[[i]]$RPD <- fun_homo_test(X1, X2, method = "ddplot", depth = "RPD", B = B)
  # end_time <- Sys.time()
  # print(end_time - start_time)
  
  # RT - 2 secs for B=10
  print("RT")
  start_time <- Sys.time()
  test_obj[[i]]$RT <- fun_homo_test(X1, X2, method = "ddplot", depth = "RT", B = B)
  end_time <- Sys.time()
  print(end_time - start_time)
  
  # # ABD - 23 secs for B=10
  # print("ABD")
  # start_time <- Sys.time()
  # test_obj[[i]]$ABD <- fun_homo_test(X1, X2, method = "ddplot", depth = "ABD", B = B)
  # end_time <- Sys.time()
  # print(end_time - start_time)
  
  # hM - 0.5 secs for B=10
  print("hM")
  start_time <- Sys.time()
  test_obj[[i]]$hM <- fun_homo_test(X1, X2, method = "ddplot", depth = "hM", B = B)
  end_time <- Sys.time()
  print(end_time - start_time)
  
  # fd2 - 6.5 secs for B=10
  print("fd2")
  start_time <- Sys.time()
  test_obj[[i]]$fd2 <- fun_homo_test(X1, X2, method = "ddplot", depth = "fd2", B = B)
  end_time <- Sys.time()
  print(end_time - start_time)
  
  # print(test_obj[[i]]$p_value)
  # print(end_time - start_time)
}
stopCluster(cl)
unregister()
# save(test_obj, file = "RData/homo_test_ddplot_eeg.RData")



### Variable selection using mirror statistics controlling FDR

# p-values
p_value <- data.frame(
  FM = sapply(test_obj, function(x){ x$FM$p_value }),
  RT = sapply(test_obj, function(x){ x$RT$p_value }),
  hM = sapply(test_obj, function(x){ x$hM$p_value }),
  fd2 = sapply(test_obj, function(x){ x$fd2$p_value })
) %>% as.matrix()

# Transform using symmetric distribution
transf <- "norm"
transf <- "cauchy"

if (transf == "norm") {
  mirror_stat <- qnorm(1 - p_value)   # probit transform  
} else if (transf == "cauchy") {
  mirror_stat <- qcauchy(1 - p_value)   # Inverse Cauchy transform  
}

# Histogram
par(mfrow = c(2, 2))
for (i in 1:4) {
  hist(mirror_stat[, i], breaks = 20, probability = T, 
       xlab = latex2exp::TeX("$F^{-1}(1-p)$"),
       main = colnames(mirror_stat)[i])   
  lines(density(mirror_stat[, i]), col = 2, lwd = 2)
  abline(v = 0, col = 1, lwd = 3)
  
  if (transf == "norm") {
    lines(seq(-3, 3, by = 0.1), dnorm(seq(-3, 3, by = 0.1)), col = 3, lwd = 2)
  } else if (transf == "cauchy") {
    lines(seq(-30, 30, by = 0.1), dnorm(seq(-30, 30, by = 0.1)), col = 3, lwd = 2)
  }
}
legend("topright",
       c("Density", "N(0,1)"),
       lty = c(1, 1), lwd = c(2, 2), col = c(2, 3))



mirror_stat_cutoff(mirror_stat[, 1], fdr = 0.1)
mirror_stat_cutoff(mirror_stat[, 2], fdr = 0.1)
mirror_stat_cutoff(mirror_stat[, 3], fdr = 0.1)
mirror_stat_cutoff(mirror_stat[, 4], fdr = 0.1)

mirror_stat_cutoff(mirror_stat[, 3], fdr = 0.2)
mirror_stat_cutoff(mirror_stat[, 4], fdr = 0.2)

# Selected variables
cutoff <- mirror_stat_cutoff(mirror_stat[, 3], fdr = 0.1)
sum(mirror_stat[, 3] > cutoff)

cutoff <- mirror_stat_cutoff(mirror_stat[, 3], fdr = 0.2)
sum(mirror_stat[, 3] > cutoff)

cutoff <- mirror_stat_cutoff(mirror_stat[, 4], fdr = 0.1)
sum(mirror_stat[, 4] > cutoff)

cutoff <- mirror_stat_cutoff(mirror_stat[, 4], fdr = 0.2)
sum(mirror_stat[, 4] > cutoff)


# Histogram with FDR cutoff
par(mfrow = c(1, 2))
i <- 3
hist(mirror_stat[, i], breaks = 15, probability = T, 
     xlab = latex2exp::TeX("$F^{-1} (p)$"),
     main = colnames(mirror_stat)[i])   
lines(density(mirror_stat[, i]), col = 2, lwd = 2)
abline(v = 0, col = 1, lwd = 3)
abline(v = mirror_stat_cutoff(mirror_stat[, i], fdr = 0.1),
       col = 3, lwd = 3)
abline(v = mirror_stat_cutoff(mirror_stat[, i], fdr = 0.2),
       col = 4, lwd = 3)
i <- 4
hist(mirror_stat[, i], breaks = 15, probability = T, 
     xlab = latex2exp::TeX("$F^{-1} (p)$"),
     main = colnames(mirror_stat)[i])   
lines(density(mirror_stat[, i]), col = 2, lwd = 2)
abline(v = 0, col = 1, lwd = 3)
abline(v = mirror_stat_cutoff(mirror_stat[, i], fdr = 0.1),
       col = 3, lwd = 3)
abline(v = mirror_stat_cutoff(mirror_stat[, i], fdr = 0.2),
       col = 4, lwd = 3)
legend("topright",
       c("FDR = 0.1", "FDR = 0.2"),
       lty = c(1, 1), lwd = c(3, 3), col = c(3, 4))
