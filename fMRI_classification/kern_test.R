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

# Remove 2 outlying curves
X <- X[-c(33, 123), , ]
y <- y[-c(33, 123)]

n <- dim(X)[1]   # number of curves
m <- dim(X)[2]   # number of timepoints
p <- dim(X)[3]  # number of functional variables
gr <- seq(0, 1, length.out = m)



##################################################
### Homogeneity test for each region
### - Wynne, G., & Duncan, A. B. (2022). A kernel two-sample test for functional data. Journal of Machine Learning Research, 23(73), 1-51.
##################################################
# Load python functions
library(reticulate)
use_condaenv("reticulate")   # Use virtual environment
py_run_file("py/kern_test.py")
np <- import("numpy")

source("R/homo_test.R")

gamma <- -1  # bandwidth; Set bandwidth parameter, -1 causes median heuristic to be used
n_perms <- 1000   # number of permutations
test_obj <- list()
for (i in 1:p) {
  print(paste(i, "th region"))
  
  set.seed(i)
  X1 <- X[y == 0, , i]
  X2 <- X[y == 1, , i]
  
  test_obj[[i]] <- list()
  
  # print("SE-T_ID")
  test_obj[[i]]$SET_ID <- py$two_sample_test(X1, X2, gamma, np$int16(n_perms), 
                                             z_alpha = 0.05, 
                                             make_K = py$K_ID, 
                                             return_p = T, seed = np$int16(i))
  
  # print("Cov")
  test_obj[[i]]$COV <- py$two_sample_test(X1, X2, gamma, np$int16(n_perms), 
                                          z_alpha = 0.05, 
                                          make_K = py$K_COV, 
                                          return_p = T, seed = np$int16(i))
  
  # print("SE-T_FPCA")
  test_obj[[i]]$SET_FPCA <- py$two_sample_test(X1, X2, gamma, np$int16(n_perms), 
                                               z_alpha = 0.05, 
                                               make_K = py$K_FPCA, 
                                               return_p = T, seed = np$int16(i))
  
  # print("SE-T_SQR")
  test_obj[[i]]$SET_SQR <- py$two_sample_test(X1, X2, gamma, np$int16(n_perms), 
                                              z_alpha = 0.05, 
                                              make_K = py$K_SQR, 
                                              return_p = T, seed = np$int16(i))
}
save(test_obj, file = "RData/homo_test_kern.RData")


###############################################################
### Variable selection using mirror statistics controlling FDR
###############################################################

# p-values
p_value <- data.frame(
  SET_ID = sapply(test_obj, function(x){ x$SET_ID[[2]] }),
  COV = sapply(test_obj, function(x){ x$COV[[2]] }),
  SET_FPCA = sapply(test_obj, function(x){ x$SET_FPCA[[2]] }),
  SET_SQR = sapply(test_obj, function(x){ x$SET_SQR[[2]] })
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
# cutoff <- mirror_stat_cutoff(mirror_stat[, 3], q = 0.2)
# sum(mirror_stat[, 3] > cutoff)

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

