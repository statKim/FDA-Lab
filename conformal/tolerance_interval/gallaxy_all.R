library(tidyverse)
source("R/RCPS.R")
source("R/npreg_tol.R")

# Load python sources
library(reticulate)
use_condaenv("torch")   # Use virtual environment
py_run_file("py/gaussian_mixture_density.py")


# Tolerance level
alpha <- 0.1    # 1 - proportion of the sampled population (guarantee level)
delta <- 0.05   # 1 - confidence level (coverage level)

# Gallaxy data
data <- read.table("data/photoz_catalogues/Happy/happy_A", skip = 1, header = F)
dim(data)
colnames(data) <- c("id", "mag_r", "u-g", "g-r", "r-i", "i-z", "z_spec", 
                    "feat1", "feat2", "feat3", "feat4", "feat5")
head(data)
y <- data$z_spec   # spectroscopic redshift
x <- data[, c("mag_r", "feat1", "feat2", "feat3", "feat4", "feat5")] %>% 
  as.matrix()

n <- nrow(x)

num_sim <- 30
res_rcps_flexcode <- array(NA, dim = c(2, 3, num_sim))
res_rcps_mdn <- array(NA, dim = c(2, 3, num_sim))
res_npreg_tol <- array(NA, dim = c(2, 3, num_sim))
res_npreg_tol_pt <- array(NA, dim = c(2, 3, num_sim))
res_npreg_tol_scale <- array(NA, dim = c(2, 3, num_sim))
res_npreg_tol_pt_scale <- array(NA, dim = c(2, 3, num_sim))
for (i in 1:num_sim) {
  print(i)
  seed <- i

  # Split data
  set.seed(seed)
  idx_test <- sample(1:n, 10000)
  idx_train_calib <- sample(setdiff(1:n, idx_test), 5000)
  idx_train <- idx_train_calib[1:2500]
  idx_calib <- setdiff(idx_train_calib, idx_train)
  
  x_train <- x[idx_train, ]
  y_train <- y[idx_train]
  
  x_calib <- x[idx_calib, ]
  y_calib <- y[idx_calib]
  
  x_test <- x[idx_test, ]
  y_test <- y[idx_test]
  
  
  # Bright and faint for red-shift
  idx_bright <- which(x_test[, 1] < median(x_test[, 1]))
  idx_faint <- which(x_test[, 1] >= median(x_test[, 1]))
  
  
  # Risk Controlling Prediction Set (RCPS) + CD-split
  # # Flexcode
  # print("RCPS-Flexcode")
  # start_time <- Sys.time()
  # fit_rcps_flexcode <- rcps_partition(x_train, y_train, x_calib, y_calib, x_test,
  #                                     cde = "flexcode", alpha = alpha, delta = delta, seed,
  #                                     regressionFunction = FlexCoDE::regressionFunction.Forest)
  # end_time <- Sys.time()
  # print(end_time - start_time)
  # res <- summary_rcps(fit_rcps_flexcode$rcps, y_test, idx_bright, idx_faint)
  # res_rcps_flexcode[, , i] <- as.matrix(res)
  # print(res)
  
  # # Mixture density network
  # print("RCPS-MDN")
  # start_time <- Sys.time()
  # fit_rcps_mdn <- rcps_partition(x_train, y_train, x_calib, y_calib, x_test,
  #                                cde = "mdn", alpha = alpha, delta = delta, seed,
  #                                n_components = 4, hidden_dim = 50)
  # end_time <- Sys.time()
  # print(end_time - start_time)
  # res <- summary_rcps(fit_rcps_mdn$rcps, y_test, idx_bright, idx_faint)
  # res_rcps_mdn[, , i] <- as.matrix(res)
  # print(res)
  # print(table(sapply(fit_rcps_mdn$rcps, length)))
  
  # # Nonparametric regression tolerance interval without scaling
  # print("npreg-tol")
  # x_sim <- rbind(x_train, x_calib)
  # y_sim <- c(y_train, y_calib)
  # set.seed(seed)
  # start_time <- Sys.time()
  # fit_npreg_tol <- npreg_tol(x_sim, y_sim, newdata = x_test, B = 1000, alpha = alpha, delta = delta)
  # end_time <- Sys.time()
  # print(end_time - start_time)
  # # Simultaneous
  # res <- summary_npreg_tol(fit_npreg_tol, y_test, type = "sim")
  # res_npreg_tol[, , i] <- as.matrix(res)
  # print(res)
  # # Pointwise
  # res <- summary_npreg_tol(fit_npreg_tol, y_test, type = "pt")
  # res_npreg_tol_pt[, , i] <- as.matrix(res)
  # print(res)
  
  # Nonparametric regression tolerance interval with scaling
  print("npreg-tol-scaling")
  x_sim <- rbind(x_train, x_calib)
  y_sim <- c(y_train, y_calib)
  set.seed(seed)
  start_time <- Sys.time()
  fit_npreg_tol <- npreg_tol(x_sim, y_sim, newdata = x_test, B = 1000,
                             alpha = alpha, delta = delta,
                             scale = T)
  end_time <- Sys.time()
  print(end_time - start_time)
  # Simultaneous
  res <- summary_npreg_tol(fit_npreg_tol, y_test, type = "sim")
  res_npreg_tol_scale[, , i] <- as.matrix(res)
  print(res)
  # Pointwise
  res <- summary_npreg_tol(fit_npreg_tol, y_test, type = "pt")
  res_npreg_tol_pt_scale[, , i] <- as.matrix(res)
  print(res)
}
# save(res_rcps_flexcode, file = "RData/gallaxy_rcps_flexcode_RF.RData")
# save(res_rcps_flexcode, file = "RData/gallaxy_rcps_flexcode_NW.RData")
# save(res_rcps_mdn, file = "RData/gallaxy_rcps_mdn.RData")
# save(res_npreg_tol, res_npreg_tol_pt, file = "RData/gallaxy_npreg_tol.RData")
# save(res_npreg_tol_scale, res_npreg_tol_pt_scale, file = "RData/gallaxy_npreg_tol_scale.RData")
# save(res_rcps_flexcode, res_rcps_mdn, res_npreg_tol, file = "RData/gallaxy_all.RData")

load("RData/gallaxy_rcps_flexcode_NW.RData")
load("RData/gallaxy_rcps_flexcode_RF.RData")
df <- paste0(
  paste(
    format(round(apply(res_rcps_flexcode, c(1,2), mean), 3), 3), 
    format(round(apply(res_rcps_flexcode, c(1,2), sd), 3), 3), 
    sep = " ("
  ), 
  ")") %>% 
  matrix(2, 3) %>% 
  as.data.frame()
colnames(df) <- c("Total","Bright","Faint")
rownames(df) <- c("Coverage","Size")
df
# NW
#                  Total        Bright         Faint
# Coverage 0.931 (0.007) 0.975 (0.011) 0.886 (0.012)
# Size     0.320 (0.026) 0.321 (0.028) 0.320 (0.026)
# RF
# Coverage 0.928 (0.007) 0.968 (0.004) 0.888 (0.012)
# Size     0.159 (0.009) 0.127 (0.005) 0.192 (0.014)


load("RData/gallaxy_rcps_mdn.RData")
df <- paste0(
  paste(
    format(round(apply(res_rcps_mdn, c(1,2), mean), 3), 3), 
    format(round(apply(res_rcps_mdn, c(1,2), sd), 3), 3),
    sep = " ("
  ), 
")") %>% 
  matrix(2, 3) %>% 
  as.data.frame()
colnames(df) <- c("Total","Bright","Faint")
rownames(df) <- c("Coverage","Size")
df
# hidden_dim = 30
#                    Total       Bright        Faint
# Coverage 0.968 (0.029) 0.980 (0.018) 0.957 (0.040)
# Size     0.292 (0.142) 0.215 (0.110) 0.369 (0.180)
# hidden_dim = 50
# Coverage 0.963 (0.030) 0.977 (0.018) 0.948 (0.041)
# Size     0.263 (0.143) 0.198 (0.110) 0.328 (0.177)
# 2 hidden layer; hidden_dim = 50
# Coverage 0.970 (0.025) 0.981 (0.016) 0.959 (0.035)
# Size     0.301 (0.163) 0.243 (0.156) 0.358 (0.179)

load("RData/gallaxy_npreg_tol.RData")
df <- paste0(
  paste(
    format(round(apply(res_npreg_tol, c(1,2), mean), 3), 3),
    format(round(apply(res_npreg_tol, c(1,2), sd), 3), 3),
    sep = " ("
  ), 
  ")") %>% 
  matrix(2, 3) %>% 
  as.data.frame()
colnames(df) <- c("Total","Bright","Faint")
rownames(df) <- c("Coverage","Size")
df
#                  Total        Bright         Faint
# Coverage 0.993 (0.002) 0.999 (0.001) 0.988 (0.004)
# Size     0.548 (0.050) 0.548 (0.050) 0.548 (0.050)
#                  Total        Bright         Faint
# Coverage 0.995 (0.001) 0.999 (0.001) 0.991 (0.002)
# Size     0.600 (0.028) 0.600 (0.028) 0.600 (0.028)
# Coverage 0.995 (0.001) 0.999 (0.001) 0.991 (0.002)
# Size     0.601 (0.028) 0.601 (0.028) 0.601 (0.028)

df <- paste0(
  paste(
    format(round(apply(res_npreg_tol_pt, c(1,2), mean), 3), 3),
    format(round(apply(res_npreg_tol_pt, c(1,2), sd), 3), 3),
    sep = " ("
  ), 
  ")") %>% 
  matrix(2, 3) %>% 
  as.data.frame()
colnames(df) <- c("Total","Bright","Faint")
rownames(df) <- c("Coverage","Size")
df
#                  Total        Bright         Faint
# Coverage 0.935 (0.005) 0.977 (0.004) 0.893 (0.007)
# Size     0.204 (0.007) 0.203 (0.007) 0.204 (0.007)
# Coverage 0.935 (0.005) 0.977 (0.004) 0.893 (0.007)
# Size     0.204 (0.007) 0.203 (0.007) 0.204 (0.007)


load("RData/gallaxy_npreg_tol_scale.RData")
df <- paste0(
  paste(
    format(round(apply(res_npreg_tol_scale, c(1,2), mean), 3), 3),
    format(round(apply(res_npreg_tol_scale, c(1,2), sd), 3), 3),
    sep = " ("
  ), 
  ")") %>% 
  matrix(2, 3) %>% 
  as.data.frame()
colnames(df) <- c("Total","Bright","Faint")
rownames(df) <- c("Coverage","Size")
df
#                  Total        Bright         Faint
# Coverage 0.996 (0.001) 0.998 (0.001) 0.994 (0.002)
# Size     0.645 (0.075) 0.518 (0.059) 0.773 (0.092)

df <- paste0(
  paste(
    format(round(apply(res_npreg_tol_pt_scale, c(1,2), mean), 3), 3),
    format(round(apply(res_npreg_tol_pt_scale, c(1,2), sd), 3), 3),
    sep = " ("
  ), 
  ")") %>% 
  matrix(2, 3) %>% 
  as.data.frame()
colnames(df) <- c("Total","Bright","Faint")
rownames(df) <- c("Coverage","Size")
df
#                  Total        Bright         Faint
# Coverage 0.927 (0.007) 0.961 (0.006) 0.893 (0.012)
# Size     0.219 (0.024) 0.176 (0.019) 0.263 (0.030)