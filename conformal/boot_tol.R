library(tidyverse)
library(FlexCoDE)
source("R/utils.R")
# Parallel computation
library(doSNOW)
library(doRNG)
library(foreach)
library(progress)  # show progress bar in `foreach`


# Find the range of the interval with dis-connection
find_split_interval <- function(interval_value, min_grid_size = NULL) {
  if (is.null(min_grid_size)) {
    # It is used to find the positions of the dis-connected interval
    min_grid_size <- median(diff(interval_value))
  }
  
  if (length(interval_value) == 0) {
    # No covered interval
    return( list(c(0, 0)) )
  } else if (length(interval_value) == 1) {
    # Only 1 data contained => make region
    return( list(c(interval_value - min_grid_size, 
                   interval_value + min_grid_size)) )
  }
  
  # Find connected components of bands
  dif <- round(diff(interval_value), 5)
  split_idx <- which(dif > min_grid_size)
  if (length(split_idx) > 0) {
    interval_split <- list()
    for (i in 1:length(split_idx)) {
      if (i == 1) {
        interval_split[[1]] <- interval_value[1:split_idx[1]]
      } else {
        interval_split[[i]] <- interval_value[(split_idx[i-1]+1):split_idx[i]]
      }
    }
    interval_split[[length(split_idx)+1]] <- interval_value[(split_idx[length(split_idx)]+1):length(interval_value)]
  } else {
    interval_split <- list(interval_value)
  }
  
  # Obtain the range of each band
  interval_split <- lapply(interval_split, range)
  
  return(interval_split)
}


#' Prediction Interval using Bootstrap
#' 
#' @param x_train n x p matrix
#' @param y_train n size vetor
#' @param x_test m x p matrix which is the test data
#' @param cde the type of the conditional density estimation. Defualt is "flexcode"
#' @param alpha 1 - proportion of the sampled population
#' @param delta 1 - confidence level (coverage level)
#' @param B the number of bootstrap
#' @param n_cores the number of cores for `foreach`
#' @param ... additional options for `FlexCoDE::fitFlexCoDE`
boot_cde_tol <- function(x_train, y_train, x_test, 
                         cde = "flexcode", alpha = 0.1, delta = 0.05, 
                         B = 100, n_cores = 1, seed = NULL, ...) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Validation set for FlexCode
  valid_prop <- 0.2  # validation set proportion
  idx_valid <- sample(1:nrow(x_train), size = round(nrow(x_train)*valid_prop))
  
  # Estimate conditional density function of training data (\hat{f}(y|x))
  fit <- FlexCoDE::fitFlexCoDE(
    xTrain = x_train[-idx_valid, ],
    zTrain = y_train[-idx_valid],
    xValidation = x_train[idx_valid, ],
    zValidation = y_train[idx_valid],
    ...
  )
  cde_train <- FlexCoDE::predict.FlexCoDE(fit, x_train)
  
  # CDF of conditional density for sampling
  cdf_train <- apply(cde_train$CDE, 1, function(x){ cumsum(x) / sum(x) })
  cdf_train <- t(cdf_train)
  
  # CDE and CDF for test data
  cde_test <- FlexCoDE::predict.FlexCoDE(fit, x_test)
  cdf_test <- apply(cde_test$CDE, 1, function(x){ cumsum(x) / sum(x) })
  cdf_test <- t(cdf_test)
  
  # par(mfrow = c(1, 2))
  # i <- 10
  # plot(cde_train$z, cdf_train[i, ], type = "l", main = "CDF")
  # plot(cde_train$z, cde_train$CDE[i, ], type = "l", main = "pdf (cde)")
  # lines(density(replicate(1000, cde_train$z[findInterval(runif(1), cdf_train[i, ])+1])), col = 2, lwd = 2)
  
  
  ### Bootstrap
  # Parallel computation
  cl <- makeCluster(n_cores)
  registerDoSNOW(cl)
  
  # Progress bar for `foreach`
  pb <- progress_bar$new(
    format = "[:bar] :percent [Execution time :elapsedfull || Estimated time remaining: :eta]",
    total = B,   # total number of ticks to complete (default 100)
    clear = FALSE,   # whether to clear the progress bar on completion (default TRUE)
    width = 80   # width of the progress bar
  )
  progress <- function(n){
    pb$tick()
  } 
  opts <- list(progress = progress)
  
  # For replication of `foreach`
  if (!is.null(seed)) {
    registerDoRNG(seed)
  }
  
  # Bootstrap using `foreach`
  n_test <- nrow(x_test)
  delta <- foreach(i = 1:B, .combine = cbind, .options.snow = opts) %dopar% {
    # Sample from CDE \hat{f}(y|x)
    y_star <- apply(cdf_train, 1, function(cdf_i){
      replicate(1, cde_train$z[findInterval(runif(1), cdf_i)+1])
    })
    y_star_test <- apply(cdf_test, 1, function(cdf_i){
      replicate(1, cde_test$z[findInterval(runif(1), cdf_i)+1])
    })
    
    # Re-estimate CDE (\hat{f}^*(y|x)) for y_star_test
    fit_boot <- FlexCoDE::fitFlexCoDE(
      xTrain = x_train[-idx_valid, ],
      zTrain = y_star[-idx_valid],
      xValidation = x_train[idx_valid, ],
      zValidation = y_star[idx_valid],
      ...
    )
    cde_test_boot <- FlexCoDE::predict.FlexCoDE(fit_boot, x_test)
    
    # Compute log likelihood ratio
    delta_b <- sapply(1:n_test, function(i) {
      # Obtain CDE of y_star_test
      f_hat_y <- mean( cde_test_boot$CDE[i, order(abs(cde_test_boot$z - y_star_test[i]))[1:2]] )
      f_hat_y <- ifelse(f_hat_y == 0, 1e-6, f_hat_y)
      
      return( log( max(cde_test_boot$CDE[i, ]) / f_hat_y) )
    })
    
    return(delta_b)
  }
  
  # End parallel backend
  stopCluster(cl) 
  
  # # Bootstrap using for loop
  # delta <- matrix(NA, n_test, B)
  # for (b in 1:B) {
  #   # Show the progress bar
  #   if (b %% 10 == 0) {
  #     progress(b, max = B)
  #   }
  #   
  #   # Sample from CDE \hat{f}(y|x)
  #   y_star <- apply(cdf_train, 1, function(cdf_i){
  #     replicate(1, cde_train$z[findInterval(runif(1), cdf_i)+1])
  #   })
  #   y_star_test <- apply(cdf_test, 1, function(cdf_i){
  #     replicate(1, cde_test$z[findInterval(runif(1), cdf_i)+1])
  #   })
  #   
  #   # Re-estimate CDE (\hat{f}^*(y|x)) for y_star_test
  #   fit_boot <- FlexCoDE::fitFlexCoDE(
  #     xTrain = x_train[-idx_valid, ],
  #     zTrain = y_star[-idx_valid],
  #     xValidation = x_train[idx_valid, ],
  #     zValidation = y_star[idx_valid],
  #     ...
  #   )
  #   cde_test_boot <- FlexCoDE::predict.FlexCoDE(fit_boot, x_test)
  #   
  #   # Compute log likelihood ratio
  #   delta[, b] <- sapply(1:n_test, function(i) {
  #     # Obtain CDE of y_star_test
  #     f_hat_y <- mean( cde_test_boot$CDE[i, order(abs(cde_test_boot$z - y_star_test[i]))[1:2]] )
  #     f_hat_y <- ifelse(f_hat_y == 0, 1e-6, f_hat_y)
  #     
  #     return( log( max(cde_test_boot$CDE[i, ]) / f_hat_y) )
  #   })
  # }
  
  
  # 1-alpha quantiles of \delta_b
  c_star <- apply(delta, 1, quantile, 1-alpha)
  
  # Cufoffs of CDE for test data
  cutoff <- apply(cde_test$CDE, 1, max) / exp(c_star)
  
  # i <- 2
  # plot(cde_test$z, cde_test$CDE[i, ], type = "l")
  # abline(h = cutoff[i], col = 2)
  
  
  # Split the dis-connected region
  grid_size_density <- round(diff(cde_test$z)[1], 5)   # equal grid size of density estimate
  tol_band <- lapply(1:n_test, function(i){
    tol_band_grid <- cde_test$z[which(cde_test$CDE[i, ] > cutoff[i])]
    find_split_interval(tol_band_grid, grid_size_density)
  })
  
  out <- list(
    cde_fit = fit,
    tol_band = tol_band,
    cde_train = cde_train,
    cde_test = cde_test
  )
  
  return(out)  
}


# Result summary function for gallaxy red-shift data
summary_tol_band <- function(tol_band, y_test, idx_bright, idx_faint) {
  df <- data.frame(
    Total = c(
      Coverage = mapply(function(tol_band_i, y_i){
        sapply(tol_band_i, function(interval){
          y_i >= interval[1] & y_i <= interval[2]
        }) %>% 
          sum()
      }, tol_band, y_test) %>% 
        mean(),
      Size = sapply(tol_band, function(tol_band_i){
        sapply(tol_band_i, function(interval){
          diff(interval)
        }) %>% 
          sum()
      }) %>%
        mean()
    ),
    Bright = c(
      Coverage = mapply(function(tol_band_i, y_i){
        sapply(tol_band_i, function(interval){
          y_i >= interval[1] & y_i <= interval[2]
        }) %>% 
          sum()
      }, tol_band[idx_bright], y_test[idx_bright]) %>% 
        mean(),
      Size = sapply(tol_band[idx_bright], function(tol_band_i){
        sapply(tol_band_i, function(interval){
          diff(interval)
        }) %>% 
          sum()
      }) %>%
        mean()
    ),
    Faint = c(
      Coverage = mapply(function(tol_band_i, y_i){
        sapply(tol_band_i, function(interval){
          y_i >= interval[1] & y_i <= interval[2]
        }) %>% 
          sum()
      }, tol_band[idx_faint], y_test[idx_faint]) %>% 
        mean(),
      Size = sapply(tol_band[idx_faint], function(tol_band_i){
        sapply(tol_band_i, function(interval){
          diff(interval)
        }) %>% 
          sum()
      }) %>%
        mean()
    )
  )
  return(df)
}





#######################################
### Gallaxy data
#######################################
data <- read.table("data/photoz_catalogues/Happy/happy_A", skip = 1, header = F)
dim(data)
colnames(data) <- c("id", "mag_r", "u-g", "g-r", "r-i", "i-z", "z_spec",
                    "feat1", "feat2", "feat3", "feat4", "feat5")
head(data)
y <- data$z_spec   # spectroscopic redshift
x <- data[, c("mag_r", "feat1", "feat2", "feat3", "feat4", "feat5")] %>%
  as.matrix()
n <- nrow(x)

# Split data
set.seed(1000)
n_test <- 5000
n_train <- 500
n_calib <- 500
# n_test <- round(n/2)
# n_train <- round((n - n_test)/2)
# n_calib <- n - n_test - n_train
idx_test <- sample(1:n, n_test)
idx_calib <- sample(setdiff(1:n, idx_test), n_calib)
# idx_train <- setdiff(1:n, c(idx_test, idx_calib))
idx_train <- sample(setdiff(1:n, c(idx_test, idx_calib)), n_train)

x_train <- x[idx_train, ]
y_train <- y[idx_train]

x_calib <- x[idx_calib, ]
y_calib <- y[idx_calib]

x_test <- x[idx_test, ]
y_test <- y[idx_test]

x <- rbind(x_train, x_calib)
y <- c(y_train, y_calib)

n <- nrow(x)

# Bright and faint for red-shift
idx_bright <- which(x_test[, 1] < median(x_test[, 1]))
idx_faint <- which(x_test[, 1] >= median(x_test[, 1]))


#######################################
### Triceps data
#######################################
triceps <- MultiKink::triceps
dim(triceps)
head(triceps)

# with(triceps, plot(age, lntriceps, ylim = c(0, 5)))
# lines(smooth.spline(triceps$age, triceps$lntriceps, cv = T), col = 2, lwd = 3)


# Training + validation set
x <- matrix(triceps$age, ncol = 1)
y <- triceps$lntriceps
n <- nrow(x)

# Split data
set.seed(1000)
n_test <- round(n/2)
n_train <- round((n - n_test)/2)
n_calib <- n - n_test - n_train
idx_test <- sample(1:n, n_test)
idx_calib <- sample(setdiff(1:n, idx_test), n_calib)
idx_train <- setdiff(1:n, c(idx_test, idx_calib))

x_train <- matrix(x[idx_train, ], ncol = 1)
y_train <- y[idx_train]

x_calib <- matrix(x[idx_calib, ], ncol = 1)
y_calib <- y[idx_calib]

x_test <- matrix(x[idx_test, ], ncol = 1)
y_test <- y[idx_test]

x <- rbind(x_train, x_calib)
y <- c(y_train, y_calib)

n <- nrow(x)

# # It uses for visualize tolerance region
# x_test <- matrix(seq(min(x), max(x), length.out = 200), ncol = 1)



#######################################
### Estimate Tolerance Interval
#######################################

# Tolerance level
alpha <- 0.1    # 1 - proportion of the sampled population (guarantee level)
delta <- 0.05   # 1 - confidence level (coverage level)

# Tolerance interval via Bootstrap for CDE
seed <- 1
fit_boot_cde_tol <- boot_cde_tol(x_train, y_train, x_test, 
                                 cde = "flexcode", alpha = 0.1, delta = 0.05, 
                                 B = 100,
                                 n_cores = 10,
                                 seed = seed,
                                 # Additional `FlexCode` options
                                 nIMax = 30,
                                 regressionFunction = FlexCoDE::regressionFunction.NW)
tol_band <- fit_boot_cde_tol$tol_band
table(sapply(tol_band, length))

# Result summary for gallaxy red-shift data
summary_tol_band(tol_band, y_test, idx_bright, idx_faint)
#             Total    Bright     Faint
# NN
# Coverage 0.8968000 0.9276000 0.8660000
# Size     0.1765288 0.1363154 0.2167421
# NW
# Coverage 0.9486000 0.9864000 0.9108000
# Size     0.3327457 0.3083176 0.3571738
# SpAM
# Coverage 0.8716000 0.929600 0.813600
# Size     0.1560755 0.139696 0.172455
# Coverage 0.8920000 0.9296000 0.8544000
# Size     0.1697641 0.1501049 0.1894233
# XGBoost
# Coverage 0.8082000 0.8796000 0.7368000
# Size     0.1259231 0.1011171 0.1507291


system.time({
  # Validation set for FlexCode
  valid_prop <- 0.2  # validation set proportion
  idx_valid <- sample(1:nrow(x_train), size = round(nrow(x_train)*valid_prop))
  
  fit_boot <- fitFlexCoDE(
    xTrain = x_train[-idx_valid, ],
    zTrain = y_train[-idx_valid],
    xValidation = x_train[idx_valid, ],
    zValidation = y_train[idx_valid],
    regressionFunction = FlexCoDE::regressionFunction.NNKernel
  )
  # Forest, XGBoost, SpAM, SDMKernel만 parallel 할 수 있음
  # 근데 core 수 늘려도 속도가 오히려 느려짐
})
# n_train = 500 기준
#   user  system elapsed 
# NN
# 0.653   0.121   0.800 
# NW
# 1.145   0.056   1.199 
# SpAM
# 0.672   0.061   3.082 
# Series
# 18.099   4.560  22.714 
# Forest
# 0.831   0.353  35.890 
# Lasso
# 1.279   0.074   1.422
# XGBoost
# 0.690   0.190   7.062 



# Visualization
y_eval <- seq(min(y), max(y), length.out = 500)
df1 <- data.frame(x = x,
                  y = y)
for (i in 1:length(tol_band)) {
  sub <- cbind(
    x_test[i],
    unlist(
      lapply(tol_band[[i]], function(tol_band_ij){
        y_eval[y_eval >= tol_band_ij[1] & y_eval <= tol_band_ij[2]]
      }) 
    )
  )
  if (i == 1) {
    df2 <- sub
  } else {
    df2 <- rbind(df2, sub)
  }
}
df2 <- data.frame(df2)
colnames(df2) <- c("x","y")
# head(df2)

ggplot() +
  geom_point(data = df1, aes(x = x, y = y)) +
  geom_point(data = df2, aes(x = x, y = y), alpha = 0.02, color = "red") +
  theme_bw() + 
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        legend.title = element_blank(),
        plot.title = element_text(size=25, face="bold"))








#######################################
### Simulation
#######################################
y <- data$z_spec   # spectroscopic redshift
x <- data[, c("mag_r", "feat1", "feat2", "feat3", "feat4", "feat5")] %>% 
  as.matrix()

n <- nrow(x)

# Number of cores
n_cores <- 10

num_sim <- 30
res_boot_tol <- array(NA, dim = c(2, 3, num_sim))
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
  
  # Tolerance interval via Bootstrap for CDE
  x_sim <- rbind(x_train, x_calib)
  y_sim <- c(y_train, y_calib)
  start_time <- Sys.time()
  fit_boot_cde_tol <- boot_cde_tol(x_sim, y_sim, x_test, 
                                   cde = "flexcode", alpha = 0.1, delta = 0.05, 
                                   B = 100,
                                   n_cores = n_cores,
                                   seed = seed,
                                   # Additional `FlexCode` options
                                   nIMax = 30,
                                   regressionFunction = FlexCoDE::regressionFunction.NW)
  end_time <- Sys.time()
  print(end_time - start_time)
  res <- summary_tol_band(fit_boot_cde_tol$tol_band, y_test, idx_bright, idx_faint)
  res_boot_tol[, , i] <- as.matrix(res)
  print(res)
}
# save(res_rcps_flexcode, res_rcps_mdn, res_npreg_tol, file = "RData/gallaxy_all.RData")

# load("RData/gallaxy_rcps_flexcode_NW.RData")
df <- paste0(
  paste(
    format(round(apply(res_boot_tol, c(1,2), mean), 3), 3), 
    format(round(apply(res_boot_tol, c(1,2), sd), 3), 3), 
    sep = " ("
  ), 
  ")") %>% 
  matrix(2, 3) %>% 
  as.data.frame()
colnames(df) <- c("Total","Bright","Faint")
rownames(df) <- c("Coverage","Size")
df