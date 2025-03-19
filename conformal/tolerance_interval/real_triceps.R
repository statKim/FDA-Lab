library(tidyverse)
source("R/RCPS.R")


# # Simulated data
# alpha <- 0.2
# delta <- 0.1
# set.seed(100)
# n <- 500
# x <- matrix(runif(n, -5, 5), ncol = 1)
# y <- 5*x + rnorm(n, 0, 1+abs(x))
# # plot(x, y)
# 
# # # bi-modal data
# # x <- matrix(runif(1000), ncol = 1)
# # y <- c(x[1:500]^2, -x[501:1000]^2) + rnorm(1000, 0, 0.2)
# # plot(x, y)
# # x_test <- matrix(seq(0, 1, length.out = 500), ncol = 1)
# 
# n_test <- 500
# x_test <- matrix(seq(-5, 5, length.out = n_test), ncol = 1)
# # x_test <- matrix(runif(n_test, -5, 5), ncol = 1)
# # # y_test <- c(x[1:500]^2, -x[501:1000]^2) + rnorm(1000, 0, 0.2)
# # y_test <- 5*x_test + rnorm(n_test, 0, 1+abs(x_test))


# Triceps data
triceps <- MultiKink::triceps
dim(triceps)
head(triceps)

with(triceps, plot(age, lntriceps, ylim = c(0, 5)))
lines(smooth.spline(triceps$age, triceps$lntriceps, cv = T), col = 2, lwd = 3)

n <- nrow(triceps)

# Tolerance level
alpha <- 0.1    # 1 - proportion of the sampled population
delta <- 0.05   # 1 - confidence level (coverage level)

# Training + validation set
set.seed(100)
x <- matrix(triceps$age, ncol = 1)
y <- triceps$lntriceps

# Test set for visualize tolerance region
n_test <- 300
x_test <- matrix(seq(min(x), max(x), length.out = n_test), ncol = 1)


# Conditional density estimates for the each neighborhood A_i in training data
# It also contains the calibration estimate.
cd_split_fit <- fit_predictionBands(x, y, 
                                    per_train = 0.4,
                                    per_val = 0.1,
                                    per_ths = 0.5, 
                                    regressionFunction = FlexCoDE::regressionFunction.NW)
# pred <- predict(cd_split_fit, x_test, type = "cd")
# # plot(pred$y_grid, pred$densities[80, ])

# Predict CDE for test data
pred_test <- FlexCoDE::predict.FlexCoDE(cd_split_fit$density_fit, x_test)

# Profile density estimate for test data
g_test <- matrix(NA, nrow(x_test), length(cd_split_fit$t_grid))
for(i in 1:nrow(x_test)){
  g_test[i, ] <- profile_density(cd_split_fit$t_grid,
                                 pred_test$z,
                                 pred_test$CDE[i, ])
}
# Find A_i which belongs to new data X_{n+1}
# It uses k nearest neighbors in data of each query
which_partition_test <- FNN::get.knnx(data = cd_split_fit$centers_kmeans,
                                      query = g_test,
                                      k = 1)$nn.index
which_partition_train <- FNN::get.knnx(data = cd_split_fit$centers_kmeans,
                                       query = cd_split_fit$g_train,
                                       k = 1)$nn.index
# length(which_partition_test)
# length(which_partition_train)

# par(mfrow = c(1, 2))
# plot(x_test, y_test, col = which_partition_test)
plot(x[cd_split_fit$splits == "Threshold"], 
     y[cd_split_fit$splits == "Threshold"], col = which_partition_train)

triceps[cd_split_fit$splits == "Threshold", c("age","lntriceps")] %>% 
  mutate(group = factor(which_partition_train)) %>% 
  ggplot(aes(x = age, y = lntriceps, color = group)) + 
  geom_point() +
  scale_y_continuous(limits = c(0, 4.5)) +
  theme_bw() +
  labs(x = "Age", y = "Triceps") +
  theme(legend.position = "none")


# Calibration set - "Threshold" indicates calibration set
idx_calib <- which(cd_split_fit$splits == "Threshold")
x_calib <- x[idx_calib]
y_calib <- y[idx_calib]


# Predict CDE for calibration data
pred_calib <- FlexCoDE::predict.FlexCoDE(cd_split_fit$density_fit, x_calib)


# Find optimal lambda for nested set-valued predcitor from calibration set
lambda_list <- -seq(1e-8, min(apply(pred_calib$CDE, 1, max))-(1e-8), length.out = 20)  # candidates of lambda
R_hat_lambda <- rep(0, length(lambda_list))
for (j in 1:length(which_partition_train)) {
  # Calibration set belongs to A_i
  idx_calib_in_A_i <- which(which_partition_train == which_partition_train[j])
  y_calib_A_i <- y_calib[idx_calib_in_A_i]
  
  # CDE for i-th test data
  f_hat <- data.frame(x = pred_calib$z,
                      y = pred_calib$CDE[j, ])
  
  # Nested set-valued predictors
  # Do not need to use Greedy algorithm for tolerance interval
  T_lambda <- lapply(lambda_list, function(lambda){
    sort(f_hat$x[f_hat$y > -lambda])
  })
  
  # # Nested set-valued predictors using Greedy algorithm
  # T_lambda <- list()
  # d <- 0.1   # step size rate
  # zeta <- 3  # bound of random variables in Hoeffding's inequality
  # for (i in 1:length(lambda_list)) {
  #   lambda <- lambda_list[i]
  #   T_lambda_i <- c()
  #   
  #   # iter <- 0
  #   zeta_update <- zeta
  #   while (zeta_update > -lambda) {
  #     # iter <- iter + 1
  #     zeta_update <- zeta_update - d*zeta_update
  #     T_lambda_i <- c(T_lambda_i, 
  #                     f_hat$x[(f_hat$y > zeta_update) & !(f_hat$x %in% T_lambda_i)])
  #   }
  #   T_lambda[[i]] <- sort(T_lambda_i)
  #   # print(iter)
  # }
  # # sapply(T_lambda, length)
  
  
  # Split the dis-connected region
  grid_size_density <- round(diff(f_hat$x)[1], 5)   # equal grid size of density estimate
  T_lambda_split <- lapply(T_lambda, function(tol_band){
    # Only 1 data contained => make region
    if (length(tol_band) < 2) {
      return( list(c(tol_band - grid_size_density, 
                     tol_band + grid_size_density)) )
    }
    
    # Find connected components of bands
    dif <- round(diff(tol_band), 5)
    split_idx <- which(dif > grid_size_density)
    if (length(split_idx) > 0) {
      tol_band_split <- list()
      for (i in 1:length(split_idx)) {
        if (i == 1) {
          tol_band_split[[1]] <- tol_band[1:split_idx[1]]
        } else {
          tol_band_split[[i]] <- tol_band[(split_idx[i-1]+1):split_idx[i]]
        }
      }
      tol_band_split[[length(split_idx)+1]] <- tol_band[(split_idx[length(split_idx)]+1):length(tol_band)]
    } else {
      tol_band_split <- list(tol_band)
    }
    
    # Obtain the range of each band
    tol_band_split <- lapply(tol_band_split, range)
    
    return(tol_band_split)
  })
  
  # # Visualization
  # par(mfrow = c(1, 2))
  # plot(f_hat, type = "l")
  # for (i in 1:length(lambda_list)) {
  #   tol_band_split <- T_lambda_split[[i]]
  #   for (j in 1:length(tol_band_split)) {
  #     lines(tol_band_split[[j]], rep(-lambda_list[[i]], 2), col = i, lwd = 3)
  #   }
  # }
  # plot(f_hat[which(f_hat$y > 0), ], type = "l")
  
  
  # Empirical risk
  R_hat_lambda <- R_hat_lambda + sapply(1:length(T_lambda_split), function(i){
    number <- sapply(T_lambda_split[[i]], function(interval){
      (y_calib[j] >= interval[1] & y_calib[j] <= interval[2])
    })
    # (1 - sum(number)) / length(y_calib_A_i)
    (1 - sum(number)) / length(y_calib)
  })
}
R_hat_lambda

## Find UCB using the concentration inequality
# Bound of concentration inequality
bound <- sqrt(log(1/delta)/(2*length(y_calib)))  # Hoeffding's inequality

# Find lambda_hat (We use max since lambda has decreasing order!)
ucb_rule <- which(R_hat_lambda < alpha - bound)
if (length(ucb_rule) == 0) {
  lambda_hat_order <- NA
} else {
  lambda_hat_order <- max(ucb_rule)
}
lambda_hat <- lambda_list[lambda_hat_order]
lambda_hat


# Risk-controliling Prediction Set (RCPS)
rcps <- list()
for (j in 1:length(which_partition_test)) {
  # Calibration set belongs to A_i
  idx_calib_in_A_i <- which(which_partition_train == which_partition_test[j])
  y_calib_A_i <- y_calib[idx_calib_in_A_i]
  
  # CDE for i-th test data
  f_hat <- data.frame(x = pred_test$z,
                      y = pred_test$CDE[j, ])
  
  # Nested set-valued predictors with lambda_hat
  # Do not need to use Greedy algorithm for tolerance interval
  T_lambda <- sort(f_hat$x[f_hat$y > -lambda_hat])
  
  # # Nested set-valued predictors using Greedy algorithm with lambda_hat
  # d <- 0.1   # step size rate
  # zeta <- 3  # bound of random variables in Hoeffding's inequality
  # T_lambda <- c()
  # # iter <- 0
  # zeta_update <- zeta
  # while (zeta_update > -lambda_hat) {
  #   # iter <- iter + 1
  #   zeta_update <- zeta_update - d*zeta_update
  #   T_lambda <- c(T_lambda, 
  #                 f_hat$x[(f_hat$y > zeta_update) & !(f_hat$x %in% T_lambda)])
  # }
  # T_lambda <- sort(T_lambda)
  # # print(iter)
  # # length(T_lambda)
  
  # Find connected components of bands
  dif <- round(diff(T_lambda), 5)
  split_idx <- which(dif > min(dif))
  if (length(split_idx) > 0) {
    T_lambda_split <- list()
    for (i in 1:length(split_idx)) {
      if (i == 1) {
        T_lambda_split[[1]] <- T_lambda[1:split_idx[1]]
      } else {
        T_lambda_split[[i]] <- T_lambda[(split_idx[i-1]+1):split_idx[i]]
      }
    }
    T_lambda_split[[length(split_idx)+1]] <- T_lambda[(split_idx[length(split_idx)]+1):length(T_lambda)]
  } else {
    T_lambda_split <- list(T_lambda)
  }
  # Obtain the range of each band
  T_lambda_split <- lapply(T_lambda_split, range)  
  
  
  # RCPS with lambda_hat
  rcps[[j]] <- T_lambda_split
}
# rcps

# names(rcps) <- round(x_test, 3)
# rcps



# Visualization
df1 <- data.frame(x = x_calib,
                  y = y_calib)
for (i in 1:length(rcps)) {
  sub <- cbind(
    x_test[i],
    as.numeric(
      sapply(rcps[[i]], function(rcps_ij){  
        seq(rcps_ij[1], rcps_ij[2], length.out = 100)
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
  # scale_y_continuous(limits = c(0, 4.5)) +
  # labs(x = "Age", y = "Triceps") +
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        legend.title = element_blank(),
        plot.title = element_text(size=25, face="bold"))

