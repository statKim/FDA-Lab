
library(tidyverse)
source("R/RCPS.R")
library(reticulate)
use_condaenv("torch")   # Use virtual environment

# Load python sources
py_run_file("py/gaussian_mixture_density.py")


# Tolerance level
alpha <- 0.1    # 1 - proportion of the sampled population
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

# Split data
set.seed(1000)
idx_test <- sample(1:nrow(x), 5000)
idx_calib <- sample(setdiff(1:nrow(x), idx_test), 5000)
idx_train <- setdiff(1:nrow(x), c(idx_test, idx_calib))

x_train <- x[idx_train, ]
y_train <- y[idx_train]

x_calib <- x[idx_calib, ]
y_calib <- y[idx_calib]

x_test <- x[idx_test, ]
y_test <- y[idx_test]

x <- rbind(x_train, x_calib) 
y <- c(y_train, y_calib)



# Mixture density network with 3 hidden layers
system.time({
  cd_split_fit <- mdn_pred_band(x_train, y_train, x_calib, y_calib, seed = 5,
                                n_components = 4, hidden_dim = 30)
})
# user   system  elapsed 
# 1004.735   34.873  148.326 
# Warning message:
# Quick-TRANSfer stage steps exceeded maximum (= 3247500) 


# Predict CDE for test data
pred_test <- list(
  z = cd_split_fit$pred_calib$z,
  CDE = matrix(NA, length(y_test), 1000)
)
# pred_test <- list()
for (i in 1:length(y_test)) {
  x_temp <- matrix(rep(x_test[i, ], 1000), ncol = ncol(x_test), byrow = T)
  pred_test$CDE[i, ] <- py$predict_mixture_density_network(cd_split_fit$density_fit,
                                                           x_temp, 
                                                           matrix(pred_test$z, ncol = 1))
}

save(x_train, x_calib, x_test, y_train, y_calib, y_test, cd_split_fit, pred_test,
     file = "RData/gallaxy_density_MDN.RData")

load("RData/gallaxy_density_MDN.RData")

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


# Predict CDE for calibration data
pred_calib <- cd_split_fit$pred_calib

# Find optimal lambda for nested set-valued predcitor from calibration set
lambda_max <- apply(pred_calib$CDE, 1, max)
print(lambda_max %>% sort() %>% head())
lambda_max <- min(lambda_max[lambda_max > 1e-4])
lambda_list <- -seq(1e-8, lambda_max-(1e-8), length.out = 20)  # candidates of lambda
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
  
  # Split the dis-connected region
  grid_size_density <- round(diff(f_hat$x)[1], 5)   # equal grid size of density estimate
  T_lambda_split <- lapply(T_lambda, function(tol_band){
    if (length(tol_band) == 0) {
      # No covered interval
      return( list(c(0, 0)) )
    } else if (length(tol_band) < 2) {
      # Only 1 data contained => make region
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
  # j <- 3071
  
  # Calibration set belongs to A_i
  idx_calib_in_A_i <- which(which_partition_train == which_partition_test[j])
  y_calib_A_i <- y_calib[idx_calib_in_A_i]
  
  # CDE for i-th test data
  f_hat <- data.frame(x = pred_test$z,
                      y = pred_test$CDE[j, ])
  
  # j <- 9
  # f_hat <- data.frame(x = pred_test$z,
  #                     y = pred_test$CDE[j, ])
  # plot(f_hat)
  # abline(h = -lambda_hat)
  
  # Nested set-valued predictors with lambda_hat
  # Do not need to use Greedy algorithm for tolerance interval
  T_lambda <- sort(f_hat$x[f_hat$y > -lambda_hat])
  if (length(T_lambda) < 1) {
    T_lambda <- c(0, 0)
  }
  
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
sapply(rcps, length) %>% table


# Result summaries
idx_bright <- which(x_test[, 1] < median(x_test[, 1]))
idx_faint <- which(x_test[, 1] >= median(x_test[, 1]))
data.frame(
  Total = c(
    Coverage = mapply(function(rcps_i, y_i){
      sapply(rcps_i, function(interval){
        y_i >= interval[1] & y_i <= interval[2]
      }) %>% 
        sum()
    }, rcps, y_test) %>% 
      mean(),
    Size = sapply(rcps, function(rcps_i){
      sapply(rcps_i, function(interval){
        diff(interval)
      }) %>% 
        sum()
    }) %>%
      mean()
  ),
  Bright = c(
    Coverage = mapply(function(rcps_i, y_i){
      sapply(rcps_i, function(interval){
        y_i >= interval[1] & y_i <= interval[2]
      }) %>% 
        sum()
    }, rcps[idx_bright], y_test[idx_bright]) %>% 
      mean(),
    Size = sapply(rcps[idx_bright], function(rcps_i){
      sapply(rcps_i, function(interval){
        diff(interval)
      }) %>% 
        sum()
    }) %>%
      mean()
  ),
  Faint = c(
    Coverage = mapply(function(rcps_i, y_i){
      sapply(rcps_i, function(interval){
        y_i >= interval[1] & y_i <= interval[2]
      }) %>% 
        sum()
    }, rcps[idx_faint], y_test[idx_faint]) %>% 
      mean(),
    Size = sapply(rcps[idx_faint], function(rcps_i){
      sapply(rcps_i, function(interval){
        diff(interval)
      }) %>% 
        sum()
    }) %>%
      mean()
  )
) %>% 
  round(4)

# seed = 10
#           Total Bright  Faint
# Coverage 0.9918 0.9948 0.9888
# Size     0.2264 0.1597 0.2932

# seed = 5
#           Total Bright  Faint
# Coverage 0.9436 0.9652 0.9220
# Size     0.1448 0.1050 0.1846

# ArXiv paper - Conformal result
#           Total Bright  Faint
# Coverage 0.900  0.906   0.895
# Size     0.131  0.082   0.181


### Different coverage levels
alpha_list <- c(0.05, 0.1, 0.15, 0.2)  # 1 - proportion of the sampled population
delta_list <- c(0.05, 0.1, 0.15, 0.2)  # 1 - confidence level (coverage level)
cand_level <- expand.grid(alpha = alpha_list,
                          delta = delta_list)
cand_level$coverage_total <- rep(NA, nrow(cand_level))
cand_level$coverage_bright <- rep(NA, nrow(cand_level))
cand_level$coverage_faint <- rep(NA, nrow(cand_level))
cand_level$size_total <- rep(NA, nrow(cand_level))
cand_level$size_bright <- rep(NA, nrow(cand_level))
cand_level$size_faint <- rep(NA, nrow(cand_level))
cand_level$lambda <- rep(NA, nrow(cand_level))

res_multimodal <- list()

for (cand in 1:nrow(cand_level)) {
  ## Find UCB using the concentration inequality
  # Bound of concentration inequality
  bound <- sqrt(log(1/cand_level$delta[cand])/(2*length(y_calib)))  # Hoeffding's inequality
  
  # Find lambda_hat (We use max since lambda has decreasing order!)
  ucb_rule <- which(R_hat_lambda < cand_level$alpha[cand] - bound)
  if (length(ucb_rule) == 0) {
    lambda_hat_order <- NA
  } else {
    lambda_hat_order <- max(ucb_rule)
  }
  lambda_hat <- lambda_list[lambda_hat_order]
  cand_level$lambda[cand] <- lambda_hat
  
  
  # Risk-controliling Prediction Set (RCPS)
  rcps <- list()
  for (j in 1:length(which_partition_test)) {
    # j <- 3071
    
    # Calibration set belongs to A_i
    idx_calib_in_A_i <- which(which_partition_train == which_partition_test[j])
    y_calib_A_i <- y_calib[idx_calib_in_A_i]
    
    # CDE for i-th test data
    f_hat <- data.frame(x = pred_test$z,
                        y = pred_test$CDE[j, ])
    
    # Nested set-valued predictors with lambda_hat
    # Do not need to use Greedy algorithm for tolerance interval
    T_lambda <- sort(f_hat$x[f_hat$y > -lambda_hat])
    if (length(T_lambda) < 1) {
      T_lambda <- c(0, 0)
    }
    
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
  res_multimodal[[cand]] <- sapply(rcps, length) %>% table
  
  
  # Result summaries
  idx_bright <- which(x_test[, 1] < median(x_test[, 1]))
  idx_faint <- which(x_test[, 1] >= median(x_test[, 1]))
  
  # Coverage
  cand_level$coverage_total[cand] <- mapply(function(rcps_i, y_i){
      sapply(rcps_i, function(interval){
        y_i >= interval[1] & y_i <= interval[2]
      }) %>% 
        sum()
    }, rcps, y_test) %>% 
      mean()
  cand_level$coverage_bright[cand] <- mapply(function(rcps_i, y_i){
    sapply(rcps_i, function(interval){
      y_i >= interval[1] & y_i <= interval[2]
    }) %>% 
      sum()
  }, rcps[idx_bright], y_test[idx_bright]) %>% 
    mean()
  cand_level$coverage_faint[cand] <- mapply(function(rcps_i, y_i){
    sapply(rcps_i, function(interval){
      y_i >= interval[1] & y_i <= interval[2]
    }) %>% 
      sum()
  }, rcps[idx_faint], y_test[idx_faint]) %>% 
    mean()
  
  # Set size
  cand_level$size_total[cand] <- sapply(rcps, function(rcps_i){
    sapply(rcps_i, function(interval){
      diff(interval)
    }) %>% 
      sum()
  }) %>%
    mean()
  cand_level$size_bright[cand] <- sapply(rcps[idx_bright], function(rcps_i){
    sapply(rcps_i, function(interval){
      diff(interval)
    }) %>% 
      sum()
  }) %>%
    mean()
  cand_level$size_faint[cand] <- sapply(rcps[idx_faint], function(rcps_i){
    sapply(rcps_i, function(interval){
      diff(interval)
    }) %>% 
      sum()
  }) %>%
    mean()
}
cand_level
res_multimodal


# Visualization
fig_list <- list()
for (i in 1:length(delta_list)) {
  # Coverage
  df <- cand_level %>% 
    filter(delta == delta_list[i]) %>% 
    select(alpha, coverage_total, coverage_bright, coverage_faint) %>% 
    gather(type, coverage, -alpha) %>% 
    mutate(type = ifelse(type == "coverage_total", "Total",
                         ifelse(type == "coverage_bright", "Bright", "Faint"))) %>% 
    mutate(type = factor(type, levels = c("Total","Bright","Faint")))
  fig_list[[i]] <- df %>% 
    ggplot(aes(alpha, coverage, color = type)) +
    geom_point(size = 2) +
    geom_line(linewidth = 0.6) +
    scale_y_continuous(limits = range(df$coverage) + c(-0.03, 0.03)) +
    theme_minimal() +
    labs(x = latex2exp::TeX("$\\alpha$"), y = "Coverage", 
         title = latex2exp::TeX(paste0("$\\delta = ", delta_list[i], "$"))) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.title = element_blank(),
          legend.position = "bottom") +
    geom_hline(yintercept = 1-delta_list[i], lty = "dashed", linewidth = 0.6)
}

for (i in 1:length(delta_list)) {
  # Set size
  df <- cand_level %>% 
    filter(delta == delta_list[i]) %>% 
    select(alpha, size_total, size_bright, size_faint) %>% 
    gather(type, size, -alpha) %>% 
    mutate(type = ifelse(type == "size_total", "Total",
                         ifelse(type == "size_bright", "Bright", "Faint"))) %>% 
    mutate(type = factor(type, levels = c("Total","Bright","Faint")))
  fig_list[[i + length(delta_list)]] <- df %>% 
    ggplot(aes(alpha, size, color = type)) +
    geom_point(size = 2) +
    geom_line(linewidth = 0.6) +
    scale_y_continuous(limits = range(df$size) + c(-0.05, 0.05)) +
    theme_minimal() +
    labs(x = latex2exp::TeX("$\\alpha$"), y = "Avg. Set Size", 
         title = latex2exp::TeX(paste0("$\\delta = ", alpha_list[i], "$"))) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.title = element_blank(),
          legend.position = "bottom")
}
ggpubr::ggarrange(plotlist = fig_list, 
                  nrow = 2, ncol = 4,
                  common.legend = TRUE, legend = "bottom")



# # Visualization
# fig_list <- list()
# for (i in 1:length(alpha_list)) {
#   # Coverage
#   df <- cand_level %>% 
#     filter(alpha == alpha_list[i]) %>% 
#     select(delta, coverage_total, coverage_bright, coverage_faint) %>% 
#     gather(type, coverage, -delta) %>% 
#     mutate(type = ifelse(type == "coverage_total", "Total",
#                          ifelse(type == "coverage_bright", "Bright", "Faint"))) %>% 
#     mutate(type = factor(type, levels = c("Total","Bright","Faint")))
#   fig_list[[i]] <- df %>% 
#     ggplot(aes(delta, coverage, color = type)) +
#     geom_point(size = 2) +
#     geom_line(linewidth = 0.6) +
#     scale_y_continuous(limits = c(0.91, 0.99)) +
#     theme_minimal() +
#     labs(x = latex2exp::TeX("$\\delta$"), y = "Coverage", 
#          title = latex2exp::TeX(paste0("$\\alpha = ", alpha_list[i], "$"))) +
#     theme(plot.title = element_text(hjust = 0.5),
#           legend.title = element_blank(),
#           legend.position = "bottom")
# }
# 
# for (i in 1:length(alpha_list)) {
#   # Set size
#   df <- cand_level %>% 
#     filter(alpha == alpha_list[i]) %>% 
#     select(delta, size_total, size_bright, size_faint) %>% 
#     gather(type, size, -delta) %>% 
#     mutate(type = ifelse(type == "size_total", "Total",
#                          ifelse(type == "size_bright", "Bright", "Faint"))) %>% 
#     mutate(type = factor(type, levels = c("Total","Bright","Faint")))
#   fig_list[[i+length(alpha_list)]] <- df %>% 
#     ggplot(aes(delta, size, color = type)) +
#     geom_point(size = 2) +
#     geom_line(linewidth = 0.6) +
#     scale_y_continuous(limits = c(0.08, 0.25)) +
#     theme_minimal() +
#     labs(x = latex2exp::TeX("$\\delta$"), y = "Avg. Set Size", 
#          title = latex2exp::TeX(paste0("$\\alpha = ", alpha_list[i], "$"))) +
#     theme(plot.title = element_text(hjust = 0.5),
#           legend.title = element_blank(),
#           legend.position = "bottom")
# }
# ggpubr::ggarrange(plotlist = fig_list, 
#                   nrow = 2, ncol = 4,
#                   common.legend = TRUE, legend = "bottom")