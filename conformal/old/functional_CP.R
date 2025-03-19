library(tidyverse)

########################################
### Berkerley Growth Data
########################################
# Load the growth velocity data
data <- fdasrvf::growth_vel
growth <- data$f[which(data$time >= 4), ]
time <- data$time[data$time >= 4]
label <- factor(c(rep("boy", 39), rep("girl", 54)))

par(mfrow = c(1, 2))
matplot(time, growth[, label == "boy"], type = "l", col = 1, ylim = c(0, 15))
matplot(time, growth[, label == "girl"], type = "l", col = 1, ylim = c(0, 15))

# i <- 10
# plot(time, growth[, i], type = "l")
# lines(x, predict(smooth.spline(time, growth[, i]), x = x)$y, col = 2)

# Smoothing with 101 timepoints
x <- seq(4, 18, length.out = 101)   # narrow timepoints
growth <- apply(growth, 2, function(x_i){ predict(smooth.spline(time, x_i, cv = T), x = x)$y })
# range(growth)


### Split conformal prediction
data <- t(growth)[label == "girl", ]
n <- nrow(data)

alpha <- 0.5  # coverage level
rho <- 0.5  # split proportion

m <- round(n * rho)   # size of training set
l <- n - m   # size of calibration set

# Split data
set.seed(1000)
idx_train <- sample(1:n, m)
idx_calib <- setdiff(1:n, idx_train)

data_train <- data[idx_train, ]
data_calib <- data[idx_calib, ]

# Point predictor
pred <- colMeans(data_train)

# Modulation function (t-function)
abs_resid_train <- apply(data_train, 1, function(x){ abs(x - pred) })  # timepoints x observation
score_H <- apply(abs_resid_train, 2, max)   # esssup_t resid
gamma <- sort(score_H)[ ceiling((1 - alpha) * (m + 1)) ]
idx_H <- which(score_H <= gamma)   # index of H_1
s_ftn <- apply(abs_resid_train[, idx_H], 1, max)
# s_ftn <- conformalInference.fd::computing_s_regression(
#   t(apply(data_train, 1, function(x) { x - pred })),
#   type = "alpha-max",
#   alpha = 0.5,
#   grid_size = length(x), 
#   tau = 1
# )$s_1
# s_ftn <- s_ftn / fdapace::trapzRcpp(x, s_ftn)   # constant product gives same result
# s_ftn <- 1   # not modulating
# s_ftn <- apply(data_train, 2, sd)   # modulate with pointwise sd


# Non-conformity score with modulation
nonconform_score <- apply(data_calib, 1, function(x){ max(abs(x - pred) / s_ftn) })
idx_cutoff <- ceiling((1 - alpha) * (l + 1))
k_s <- sort(nonconform_score)[idx_cutoff]

# Coverage check
sum(nonconform_score <= k_s) / l

# Conformal prediction band
lb <- pred - k_s*s_ftn
ub <- pred + k_s*s_ftn

# Visualize conformal prediction band
# matplot(x, t(data_calib), col = 1, type = "l", lty = 1)
# polygon(c(x, rev(x)),
#         c(pred + k_s*s_ftn, rev(pred - k_s*s_ftn)),
#         border = F, col = "gray")
# matlines(x, t(data_calib), col = 1, type = "l", lty = 1)
# lines(x, pred, col = 2, lwd = 3)
df <- data.frame(
  x = x,
  pred = pred,
  lb = lb,
  ub = ub
)
df2 <- data.frame(
  x = rep(x, l),
  y = as.numeric( t(data_calib) ),
  group = rep(1:l, each = length(x))
)
ggplot() +
  geom_line(data = df2, aes(x = x, y = y, group = group), color = "gray") +
  geom_line(data = df, aes(x = x, y = pred), color = "blue", linewidth = 1) +
  geom_ribbon(data = df, aes(x = x, ymin = lb, ymax = ub), 
              color = "lightblue", fill = "lightblue", alpha = 0.2, linewidth = 1) +
  theme_light() +
  labs(x = "", y = "")









library(conformalInference.fd)
library(roahd)
mfd <- mfData(time, list(t(growth)[label == "girl", ]))

set.seed(1000)
conformal.obj <- conformal.fun.split(
  x = NULL, t_x = NULL,
  y = mfd, t_y = NULL,
  x0 = list(as.list(time)),
  train.fun = mean_lists()$train.fun,
  predict.fun = mean_lists()$predict.fun,
  alpha = 0.5,
  # split = NULL,
  # seed = FALSE,
  # randomized = FALSE,
  # seed.rand = FALSE,
  rho = 0.5,
  verbose = TRUE,
  s.type = "alpha-max"
)
plot_fun(conformal.obj)


# Bike data
gr <- seq(0, 1, length.out = length(bike_log[[1]]$start))
mfd_y <- mfData(gr, list(
  sapply(bike_log, function(x){ x$start }) %>% t(),
  sapply(bike_log, function(x){ x$end }) %>% t()
))
mfd_x <- mfData(gr, list(
  sapply(bike_regressors, function(x){ x$days_week }) %>% t(),
  sapply(bike_regressors, function(x){ x$rain }) %>% t(),
  sapply(bike_regressors, function(x){ x$temp }) %>% t(),
  sapply(bike_regressors, function(x){ x$rain_week }) %>% t()
))

set.seed(1000)
conformal.obj <- conformal.fun.split(
  x = mfd_x, t_x = NULL,
  y = mfd_y, t_y = NULL,
  x0 = list(as.list(gr)),
  train.fun = mean_lists()$train.fun,
  predict.fun = mean_lists()$predict.fun,
  alpha = 0.5,
  # split = NULL,
  # seed = FALSE,
  # randomized = FALSE,
  # seed.rand = FALSE,
  rho = 0.5,
  verbose = TRUE,
  s.type = "alpha-max"
)
plot_fun(conformal.obj)






########################################
### Bike Mobility in Milan Data
########################################
# Bike data
library(conformalInference.fd)
gr <- seq(0, 1, length.out = length(bike_log[[1]]$start))
data <- list(
  end = sapply(bike_log, function(x){ x$end }) %>% t() %>% exp(),
  start = sapply(bike_log, function(x){ x$start }) %>% t() %>% exp()
)

par(mfrow = c(1, 2))
matplot(t(data[[1]]), type = "l")
matplot(t(data[[2]]), type = "l")

# # Smoothing with 101 timepoints
# x <- seq(4, 18, length.out = 101)   # narrow timepoints
# growth <- apply(growth, 2, function(x_i){ predict(smooth.spline(time, x_i, cv = T), x = x)$y })
# # range(growth)


### Split conformal prediction
n <- nrow(data[[1]])
p <- length(data)

alpha <- 0.5  # coverage level
rho <- 0.5  # split proportion

m <- round(n * rho)   # size of training set
l <- n - m   # size of calibration set

# Split data
set.seed(1000)
idx_train <- sample(1:n, m)
idx_calib <- setdiff(1:n, idx_train)

data_train <- lapply(data, function(x){ x[idx_train, ] })
data_calib <- lapply(data, function(x){ x[idx_calib, ] })

# Point predictor
pred <- lapply(data_train, function(x){ colMeans(x) })

# Modulation function (t-function)
abs_resid_train <- mapply(function(X_p, pred_p) {
  apply(X_p, 1, function(x){ abs(x - pred_p) })  # timepoints x observation
}, data_train, pred, SIMPLIFY = F)
score_H <- sapply(abs_resid_train, function(resid_p) { 
  apply(resid_p, 2, max)   # esssup_t resid
}) %>% 
  apply(1, max)
gamma <- sort(score_H)[ ceiling((1 - alpha) * (m + 1)) ]
idx_H <- which(score_H <= gamma)   # index of H_1
s_ftn <- lapply(abs_resid_train, function(resid_p) {
  apply(resid_p[, idx_H], 1, max)
}) 


# Non-conformity score with modulation
nonconform_score <- mapply(function(X_p, pred_p, s_ftn_p){
  apply(X_p, 1, function(x){ max(abs(x - pred_p) / s_ftn_p) })
}, data_calib, pred, s_ftn) %>% 
  apply(1, max)
idx_cutoff <- ceiling((1 - alpha) * (l + 1))
k_s <- sort(nonconform_score)[idx_cutoff]

# Coverage check
sum(nonconform_score <= k_s) / l

# Conformal prediction band
lb <- mapply(function(pred_p, s_ftn_p){
  pred_p - k_s*s_ftn_p
}, pred, s_ftn, SIMPLIFY = F)
ub <- mapply(function(pred_p, s_ftn_p){
  pred_p + k_s*s_ftn_p
}, pred, s_ftn, SIMPLIFY = F)

lb <- lapply(lb, function(x){ ifelse(x < 0, 0, x) })
ub <- lapply(ub, function(x){ ifelse(x < 0, 0, x) })

# Visualize conformal prediction band
fig_list <- list()
for (i in 1:p) {
  df <- data.frame(
    x = gr,
    pred = pred[[i]],
    lb = lb[[i]],
    ub = ub[[i]]
  )
  df2 <- data.frame(
    x = rep(gr, l),
    y = as.numeric( t(data_calib[[i]]) ),
    group = rep(1:l, each = length(gr))
  )
  fig_list[[i]] <- ggplot() +
    geom_line(data = df2, aes(x = x, y = y, group = group), color = "gray") +
    geom_line(data = df, aes(x = x, y = pred), color = "blue", linewidth = 1) +
    geom_ribbon(data = df, aes(x = x, ymin = lb, ymax = ub), 
                color = "lightblue", fill = "lightblue", alpha = 0.2, linewidth = 1) +
    theme_light() +
    labs(x = "", y = "")
}
gridExtra::grid.arrange(grobs = fig_list, ncol = 1)




########################################
### fMRI Data
########################################
# Class label of 191 subjects
y <- read.csv("../fMRI_classification/fMRI_Classification/class_label.csv", header = T)
y <- y$y

# Functional covariates from 82 retions
file_list <- list.files("../fMRI_classification/fMRI_Classification/AlignedSubject/")
ord <- sapply(strsplit(file_list, ".csv"), function(x){ strsplit(x, "AlignedSub")[[1]][2] }) %>% 
  as.integer() %>% 
  order()
file_list <- file_list[ord]
for (i in 1:length(file_list)) {
  df <- read.csv(paste0("../fMRI_classification/fMRI_Classification/AlignedSubject/", file_list[i]), header = T)
  df <- df[, -1] %>% as.matrix()
  
  if (i == 1) {
    X <- array(0, c(length(file_list), ncol(df), nrow(df)))
  }
  X[i, , ] <- t(df)
}
dim(X)

# # Remove 2 outlying curves
# X <- X[-c(33, 123), , ]
# y <- y[-c(33, 123)]

n <- dim(X)[1]   # number of curves
m <- dim(X)[2]   # number of timepoints
p <- dim(X)[3]  # number of functional variables
gr <- seq(0, 1, length.out = m)
data <- lapply(1:p, function(i){ X[y == 0, , i] })   # Control group
# data <- lapply(1:p, function(i){ X[, , i] })   # Overall data

par(mfrow = c(1, 2))
matplot(t(data[[1]]), type = "l")
matplot(t(data[[2]]), type = "l")

# # Smoothing with 101 timepoints
# x <- seq(4, 18, length.out = 101)   # narrow timepoints
# growth <- apply(growth, 2, function(x_i){ predict(smooth.spline(time, x_i, cv = T), x = x)$y })
# # range(growth)


# True outlier index - c(17, 70)
df <- data.frame(
  id = 1:nrow(X),
  y = y
) %>% 
  filter(y == 0) 
which(df$id %in% c(33, 123))


### Split conformal prediction
n <- nrow(data[[1]])
p <- length(data)

alpha <- 0.05  # coverage level
rho <- 0.8  # split proportion (train + calib)

m <- round(n * rho)   # size of training + calib
# l <- n - m   # size of test set

# Split data
set.seed(1000)
idx_train_calib <- sample(setdiff(1:n, c(17, 70)), m)
idx_train <- sample(idx_train_calib, round(m/2))
idx_calib <- setdiff(idx_train_calib, idx_train)
idx_test <- setdiff(1:n, c(idx_train, idx_calib))

data_train <- lapply(data, function(x){ x[idx_train, ] })
data_calib <- lapply(data, function(x){ x[idx_calib, ] })
data_test <- lapply(data, function(x){ x[idx_test, ] })

# Point predictor
pred <- lapply(data_train, function(x){ colMeans(x) })

# Modulation function (t-function)
abs_resid_train <- mapply(function(X_p, pred_p) {
  apply(X_p, 1, function(x){ abs(x - pred_p) })  # timepoints x observation
}, data_train, pred, SIMPLIFY = F)
score_H <- sapply(abs_resid_train, function(resid_p) { 
  apply(resid_p, 2, max)   # esssup_t resid
}) %>% 
  apply(1, max)
gamma <- sort(score_H)[ ceiling((1 - alpha) * (length(idx_train) + 1)) ]
idx_H <- which(score_H <= gamma)   # index of H_1
s_ftn <- lapply(abs_resid_train, function(resid_p) {
  apply(resid_p[, idx_H], 1, max)
}) 


# Non-conformity score with modulation
nonconform_score <- mapply(function(X_p, pred_p, s_ftn_p){
  apply(X_p, 1, function(x){ max(abs(x - pred_p) / s_ftn_p) })
}, data_calib, pred, s_ftn) %>% 
  apply(1, max)
idx_cutoff <- ceiling((1 - alpha) * (length(idx_calib) + 1))
k_s <- sort(nonconform_score)[idx_cutoff]

# Coverage check
sum(nonconform_score <= k_s) / length(idx_calib)

# Conformal prediction band
lb <- mapply(function(pred_p, s_ftn_p){
  pred_p - k_s*s_ftn_p
}, pred, s_ftn, SIMPLIFY = F)
ub <- mapply(function(pred_p, s_ftn_p){
  pred_p + k_s*s_ftn_p
}, pred, s_ftn, SIMPLIFY = F)


### Outlier Detection - True: c(17, 70)
nonconform_score_test <- mapply(function(X_p, pred_p, s_ftn_p){
  apply(X_p, 1, function(x){ max(abs(x - pred_p) / s_ftn_p) })
}, data_test, pred, s_ftn) %>% 
  apply(1, max)
idx_test[which(nonconform_score_test > k_s)]
# alpha = 0.1: 17 70 89 91
# alpha = 0.05: 17 70 89 91
# alpha = 0.03: 17 70 89 91 (0.02부터는 계산 안됨)

sort(sapply(data, function(x){ max(x[c(17, 70), ]) }), decreasing = T)
order(sapply(data, function(x){ max(x[c(17, 70), ]) }), decreasing = T)

# Visualize conformal prediction band
idx_list <- order(sapply(data, function(x){ max(x[idx_test[which(nonconform_score_test > k_s)], ]) }), decreasing = T)[1:4]
# idx_list <- sample(1:82, 4)
fig_list <- list()
for (j in 1:4) {
  i <- idx_list[j]
  df <- data.frame(
    x = gr,
    pred = pred[[i]],
    lb = lb[[i]],
    ub = ub[[i]]
  )
  df2 <- data.frame(
    x = rep(gr, length(idx_test)),
    y = as.numeric( t(data_test[[i]]) ),
    group = rep(idx_test, each = length(gr))
  )
  fig_list[[j]] <- ggplot() +
    geom_line(data = df2, aes(x = x, y = y, group = group), color = "gray") +
    geom_line(data = df2[df2$group %in% idx_test[which(nonconform_score_test > k_s)], ], aes(x = x, y = y, group = group, color = factor(group))) +
    geom_line(data = df, aes(x = x, y = pred), color = "blue", linewidth = 1) +
    geom_ribbon(data = df, aes(x = x, ymin = lb, ymax = ub), 
                color = "lightblue", fill = "lightblue", alpha = 0.2, linewidth = 0.7) +
    labs(x = "", y = "", title = paste(i, "th region")) +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.title = element_blank(),
          legend.position = "bottom")
    
}
ggpubr::ggarrange(plotlist = fig_list, 
                  nrow = 2, ncol = 2,
                  common.legend = TRUE, 
                  legend = "bottom")


### Alternative method - Functional boxplot using whole data
### - Use nonconformity measure instead of band depth
### - See J. Diquigiovanni, M. Fontana and S. Vantini (2022)
# Point predictor
data_train_calib <- lapply(data, function(x){ x[idx_train_calib, ] })
pred <- lapply(data_train_calib, function(x){ colMeans(x) })

# Modulation function (t-function)
abs_resid <- mapply(function(X_p, pred_p) {
  apply(X_p, 1, function(x){ abs(x - pred_p) })  # timepoints x observation
}, data_train_calib, pred, SIMPLIFY = F)
score_H <- sapply(abs_resid, function(resid_p) { 
  apply(resid_p, 2, max)   # esssup_t resid
}) %>% 
  apply(1, max)
gamma <- sort(score_H)[ ceiling((1 - alpha) * (length(idx_train_calib) + 1)) ]
idx_H <- which(score_H <= gamma)   # index of H_1
s_ftn <- lapply(abs_resid, function(resid_p) {
  apply(resid_p[, idx_H], 1, max)
}) 

# Non-conformity score with modulation
nonconform_score <- mapply(function(X_p, pred_p, s_ftn_p){
  apply(X_p, 1, function(x){ max(abs(x - pred_p) / s_ftn_p) })
}, data_train_calib, pred, s_ftn) %>% 
  apply(1, max)
idx_band <- order(nonconform_score)[1:ceiling((1 - alpha) * length(idx_train_calib))]
# idx_outlier <- setdiff(1:n, idx_band)

# Band
lb <- lapply(data_train_calib, function(X_p) {
  apply(X_p[idx_band, ], 2, min)
})
ub <- lapply(data_train_calib, function(X_p) {
  apply(X_p[idx_band, ], 2, max)
})

# Outlier detection for test data
idx_outlier <- mapply(function(X_p, lb_p, ub_p){
  apply(X_p, 1, function(x){ sum( (x < lb_p) | (x > ub_p) ) })
}, data_test, lb, ub) %>% 
  rowSums() %>% 
  sign()
idx_outlier <- idx_test[which(idx_outlier == 1)]
idx_outlier


# Visualize band
fig_list <- list()
for (j in 1:4) {
  i <- idx_list[j]
  df <- data.frame(
    x = gr,
    pred = pred[[i]],
    lb = lb[[i]],
    ub = ub[[i]]
  )
  df2 <- data.frame(
    x = rep(gr, length(idx_test)),
    y = as.numeric( t(data_test[[i]]) ),
    group = rep(idx_test, each = length(gr))
  )
  fig_list[[j]] <- ggplot() +
    geom_line(data = df2, aes(x = x, y = y, group = group), color = "gray") +
    geom_line(data = df2[df2$group %in% idx_outlier, ], aes(x = x, y = y, group = group, color = factor(group))) +
    geom_line(data = df, aes(x = x, y = pred), color = "blue", linewidth = 1) +
    geom_ribbon(data = df, aes(x = x, ymin = lb, ymax = ub), 
                color = "lightblue", fill = "lightblue", alpha = 0.2, linewidth = 0.7) +
    labs(x = "", y = "", title = i) +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.title = element_blank(),
          legend.position = "bottom")
  
}
ggpubr::ggarrange(plotlist = fig_list, 
                  nrow = 2, ncol = 2,
                  common.legend = TRUE, 
                  legend = "bottom")

