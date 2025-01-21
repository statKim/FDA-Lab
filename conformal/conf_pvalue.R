########################################
### fMRI Data
########################################
library(tidyverse)
# Class label of 191 subjects
y <- read.csv("../fMRI_classification/fMRI_data/PekingUniv/class_label.csv", header = T)
y <- y$y

# Functional covariates from 82 retions
file_list <- list.files("../fMRI_classification/fMRI_data/PekingUniv/AlignedSubject/")
ord <- sapply(strsplit(file_list, ".csv"), function(x){ strsplit(x, "AlignedSub")[[1]][2] }) %>% 
  as.integer() %>% 
  order()
file_list <- file_list[ord]
for (i in 1:length(file_list)) {
  df <- read.csv(paste0("../fMRI_classification/fMRI_data/PekingUniv/AlignedSubject/", file_list[i]), header = T)
  df <- df[, -1] %>% as.matrix()
  
  if (i == 1) {
    X <- array(0, c(length(file_list), ncol(df), nrow(df)))
  }
  X[i, , ] <- t(df)
}
dim(X)

# Fast Fourier Transform with smoothing splines
# Guo, X., Li, Y., & Hsing, T. (2023). An RKHS Approach for Variable Selection in High-dimensional Functional Linear Models. arXiv preprint arXiv:2310.14419.
X_fft <- X
X_fft_sm <- X
for (i in 1:p) {
  print(i)
  X_i_fft <- apply(X[, , i], 1, function(x) {
    Mod(fft(x)) * (2/m)
    # smooth.spline(Mod(fft(x)) * (2/m))$y
  })
  
  X_i_fft_sm <- apply(X[, , i], 1, function(x) {
    # Mod(fft(x)) * (2/m)
    smooth.spline(Mod(fft(x)) * (2/m))$y
  })
  
  X_fft[, , i] <- t(X_i_fft)
  X_fft_sm[, , i] <- t(X_i_fft_sm)
}
X <- X_fft
X <- X_fft_sm


n <- dim(X)[1]   # number of curves
m <- dim(X)[2]   # number of timepoints
p <- dim(X)[3]  # number of functional variables
gr <- seq(0, 1, length.out = m)


# True outlier index - c(17, 70)
df <- data.frame(
  id = 1:nrow(X),
  y = y
) %>% 
  filter(y == 0)
idx_outliers <- which(df$id %in% c(33, 123))
data <- lapply(1:p, function(i){ X[y == 0, , i] })   # Control group

# Add the 10 ADHD labels as outliers
# - Choose the 20 candidates having the higher value
set.seed(100)
idx_adhd_cand <- which(y == 1)[ order(sapply(which(y == 1), function(i){ max(X[i, , ]) }), decreasing = T)[1:20] ]
# which(y == 1)[ order(sapply(which(y == 1), function(i){ var(apply(X[i, , ], 2, var)) }), decreasing = T)[1:20] ]
# which(y == 1)[ order(sapply(which(y == 1), function(i){ max(apply(X[i, , ], 2, var)) }), decreasing = T)[1:20] ]
idx_adhd <- sample(idx_adhd_cand, 10)   # 10 sampled ADHD curves
data <- lapply(1:p, function(i){ rbind(data[[i]], X[idx_adhd, , i]) })
idx_outliers <- c(idx_outliers, 
                  (nrow(data[[1]])-9):nrow(data[[1]]))
idx_outliers  # 12 outliers

# idx_adhd <- sample(idx_adhd_cand, 10)   # 10 sampled ADHD curves
# df <- data.frame(
#   id = 1:nrow(X),
#   y = y
# ) %>% 
#   filter(y == 0 | id %in% idx_adhd) 
# idx_outliers <- which(df$id %in% c(33, 123, idx_adhd))
# df$y[idx_outliers]
# data <- lapply(1:p, function(i){ X[df$id, , i] })   # data with 12 outliers



source("R/foutlier_cp.R")

### Split conformal prediction
# set.seed(10)
set.seed(1000)
n <- nrow(data[[1]])
p <- length(data)

alpha <- 0.2  # coverage level
prop_train <- 0.8  # proportion of training set

n_train <- round(n * prop_train)   # size of training set
n_test <- n - n_train   # size of test set

# Split training and test data
idx_train <- sample(setdiff(1:n, idx_outliers), n_train)
idx_test <- setdiff(1:n, idx_train)

data_train <- lapply(data, function(x){ x[idx_train, ] })
data_test <- lapply(data, function(x){ x[idx_test, ] })

# Marginal and CCV conformal p-value
cp_obj <- split_conformal_fd(X = data_train, X_test = data_test,
                             type = "depth_transform", type_depth = "projdepth", 
                             depthOptions = list(type = "Rotation"))
conf_pvalue <- cp_obj$conf_pvalue
# conf_pvalue
plot(cp_obj$nonconform_score_calib, ylim = c(-0.4, 0.05))
points(cp_obj$nonconform_score_test, col = 2)



# remove.packages("mrfDepth")
# devtools::install_github("statKim/mrfDepth")
# library(mrfDepth)
# 
# xx <- sapply(data_train, function(x){x[, 100]})
# dim(xx)
# svd(scale(xx))$d
# scale(xx) %*% svd(scale(xx))$v %>% dim
# outlyingness(xx, options = list(type = "Rotation"))
# projdepth(xx, 
#           sapply(data_test, function(x){x[, 100]}),
#           options = list(type = "Rotation"))



### Outlier detection - True: c(17, 70) + 10 ADHDs (115 ~ 124)
# # Single test
# idx_single <- apply(conf_pvalue, 2, function(x){ idx_test[which(x < alpha)] }, simplify = F)
# idx_single
# sapply(idx_single, function(x){
#   get_fdr(x, idx_outliers)
# })
# 
# # Bonfferonni correction - Not reject
# apply(conf_pvalue, 2, function(x){ idx_test[which(x < alpha/n_test)] }, simplify = F)

# BH procedure - Not reject
idx_bh <- apply(conf_pvalue, 2, function(x){ 
  if (sum(sort(x) < (1:n_test)/n_test * alpha) == 0) {
    return(integer(0))
  } else {
    order(x)[1:max(which(sort(x) < (1:n_test)/n_test * alpha))]
  }
  }, simplify = F)
idx_bh <- lapply(idx_bh, function(x){ idx_test[x] })
idx_bh
sapply(idx_bh, function(x){
  get_fdr(x, idx_outliers)
})
sapply(idx_bh, function(x){
  get_tpr(x, idx_outliers)
})
conf_pvalue[idx_test %in% idx_outliers, ]

# # Fisher's combination test (Global test) - Not reject
# apply(conf_pvalue, 2, function(x){ -2*sum(log(x)) >= qchisq(1-alpha, df = 2*n_test) }, simplify = F)


### Existing functional outlier detection (Coverage guarantee X)
library(fdaoutlier)
### Only use test set
df <- abind::abind(data_test, along = 3)

# # MS plot 
# # - Error occurs because of high-dimensionality
# idx_test[msplot(dts = df, plot = F)$outliers]

# Sequential method
seqobj <- seq_transform(df, sequence = "O", depth_method = "erld",
                        erld_type = "one_sided_right", save_data = F)
idx_test[seqobj$outliers$O]
get_fdr(idx_test[seqobj$outliers$O], idx_outliers)
get_tpr(idx_test[seqobj$outliers$O], idx_outliers)



## Use all data set (train + test)
df <- abind::abind(abind::abind(data_train, along = 3),
                   abind::abind(data_test, along = 3), along = 1)
# MS plot
msplot(dts = df, plot = F)$outliers
# Sequential transformation
seqobj <- seq_transform(df, sequence = "O", depth_method = "erld",
                        erld_type = "one_sided_right", save_data = F)
seqobj$outliers$O



# Visualize conformal prediction band
outlier_idx <- idx_test[conf_pvalue$marg <= alpha]
idx_list <- order(sapply(data, function(x){ max(x[outlier_idx, ]) }), decreasing = T)[1:4]
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
    geom_line(data = df2[df2$group %in% outlier_idx, ], 
              aes(x = x, y = y, group = group, color = factor(group))) +
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




##################################################
### EEG data
##################################################
library(tidyverse)

# Load EEG dataset
X1 <- data.table::fread("../fMRI_classification/eeg/alcoholic_data.txt")
X2 <- data.table::fread("../fMRI_classification/eeg/control_data.txt")
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
dim(X)
par(mfrow = c(2, 2))
for (p in c(1,2,20,30)) {
  matplot(t(X[78:122, , p]), type = "l", col = 1)
  matlines(t(X[1:77, , p]), type = "l", col = 2)
}
# par(mfrow = c(2, 2))
# for (p in c(1,2,20,30)) {
#   matplot(t(X[1:77, , p]), type = "l", col = 1)
# }


# Class labels
y <- c(rep(0, 77), rep(1, 45))


data <- lapply(1:p, function(i){ X[, , i] })

### Split conformal prediction
set.seed(1)
n <- nrow(data[[1]])
p <- length(data)

alpha <- 0.2  # coverage level
# prop_train <- 0.7  # proportion of training set
# 
# n_train <- round(n * prop_train)   # size of training set
# n_test <- n - n_train   # size of test set

# Split training and test data
n_train <- 70
idx_outliers <- which(y == 1)
idx_train <- sample(setdiff(1:n, idx_outliers), n_train)
idx_test <- setdiff(1:n, idx_train)
# idx_train <- 1:70
# idx_test <- 71:122
n_test <- length(idx_test)
# idx_train <- 1:77
# idx_test <- 78:122

data_train <- lapply(data, function(x){ x[idx_train, ] })
data_test <- lapply(data, function(x){ x[idx_test, ] })

# Marginal and CCV conformal p-value
cp_obj <- split_conformal_fd(X = data_train, X_test = data_test,
                             type = "depth_transform", type_depth = "projdepth", 
                             depthOptions = list(type = "Rotation"))
conf_pvalue <- cp_obj$conf_pvalue
conf_pvalue

# BH procedure
idx_bh <- apply(conf_pvalue, 2, function(x){ 
  if (sum(sort(x) < (1:n_test)/n_test * alpha) == 0) {
    return(integer(0))
  } else {
    order(x)[1:max(which(sort(x) < (1:n_test)/n_test * alpha))]
  }
}, simplify = F)
idx_bh <- lapply(idx_bh, function(x){ idx_test[x] })
sapply(idx_bh, function(x){
  get_fdr(x, idx_outliers)
})
sapply(idx_bh, function(x){
  get_tpr(x, idx_outliers)
})


### Existing functional outlier detection (Coverage guarantee X)
library(fdaoutlier)
### Only use test set
df <- abind::abind(data_test, along = 3)

# MS plot
idx_test[msplot(dts = df, plot = F)$outliers]
# Sequential method
seqobj <- seq_transform(df, sequence = "O", depth_method = "erld",
                        erld_type = "one_sided_right", save_data = F)
idx_test[seqobj$outliers$O]


## Use all data set (train + test)
df <- abind::abind(abind::abind(data_train, along = 3),
                   abind::abind(data_test, along = 3), along = 1)
# MS plot
msplot(dts = df, plot = F)$outliers
# Sequential transformation
seqobj <- seq_transform(df, sequence = "O", depth_method = "erld",
                        erld_type = "one_sided_right", save_data = F)
seqobj$outliers$O
