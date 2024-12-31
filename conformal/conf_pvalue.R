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
set.seed(1)
n <- nrow(data[[1]])
p <- length(data)

alpha <- 0.2  # coverage level
prop_train <- 0.7  # proportion of training set

n_train <- round(n * prop_train)   # size of training set
n_test <- n - n_train   # size of test set

# Split training and test data
idx_train <- sample(setdiff(1:n, idx_outliers), n_train)
idx_test <- setdiff(1:n, idx_train)

data_train <- lapply(data, function(x){ x[idx_train, ] })
data_test <- lapply(data, function(x){ x[idx_test, ] })

# Marginal and CCV conformal p-value
cp_obj <- split_conformal_fd(X = data_train[1:30], X_test = data_test[1:30],
                             type = "depth", type_depth = "projdepth")
conf_pvalue <- cp_obj$conf_pvalue
conf_pvalue

### Outlier detection - True: c(17, 70) + 10 ADHDs (115 ~ 124)
# Single test
idx_single <- apply(conf_pvalue, 2, function(x){ idx_test[which(x < alpha)] }, simplify = F)
idx_single
sapply(idx_single, function(x){
  get_fdr(x, idx_outliers)
})

# Bonfferonni correction - Not reject
apply(conf_pvalue, 2, function(x){ idx_test[which(x < alpha/n_test)] }, simplify = F)

# BH procedure - Not reject
idx_bh <- apply(conf_pvalue, 2, function(x){ 
  if (sum(sort(x) < (1:n_test)/n_test * alpha) == 0) {
    return(NA)
  } else {
    order(x)[1:max(which(sort(x) < (1:n_test)/n_test * alpha))]
  }
  }, simplify = F)
idx_bh
sapply(idx_bh, function(x){
  get_fdr(x, idx_outliers)
})

# Fisher's combination test (Global test) - Not reject
apply(conf_pvalue, 2, function(x){ -2*sum(log(x)) >= qchisq(1-alpha, df = 2*n_test) }, simplify = F)


conf_pvalue[idx_test %in% idx_outliers, ]


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
