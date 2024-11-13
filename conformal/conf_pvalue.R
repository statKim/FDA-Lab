########################################
### fMRI Data
########################################
library(tidyverse)
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



### Split conformal prediction
n <- nrow(data[[1]])
p <- length(data)

alpha <- 0.1  # coverage level
rho <- 0.7  # split proportion (train + calib)

m <- round(n * rho)   # size of training + calib
# l <- n - m   # size of test set

# Split data
idx_train_calib <- sample(setdiff(1:n, idx_outliers), m)
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
nonconform_score_calib <- mapply(function(X_p, pred_p, s_ftn_p){
  apply(X_p, 1, function(x){ max(abs(x - pred_p) / s_ftn_p) })
}, data_calib, pred, s_ftn) %>% 
  apply(1, max)
idx_cutoff <- ceiling((1 - alpha) * (length(idx_calib) + 1))
k_s <- sort(nonconform_score_calib)[idx_cutoff]

# Coverage check
sum(nonconform_score_calib <= k_s) / length(idx_calib)

# Conformal prediction band
lb <- mapply(function(pred_p, s_ftn_p){
  pred_p - k_s*s_ftn_p
}, pred, s_ftn, SIMPLIFY = F)
ub <- mapply(function(pred_p, s_ftn_p){
  pred_p + k_s*s_ftn_p
}, pred, s_ftn, SIMPLIFY = F)


### Conformal p-value (marginal)
m <- length(idx_test)
nonconform_score_test <- mapply(function(X_p, pred_p, s_ftn_p){
  apply(X_p, 1, function(x){ max(abs(x - pred_p) / s_ftn_p) })
}, data_test, pred, s_ftn) %>% 
  apply(1, max)
conf_pvalue_marg <- sapply(nonconform_score_test, function(s){
  (1 + sum(nonconform_score_calib >= s)) / (length(idx_calib) + 1)
})


### Conformal p-value (CCV)
ccv_conf_pvalue <- function(marg_conf_pvalue, method = "simes", delta = 0.1, n_cal, k = NULL) {
  n <- n_cal
  
  if (method == "simes") {
    # Simes adjustment when k = n/2
    if (is.null(k)) {
      k <- ceiling(n/2)
    }
    
    b <- rep(1, n)
    b[1] <- 1 - delta^(1/k)
    sub <- delta
    for (i in 2:(k+1)) {
      sub <- sub * (n-k+2-i)/(n+2-i)
      b[i] <- 1 - sub^(1/k)
    }
    
    h <- function(t) { 
      idx <- ceiling((n+1)*t)
      out <- ifelse(idx == 0, 0,
                    ifelse(idx == n+1, 1, b[idx]))
      return(out)
    }
  } else if (method == "asymp") {
    # Asymptotic adjustment
    c_n <- (-log(-log(1-delta)) + 2*log(log(n)) + 1/2*log(log(log(n))) - 1/2*log(pi)) / sqrt(2*log(log(n)))
    
    b <- sapply(1:n, function(i){
      min(1, i/n + c_n*sqrt( i*(n-i)/(n^3) ))
    })
  } else if (method == "mc") {
    
  }
  
  # Adjusted function for marginal conformal p-value
  h <- function(t) { 
    idx <- ceiling((n+1)*t)
    out <- ifelse(idx == 0, 0,
                  ifelse(idx == n+1, 1, b[idx]))
    return(out)
  }
  
  # Calibration-conditional valid p-value
  conf_pvalue_ccv <- h(marg_conf_pvalue)
  
  return(conf_pvalue_ccv)
}
conf_pvalue_simes <- ccv_conf_pvalue(conf_pvalue_marg, method = "simes", delta = 0.1, n_cal = length(idx_calib))
conf_pvalue_asymp <- ccv_conf_pvalue(conf_pvalue_marg, method = "asymp", delta = 0.1, n_cal = length(idx_calib))

# Conformal p-values
conf_pvalue <- data.frame(
  marg = conf_pvalue_marg,
  simes = conf_pvalue_simes,
  asymp = conf_pvalue_asymp
)

### Outlier detection - True: c(17, 70) + 10 ADHDs (115 ~ 124)
# Single test
idx_single <- apply(conf_pvalue, 2, function(x){ idx_test[which(x < alpha)] }, simplify = F)
# idx_single <- list(
#   marg = idx_test[which(conf_pvalue_marg < alpha)],
#   simes = idx_test[which(conf_pvalue_simes < alpha)],
#   asymp = idx_test[which(conf_pvalue_asymp < alpha)]
# )
idx_single
fdr_single <- sapply(idx_single, function(x){
  sum(!(x %in% idx_outliers)) / max(1, length(x))
})
fdr_single

# Bonfferonni correction - Not reject
apply(conf_pvalue, 2, function(x){ idx_test[which(x < alpha/m)] }, simplify = F)

# BH procedure - Not reject
idx_bh <- apply(conf_pvalue, 2, function(x){ 
  if (sum(sort(x) < (1:m)/m * alpha) == 0) {
    return(NA)
  } else {
    order(x)[1:max(which(sort(x) < (1:m)/m * alpha))]
  }
  }, simplify = F)
idx_bh
fdr_bh <- sapply(idx_bh, function(x){
  sum(!(x %in% idx_outliers)) / max(1, length(x))
})
fdr_bh

# Fisher's combination test (Global test) - Not reject
apply(conf_pvalue, 2, function(x){ -2*sum(log(x)) >= qchisq(1-alpha, df = 2*m) }, simplify = F)


conf_pvalue[idx_test %in% idx_outliers, ]


# Visualize conformal prediction band
outlier_idx <- idx_test[conf_pvalue_marg <= alpha]
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
