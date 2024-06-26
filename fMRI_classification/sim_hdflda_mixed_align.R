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

# Scalar variabls of 191 subjects
X2 <- read.csv("./fMRI_Classification/X2.csv", header = T)
colnames(X2) <- c("Gender","Age")
X2$Gender <- ifelse(X2$Gender == "Male", 1, 0)

# Remove 2 outlying curves
X <- X[-c(33, 123), , ]
y <- y[-c(33, 123)]
X2 <- X2[-c(33, 123), ]

n <- dim(X)[1]   # number of curves
m <- dim(X)[2]   # number of timepoints
p <- dim(X)[3]  # number of functional variables
gr <- seq(0, 1, length.out = m)

# Curve length
X_length <- matrix(0, n, p)
for (i in 1:p) {
  X_length[, i] <- apply(X[, , i], 1, function(x){ sum(abs(diff(x))) })
}
X_length <- as.data.frame(X_length)
colnames(X_length) <- paste0("length_", 1:p)


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

# i <- 82
# par(mfrow = c(3, 1))
# matplot(gr, t(X[(y == 0)[1:15], , i]), type = "l", col = 1, lty = 1, lwd = 0.5,
#         xlab = "Time", ylab = "", main = "Raw Furves")
# matlines(gr, t(X[(y == 1)[1:15], , i]), col = 2, lty = 1, lwd = 0.5)
# grid()
# matplot(t(X_fft[(y == 0)[1:15], , i]), type = "l", col = 1, lty = 1, lwd = 0.5,
#         xlab = "Hz", ylab = "", main = "Periodogram using FFT")
# matlines(t(X_fft[(y == 0)[1:15], , i]), col = 2, lty = 1, lwd = 0.5)
# grid()
# matplot(t(X_fft_sm[(y == 0)[1:15], , i]), type = "l", col = 1, lty = 1, lwd = 0.5,
#         xlab = "Hz", ylab = "", main = "Smoothed Periodogram using FFT")
# matlines(t(X_fft_sm[(y == 0)[1:15], , i]), col = 2, lty = 1, lwd = 0.5)
# grid()


# Add curve length into scalar covariates
X2 <- cbind(X2, X_length)

# Scaling the scalar covariates
X2[, -1] <- scale(X2[, -1])


##################################################
### Classification
##################################################
library(hdfda)
library(doParallel)
basis <- "bspline"
n_basis_list <- c(4:9, seq(10, min(40, round(m/2)), by = 1))

# Parallel computing
# cl <- makePSOCKcluster(detectCores()/2)
cl <- makePSOCKcluster(40)
registerDoParallel(cl)

compare_methods <- c("hdflda","hdflda_fft","hdflda_fft_sm",
                     "hdflda_mix","hdflda_fft_mix","hdflda_fft_sm_mix")
num_sim <- 100
model_obj <- list()
pred_list <- list()
prop1_sim <- data.frame(matrix(0, num_sim, 2))
colnames(prop1_sim) <- c("P(Y_train=1)", "P(Y_test=1)")
acc_sim <- data.frame(matrix(0, num_sim, length(compare_methods)))
colnames(acc_sim) <- compare_methods
error_TF <- rep(FALSE, num_sim)
for (sim in 1:num_sim) {
  set.seed(sim)
  print(sim)
  idx_train <- sample(1:n, round(n*0.7))
  
  X_train <- X[idx_train, , ]
  X_test <- X[-idx_train, , ]
  
  y_train <- y[idx_train]
  y_test <- y[-idx_train]
  
  X_fft_train <- X_fft[idx_train, , ]
  X_fft_test <- X_fft[-idx_train, , ]
  
  X_fft_sm_train <- X_fft_sm[idx_train, , ]
  X_fft_sm_test <- X_fft_sm[-idx_train, , ]
  
  X2_train <- X2[idx_train, ]
  X2_test <- X2[-idx_train, ]
  
  # Proportion of ADHD
  prop1_sim[sim, ] <- c(mean(y_train), mean(y_test))
  print(paste("P(Y_train=1) =", round(mean(y_train), 2), "; P(Y_test=1) =", round(mean(y_test), 2)))
  
  # Predictions
  pred_mat <- data.frame(y_test = y_test)
  model_fit <- list()
  
  # hdflda
  set.seed(sim)
  start_time <- Sys.time()
  tryCatch({
    cv.fit <- cv.hdflda(X_train, y_train, n_basis_list = n_basis_list)
    pred <- predict(cv.fit$opt_fit, X_test)
  },
  error = function(err){
    print(err)
    error_TF[sim] <<- TRUE
    cv.fit <<- NULL
    pred <<- rep(0, length(y_test))
  })
  end_time <- Sys.time()
  print(end_time - start_time)
  pred_mat$hdflda <- pred
  model_fit$hdflda <- cv.fit
  
  
  # hdflda + FFT
  set.seed(sim)
  start_time <- Sys.time()
  tryCatch({
    cv.fit <- cv.hdflda(X_fft_train, y_train, n_basis_list = n_basis_list)
    pred <- predict(cv.fit$opt_fit, X_fft_test)
  },
  error = function(err){
    print(err)
    error_TF[sim] <<- TRUE
    cv.fit <<- NULL
    pred <<- rep(0, length(y_test))
  })
  end_time <- Sys.time()
  print(end_time - start_time)
  pred_mat$hdflda_fft <- pred
  model_fit$hdflda_fft <- cv.fit
  
  
  # hdflda + FFT + smoothing
  set.seed(sim)
  start_time <- Sys.time()
  tryCatch({
    cv.fit <- cv.hdflda(X_fft_sm_train, y_train, n_basis_list = n_basis_list)
    pred <- predict(cv.fit$opt_fit, X_fft_sm_test)
  },
  error = function(err){
    print(err)
    error_TF[sim] <<- TRUE
    cv.fit <<- NULL
    pred <<- rep(0, length(y_test))
  })
  end_time <- Sys.time()
  print(end_time - start_time)
  pred_mat$hdflda_fft_sm <- pred
  model_fit$hdflda_fft_sm <- cv.fit
  
  
  # hdflda + scalar covariates (include curve length)
  set.seed(sim)
  start_time <- Sys.time()
  tryCatch({
    cv.fit <- cv.hdflda(X_train, y_train, X2_train, n_basis_list = n_basis_list)
    pred <- predict(cv.fit$opt_fit, X_test, X2_test)
  },
  error = function(err){
    print(err)
    error_TF[sim] <<- TRUE
    cv.fit <<- NULL
    pred <<- rep(0, length(y_test))
  })
  end_time <- Sys.time()
  print(end_time - start_time)
  pred_mat$hdflda_mix <- pred
  model_fit$hdflda_mix <- cv.fit
  
  
  # hdflda + FFT + scalar covariates (include curve length)
  set.seed(sim)
  start_time <- Sys.time()
  tryCatch({
    cv.fit <- cv.hdflda(X_fft_train, y_train, X2_train, n_basis_list = n_basis_list)
    pred <- predict(cv.fit$opt_fit, X_fft_test, X2_test)
  },
  error = function(err){
    print(err)
    error_TF[sim] <<- TRUE
    cv.fit <<- NULL
    pred <<- rep(0, length(y_test))
  })
  end_time <- Sys.time()
  print(end_time - start_time)
  pred_mat$hdflda_fft_mix <- pred
  model_fit$hdflda_fft_mix <- cv.fit
  
  
  # hdflda + FFT + smoothing + scalar covariates (include curve length)
  set.seed(sim)
  start_time <- Sys.time()
  tryCatch({
    cv.fit <- cv.hdflda(X_fft_sm_train, y_train, X2_train, n_basis_list = n_basis_list)
    pred <- predict(cv.fit$opt_fit, X_fft_sm_test, X2_test)
  },
  error = function(err){
    print(err)
    error_TF[sim] <<- TRUE
    cv.fit <<- NULL
    pred <<- rep(0, length(y_test))
  })
  end_time <- Sys.time()
  print(end_time - start_time)
  pred_mat$hdflda_fft_sm_mix <- pred
  model_fit$hdflda_fft_sm_mix <- cv.fit
  
  
  pred_list[[sim]] <- pred_mat
  acc_sim[sim, ] <- apply(pred_mat, 2, function(pred){ mean(pred == y_test) })[-1]
  model_obj[[sim]] <- model_fit
  
  # print(cv.fit$opt_params)
  print(round(acc_sim[sim, ], 3))
  print(round(colMeans(acc_sim[1:sim, ]), 3))
  
  # if (sim %% 2 == 0) {
  #   save(model_obj, pred_list, acc_sim, error_TF, file = "RData/sim_hdflda_align.RData")
  # }
  save(model_obj, pred_list, acc_sim, error_TF, file = "RData/sim_hdflda_align.RData")
}
stopCluster(cl)
unregister()
save(model_obj, pred_list, acc_sim, error_TF, file = "RData/sim_hdflda_align.RData")

print(colMeans(acc_sim))
print(apply(acc_sim, 2, sd))

# sapply(pred_list, function(i){ sum(i[, 2]) }) %>% print()

# accuracy, sensitivity and specificity
rbind(
  acc = colMeans(acc_sim),
  sens = sapply(pred_list, function(p_mat) {
    p_mat %>% 
      filter(y_test == 1) %>% 
      dplyr::select(-y_test) %>% 
      colMeans()
  }) %>% 
    rowMeans(),
  spec = 1 - sapply(pred_list, function(p_mat) {
    p_mat %>% 
      filter(y_test == 0) %>% 
      dplyr::select(-y_test) %>% 
      colMeans()
  }) %>% 
    rowMeans()
) %>% 
  round(3) %>% 
  print()

sum(error_TF)



