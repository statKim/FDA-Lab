
##################################################
### Load fMRI data
##################################################
library(tidyverse)

# Class label of 191 subjects
y <- read.csv("./fMRI_Classification/class_label.csv", header = T)
y <- y$y

# Functional covariates from 82 retions
file_list <- list.files("./fMRI_Classification/AlignedSubject/")
for (i in 1:length(file_list)) {
  df <- read.csv(paste0("./fMRI_Classification/AlignedSubject/", file_list[i]), header = T)
  df <- df[, -1] %>% as.matrix()
  
  if (i == 1) {
    X <- array(0, c(length(file_list), ncol(df), nrow(df)))
  }
  X[i, , ] <- t(df)
}
dim(X)

n <- dim(X)[1]   # number of curves
m <- dim(X)[2]   # number of timepoints
p <- dim(X)[3]  # number of functional variables
gr <- seq(0, 1, length.out = m)



load("fMRI_Classification/fMRI.RData")
source("R/make_basis_mf.R")
source("R/hdfd_lda.R")

n <- dim(X)[1]   # number of curves
m <- dim(X)[2]   # number of timepoints
p <- dim(X)[3]  # number of functional variables
gr <- seq(0, 1, length.out = m)

# Fast Fourier Transform
for (i in 1:p) {
  print(i)
  X_fft <- apply(X[, , i], 1, function(x) {
    # Mod(fft(x)) * (2/m)
    smooth.spline(Mod(fft(x)) * (2/m))$y
  })
  X_fft <- t(X_fft)
  
  X[, , i] <- X_fft
}

##################################################
### Classification
##################################################
basis <- "bspline"

# Parallel computing
# cl <- makePSOCKcluster(detectCores()/2)
cl <- makePSOCKcluster(40)
registerDoParallel(cl)

compare_methods <- c("hdfd_lda")

num_sim <- 100
model_obj <- list()
pred_list <- list()
prop1_sim <- data.frame(matrix(0, num_sim, 2))
colnames(prop1_sim) <- c("P(Y_train=1)", "P(Y_test=1)")
acc_sim <- data.frame(matrix(0, num_sim, length(compare_methods)))
colnames(acc_sim) <- compare_methods
for (sim in 1:num_sim) {
  set.seed(sim)
  print(sim)
  idx_train <- sample(1:n, round(n*0.7))
  
  X_train <- X[idx_train, , ]
  y_train <- y[idx_train]
  X_test <- X[-idx_train, , ]
  y_test <- y[-idx_train]
  
  # Proportion of ADHD
  prop1_sim[sim, ] <- c(mean(y_train), mean(y_test))
  print(paste("P(Y_train=1) =", round(mean(y_train), 2), "; P(Y_test=1) =", round(mean(y_test), 2)))
  
  # Predictions
  pred_mat <- data.frame(y_test = y_test)
  
  # High-dimensional functional LDA
  cv.fit <- cv.hdfd_lda(X_train, y_train)
  pred <- predict(cv.fit$opt_fit, X_test)
  pred_mat$hdfd_lda <- pred
  
  pred_list[[sim]] <- pred_mat
  acc_sim[sim, ] <- apply(pred_mat, 2, function(pred){ mean(pred == y_test) })[-1]
  model_obj[[sim]] <- cv.fit
  
  print(cv.fit$opt_params)
  print(round(acc_sim[sim, ], 3))
  print(round(mean(acc_sim[1:sim, ]), 3))
  
  if (sim %% 5 == 0) {
    save(model_obj, pred_list, acc_sim, file = "RData/sim_hdfd_lda.RData")
  }
}
stopCluster(cl)
unregister()
save(model_obj, pred_list, acc_sim, file = "RData/sim_hdfd_lda.RData")

print(colMeans(acc_sim))
print(apply(acc_sim, 2, sd))

# sapply(pred_list, function(i){ sum(i[, 2]) }) %>% print()

# accuracy, sensitivity and specificity
rbind(
  acc = mean(acc_sim[1:30,]),
  sens = sapply(pred_list, function(p_mat) {
    p_mat %>% 
      filter(y_test == 1) %>% 
      dplyr::select(-y_test) %>% 
      colMeans()
  }) %>% 
    mean(),
  spec = 1 - sapply(pred_list, function(p_mat) {
    p_mat %>% 
      filter(y_test == 0) %>% 
      dplyr::select(-y_test) %>% 
      colMeans()
  }) %>% 
    mean()
) %>% 
  round(3) %>% 
  print()









##################################################
### Load EEG data
##################################################
library(tidyverse)

# Load EEG dataset
X1 <- data.table::fread("eeg/alcoholic_data.txt")
X2 <- data.table::fread("eeg/control_data.txt")
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
# dim(X)
# par(mfrow = c(2, 2))
# for (p in c(1,2,20,30)) {
#   matplot(t(X[90:120, , p]), type = "l", col = 1)
#   matlines(t(X[1:30, , p]), type = "l", col = 2)
# }

# Class labels
y <- c(rep(0, 77), rep(1, 45))


##################################################
### Classification
##################################################

basis <- "bspline"
n_basis_list <- 4:8
lambda_list <- seq(0.1, 0.5, length.out = 10)

# Parallel computing
# cl <- makePSOCKcluster(detectCores()/2)
cl <- makePSOCKcluster(10)
registerDoParallel(cl)

compare_methods <- c("hdfd_lda")

num_sim <- 100
model_obj <- list()
pred_list <- list()
prop1_sim <- data.frame(matrix(0, num_sim, 2))
colnames(prop1_sim) <- c("P(Y_train=1)", "P(Y_test=1)")
acc_sim <- data.frame(matrix(0, num_sim, length(compare_methods)))
colnames(acc_sim) <- compare_methods
# num_nonzero <- data.frame(matrix(0, num_sim, length(compare_methods)))
# colnames(num_nonzero) <- compare_methods
for (sim in 1:num_sim) {
  set.seed(sim)
  print(sim)
  idx_train <- sample(1:n, round(n*0.7))
  
  X_train <- X[idx_train, , ]
  y_train <- y[idx_train]
  X_test <- X[-idx_train, , ]
  y_test <- y[-idx_train]
  
  # Proportion of ADHD
  prop1_sim[sim, ] <- c(mean(y_train), mean(y_test))
  print(paste("P(Y_train=1) =", round(mean(y_train), 2), "; P(Y_test=1) =", round(mean(y_test), 2)))
  
  # Predictions
  pred_mat <- data.frame(y_test = y_test)
  
  # High-dimensional functional LDA
  cv.fit <- cv.hdfd_lda(X_train, y_train, 
                        n_basis_list = n_basis_list, lambda_list = lambda_list)
  pred <- predict(cv.fit$opt_fit, X_test)
  pred_mat$hdfd_lda <- pred
  
  pred_list[[sim]] <- pred_mat
  acc_sim[sim, ] <- apply(pred_mat, 2, function(pred){ mean(pred == y_test) })[-1]
  model_obj[[sim]] <- cv.fit
  
  print(cv.fit$opt_params)
  print(round(acc_sim[sim, ], 3))
  
  # if (sim %% 5 == 0) {
  #   save(model_obj, pred_list, acc_sim, file = "RData/sim_hdfd_lda.RData")
  # }
}
stopCluster(cl)
unregister()
save(model_obj, pred_list, acc_sim, file = "RData/EEG_sim_hdfd_lda.RData")

print(colMeans(acc_sim))
print(apply(acc_sim, 2, sd))

sapply(pred_list, function(i){ sum(i[, 2]) }) %>% print()

# accuracy, sensitivity and specificity
rbind(
  acc = colMeans(acc_sim),
  sens = sapply(pred_list, function(p_mat) {
    p_mat %>% 
      filter(y_test == 1) %>% 
      dplyr::select(-y_test) %>% 
      colMeans()
  }) %>% 
    mean(),
  spec = 1 - sapply(pred_list, function(p_mat) {
    p_mat %>% 
      filter(y_test == 0) %>% 
      dplyr::select(-y_test) %>% 
      colMeans()
  }) %>% 
    mean()
) %>% 
  round(3) %>% 
  print()


sapply(pred_list, function(i){ mean(i[, 1] == i[, 2]) })








##################################################
### K-fold Replications similar to the original paper
##################################################
K <- 5   # 5-fold CV
basis <- "bspline"
n_basis_list <- 4:8
lambda_list <- seq(0.1, 0.5, length.out = 10)

# Parallel computing
# cl <- makePSOCKcluster(detectCores()/2)
cl <- makePSOCKcluster(30)
registerDoParallel(cl)

compare_methods <- c("hdfd_lda")

num_sim <- 100
model_obj <- list()
# pred_list <- list()
# prop1_sim <- data.frame(matrix(0, num_sim, 2))
# colnames(prop1_sim) <- c("P(Y_train=1)", "P(Y_test=1)")
acc_sim <- data.frame(matrix(NA, num_sim, length(compare_methods)))
colnames(acc_sim) <- compare_methods
# num_nonzero <- data.frame(matrix(0, num_sim, length(compare_methods)))
# colnames(num_nonzero) <- compare_methods
for (sim in 1:num_sim) {
  set.seed(sim)
  print(sim)
  
  # Folds
  fold_list <- sample(1:K, n, replace = T)
  
  # Predictions
  # pred_mat <- data.frame(y = y,
  #                        hdfd_lda = NA)
  model_obj[[sim]] <- list()
  
  # K-fold prediction
  acc <- rep(0, K)
  for (j in 1:K) {
    print(paste("Fold", j))
    # Split data
    X_train <- X[fold_list != j, , ]
    X_test <- X[fold_list == j, , ]
    y_train <- y[fold_list != j]
    y_test <- y[fold_list == j]
    
    # # Proportion of alcoholic
    # prop1_sim[sim, ] <- c(mean(y_train), mean(y_test))
    # print(paste("P(Y_train=1) =", round(mean(y_train), 2), "; P(Y_test=1) =", round(mean(y_test), 2)))
    
    # High-dimensional functional LDA
    cv.fit <- cv.hdfd_lda(X_train, y_train, 
                          n_basis_list = n_basis_list, lambda_list = lambda_list)
    pred <- predict(cv.fit$opt_fit, X_test)
    # pred_mat$hdfd_lda[fold_list == j] <- pred
    acc[j] <- mean(pred == y_test)
  }
  # pred_list[[sim]] <- pred_mat
  # acc_sim[sim, ] <- apply(pred_mat, 2, function(pred){ mean(pred == y) })[-1]
  acc_sim[sim, ] <- mean(acc)
  model_obj[[sim]][[j]] <- cv.fit
  
  print(paste(sim, "th accuracy:", round(acc_sim[sim, ], 3)))
  print(paste("Avg accuracy:", round(colMeans(acc_sim, na.rm = T), 3)))
}
stopCluster(cl)
unregister()
save(model_obj, acc_sim, file = "RData/EEG_sim_hdfd_lda_CV.RData")

print(colMeans(acc_sim))
print(apply(acc_sim, 2, sd))

# # accuracy, sensitivity and specificity
# rbind(
#   acc = colMeans(acc_sim),
#   sens = sapply(pred_list, function(p_mat) {
#     p_mat %>% 
#       filter(y_test == 1) %>% 
#       dplyr::select(-y_test) %>% 
#       colMeans()
#   }) %>% 
#     mean(),
#   spec = 1 - sapply(pred_list, function(p_mat) {
#     p_mat %>% 
#       filter(y_test == 0) %>% 
#       dplyr::select(-y_test) %>% 
#       colMeans()
#   }) %>% 
#     mean()
# ) %>% 
#   round(3) %>% 
#   print()




