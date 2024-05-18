##################################################
### Load fMRI data
##################################################

library(tidyverse)
library(fda)
library(fda.usc)

# Class label of 191 subjects
y <- read.csv("./fMRI_Classification/class_label.csv", header = T)
y <- y$y

# ## 82 regions with 232 timepoints per 1 subject
# df <- read.csv("./fMRI_Classification/AlignedSubject/AlignedSub0.csv", header = T)
# df <- df[, -1]
# dim(df)   # 82 232
# 
# par(mfrow = c(2, 2))
# matplot(df[, 1:10], type = "l")
# lines(rowMeans(df), lwd = 3)
# 
# matplot(t(df)[, 1:10], type = "l")
# lines(colMeans(df), lwd = 3)
# 
# plot(rowMeans(df), type = "l")
# plot(colMeans(df), type = "l")


# 191 curves from 1st region
idx <- 75
file_list <- list.files("./fMRI_Classification/AlignedSubject/")
for (i in 1:length(file_list)) {
  df <- read.csv(paste0("./fMRI_Classification/AlignedSubject/", file_list[i]), header = T)
  df <- df[idx, -1] %>% unlist() %>% as.numeric()
  
  if (i == 1) {
    X <- df
  } else {
    X <- rbind(X, df)
  }
}
dim(X)  # 191 232
class(X)

matplot(t(X), type = "l", col = y+1)

par(mfrow = c(1, 2))
matplot(t(X[y == 0, ]), type = "l", ylim = c(-30, 30))
matplot(t(X[y == 1, ]), type = "l", ylim = c(-30, 30))

par(mfrow = c(1, 2))
matplot(t(X[y == 0, ][1:5, ]), type = "l", lty = 1, ylim = c(-30, 30))
matplot(t(X[y == 1, ][1:5, ]), type = "l", lty = 1, ylim = c(-30, 30))



n <- nrow(X)   # number of curves
p <- ncol(X)   # number of timepoints
gr <- seq(0, 1, length.out = p)


##################################################
### Classification
##################################################
num_sim <- 100
acc_sim <- data.frame(matrix(0, num_sim, 4))
colnames(acc_sim) <- c("fused", "vpc", "FGLM", "Eig")
for (sim in 1:num_sim) {
  set.seed(sim)
  print(sim)
  idx_train <- sample(1:n, round(n*0.8))
  
  X_train <- X[idx_train, ]
  y_train <- y[idx_train]
  X_test <- X[-idx_train, ]
  y_test <- y[-idx_train]
  
  print(paste("P(Y_train=1) =", round(mean(y_train), 2), "; P(Y_test=1) =", round(mean(y_test), 2)))
  
  # Freg with fused lasso
  fit_fused <- freg_fused(X_train, y_train, lambda = c(0.1, 0.1), n_knots = 60)
  acc_sim[sim, 1] <- mean(predict(fit_fused, X_test) == y_test)
  # mean(predict(fit_fused, X_train) == y_train)
  
  # VPC
  fit_vpc <- VPC(X_train, y_train, X_test, pve = 0.95)
  acc_sim[sim, 2] <- mean(fit_vpc$pred == y_test)
  
  # Functional GLM - logistic reg
  ldat <- ldata("df" = data.frame(y = y_train),
                "x" = fdata(X_train, argvals = gr))
  a1 <- classif.glm(y ~ x, data = ldat)
  newldat <- ldata("df" = data.frame(y = y_test),
                   "x" = fdata(X_test, argvals = gr))
  p1 <- predict(a1, newldat)
  acc_sim[sim, 3] <- mean(p1 == y_test)
  
  # Sum of first d eigenvalues
  fit_eig <- classif_eigvalue(X_train, y_train, X_test, d = 5)
  acc_sim[sim, 4] <- mean(fit_eig == y_test)
  
  print( round(acc_sim[sim, ], 3) )
}

colMeans(acc_sim)
apply(acc_sim, 2, sd)

### 73 region
# fused       vpc      FGLM       Eig 
# 0.4797368 0.4152632 0.5728947 0.5873684 


### Patterns between the number of eigenvalues
num_sim <- 100
max_num_d <- 15
acc_sim <- data.frame(matrix(0, num_sim, max_num_d))
colnames(acc_sim) <- c(1:max_num_d)
for (sim in 1:num_sim) {
  set.seed(sim)
  print(sim)
  idx_train <- sample(1:n, round(n*0.8))
  
  X_train <- X[idx_train, ]
  y_train <- y[idx_train]
  X_test <- X[-idx_train, ]
  y_test <- y[-idx_train]
  
  # Sum of first d eigenvalues
  for (j in 1:max_num_d) {
    fit_eig <- classif_eigvalue(X_train, y_train, X_test, d = j)
    acc_sim[sim, j] <- mean(fit_eig == y_test)
  }
  
  print(round(acc_sim[sim, ], 2))
}
colMeans(acc_sim)
apply(acc_sim, 2, sd)

plot(colMeans(acc_sim), type = "o")



### Classification using curve_length data
# Curve length data
df <- read.csv("./fMRI_Classification/curve_length.csv", header = T)
df <- t(df)
dim(df)
matplot(t(df), type = "l", col = y+1)

par(mfrow = c(1, 2))
matplot(t(df[y == 0, ]), type = "l", ylim = c(0, 1600), col = "gray")
matplot(t(df[y == 1, ]), type = "l", ylim = c(0, 1600), col = "gray")

X <- df
n <- nrow(X)

num_sim <- 100
acc_sim <- data.frame(matrix(0, num_sim, 2))
colnames(acc_sim) <- c("fused", "vpc")
for (sim in 1:num_sim) {
  set.seed(sim)
  print(sim)
  idx_train <- sample(1:n, round(n*0.8))
  
  X_train <- X[idx_train, ]
  y_train <- y[idx_train]
  X_test <- X[-idx_train, ]
  y_test <- y[-idx_train]
  
  # Freg with fused lasso
  fit_fused <- freg_fused(X_train, y_train, lambda = c(0.1, 0.1), n_knots = 50)
  acc_sim[sim, 1] <- mean(predict(fit_fused, X_test) == y_test)
  
  # VPC
  fit_vpc <- VPC(X_train, y_train, X_test, pve = 0.9)
  acc_sim[sim, 2] <- mean(fit_vpc$pred == y_test)
  
  print(acc_sim[sim, ])
}

colMeans(acc_sim)
apply(acc_sim, 2, sd)
# > colMeans(acc_sim)
# fused       vpc 
# 0.4618421 0.4276316 
# > apply(acc_sim, 2, sd)
# fused        vpc 
# 0.08589613 0.08466575 





##################################################
### Multivariate functional GLM
##################################################

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

### Classification
num_sim <- 100
# acc_sim <- data.frame(matrix(0, num_sim, 3))
acc_sim <- c()
# colnames(acc_sim) <- c("fused", "vpc", "FGLM")
for (sim in 1:num_sim) {
  set.seed(sim)
  print(sim)
  idx_train <- sample(1:n, round(n*0.8))
  
  X_train <- X[idx_train, , ]
  y_train <- y[idx_train]
  X_test <- X[-idx_train, , ]
  y_test <- y[-idx_train]
  
  # # Freg with fused lasso
  # fit_fused <- freg_fused(X_train, y_train, lambda = c(0.1, 0.1), n_knots = 50)
  # acc_sim[sim, 1] <- mean(predict(fit_fused, X_test) == y_test)
  # 
  # # VPC
  # fit_vpc <- VPC(X_train, y_train, X_test, pve = 0.9)
  # acc_sim[sim, 2] <- mean(fit_vpc$pred == y_test)
  
  
  X_train <- lapply(1:p, function(i){ fdata(X_train[, , i], argvals = gr) })
  X_test <- lapply(1:p, function(i){ fdata(X_test[, , i], argvals = gr) })
  names(X_train) <- paste0("x", 1:p)
  names(X_test) <- paste0("x", 1:p)
  
  # Functional GLM - logistic reg
  ldat <- ldata(
    df = data.frame(y = y_train),
    mfdata = X_train
  )
  
  formula_fglm <- as.formula( paste("y ~", paste(names(X_train), collapse = "+")) )
  a1 <- classif.glm(formula_fglm, data = ldat)
  # a1 <- classif.gsam(formula_fglm, data = ldat)
  # a1 <- classif.gkam(formula_fglm, data = ldat)  # Time-consuming
  newldat <- ldata(
    df = data.frame(y = y_test),
    mfdata = X_test
  )
  p1 <- predict(a1, newldat)
  # Warning message:
  #   In predict.lm(object, newdata, se.fit, scale = 1, type = if (type ==  : prediction from rank-deficient fit; attr(*, "non-estim") has doubtful cases
  acc_sim[sim] <- mean(p1 == y_test)
  
  print(acc_sim[sim])
}

mean(acc_sim)
sd(acc_sim)





### Classification using functional GLM with lasso
FVE <- 0.7
K <- NULL
n_knots <- 30
alpha <- 0.2

num_sim <- 100
pred_list <- list()
prop1_sim <- data.frame(matrix(0, num_sim, 2))
colnames(prop1_sim) <- c("P(Y_train=1)", "P(Y_test=1)")
acc_sim <- data.frame(matrix(0, num_sim, 4))
colnames(acc_sim) <- c("fglm-bspline", "fglm-fpca", "flasso-bspline", "flasso-fpca")
num_nonzero <- data.frame(matrix(0, num_sim, 2))
colnames(num_nonzero) <- c("flasso-bspline", "flasso-fpca")
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
  
  # Functional GLM with B-spline
  fit <- fglm(X_train, y_train, 
              family = "binomial",
              grid = gr, 
              basis = "bspline",
              n_knots = n_knots)
  pred <- predict(fit, newdata = X_test)
  pred_mat$fglm_bspline <- pred
  
  
  # Functional GLM with FPCA using 5 components
  fit <- fglm(X_train, y_train, 
              family = "binomial",
              grid = gr, 
              basis = "fpca",
              FVE = FVE,
              K = K)
  pred <- predict(fit, newdata = X_test)
  pred_mat$fglm_fpca <- pred
  
  
  # Functional lasso with B-spline basis 
  fit_flasso <- flasso(X_train, y_train, 
                       grid = gr, 
                       basis = "bspline",
                       # lambda = 0.1,
                       alpha = alpha,
                       # FVE = 0.90,
                       # K = NULL, 
                       n_knots = n_knots,
                       cv = TRUE)
  pred <- predict(fit_flasso, X_test)
  pred_mat$flasso_bspline <- pred
  
  beta <- coef(fit_flasso$model.obj, s = "lambda.min") %>% as.numeric()
  num_nonzero[sim, 1] <- sum(abs(beta) > 0)
  print( paste("# of non-zero coef of B-spline =", num_nonzero[sim, 1]) )
  
  # Functional lasso with FPCA
  fit_flasso <- flasso(X_train, y_train, 
                       grid = gr, 
                       basis = "fpca",
                       # lambda = 0.1,
                       alpha = alpha,
                       FVE = FVE,
                       K = K,
                       # n_knots = 20, 
                       cv = TRUE)
  pred <- predict(fit_flasso, X_test)
  pred_mat$flasso_fpca <- pred
  
  beta <- coef(fit_flasso$model.obj, s = "lambda.min") %>% as.numeric()
  num_nonzero[sim, 2] <- sum(abs(beta) > 0)
  print( paste("# of non-zero coef of FPCA =", num_nonzero[sim, 2]) )
  
  
  pred_list[[sim]] <- pred_mat
  acc_sim[sim, ] <- apply(pred_mat[, -1], 2, function(pred){ mean(pred == y_test) })
  
  print(acc_sim[sim, ])
}

colMeans(acc_sim)
apply(acc_sim, 2, sd)

# num_nonzero
colMeans(num_nonzero)

sapply(pred_list, function(i){ sum(i[, 4]) })
sapply(pred_list, function(i){ sum(i[, 5]) })


# accuracy, sensitivity and specificity
rbind(
  acc = colMeans(acc_sim),
  sens = sapply(pred_list, function(p_mat) {
    p_mat %>% 
      filter(y_test == 1) %>% 
      dplyr::select(-y_test) %>% 
      colMeans()
  }) %>% 
    rowmeans(),
  spec = 1 - sapply(pred_list, function(p_mat) {
    p_mat %>% 
      filter(y_test == 0) %>% 
      dplyr::select(-y_test) %>% 
      colMeans()
  }) %>% 
    rowmeans()
) %>% 
  round(3)




