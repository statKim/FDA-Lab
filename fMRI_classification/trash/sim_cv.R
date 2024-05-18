source("functions.R")

##################################################
### Load fMRI data
##################################################

library(tidyverse)
library(fda)
library(CVXR)
library(fda.usc)

# Class label of 191 subjects
y <- read.csv("./fMRI_Classification/class_label.csv", header = T)
y <- y$y

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


n <- nrow(X)   # number of curves
p <- ncol(X)   # number of timepoints
gr <- seq(0, 1, length.out = p)



# Parallel computing
# cl <- makePSOCKcluster(detectCores()/2)
cl <- makePSOCKcluster(10)
registerDoParallel(cl)

num_sim <- 100
acc_sim <- data.frame(matrix(0, num_sim, 4))
colnames(acc_sim) <- c("fused", "vpc", "FGLM", "Eig")
for (sim in 1:num_sim) {
  print(sim)
  set.seed(sim)
  idx_train <- sample(1:n, round(n*0.8))
  
  X_train <- X[idx_train, ]
  y_train <- y[idx_train]
  X_test <- X[-idx_train, ]
  y_test <- y[-idx_train]
  
  print(paste("P(Y_train=1) =", round(mean(y_train), 2), "; P(Y_test=1) =", round(mean(y_test), 2)))
  
  # Freg with fused lasso
  fit_fused <- tryCatch({
    freg_fused(X_train, y_train, cv = TRUE, n_knots = 30)
  },
  error = function(err){
    print(err)
    idx_train <<- sample(1:n, round(n*0.8))
    
    X_train <<- X[idx_train, ]
    y_train <<- y[idx_train]
    X_test <<- X[-idx_train, ]
    y_test <<- y[-idx_train]
    
    fit_fused <<- freg_fused(X_train, y_train, cv = TRUE, n_knots = 30)
  })
  # fit_fused <- freg_fused(X_train, y_train, cv = TRUE, n_knots = 30)
  acc_sim[sim, 1] <- mean(predict(fit_fused, X_test) == y_test)
  # mean(predict(fit_fused, X_train) == y_train)
  print(fit_fused$lambda)
  
  # VPC
  fit_vpc <- VPC(X_train, y_train, X_test, pve = 0.95)
  acc_sim[sim, 2] <- mean(fit_vpc$pred == y_test)
  
  # Functional GLM - logistic reg
  ldat <- ldata("df" = data.frame(y = y_train),
                "x" = fdata(X_train, argvals = gr))
  a1 <- classif.glm(y ~ x, data = ldat, CV = T)  # CV 안돌고 있음; 애초에 CV 안돌게 되어있음
  newldat <- ldata("df" = data.frame(y = y_test),
                   "x" = fdata(X_test, argvals = gr))
  p1 <- predict(a1, newldat)
  acc_sim[sim, 3] <- mean(p1 == y_test)
  mean(predict(a1, ldat) == y_train)
  
  # Sum of first d eigenvalues
  fit_eig <- classif_eigvalue(X_train, y_train, X_test, d = 5)
  acc_sim[sim, 4] <- mean(fit_eig == y_test)
  
  print( round(acc_sim[sim, ], 3) )
}
stopCluster(cl)  # End parallel

colMeans(acc_sim)
# 30 knots
# fused       vpc      FGLM       Eig 
# 0.5434211 0.4468421 0.5923684 0.5563158 

print(colMeans(acc_sim))
print(apply(acc_sim, 2, sd))

