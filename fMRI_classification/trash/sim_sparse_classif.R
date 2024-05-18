
##################################################
### Load fMRI data
##################################################
library(tidyverse)
source("functions.R")

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


##################################################
### Multivariate functional GLM with Group Lasso
##################################################

basis <- "bspline"
FVE <- 0.95
K <- NULL
n_knots <- 50

library(sparseLDA)
library(HiDimDA)
library(HDclassif)

compare_methods <- c("sparseLDA","Dlda","Mlda","Slda","RFlda","hdda")

num_sim <- 100
pred_list <- list()
prop1_sim <- data.frame(matrix(0, num_sim, 2))
colnames(prop1_sim) <- c("P(Y_train=1)", "P(Y_test=1)")
acc_sim <- data.frame(matrix(0, num_sim, length(compare_methods)))
colnames(acc_sim) <- compare_methods
num_nonzero <- data.frame(matrix(0, num_sim, length(compare_methods)))
colnames(num_nonzero) <- compare_methods
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
  
  # Basis representation for each functional covariate
  basis_obj <- make_basis_mf(X_train, grid = gr, 
                             basis = basis,
                             FVE = FVE,
                             K = K, 
                             n_knots = n_knots)
  X_coef <- basis_obj$X_coef
  # Group indicator for each functional covariate
  groups <- basis_obj$groups
  
  # Make basis coefficient matrix
  X_coef_test <- predict.make_basis_mf(basis_obj, X_test)
  
  
  # Clemmensen, L., Hastie, T. Witten, D. and Ersboell, K. (2011) "Sparse discriminant analysis", Technometrics, To appear
  # https://cran.r-project.org/web/packages/sparseLDA/sparseLDA.pdf
  # library(sparseLDA)
  y_sub <- cbind(y_train, 1-y_train)
  colnames(y_sub) <- c("1", "0")
  fit <- sda(X_coef, 
             y_sub,
             lambda = 1e-6,
             stop = 0.001,
             maxIte = 25,
             trace = F)
  print(paste("sparseLDA:", nrow(fit$beta), "/", ncol(X_coef)))
  pred <- as.integer(predict(fit, X_coef_test)$class) - 1
  pred_mat$sparseLDA <- pred
  num_nonzero$sparseLDA <- nrow(fit$beta)
  
  # https://cran.r-project.org/web/packages/HiDimDA/HiDimDA.pdf
  # Dlda: Diagonal Linear Discriminant Analysis
  # library(HiDimDA)
  dlda_fit <- Dlda(X_coef, factor(y_train))
  print(paste("Dlda:", dlda_fit$nvkpt, "/", ncol(X_coef)))
  pred <- as.integer(predict(dlda_fit, X_coef_test, grpcodes = levels(factor(y_train)))$class) - 1
  pred_mat$Dlda <- pred
  num_nonzero$Dlda <- dlda_fit$nvkpt
  
  # Mlda: Maximum-uncertainty Linear Discriminant Analysis
  # Thomaz, Kitani and Gillies (2006) “A maximum uncertainty LDA-based approach for limited sample size problems - with application to face recognition”, Journal of the Brazilian Computer Society, 12 (2), 7-18.
  mlda_fit <- Mlda(X_coef, factor(y_train))
  print(paste("Mlda:", mlda_fit$nvkpt, "/", ncol(X_coef)))
  pred <- as.integer(predict(mlda_fit, X_coef_test, grpcodes = levels(factor(y_train)))$class) - 1
  pred_mat$Mlda <- pred
  num_nonzero$Mlda <- mlda_fit$nvkpt
  
  # Slda: Shrunken Linear Discriminant Analysis
  # Ledoit, O. and Wolf, M. (2004) “A well-conditioned estimator for large-dimensional covariance matrices.”, Journal of Multivariate Analysis, 88 (2), 365-411.
  slda_fit <- Slda(X_coef, factor(y_train))
  print(paste("Slda:", slda_fit$nvkpt, "/", ncol(X_coef)))
  pred <- as.integer(predict(slda_fit, X_coef_test, grpcodes = levels(factor(y_train)))$class) - 1
  pred_mat$Slda <- pred
  num_nonzero$Slda <- slda_fit$nvkpt
  
  # RFlda: Factor-model Linear Discriminant Analysis
  # Pedro Duarte Silva, A. (2011) “Two Group Classification with High-Dimensional Correlated Data: A Factor Model Approach”, Computational Statistics and Data Analysis, 55 (1), 2975-2990.
  rflda_fit <- RFlda(X_coef, factor(y_train))
  print(paste("RFlda:", rflda_fit$nvkpt, "/", ncol(X_coef)))
  pred <- as.integer(predict(rflda_fit, X_coef_test, grpcodes = levels(factor(y_train)))$class) - 1
  pred_mat$RFlda <- pred
  num_nonzero$RFlda <- rflda_fit$nvkpt
  
  # https://cran.r-project.org/web/packages/HDclassif/HDclassif.pdf
  # Berge, L. Bouveyron, C. and Girard, S. (2012) “HDclassif: An R Package for Model-Based Clustering and Discriminant Analysis of High-Dimensional Data”, Journal of Statistical Software, 46(6), 1–29, url: doi:10.18637/jss.v046.i06
  # library(HDclassif)
  fit <- hdda(X_coef, y_train)
  print(paste("hdda -", 0:1, ":", fit$d, "/", ncol(X_coef)))
  pred <- as.integer(predict(fit, X_coef_test)$class) - 1
  pred_mat$hdda <- pred
  num_nonzero$hdda <- max(fit$d)
  
  
  
  pred_list[[sim]] <- pred_mat
  acc_sim[sim, ] <- apply(pred_mat, 2, function(pred){ mean(pred == y_test) })[-1]
  
  print(acc_sim[sim, ])
}


print(colMeans(acc_sim))
print(apply(acc_sim, 2, sd))

# num_nonzero
colMeans(num_nonzero) %>% print()
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



