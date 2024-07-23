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

# Remove 2 outlying curves
X <- X[-c(33, 123), , ]
y <- y[-c(33, 123)]

n <- dim(X)[1]   # number of curves
m <- dim(X)[2]   # number of timepoints
p <- dim(X)[3]  # number of functional variables
gr <- seq(0, 1, length.out = m)


fit1 <- uFPCA(X[y == 1, , 1], K = 50)
fit1$eig.val %>% round(3)
cumsum(fit1$eig.val) / sum(fit1$eig.val)
fit2 <- uFPCA(X[y == 0, , 1], K = 50)
fit2$eig.val %>% round(3)

8^(-(1:50)) %>% round(3)

2^(-(1:50)) * 1:50 %>% round(3)

##################################################
### Classification
##################################################
library(hdfda)
library(doParallel)
basis <- "bspline"
# n_basis_list <- c(4:9, seq(10, min(40, round(m/2)), by = 1))
n_basis_list <- 4:40

# Parallel computing
# cl <- makePSOCKcluster(detectCores()/2)
cl <- makePSOCKcluster(30)
registerDoParallel(cl)

compare_methods <- c("flasso","hdflda")

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
  
  
  # Proportion of ADHD
  prop1_sim[sim, ] <- c(mean(y_train), mean(y_test))
  print(paste("P(Y_train=1) =", round(mean(y_train), 2), "; P(Y_test=1) =", round(mean(y_test), 2)))
  
  # Predictions
  pred_mat <- data.frame(y_test = y_test)
  model_fit <- list()
  
  
  # flasso
  set.seed(sim)
  start_time <- Sys.time()
  tryCatch({
    cv.fit <- flasso(X_train, y_train, 
                     family = "binomial",
                     basis = "bspline",
                     lambda = NULL,
                     n_basis = n_basis_list)
    pred <- predict(cv.fit, X_test)
  },
  error = function(err){
    print(err)
    error_TF[sim] <<- TRUE
    cv.fit <<- NULL
    pred <<- rep(0, length(y_test))
  })
  end_time <- Sys.time()
  print(end_time - start_time)
  pred_mat$flasso <- pred
  model_fit$flasso <- cv.fit
  
  # mean(pred == y_test)
  
  
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
  
  
  pred_list[[sim]] <- pred_mat
  acc_sim[sim, ] <- apply(pred_mat, 2, function(pred){ mean(pred == y_test) })[-1]
  model_obj[[sim]] <- model_fit
  
  # print(cv.fit$opt_params)
  print(round(acc_sim[sim, ], 3))
  print(round(colMeans(acc_sim[1:sim, ]), 3))
  
  # if (sim %% 2 == 0) {
  #   save(model_obj, pred_list, acc_sim, error_TF, file = "RData/sim_hdflda_align.RData")
  # }
  # save(model_obj, pred_list, acc_sim, error_TF, file = "RData/sim_hdflda_align.RData")
}
stopCluster(cl)
unregister()
save(model_obj, pred_list, acc_sim, error_TF, file = "RData/fmri_yeji_thesis.RData")


### Summary results
load("RData/fmri_yeji_thesis.RData")
load("RData/fmri_yeji_thesis_with_outlier.RData")

print(colMeans(acc_sim))
print(apply(acc_sim, 2, sd))

sapply(pred_list, function(i){ sum(i[, 2]) }) %>% print()
sapply(pred_list, function(i){ sum(i[, 3]) }) %>% print()

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

# Accuracy, specificity, sensitivity with its sd
res_all <- list(
  flasso = sapply(pred_list, function(p_mat) {
    p_mat %>% 
      summarise(
        Accuracy = mean(y_test == flasso),
        Sensitivity = mean(flasso[y_test == 1]),
        Specificity = 1 - mean(flasso[y_test == 0])
      ) %>% 
      unlist()
  }),
  hdflda = sapply(pred_list, function(p_mat) {
    p_mat %>% 
      summarise(
        Accuracy = mean(y_test == hdflda),
        Sensitivity = mean(hdflda[y_test == 1]),
        Specificity = 1 - mean(hdflda[y_test == 0])
      ) %>% 
      unlist()
  }) 
)
res <- data.frame(
  flasso = paste0(
    format(round(rowMeans(res_all$flasso), 3), nsmall = 3),
    " (",
    format(round(apply(res_all$flasso, 1, sd), 3), nsmall = 3),
    sep = ")"
  ),
  hdflda = paste0(
    format(round(rowMeans(res_all$hdflda), 3), nsmall = 3),
    " (",
    format(round(apply(res_all$hdflda, 1, sd), 3), nsmall = 3),
    sep = ")"
  )
)
rownames(res) <- rownames(res_all$flasso)
res
res2 <- res

# Python에서의 summary scores
library(reticulate)
use_condaenv("reticulate")   # Use virtual environment
sk <- import("sklearn.metrics")
res_all <- list(
  # flasso = sapply(pred_list, function(p_mat) {
  #   p_mat %>% 
  #     summarise(
  #       Accuracy = Accuracy(y_test, flasso),
  #       Precision = Precision(y_test, flasso),
  #       Recall = Recall(y_test, flasso),
  #       F1_score = F1_Score(y_test, flasso)
  #     ) %>% 
  #     unlist()
  # }),
  hdflda = sapply(pred_list, function(p_mat) {
    p_mat %>% 
      summarise(
        # Accuracy = sk$accuracy_score(y_test, hdflda),
        Precision = sk$precision_score(y_test, hdflda),
        Recall = sk$recall_score(y_test, hdflda),
        F1_score = sk$f1_score(y_test, hdflda)
      ) %>% 
      unlist()
  }) 
)
res <- data.frame(
  # flasso = paste0(
  #   format(round(rowMeans(res_all$flasso), 3), nsmall = 3),
  #   " (",
  #   format(round(apply(res_all$flasso, 1, sd), 3), nsmall = 3),
  #   sep = ")"
  # ),
  hdflda = paste0(
    format(round(rowMeans(res_all$hdflda), 3), nsmall = 3),
    " (",
    format(round(apply(res_all$hdflda, 1, sd), 3), nsmall = 3),
    sep = ")"
  )
)
rownames(res) <- rownames(res_all$hdflda)
res

rbind(res2, cbind(flasso = "NA", res))


### Python 결과랑 조금 다름
library(MLmetrics)
res_all <- list(
  flasso = sapply(pred_list, function(p_mat) {
    p_mat %>% 
      summarise(
        Accuracy = Accuracy(y_test, flasso),
        Precision = Precision(y_test, flasso),
        Recall = Recall(y_test, flasso),
        F1_score = F1_Score(y_test, flasso)
      ) %>% 
      unlist()
  }),
  hdflda = sapply(pred_list, function(p_mat) {
    p_mat %>% 
      summarise(
        ccuracy = Accuracy(y_test, hdflda),
        Precision = Precision(y_test, hdflda),
        Recall = Recall(y_test, hdflda),
        F1_score = F1_Score(y_test, hdflda)
      ) %>% 
      unlist()
  }) 
)
res <- data.frame(
  flasso = paste0(
    format(round(rowMeans(res_all$flasso), 3), nsmall = 3),
    " (",
    format(round(apply(res_all$flasso, 1, sd), 3), nsmall = 3),
    sep = ")"
  ),
  hdflda = paste0(
    format(round(rowMeans(res_all$hdflda), 3), nsmall = 3),
    " (",
    format(round(apply(res_all$hdflda, 1, sd), 3), nsmall = 3),
    sep = ")"
  )
)
rownames(res) <- rownames(res_all$flasso)
res








