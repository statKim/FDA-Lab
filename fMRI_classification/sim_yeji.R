##################################################
### Simulation of Yeji's thesis
### - functional data classification methods
##################################################
library(tidyverse)
library(hdfda)
library(doParallel)
library(reticulate)
use_condaenv("reticulate")   # Use virtual environment
py_run_file("py/generate_sim_data.py")

n <- 100   # number of curves
m <- 200   # number of timepoints
p <- 50  # number of functional variables
gr <- seq(0, 1, length.out = m)

### Classification
n_basis_list <- 4:40
# n_basis_list <- c(4:9, seq(10, min(40, round(m/2)), by = 1))

# Parallel computing
# cl <- makePSOCKcluster(detectCores()/2)
cl <- makePSOCKcluster(10)
registerDoParallel(cl)

compare_methods <- c("flasso","hdflda")

num_sim <- 100
model_obj <- list()
pred_list <- list()
acc_sim <- data.frame(matrix(0, num_sim, length(compare_methods)))
colnames(acc_sim) <- compare_methods
error_TF <- rep(FALSE, num_sim)
for (sim in 84:num_sim) {
  set.seed(sim)
  print(sim)

  # Generate data and split
  obj <- py$gen_sim_data(
    py$f_bar_star_g1,
    py$f_bar_star_g2,
    py$cov_f_star_g1,
    py$cov_f_star_g2,
    py$n1,
    py$n2,
    y = py$y_,
    seed = py$np$int16(sim)
  )
  X_train <- obj[[1]]
  X_test <- obj[[2]]
  y_train <- obj[[3]]
  y_test <- obj[[4]]
  
  plot(mvtnorm::rmvnorm(1, mean = py$f_bar_star_g1[, 1], sigma = py$cov_f_star_g1)[1, ], type = "l")
  # par(mfrow = c(1, 2))
  # matplot(t(X_train[, , 5]), type = "l", col = y_train+1)
  
  
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
  acc_sim[sim, ] <- apply(pred_mat[, -1], 2, function(pred){ mean(pred == y_test) })
  print(round(acc_sim[sim, ], 3))
}
stopCluster(cl)
# unregister()
save(model_obj, pred_list, acc_sim, error_TF, file = "RData/sim_yeji_thesis.RData")
colMeans(acc_sim)
apply(acc_sim, 2, sd)

sapply(pred_list, function(i){ sum(i[, 2]) })
sapply(pred_list, function(i){ sum(i[, 3]) })

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


