
load("fMRI_Classification/fMRI.RData")

##################################################
### Classification
##################################################
basis <- "bspline"

# Parallel computing
# cl <- makePSOCKcluster(detectCores()/2)
cl <- makePSOCKcluster(80)
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
  # print(paste("P(Y_train=1) =", round(mean(y_train), 2), "; P(Y_test=1) =", round(mean(y_test), 2)))
  
  # Predictions
  pred_mat <- data.frame(y_test = y_test)
  
  # High-dimensional functional LDA
  # cv.fit <- cv.hdfd_lda(X_train, y_train)
  cv.fit <- cv.hdfd_lda(X_train, y_train, 
                        n_basis_list = 4:40)
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
save(model_obj, pred_list, acc_sim, file = "RData/sim_hdfd_lda.RData")

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
