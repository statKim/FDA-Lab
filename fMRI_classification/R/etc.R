classif_eigvalue <- function(train_data, train_label, test_data, d = 5) {
  # group index
  gr_list <- unique(train_label)
  idx_g1 <- which(train_label == gr_list[[1]])
  idx_g2 <- which(train_label == gr_list[[2]])
  
  # data for each group
  X1 <- train_data[idx_g1, ]
  X2 <- train_data[idx_g2, ]
  N1 <- nrow(X1)
  N2 <- nrow(X2)
  N_test <- nrow(test_data)  # number of test data
  
  # sum of first d eigenvalues
  eig_sum_g1 <- sum(svd(X1)$d[1:d])
  eig_sum_g2 <- sum(svd(X2)$d[1:d])
  
  # Classifier using sum of first d eigenvalues
  pred <- rep(0, N_test)
  for (i in 1:N_test) {
    X1_new<- rbind(X1, test_data[i, ])
    X2_new<- rbind(X2, test_data[i, ])
    
    eig_sum_g1_new <- sum(svd(X1_new)$d[1:d])
    eig_sum_g2_new <- sum(svd(X2_new)$d[1:d])
    
    if (abs(eig_sum_g1_new - eig_sum_g1) < abs(eig_sum_g2_new - eig_sum_g2)) {
      pred[i] <- gr_list[[1]]
    } else {
      pred[i] <- gr_list[[2]]
    }
  }
  
  return(pred)
}


# # set.seed(1000)
# num_sim <- 100
# acc_sim <- rep(0, num_sim)
# for (sim in 1:num_sim) {
#   set.seed(sim)
#   print(sim)
#   idx_train <- sample(1:n, round(n*0.8))
#   
#   X_train <- X[idx_train, ]
#   y_train <- y[idx_train]
#   X_test <- X[-idx_train, ]
#   y_test <- y[-idx_train]
#   
#   
#   fit_eig <- classif_eigvalue(X_train, y_train, X_test, d = 10)
#   acc_sim[sim] <- mean(fit_eig == y_test)
#   
#   print(acc_sim[sim])
# }
# mean(acc_sim)

