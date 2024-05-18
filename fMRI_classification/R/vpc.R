### Jiao, S., Frostig, R. D., & Ombao, H. (2023). Variation pattern classification of functional data. Canadian Journal of Statistics, 51(4), 943-958.
# Inputs:
#   - train_data: n x p matrix; n curves from p timepoints (training data)
#   - train_label: n x 1 vector; class label of `train_data`
#   - test_data: m x p matrix; m curves from p timepoints (test data)
#   - pve: PVE (proportion of variance explained) to choose the number of the feature functions
# Outputs: predicted class label of `test_data`
VPC <- function(train_data, train_label, test_data, pve = 0.9) {
  # group index
  gr_list <- unique(train_label)
  idx_g1 <- which(train_label == gr_list[[1]])
  idx_g2 <- which(train_label == gr_list[[2]])
  
  # data for each group
  X1 <- train_data[idx_g1, ]
  X2 <- train_data[idx_g2, ]
  N1 <- nrow(X1)
  N2 <- nrow(X2)
  
  p <- ncol(train_data)  # number of time points
  N_test <- nrow(test_data)  # number of test data
  
  # difference of covariances
  C1 <- cov(X1)
  C2 <- cov(X2)
  C.diff <- C1 - C2
  C.diff2 <- C.diff %*% C.diff  # R_h
  
  # feature functions
  eig.obj <- eigen(C.diff2, symmetric = T)
  eig.val <- eig.obj$values
  d <- which(cumsum(eig.val) / sum(eig.val) > pve)[1]
  feat.ftn <- eig.obj$vectors[, 1:d] * sqrt(p)  # normalization for numerical integration
  
  # scores using train data using kappa_g
  s_g1 <- t(feat.ftn) %*% C1 %*% feat.ftn / p^2
  s_g2 <- t(feat.ftn) %*% C2 %*% feat.ftn / p^2
  
  # <y, feat.ftn>: to compute scores using test data
  y_nu_inner <- test_data %*% feat.ftn / p
  
  # compute distance and construct classifier
  # - Since we did not consider the lag, W(h) which is described in the paper does not required.
  pred <- rep(0, N_test)
  for (i in 1:N_test) {
    s_test <- y_nu_inner[i, ] %o% y_nu_inner[i, ]
    D_g1 <- sum((s_test - s_g1)^2)
    D_g2 <- sum((s_test - s_g2)^2)
    if (D_g1 < D_g2) {
      pred[i] <- gr_list[[1]]
    } else {
      pred[i] <- gr_list[[2]]
    }
  }
  
  res <- list(
    pred = pred,
    feat.ftn = feat.ftn,
    cov = list(g1 = C1,
               g2 = C2),
    num.feat.ftn = d
  )
  
  return(res)
}