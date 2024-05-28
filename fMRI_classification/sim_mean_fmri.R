library(tidyverse)
df <- data.table::fread("./fMRI_Classification/Clean_fMRI_PekingUniversity/sub7.csv")
dim(df)
colnames(df)
unique(df$Region_ID) %>% sort()

head(df[, 1:17])


data <- as.matrix(df[, which(colnames(df) %in% "t1"):which(colnames(df) %in% "t232")])
plot(data[1, ], type = "l")


matplot(t(data[1:10, ]), type = "l")
lines(colMeans(data[which(df$Region_ID == 1), ]), col = 1, lwd = 2)
lines(colMeans(data[which(df$Region_ID == 2), ]), col = 2, lwd = 2)
lines(colMeans(data[which(df$Region_ID == 3), ]), col = 3, lwd = 2)



data <- matrix(0, length(unique_region_id), 232)
for (i in 1:length(unique_region_id)) {
  data[i, ] <- colMeans(df[which(region_id == i), 18:ncol(df)])
}

matplot(t(data), type = "l")





##################################################
### Load fMRI data
##################################################
library(tidyverse)

# # Class label of 191 subjects
# y <- read.csv("./fMRI_Classification/class_label.csv", header = T)
# y <- y$y

# Functional covariates from 82 retions
n <- 191
m <- 232
p <- 82
X <- array(0, dim = c(n, m, p))
dir_name <- "./fMRI_Classification/Clean_fMRI_PekingUniversity/"
f_list <- list.files(dir_name)
# substr(f_list, 4, 6)

# subject_id <- sapply(f_list, function(x){
#   strsplit(strsplit(x, "sub")[[1]][2], ".csv")[[1]]
# }) %>% 
#   as.integer()
# subject_id[order(subject_id)]
# strsplit(f_list, "sub" ,".csv")

X2 <- data.frame(Gender = rep(0, n),
                 Age = rep(0, n))
y <- rep(0, n)
for (i in 1:length(f_list)) {
  print(i)
  df <- data.table::fread(paste0(dir_name, f_list[i]))
  
  region_id <- df$Region_ID
  unique_region_id <- sort(unique(region_id))
  
  for (j in 1:length(unique_region_id)) {
    X[i, , j] <- colMeans(df[which(region_id == unique_region_id[j]), 
                             which(colnames(df) %in% "t1"):which(colnames(df) %in% "t232")])
    if (sum(is.na(X[i, , j])) > 0) {
      print(paste(i, j))
    }
  }
  
  if (df$Phenotype_Diagnosis_L1[1] == "ADHD Diagnosed") {
    y[i] <- 1
  }
  
  X2$Gender[i] <- ifelse(df$Phenotype_Gender[1] == "Male", 1, 0)
  X2$Age[i] <- df$Phenotype_Age[1]
}
dim(X)
table(y)
head(X2)

sum(is.na(X))


save(X, y, X2, file = "fMRI_Classification/fMRI.RData")

# n <- dim(X)[1]   # number of curves
# m <- dim(X)[2]   # number of timepoints
# p <- dim(X)[3]  # number of functional variables


df_train <- cbind(y = y_train, 
                  X2[idx_train, ])
fit <- glm(y ~ ., df_train, family = "binomial")
mean(y_test == ifelse(predict(fit, X2[-idx_train, ], type = "response") > 0.5, 1, 0))


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
