##################################################
### Curve alignment
### - Align curves for randomly selected curve of each patient
### - Average it
##################################################
library(tidyverse)
library(fdasrvf)
library(foreach)
library(doParallel)
library(doSNOW)
library(progress)

# Functional covariates from 82 retions
n <- 191
m <- 232
p <- 82
X <- array(0, dim = c(n, m, p))
dir_name <- "../fMRI_classification/fMRI_data/PekingUniv/Clean_fMRI_PekingUniversity/"
f_list <- list.files(dir_name)
# substr(f_list, 4, 6)
idx_order <- sapply(strsplit(f_list, c(".csv")), function(x){ strsplit(x, "sub")[[1]][2] }) %>% 
  as.integer() %>% 
  order()

# y[c(33, 123)]
# unlist(x2)[idx_order][c(33, 123)]
# unlist(x2)[c(33, 123)]
# subject_id <- sapply(f_list, function(x){
#   strsplit(strsplit(x, "sub")[[1]][2], ".csv")[[1]]
# }) %>% 
#   as.integer()
# subject_id[order(subject_id)]
# strsplit(f_list, "sub" ,".csv")

n_cores <- 40   # number of cores for competing methods

# Progress bar
pb <- progress_bar$new(
  format = "[:bar] :percent [Execution time :elapsedfull || Estimated time remaining: :eta]",
  total = length(f_list),   # total number of ticks to complete (default 100)
  clear = FALSE,   # whether to clear the progress bar on completion (default TRUE)
  width = 80   # width of the progress bar
)
progress <- function(n){
  pb$tick()
} 

X <- array(NA, c(n, m, p))
X2 <- data.frame(Gender = rep(0, n),
                 Age = rep(0, n))
y <- rep(0, n)
# X_standard <- matrix(NA, p, m)

# Parallel computation
cl <- makeCluster(n_cores)
registerDoSNOW(cl)
pkgs <- c("fdasrvf")
# Align per each individual
for (i in 1:length(f_list)) {
  # print(i)
  
  # Show the progress bar
  progress(i)
  
  
  df <- data.table::fread(paste0(dir_name, f_list[i]))
  
  region_id <- df$Region_ID
  unique_region_id <- sort(unique(region_id))
  
  
  # # 1st curve of 1st patient (Not align)
  # if (i == 1) {
  #   for (j in 1:p) {
  #     X_standard[j, ] <- df[which(region_id == unique_region_id[j])[1], 
  #                           which(colnames(df) %in% "t1"):which(colnames(df) %in% "t232")] %>% 
  #       as.numeric()
  #   }
  # } 
  
  # Align per each region
  for (j in 1:length(unique_region_id)) {
    # X[i, , j] <- colMeans(df[which(region_id == unique_region_id[j]), 
    #                          which(colnames(df) %in% "t1"):which(colnames(df) %in% "t232")])
    
    idx <- which(region_id == unique_region_id[j])
    x <- df[idx, 
            which(colnames(df) %in% "t1"):which(colnames(df) %in% "t232")] %>% 
      as.matrix()
    
    # Remove all 0 signals
    idx <- setdiff(idx, idx[which(apply(x, 1, sum) == 0)])
    if (length(which(apply(x, 1, sum) == 0)) > 0) {
      x <- x[-which(apply(x, 1, sum) == 0), ]
    }
    
    # Randomly selected curve from i-th patient
    set.seed(i*j)
    idx_selected <- sample(1:length(idx), 1)
    
    # mu <- colMeans(x)
    # obj <- eigen(cov(x), symmetric = T)
    # k <- min(which(cumsum(obj$values) / sum(obj$values) > 0.95))
    # X[i, , j] <- mu + as.numeric(
    #   matrix(colMeans(x %*% obj$vectors[, 1:k]), nrow = 1) %*% t(obj$vectors[, 1:k])
    # )
    
    # Alginment and average it
    # x_align <- matrix(NA, length(idx), m)
    # for (k in 1:length(idx)) {
    #   if (i == 1 & k == 1) {
    #     x_align[k, ] <- x[k, ]
    #   } else {
    #     x_align[k, ] <- curve_pair_align(
    #       matrix(X_standard[j, ], nrow = 1),
    #       matrix(x[k, ], nrow = 1),
    #       # rotation = F,
    #       scale = F
    #     )$beta2n
    #   }
    # }
    tryCatch({
      x_align <- foreach(k = 1:length(idx), .packages = pkgs, .combine = "rbind") %dopar% {
        if (k == idx_selected) {
          out <- x[k, ]
        } else {
          out <- curve_pair_align(
            matrix(x[idx_selected, ], nrow = 1),
            matrix(x[k, ], nrow = 1),
            # rotation = F,
            # scale = F
            scale = T
          )$beta2n
        }
        
        # if (i == 1 & k == 1) {
        #   out <- x[k, ]
        # } else {
        #   out <- curve_pair_align(
        #     matrix(X_standard[j, ], nrow = 1),
        #     matrix(x[k, ], nrow = 1),
        #     # rotation = F,
        #     scale = F
        #   )$beta2n
        # }
        
        return(out)
      }
      X[i, , j] <- colMeans(x_align, na.rm = T)
    },
    error = function(err) {
      X[i, , j] <<- rep(NA, m)
      return(err)
    })

    # matplot(t(x_align), type = "l", col = "gray")
    # lines(colMeans(x_align), lwd = 2)
    
    
    if (sum(is.na(X[i, , j])) > 0) {
      print(paste(i, j))
    }
  }
  
  
  # arr <- array(NA, dim = c(1, ncol(x_align), nrow(x_align)))
  # arr[1, , ] <- t(x_align)
  # obj <- curve_karcher_mean(arr)
  # obj <- multivariate_karcher_mean(arr)
  
  
  if (df$Phenotype_Diagnosis_L1[1] == "ADHD Diagnosed") {
    y[i] <- 1
  }
  
  X2$Gender[i] <- ifelse(df$Phenotype_Gender[1] == "Male", 1, 0)
  X2$Age[i] <- df$Phenotype_Age[1]
}
# End parallel backend
stopCluster(cl) 

dim(X)
table(y)
# head(X2)

sum(is.na(X))


save(X, y, X2, file = "./fMRI_ex.RData")

matplot(t(X[, , 1]), type = "l", col = y+1)


# Indices of error occuring when curve alignment
idx_NA <- matrix(NA, 0, 2)
for (i in 1:n) {
  for (j in 1:p) {
    if (is.na(sum(X[i, , j]))) {
      idx_NA <- rbind(idx_NA, c(i, j))
    }
  }
}
idx_NA




