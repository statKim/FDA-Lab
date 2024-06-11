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
### Averaging for each subject
##################################################
library(tidyverse)

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



##################################################
### Alignment using FPCA for each subject
##################################################
library(tidyverse)
library(fdasrvf)
# Functional covariates from 82 regions
n <- 191
m <- 232
p <- 82
X <- array(0, dim = c(n, m, p))
dir_name <- "./fMRI_Classification/Clean_fMRI_PekingUniversity/"
f_list <- list.files(dir_name)
gr <- seq(0, 1, length.out = m)

X2 <- data.frame(Gender = rep(0, n),
                 Age = rep(0, n))
y <- rep(0, n)

n_cores <- 40   # parallet computation
for (i in 1:length(f_list)) {
  print(i)
  df <- data.table::fread(paste0(dir_name, f_list[i]))
  
  region_id <- df$Region_ID
  unique_region_id <- sort(unique(region_id))
  
  for (j in 1:length(unique_region_id)) {
    print(paste(i, "-", j))
    df2 <- df[which(region_id == unique_region_id[j]), 
              which(colnames(df) %in% "t1"):which(colnames(df) %in% "t232")]
    
    # df2_align <- multiple_align_functions(f = t(df2), time = gr, colMeans(df2))
    df2_align <- align_fPCA(f = t(df2), time = gr, num_comp = 50, showplot = F,
                            parallel = T, cores = n_cores)
    X[i, , j] <- rowMeans(df2_align$fn)
    
    # par(mfrow = c(2, 1))
    # matplot(t(df2), type = "l")
    # lines(colMeans(df2), col = 1, lwd = 3)
    # matplot(df2_align$fn, type = "l")
    # lines(rowMeans(df2_align$fn), col = 1, lwd = 3)
    
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


save(X, y, X2, file = "fMRI_Classification/fMRI_align_fpca.RData")

