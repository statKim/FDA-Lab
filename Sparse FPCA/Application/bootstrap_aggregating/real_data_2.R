#########################################
### Real data analysis
###   Berkely growth data
#########################################
setwd("C:\\Users\\user\\Desktop\\KHS\\FDA-Lab\\Sparse FPCA\\Application\\bootstrap_aggregating")

library(tidyverse)
library(doParallel)   # parallel computing
library(fdapace)
library(e1071)
library(MASS)
library(data.table)
library(reshape2)
library(fda)
library(xtable)
source("bagFPCA.R")

# parallel computing setting
ncores <- detectCores() - 2
registerDoParallel(ncores)

packages <- c("fdapace","e1071","MASS")   # foreach에서 사용할 package 정의


### load data and preprocessing
data <- rbind(t(growth$hgtm),
              t(growth$hgtf))
y <- factor(c(rep(0, ncol(growth$hgtm)), 
              rep(1, ncol(growth$hgtf))), levels = c(0, 1))
time <- growth$age
N <- length(y)

# graph of data
cbind(as.data.frame(data),
      y = ifelse(y == 1, "female", "male"),
      id = 1:nrow(data)) %>% 
  melt(id.vars = c("id", "y"), 
       variable.name = "time",
       value.name = "height") %>% 
  ggplot(aes(x=time, y=height, group=id, color=y)) +
  geom_line() +
  xlab("Time") +
  ylab("Height") +
  theme_bw() +
  theme(legend.title = element_blank())

# random sparsify
set.seed(1000)
data <- Sparsify(data, time, sparsity = c(12, 15))
id <- unlist( mapply(function(x, y){ rep(x, y) },
                     1:N,
                     sapply(data$Lt, length)) )
time <- unlist(data$Lt)
val <- unlist(data$Ly)
y <- unlist( mapply(function(x, z){ rep(x, z) },
                    y,
                    sapply(data$Lt, length)) )
data <- data.frame(id, time, val, y)


### construct classification models
result <- list()
# set.seed(1000)
seed <- sample(1:10000, 100)
for (simm in 1:100) {
  print( paste(simm, ":", seed[simm]) )
  set.seed(seed[simm])
  
  ### train, test split
  train_test <- train_test_split(data, train.prop = 2/3)
  X.train <- train_test$X.train
  X.test <- train_test$X.test
  y.train <- train_test$y.train
  y.test <- train_test$y.test
  
  ### Single classifier
  err.single <- get_single_err(X.train, X.test, y.train, y.test)
  
  ## Bagging
  # Bootstrap aggregating
  err.bag <- get_bag_err(X.train, X.test, y.train, y.test, B = 100, packages = packages, hyper.para = err.single$hyper.para)
  
  # save result
  res <- as.data.frame(rbind(err.single = err.single$err.single,
                             err.majority = err.bag$err.majority,
                             err.oob = err.bag$err.oob))
  colnames(res) <- c("Logit","SVM(Linear)","SVM(Gaussian)","LDA","QDA","NaiveBayes")
  
  result[[simm]] <- res
}

save(result, file="RData/real_data_2.RData")

## 결과 정리
result <- result[!sapply(result, is.null)]

res <- sapply(1:3, function(i){
  paste(lapply(result, function(x){ x[i, ]*100 }) %>% 
          rbindlist %>% 
          colMeans %>% 
          round(1),
        " (",
        apply(lapply(result[!sapply(result, is.null)], 
                     function(x){ x[i, ]*100 }) %>% 
                rbindlist, 2, sd) %>% 
          round(2),
        ")",
        sep="")
}) %>% 
  t() %>% 
  as.data.frame
colnames(res) <- c("Logit","SVM(Linear)","SVM(Gaussian)","LDA","QDA","NaiveBayes")
rownames(res) <- c("Single","Majority vote","OOB weight")
res

xtable(res)


