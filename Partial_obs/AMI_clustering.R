################################################
### Real data - AMI data
### - Clustering using robust FPCA
### 1. PACE
### 2. Lin & Wang (2020)
### 3. Robust method (Huber loss)
### 4. WRM
### - For all methods, 5-fold CV are performed.
################################################

### Load packages and data
library(GA)   # persp plot
library(mvtnorm)
library(fdapace)   # 1, 2
library(mcfda)   # 7
library(synfd)   # 7
library(doParallel)   # parallel computing
library(doRNG)   # set.seed for foreach
library(MASS)   # huber, rlm
library(latex2exp)
library(tidyverse)
library(robfilter)
source("R/functions.R")
source("Kraus(2015)/pred.missfd.R")
source("Kraus(2015)/simul.missfd.R")
source("sim_kraus.R")

load("AMI_data/AMI.RData")

dim(AMI)

# AMI$Time <- as.integer(AMI$Time)
AMI
length(unique(AMI$ID))   # 1000
colnames(AMI)

df <- spread(AMI, Time, Demand)
dim(df)
df <- df[, -1] %>% 
  as.data.frame

plot(1:240, df[1, 1:240], type = "l")
