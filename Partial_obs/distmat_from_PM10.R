#############################################
### Calculate distance matrix from real data
### PM10
### - Example of partially observed data
#############################################

library(tidyverse)
### Data Load
load("/Users/hyunsung/GoogleDrive/Lab/KHS/partiall_obs/real_data/PM10/pm_data_korea.RData")
head(new_data)
dim(new_data)

data <- new_data[, c(1, 5, 6)] %>% 
  drop_na() %>% 
  distinct()
dim(data)
data
# mat <- dist(data[, -1]) %>% 
#   as.matrix()
# dim(mat)
# mat

library(geosphere)
dist.mat <- distm(data[, -1])
dim(dist.mat)
dist.mat[1:5, 1:5]



# hav.dist <- function(long1, lat1, long2, lat2) {
#   R <- 6371
#   diff.long <- (long2 - long1)
#   diff.lat <- (lat2 - lat1)
#   a <- sin(diff.lat/2)^2 + cos(lat1) * cos(lat2) * sin(diff.long/2)^2
#   b <- 2 * asin(pmin(1, sqrt(a))) 
#   d = R * b
#   return(d)
# }
# hav.dist(data[1, 2], data[1, 3], data[2, 2], data[2, 3])
# 
# distm(data[1:3, -1], fun = distVincentyEllipsoid)


