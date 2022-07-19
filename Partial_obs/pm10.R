################################################
### Real data analysis - PM10
################################################
### Download required packages
# devtools::install_github("statKim/robfpca")
# devtools::install_github("statKim/mcfda.rob")

### Load packages
library(robfpca)   # proposed methods
library(mcfda.rob)   # R-Kraus
library(tidyverse)
library(fdapace)
library(univOutl)

# Codes can be obtained from Kraus(2015), JRSS-B.
source("sim_utills/pred.missfd.R")
source("sim_utills/simul.missfd.R")

# R-Kraus
source("sim_utills/robust_Kraus.R")


#########################################
### Data Load
#########################################
### PM10 데이터
load("/Users/hyunsung/GoogleDrive/Lab/KHS/partiall_obs/real_data/pm_data_korea.RData")
head(new_data)

# Unique한 측정소 id
loc <- unique(new_data[, 1]) %>% 
  as.numeric()

### 측정소 위치 데이터
full_data <- readRDS("/Users/hyunsung/GoogleDrive/Lab/KHS/partiall_obs/real_data/KoreaPM10.RDS")
class(full_data)
full_data$place

# 측정소별 위경도
temp_loc <- unique(cbind(new_data$lon, 
                         new_data$lat))

# 한반도 지도 위에 측정소 위치 찍기
plot(full_data$shape)
points(temp_loc[, 1], temp_loc[, 2], pch = 16, col = "red")

length(loc)   # 336

unique_loc <- NA
for (loc_num in 1:length(loc)) {
  unique_loc <- rbind(unique_loc, 
                      new_data[which(new_data[,1] == loc[loc_num])[1], 5:6]) 
}
unique_loc <- unique_loc[-1, ]
dim(unique_loc)   # 336 2


########################################
### 3월 자료만 가져오기
### n = location*day, T = hour(24)
#########################################
data_3 <- matrix(NA, nrow = length(loc)*31, ncol = (24+1+2))
for (loc_num in 1:length(loc)){
  data <- filter(new_data, new_data[,1] == loc[loc_num])
  data_3[c( ((loc_num-1)*31+1) : (31*loc_num) ), 1] <- loc[loc_num]
  data_3[c( ((loc_num-1)*31+1) : (31*loc_num) ), 2:3] <- t( matrix(rep(as.numeric(unique_loc[loc_num, ]), 31), nrow = 2) )
  
  if (length( which(data[, 2] == '2017030101') ) != 0) {
    time_index <- c( which(data[, 2] == '2017030101') : which(data[, 2] == '2017033124') )
    data_3[c( ((loc_num-1)*31+1) : (31*loc_num) ), -c(1:3) ] <- t( matrix(data[time_index, 3], ncol = 31) )
  }
}


par(mfrow = c(1, 3))
for (loc in c(28,81,300)) {
  temp_data_3 <- data_3[c((loc-1)*31+1):(loc*31), ]
  plot(temp_data_3[1, -c(1:3)], 
       type = 'n', ylim = c(0, 240), 
       main = '', ylab = 'PM10', xlab = 'Hour')
  for (k in 1:31){
    lines(temp_data_3[k, -c(1:3)], col = gray(0.8), lwd = 2)
  }
  abline(h = 100, lty = 2)
  print(loc)
  print( length(which(is.na(temp_data_3) == T)) / (31*24)*100 )
}


aa <- vector()
for (loc in 1:100){
  temp_data_3 <- data_3[c((loc-1)*31+1):(loc*31), ]
  aa[loc] <- length(which(is.na(temp_data_3) == T)) / (31*24)*100
}
mean(aa)


### 최종 분석 데이터 생성
x_without_grid <- data_3[, -c(1:3)]
dim(x_without_grid)
all_na_ind <- which( apply(x_without_grid, 1, function(row){ sum(is.na(row)) == ncol(x_without_grid) }) )
elimintaed_pt <- unique(data_3[all_na_ind,1])
com.x <- x_without_grid[-all_na_ind,]
dim(com.x)
com.x <- com.x[-7132, ]
dim(com.x)   # 9913 24

### PM10과 대응되는 측정소 위치 및 이름 데이터
loc.x <- data_3[-all_na_ind, 1:3]
loc.x <- data.frame(loc.x[-7132, ])
colnames(loc.x) <- c("id","lon","lat")
dim(loc.x)   # 9913 3
head(loc.x)

place <- full_data$place[, 2:3]
colnames(place) <- c("city","id")
place <- loc.x %>% 
  left_join(place, by = "id") %>% 
  mutate(city = as.character(city))
dim(place)   # 9913 4
head(place)




##########################################################
### Functional principal component analysis
##########################################################
x <- com.x   # matrix type
x.2 <- matrix2list(x)   # Lt, Ly type

n <- nrow(x)   # 9913
p <- ncol(x)
work.grid <- seq(0, 1, length.out = p)

K <- 3   # number of PCs
seed <- 100   # seed number

#############################
### Proposed FPCA
#############################
set.seed(seed)
cov.obj.cv <- cv.cov_ogk(x,
                         K = 5,
                         # bw_cand = bw_cand,
                         MM = TRUE,
                         type = 'huber')
print(cov.obj.cv$selected_bw)
cov.obj <- cov_ogk(x,
                   type = "huber",
                   MM = TRUE,
                   smooth = T,
                   bw = cov.obj.cv$selected_bw)
mu.ogk <- cov.obj$mean
cov.ogk <- cov.obj$cov
pca.ogk.obj <- funPCA(x.2$Lt, x.2$Ly,
                      mu.ogk, cov.ogk, sig2 = 0,
                      work.grid, K = K)
prop_pc <- pca.ogk.obj$eig.fun
pcs_prop <- pca.ogk.obj$pc.score

range(pcs_prop, na.rm = T)

max.grid <- which.max(apply(x, 1, max, na.rm = T))
x[max.grid, ]
pcs_prop[max.grid, ]



#############################
### Kraus
#############################
K <- 3   # number of PCs

mu.kraus <- mean.missfd(x)
cov.kraus <- var.missfd(x)

eig.kraus <- get_eigen(cov.kraus, work.grid)
pc_kraus <- eig.kraus$phi[, 1:K]
load("RData/PM10_score_kraus.RData")
dim(pcs_kraus)
# # 엄청 오래 걸림... RData로 저장함
# # user  system elapsed 
# # 319.123   7.051 326.107 
# system.time({
#   pcs_kraus <- matrix(0, n, K)
#   for (i in 1:n) {
#     if (i %% 500 == 0) {
#       print(i)
#     }
#     pcs_kraus[i, ] <- pred.score.missfd(x[i, ], 
#                                         phi = pc_kraus, 
#                                         x = x)
#   }
# })
# dim(pcs_kraus)
# save(pcs_kraus, file = "RData/PM10_score_kraus.RData")


#############################
### Robust Kraus
#############################
mu.Mest <- mean_Mest(x)	
cov.Mest <- cov_Mest(x)

eig.Mkraus <- get_eigen(cov.Mest, work.grid)
pc_Rkraus <- eig.Mkraus$phi[, 1:K]
load("RData/PM10_score_Rkraus.RData")
dim(pcs_Rkraus)
# # 엄청 오래 걸림... RData로 저장함
# # user   system  elapsed 
# # 4760.884 1147.125 5907.102 
# pcs_Rkraus <- matrix(0, n, K)
# for (i in 1:n) {
#   if (i %% 500 == 0) {
#     print(i)
#   }
#   pcs_Rkraus[i, ] <- pred.score.rob.missfd(x[i, ], 
#                                            phi = pc_Rkraus, 
#                                            x = x,
#                                            R = cov.Mest,
#                                            mu = mu.Mest)
# }
# save(pcs_Rkraus, file = "RData/PM10_score_Rkraus.RData")
if (!is_null(K)) {
  K_Mkraus <- K
  pve_Mkraus <- eig.Mkraus$PVE[K_Mkraus]
} else {
  K_Mkraus <- which(eig.Mkraus$PVE > pve)[1]
  pve_Mkraus <- eig.Mkraus$PVE[K_Mkraus]
}
pca.Mkraus.obj <- list(K = K_Mkraus,
                       PVE = pve_Mkraus,
                       mu = mu.Mest,
                       cov = cov.Mest,
                       lambda = eig.Mkraus$lambda[1:K_Mkraus],
                       eig.fun = eig.Mkraus$phi[, 1:K_Mkraus],
                       # additionally need
                       data = x.2,
                       sig2 = 0,
                       pc.score = pcs_Rkraus,
                       work.grid = work.grid,
                       eig.obj = eig.Mkraus)
class(pca.Mkraus.obj) <- "funPCA"


#############################
### Yao
#############################
set.seed(seed)
kern <- "epan"
# optns <- list(methodXi = "CE", dataType = "Sparse", kernel = kern, verbose = FALSE,
#               kFoldMuCov = 3, methodBwMu = "CV", methodBwCov = "CV", 
#               useBinnedCov = FALSE, error = FALSE)
optns <- list(methodXi = "CE", dataType = "Sparse", kernel = kern, verbose = FALSE,
              kFoldMuCov = 5, methodBwMu = "GCV", methodBwCov = "GCV", 
              error = FALSE, useBinnedData = "OFF")

mu.yao.obj <- GetMeanCurve(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)
cov.yao.obj <- GetCovSurface(Ly = x.2$Ly, Lt = x.2$Lt, optns = optns)
mu.yao <- mu.yao.obj$mu
cov.yao <- cov.yao.obj$cov

# 51 grids로 estimate된 것을 24 grids로 변환
if (length(work.grid) != 51) {
  mu.yao <- ConvertSupport(fromGrid = cov.yao.obj$workGrid, toGrid = work.grid,
                           mu = mu.yao.obj$mu)
  cov.yao <- ConvertSupport(fromGrid = cov.yao.obj$workGrid, toGrid = work.grid,
                            Cov = cov.yao.obj$cov)
}
pca.yao.obj <- funPCA(x.2$Lt, x.2$Ly, 
                      mu.yao, cov.yao, sig2 = 0, 
                      work.grid, PVE = pve, K = K)
pc_yao <- pca.yao.obj$eig.fun
pcs_yao <- pca.yao.obj$pc.score




#############################
### Plot Eigenfunction
#############################
par(mfrow = c(1, 3))
mname <- c("(a)","(b)","(c)")
for (k in 1:3){
  plot(prop_pc[, k], type = 'l', 
       main = mname[k], ylab = '', xlab = 'Hour', 
       ylim = range(prop_pc[, k], pc_yao[, k], pc_Rkraus[, k]),
       lwd = 2, cex.lab = 1.3, cex.axis = 1.3, cex.main = 1.3)
  lines(pc_Rkraus[, k], col = 2, lwd = 2, lty = 3)
  lines(pc_yao[, k], col = 3, lwd = 2, lty = 2)
  
  if (k == 1) {
    legend('topright', c('Proposed FPCA','Sparse FPCA',"Robust Kraus"), 
           col = c(1,3,2), lty = c(1,2,3), lwd = 2, cex = 1.2)     
  }
}

for (k in 1:3){
  plot(prop_pc[, k], type = 'l', 
       main = mname[k], ylab = '', xlab = 'Hour', 
       ylim = range(prop_pc[, k], pc_yao[, k], pc_kraus[, k]),
       lwd = 2, cex.lab = 1.3, cex.axis = 1.3, cex.main = 1.3)
  lines(pc_kraus[, k], col = 2, lwd = 2, lty = 3)
  lines(pc_yao[, k], col = 3, lwd = 2, lty = 2)
  
  if (k == 1) {
    legend('topright', c('Proposed FPCA','Sparse FPCA',"Kraus"), 
           col = c(1,3,2), lty = c(1,2,3), lwd = 2, cex = 1.2)     
  }
}

# for(k in 2){
#   plot( prop_pc[,k],type='l', main='(b)', ylim=range(pc_yao[,k], prop_pc[,k] ), ylab='', xlab='Hour', lwd=2, cex.lab=1.3, cex.axis=1.3, cex.main=1.3)
#   lines( pc_Rkraus[,k],col=2, lwd=2, lty=3)
#   lines( pc_yao[,k], col=3, lwd=2, lty=2)
#   # legend('topleft', c('Kraus', "Rbust Kraus",'Proposed'), col=c(1,3,2), lty=1, lwd=2) 
# }
# for(k in 3){
#   plot( prop_pc[,k],type='l', main='(c)', ylim=range(pc_yao[,k], prop_pc[,k] ), ylab='', xlab='Hour', lwd=2, cex.lab=1.3, cex.axis=1.3, cex.main=1.3)
#   lines( pc_Rkraus[,k],col=2, lwd=2, lty=3)
#   lines( pc_yao[,k], col=3, lwd=2, lty=2)
#   # legend('topleft', c('Kraus', "Rbust Kraus",'Proposed'), col=c(1,3,2), lty=1, lwd=2) 
# }




#### Compute average PC scores
num_pc <- 1

result_aa <- cbind(data_3[-all_na_ind, 2], 
                   data_3[-all_na_ind, 3])  
col.pall1 <- pcs_prop[, num_pc]
col.pall2 <- pcs_yao[, num_pc]
col.pall3 <- pcs_Rkraus[, num_pc]

unique_index <- seq(1, nrow(result_aa), by = 31)

avg_pcs_prop <- avg_pcs_Rkraus <- avg_pcs_yao <- vector()
for (k in 1:length(unique_index)) {
  temp_index <- which( 
    result_aa[, 1] == result_aa[unique_index[k], 1] & 
      result_aa[, 2] == result_aa[unique_index[k], 2] 
  )
  avg_pcs_prop[k] <- median(col.pall1[temp_index], na.rm = T)
  avg_pcs_Rkraus[k] <- median(col.pall3[temp_index], na.rm = T)
  avg_pcs_yao[k] <- median(col.pall2[temp_index], na.rm = T)
}


par(mfrow = c(1, 1))
plot(full_data$shape, main = '')
rbPal <- colorRampPalette(c("green","yellow","red"))
br <- 8
br2 <- 8
col1 <- rbPal(br)[as.numeric(cut(avg_pcs_prop, breaks = br2))] #colorRampPalette(c("red","yellow","green"))(length(col.pall))
points(result_aa[unique_index, 1], result_aa[unique_index, 2], 
       pch = 19, col = col1, cex = 2)
legend('bottomright', levels(cut(avg_pcs_prop, breaks = br)), 
       col = rbPal(br), pch = 19, cex = 1.3)
box()



plot(full_data$shape, main = 'Kraus')
col1 <- rbPal(br)[as.numeric(cut(avg_pcs_kraus, breaks = br2))]
points(result_aa[unique_index, 1], result_aa[unique_index, 2], 
       pch = 19, col = col1, cex = 2)
legend('bottomright', levels(cut(col.pall2, breaks = br)),
       col = rbPal(br), pch = 19)
box()


plot(full_data$shape, main = 'Robust Kraus')
col1 <- rbPal(br)[as.numeric(cut(avg_pcs_Rkraus, breaks = br2))]
points(result_aa[unique_index, 1], result_aa[unique_index, 2], 
       pch = 19, col = col1, cex = 2)
legend('bottomright', levels(cut(col.pall3, breaks = br)), 
       col = rbPal(br), pch = 19)
box()




##############################
### Outlier detection
##############################
num_pc <- 1
col.pall1 <- pcs_prop[, num_pc]
col.pall3 <- pcs_Rkraus[, num_pc]
col.pall2 <- pcs_yao[, num_pc]

cut_off <- 2.2
adj_boxplot <- boxB(col.pall1, method = 'asymmetric', k = cut_off)
prop_outlier <- adj_boxplot$outliers
place$city[prop_outlier]

adj_boxplot <- boxB(col.pall3, method = 'asymmetric', k = cut_off)
Rkraus_outlier <- adj_boxplot$outliers
place$city[Rkraus_outlier]

adj_boxplot <- boxB(col.pall2, method = 'asymmetric', k = cut_off)
yao_outlier <- adj_boxplot$outliers
place$city[yao_outlier]

setdiff(prop_outlier, Rkraus_outlier)
setdiff(prop_outlier, yao_outlier)
setdiff(Rkraus_outlier, prop_outlier)
setdiff(yao_outlier, prop_outlier)

### Plot
par(mfrow = c(1, 3))
plot(x[1, ], type = 'n',
     ylim = c(0, 350), main = 'Proposed FPCA', ylab = 'PM10', xlab = 'Hour')
set.seed(30)
for (k in c(sample(c(1:nrow(x)), 40), 
            prop_outlier,
            yao_outlier)) {
  lines(x[k, ], col = gray(0.8), lwd = 2)
}
abline(h = 100, lty = 2)
for (j in prop_outlier) {
  lines(x[j, ], col = 2)
}

plot(x[1, ], type = 'n', 
     ylim = c(0, 350), main = 'Sparse FPCA', ylab = 'PM10', xlab = 'Hour')
set.seed(30)
for (k in c(sample(c(1:nrow(x)), 40), 
            prop_outlier,
            yao_outlier)) {
  lines(x[k, ], col = gray(0.8), lwd = 2)
}
abline(h = 100, lty = 2)
for (j in yao_outlier) {
  lines(x[j, ], col = 2)	
}
for (j in setdiff(prop_outlier, yao_outlier)) {
  lines(x[j, ], col = "blue", lwd = 3)
}
# lines(x[2522, ], col = "blue", lwd = 3)

plot(x[1, ], type = 'n', 
     ylim = c(0, 350), main = 'Robust Kraus', ylab = 'PM10', xlab = 'Hour')
set.seed(30)
for (k in c(sample(c(1:nrow(x)), 40),
            prop_outlier,
            yao_outlier)) {
  lines(x[k, ], col = gray(0.8), lwd = 2)
}
abline(h = 100, lty = 2)
for (j in Rkraus_outlier) {
  lines(x[j, ], col = 2)	
}
for (j in setdiff(prop_outlier, Rkraus_outlier)) {
  lines(x[j, ], col = "blue", lwd = 3)
}
# lines(x[9778, ], col = "blue", lwd = 3)



### Plot of paper
par(mfrow = c(1, 3))
plot(com.x[1, ], type = 'n',
     ylim = c(0, 350), main = 'Proposed FPCA', ylab = 'PM10', xlab = 'Hour', cex.main=1.5, cex.lab=1.2, cex.axis=1.2)
set.seed(30)
for (k in c(sample(c(1:nrow(com.x)), 40)) ) {
  lines(com.x[k, ], col = gray(0.8), lwd = 1.5)
}
abline(h = 100, lty = 2)
for (j in 1:length(prop_outlier)) {
  lines(com.x[prop_outlier, ][j, ], col = 2, lwd=1.5)
}
legend('topright',  c('Observation', 'Detected outliers'), col=c(gray(0.8),2) ,lty=1, lwd=2,cex=1.6)

plot(com.x[1, ], type = 'n', 
     ylim = c(0, 350), main = 'Sparse FPCA', ylab = 'PM10', xlab = 'Hour', cex.main=1.5, cex.lab=1.2, cex.axis=1.2)
set.seed(30)
for (k in c(sample(c(1:nrow(com.x)), 40)) ) {
  lines(com.x[k, ], col = gray(0.8), lwd = 1.5)
}
abline(h = 100, lty = 2)
for (j in 1:length(yao_outlier)) {
  lines(com.x[yao_outlier, ][j, ], col = 2, lwd = 1.5)	
}
for (k in setdiff(prop_outlier,
                  intersect(prop_outlier, yao_outlier))) {
  lines(x[k, ], col = 4, lwd = 3)
}
legend('topright',  c('Observation', 'Detected outliers', 'Only by proposed FPCA'), col=c(gray(0.8),2,4) ,lty=1, lwd=2,cex=1.6)


plot(com.x[1, ], type = 'n', 
     ylim = c(0, 350), main = 'Robust Kraus', ylab = 'PM10', xlab = 'Hour', cex.main=1.5, cex.lab=1.2, cex.axis=1.2)
set.seed(30)
for (k in c(sample(c(1:nrow(com.x)), 40))) {
  lines(com.x[k, ], col = gray(0.8), lwd = 1.5)
}
abline(h = 100, lty = 2)
for (j in 1:length(Rkraus_outlier)) {
  lines(com.x[Rkraus_outlier, ][j, ], col = 2, lwd = 1.5)	
}
for (k in setdiff(prop_outlier,
                  intersect(prop_outlier, Rkraus_outlier))) {
  lines(x[k, ], col = 4, lwd = 3)
}
legend('topright',  c('Observation', 'Detected outliers', 'Only by proposed FPCA'), col=c(gray(0.8),2,4) ,lty=1, lwd=2,cex=1.6)








####################################################
### outlier detection에서 이것저것 test
####################################################
col.pall1[prop_outlier]
hist(col.pall1)
boxplot(col.pall1)
adjbox(col.pall1)
x[prop_outlier, ]

df <- cbind(col.pall1,
            col.pall2,
            col.pall3)
boxplot(df)
par(mfrow = c(1, 3))
adjbox(col.pall1)
adjbox(col.pall2)
adjbox(col.pall3)

par(mfrow = c(1, 3))
plot(col.pall1, ylim = range(as.numeric(df)))
plot(col.pall2, ylim = range(as.numeric(df)))
plot(col.pall3, ylim = range(as.numeric(df)))

summary(col.pall1)
summary(col.pall2)
summary(col.pall3)

q <- 0.8
quantile(col.pall1, q)
quantile(col.pall2, q)
quantile(col.pall3, q)

### Residual로 adjbox plot
# reconstructed curves
pred_reconstr <- list(
  predict(pca.ogk.obj, K = K),
  predict(pca.Mkraus.obj, K = K),   # reconstruction using K-L expansion
  predict(pca.yao.obj, K = K)
  # NA   # R-Kraus does not do reconstruction
)

# MISE of reconstruction
sse_reconstr <- lapply(pred_reconstr, function(method){
  rowMeans((method - x)^2, na.rm = T)
})
lapply(sse_reconstr, summary)
# q <- 0.9
# lapply(sse_reconstr, quantile, q)
# 
# cut_off <- 800
# sapply(sse_reconstr, function(method) {
#   length(which(method > cut_off))
# })

cut_off <- 3
adj_boxplot <- boxB(sse_reconstr[[1]], method = 'asymmetric', k = cut_off)
adj_boxplot <- boxB(sse_reconstr[[2]], method = 'asymmetric', k = cut_off)
adj_boxplot <- boxB(sse_reconstr[[3]], method = 'asymmetric', k = cut_off)

adj_boxplot <- boxB(sse_reconstr[[1]], method = 'adjbox')
adj_boxplot <- boxB(sse_reconstr[[2]], method = 'adjbox')
adj_boxplot <- boxB(sse_reconstr[[3]], method = 'adjbox')


cut_off <- 8
adj_boxplot <- boxB(sse_reconstr[[1]], method = 'asymmetric', k = cut_off)
prop_outlier <- adj_boxplot$outliers
place$city[prop_outlier] %>% unique() %>% sort()

adj_boxplot <- boxB(sse_reconstr[[2]], method = 'asymmetric', k = cut_off)
Rkraus_outlier <- adj_boxplot$outliers
place$city[Rkraus_outlier] %>% unique() %>% sort()

adj_boxplot <- boxB(sse_reconstr[[3]], method = 'asymmetric', k = cut_off)
yao_outlier <- adj_boxplot$outliers
place$city[yao_outlier] %>% unique() %>% sort()


### Plot
{
par(mfrow = c(1, 3))
set.seed(10)
ind_sample <- c(sample(c(1:nrow(com.x)), 40), 
                prop_outlier,
                Rkraus_outlier,
                yao_outlier)
plot(com.x[1, ], type = 'n',
     ylim = c(0, 350), main = 'Proposed FPCA', ylab = 'PM10', xlab = 'Hour')
set.seed(30)
for (k in ind_sample) {
  lines(com.x[k, ], col = gray(0.8), lwd = 2)
}
abline(h = 100, lty = 2)
for (j in 1:length(prop_outlier)) {
  lines(com.x[prop_outlier, ][j, ], col = 2)
}

plot(com.x[1, ], type = 'n', 
     ylim = c(0, 350), main = 'Sparse FPCA', ylab = 'PM10', xlab = 'Hour')
set.seed(30)
for (k in ind_sample) {
  lines(com.x[k, ], col = gray(0.8), lwd = 2)
}
abline(h = 100, lty = 2)
for (j in 1:length(yao_outlier)) {
  lines(com.x[yao_outlier, ][j, ], col = 2)	
}

plot(com.x[1, ], type = 'n', 
     ylim = c(0, 350), main = 'Robust Kraus', ylab = 'PM10', xlab = 'Hour')
set.seed(30)
for (k in ind_sample) {
  lines(com.x[k, ], col = gray(0.8), lwd = 2)
}
abline(h = 100, lty = 2)
for (j in 1:length(Rkraus_outlier)) {
  lines(com.x[Rkraus_outlier, ][j, ], col = 2)	
}
}


plot(colMeans(x[-prop_outlier, ], na.rm = T), type = "l", 
     lwd = 2, ylim = c(0, 350))
lines(colMeans(x[prop_outlier, ], na.rm = T), lwd = 2, col = 2)
plot(colMeans(x[-Rkraus_outlier, ], na.rm = T), type = "l", 
     lwd = 2, ylim = c(0, 350))
lines(colMeans(x[Rkraus_outlier, ], na.rm = T), lwd = 2, col = 2)
plot(colMeans(x[-yao_outlier, ], na.rm = T), type = "l", 
     lwd = 2, ylim = c(0, 350))
lines(colMeans(x[yao_outlier, ], na.rm = T), lwd = 2, col = 2)



par(mfrow = c(1, 3))
adjbox(sse_reconstr[[1]])
adjbox(sse_reconstr[[2]])
adjbox(sse_reconstr[[3]])
plot(density(sse_reconstr[[1]]))
plot(density(sse_reconstr[[2]]))
plot(density(sse_reconstr[[3]]))


plot(pcs_prop[, 1:2], xlim = c(-50, 130), ylim = c(-50, 50))
plot(pcs_Rkraus[, 1:2], xlim = c(-50, 130), ylim = c(-50, 50))
plot(pcs_yao[, 1:2], xlim = c(-50, 130), ylim = c(-50, 50))



source("./old_code/20211209 이전/function_outlier.R")
prop_outlier <- fun_outlier(pca.ogk.obj, method = "robMah")
prop_outlier <- which(prop_outlier == 1)
length(prop_outlier)

Rkraus_outlier <- fun_outlier(pca.Mkraus.obj, method = "robMah")
Rkraus_outlier <- which(Rkraus_outlier == 1)
length(Rkraus_outlier)

yao_outlier <- fun_outlier(pca.yao.obj, method = "robMah")
yao_outlier <- which(yao_outlier == 1)
length(yao_outlier)

par(mfrow = c(2, 4))
for (i in prop_outlier) {
  plot(com.x[i, ], type = 'l',
       ylim = c(0, 350), main = i, ylab = 'PM10', xlab = 'Hour')
}
plot(score_dist(pca.ogk.obj), orthogonal_dist(pca.ogk.obj))


setdiff(prop_outlier,
        intersect(prop_outlier, Rkraus_outlier))
setdiff(Rkraus_outlier,
        intersect(prop_outlier, Rkraus_outlier))

setdiff(prop_outlier,
        intersect(prop_outlier, yao_outlier))
setdiff(yao_outlier,
        intersect(prop_outlier, yao_outlier))

for (i in setdiff(prop_outlier,
                  intersect(prop_outlier, yao_outlier))) {
  plot(com.x[i, ], type = 'l',
       ylim = c(0, 350), main = i, ylab = 'PM10', xlab = 'Hour')
}


for (i in sample(setdiff(1:n, prop_outlier), 4)) {
  plot(com.x[i, ], type = 'l',
       ylim = c(0, 350), main = i, ylab = 'PM10', xlab = 'Hour')
}
for (i in sample(prop_outlier, 4)) {
  plot(com.x[i, ], type = 'l',
       ylim = c(0, 350), main = i, ylab = 'PM10', xlab = 'Hour')
}


par(mfrow = c(1, 3))
set.seed(10)
ind_sample <- c(
  order(apply(x[prop_outlier, ], 1, max, na.rm = T), decreasing = T)[1:20],
  order(apply(x[yao_outlier, ], 1, max, na.rm = T), decreasing = T)[1:20],
  order(apply(x[Rkraus_outlier, ], 1, max, na.rm = T), decreasing = T)[1:20],
  sample(setdiff(1:n, prop_outlier), 50)   # no outlier
)
plot(x[1, ], type = 'n',
     ylim = c(0, 350), main = 'Proposed FPCA', ylab = 'PM10', xlab = 'Hour')
for (k in ind_sample) {
  lines(x[k, ], col = gray(0.8), lwd = 2)
}
abline(h = 100, lty = 2)
for (j in order(apply(x[prop_outlier, ], 1, max, na.rm = T), decreasing = T)[1:20]) {
  lines(x[prop_outlier[j], ], col = 2)
}

plot(x[1, ], type = 'n', 
     ylim = c(0, 350), main = 'Sparse FPCA', ylab = 'PM10', xlab = 'Hour')
for (k in ind_sample) {
  lines(x[k, ], col = gray(0.8), lwd = 2)
}
abline(h = 100, lty = 2)
for (j in order(apply(x[yao_outlier, ], 1, max, na.rm = T), decreasing = T)[1:20]) {
  lines(x[yao_outlier[j], ], col = 2)
}

plot(x[1, ], type = 'n', 
     ylim = c(0, 350), main = 'Robust Kraus', ylab = 'PM10', xlab = 'Hour')
for (k in ind_sample) {
  lines(x[k, ], col = gray(0.8), lwd = 2)
}
abline(h = 100, lty = 2)
for (j in order(apply(x[Rkraus_outlier, ], 1, max, na.rm = T), decreasing = T)[1:20]) {
  lines(x[Rkraus_outlier[j], ], col = 2)
}




par(mfrow = c(2, 3))
{
set.seed(100)
ind <- sample(prop_outlier, 30)   # proposed sample
plot(x[1, ], type = 'n',
     ylim = c(0, 350), main = 'Proposed FPCA', ylab = 'PM10', xlab = 'Hour')
for (k in prop_outlier) {
  lines(x[k, ], col = k, lwd = 1.2)
}

plot(x[1, ], type = 'n', 
     ylim = c(0, 350), main = 'Proposed - Robust Kraus', ylab = 'PM10', xlab = 'Hour')
for (k in ind) {
  lines(x[k, ], col = gray(0.8), lwd = 2)
}
for (k in setdiff(prop_outlier,
                  intersect(prop_outlier, Rkraus_outlier))) {
  lines(x[k, ], col = 2, lwd = 2)
}

plot(x[1, ], type = 'n', 
     ylim = c(0, 350), main = 'Proposed - Sparse FPCA', ylab = 'PM10', xlab = 'Hour')
for (k in ind) {
  lines(x[k, ], col = gray(0.8), lwd = 2)
}
for (k in setdiff(prop_outlier,
                  intersect(prop_outlier, yao_outlier))) {
  lines(x[k, ], col = 2, lwd = 2)
}

plot(x[1, ], type = 'n', 
     ylim = c(0, 350), main = 'Robust Kraus - Proposed', ylab = 'PM10', xlab = 'Hour')
for (k in ind) {
  lines(x[k, ], col = gray(0.8), lwd = 2)
}
for (k in setdiff(Rkraus_outlier,
                  intersect(prop_outlier, Rkraus_outlier))) {
  lines(x[k, ], col = 2, lwd = 2)
}

plot(x[1, ], type = 'n', 
     ylim = c(0, 350), main = 'Sparse FPCA - Proposed', ylab = 'PM10', xlab = 'Hour')
for (k in ind) {
  lines(x[k, ], col = gray(0.8), lwd = 2)
}
for (k in setdiff(yao_outlier,
                  intersect(prop_outlier, yao_outlier))) {
  lines(x[k, ], col = 2, lwd = 2)
}
}

## 각각 그리기
par(mfrow = c(3, 4))
for (k in 1:12) {
  ind <- 10*(k-1) + (1:10)
  matplot(t(x[prop_outlier[ind], ]), type = "l", lty = 1)
  # plot(x[1, ], type = 'n',
  #      ylim = c(0, 350), main = 'Proposed FPCA', ylab = 'PM10', xlab = 'Hour')
  # lines(x[, ], col = k, lwd = 1.2)
}




plot(full_data$shape, main='Proposed')
rbPal <- colorRampPalette(c("green","yellow","red"))
col1 <- rep(1, length(col.pall1))
col1[prop_outlier]=2
points((result_aa[,1]),(result_aa[,2]), pch=19, col=col1,cex=1)
points((result_aa[prop_outlier,1]),(result_aa[prop_outlier,2]), pch=19, col=2,cex=2)
box()


plot(full_data$shape, main='Kraus')
rbPal <- colorRampPalette(c("green","yellow","red"))
col1 <- rep(1, length(col.pall2))
col1[kraus_outlier]=2
points((result_aa[,1]),(result_aa[,2]), pch=19, col=col1,cex=1)
points((result_aa[kraus_outlier,1]),(result_aa[kraus_outlier,2]), pch=19, col=2,cex=2)
box()



plot(full_data$shape, main='Robust Kraus')
rbPal <- colorRampPalette(c("green","yellow","red"))
col1 <- rep(1, length(col.pall3))
col1[kraus_outlier]=2
points((result_aa[,1]),(result_aa[,2]), pch=19, col=col1,cex=1)
points((result_aa[Rkraus_outlier,1]),(result_aa[Rkraus_outlier,2]), pch=19, col=2,cex=2)
box()








plot(full_data$shape, main='Average PM10')
col.pall= data_avg[is.na(data_avg)==F]
col1 <- rbPal(br)[as.numeric(cut(col.pall,breaks =br2))]
points(data_avg[,2], data_avg[,3], pch=19, col=col1,cex=2)
legend('bottomright',  levels(cut(col.pall,breaks = br)), col=(rbPal(br)) ,pch=19)
box()


plot(full_data$shape, main='Maximum PM10')
col.pall= data_max[is.na(data_avg)==F]
col1 <- rbPal(br)[as.numeric(cut(col.pall,breaks =br2))]
points(data_avg[,2], data_avg[,3], pch=19, col=col1,cex=2)
legend('bottomright',  levels(cut(col.pall,breaks = br)), col=(rbPal(br)) ,pch=19)
box()





full_data$shape_new <- subset(full_data$shape, full_data$shape@data$CTP_ENG_NM=="Seoul" | full_data$shape@data$CTP_ENG_NM=="Gyeonggi-do" )

plot(full_data$shape_new, main='Proposed')
br=10 
br2=10
col1 <- rbPal(br)[as.numeric(cut(avg_pcs_prop,breaks =br2))] #colorRampPalette(c("red","yellow","green"))(length(col.pall))
points((result_aa[unique_index,1]),(result_aa[unique_index,2]), pch=19, col=col1,cex=2)
legend('bottomleft',  levels(cut(col.pall1,breaks = br)), col=(rbPal(br)) ,pch=19)
box()


plot(full_data$shape_new, main='Conventinal PCA')
col1 <- rbPal(br)[as.numeric(cut(avg_pcs_conv,breaks = br2))]
points((result_aa[unique_index,1]),(result_aa[unique_index,2]), pch=19, col=col1,cex=2)
legend('bottomleft',  levels(cut(col.pall2,breaks = br)), col=(rbPal(br)) ,pch=19)
box()

plot(full_data$shape_new, main='Kraus')
col1 <- rbPal(10)[as.numeric(cut(avg_pcs_kraus,breaks =br2))]
points((result_aa[unique_index,1]),(result_aa[unique_index,2]), pch=19, col=col1,cex=2)
legend('bottomleft',  levels(cut(col.pall3,breaks = br)), col=(rbPal(br)) ,pch=19)
box()


plot(full_data$shape_new, main='Average PM10')
col.pall= data_avg[is.na(data_avg)==F]
col1 <- rbPal(br)[as.numeric(cut(col.pall,breaks =br2))]
points(data_avg[,2], data_avg[,3], pch=19, col=col1,cex=2)
legend('bottomright',  levels(cut(col.pall,breaks = br)), col=(rbPal(br)) ,pch=19)
box()


plot(full_data$shape_new, main='Maximum PM10')
col1 <- rbPal(br)[as.numeric(cut(col.pall,breaks =br2))]
points(data_avg[,2], data_avg[,3], pch=19, col=col1,cex=2)
legend('bottomright',  levels(cut(col.pall,breaks = br)), col=(rbPal(br)) ,pch=19)
box()






hist(pca.ogk.obj$pc.score[, num_pc][-7132])
hist(pcs_conv[, num_pc])
hist(pcs_kraus[, num_pc])










data_avg=data_max= matrix(nrow=length(loc), ncol=3)
for(loc_num in 1:length(loc)){
  
  data= filter(new_data, new_data[,1] == loc[loc_num])
  
  if(  length(which(data[,2]=='2017030101')) !=0){
    time_index= c( which(data[,2]=='2017030101') : which(data[,2]=='2017033124') )
    data_avg[loc_num,1  ]= mean( data[time_index,3 ], na.rm=T  ) 
    data_max[loc_num,1  ]= max( data[time_index,3 ], na.rm=T  )
    data_max[loc_num ,c(2:3) ]=data_avg[loc_num ,c(2:3) ]= as.numeric(data[time_index[1],5:6 ])
  }
  
  
}

