
data <- foutlier_sim_mfd(n = 500, m = 51, p = 20, outlier_rate = 0.1, model = 5, rho = 0.3)
# length(data$data)
# data$idx_outliers
plot(data, p = 1)



library(mvtnorm)

n <- 100
m <- 51
p <- 20

rho <- 0.3
n_comp <- 9

gr <- seq(0, 1, length.out = m)


# autoregressive correlation matrix
ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                    (1:n - 1))
  rho^exponent
}
sigma <- ar1_cor(p, rho)*0.49

# nu_k <- sqrt(exp(-(1:p)/4))

j <- 2

X <- array(NA, dim = c(n, m, p))
# u <- sapply(1:p, function(i){ 0.5*sin(i/2*pi*gr) + 0.5*cos(i/2*pi*gr) })
# u <- sapply(c(rep(4, 4), 5:20), function(i){ 0.5*sin(i/2*pi*gr) + 0.5*cos(i/2*pi*gr) })
# u <- 4*gr
u <- 0.5*sin(p/4*pi*gr) + 0.5*cos(p/4*pi*gr)
for (i in 1:n) {
  eps <- rmvnorm(m, sigma = sigma)
  # u <- runif(p, -1.1, 1.1)
  # u <- rnorm(1, sd = 0.5)
  X[i, , ] <- eps + u
}
matplot(t(X[, , j]), type = "l", lty = 1, col = "gray")


X <- array(NA, dim = c(n, m, p))
# u <- sapply(1:p, function(i){ ((i+1) %% 2)*sin(i/2*pi*gr) + (i %% 2)*cos(i/2*pi*gr) })
# u <- sapply(1:p, function(i){ (i %% 2)*sin(i*pi*gr) + ((i+1) %% 2)*cos(i*pi*gr) })
# u <- sapply(1:p, function(i){ 0.5*sin(i*pi*gr) + 0.5*cos(i*pi*gr) })
u <- 0.5*sin(p/2*pi*gr) + 0.5*cos(p/2*pi*gr)
for (i in 1:n) {
  eps <- rmvnorm(m, sigma = sigma)
  # u <- sin(4*pi*gr)
  # u <- cos(8*pi*gr)
  X[i, , ] <- eps + u
}
matlines(t(X[, , j]), type = "l", lty = 1, col = 2)





mu <- 4*gr
# mu <- 5*sin(2*pi*gr)
# mu <- 5*(gr-1)^2
mu <- 4*gr + 0.5*sin(p/4*pi*gr) + 0.5*cos(p/4*pi*gr)
n_comp <- 9
basis_obj <- create.fourier.basis(nbasis = n_comp)
basis_ftn <- eval.basis(gr, basis_obj) %>% 
  apply(2, function(x){ x / sqrt(sum(x^2)) })

# nu_k <- sqrt(exp(-(1:n_comp)/4))
nu_k <- (n_comp + 1 - (1:n_comp)) / n_comp
# sigma_eps <- diag(runif(p, 0.1, 0.3))
sigma_eps <- ar1_cor(m, rho)*0.1
# sigma <- ar1_cor(p, rho)*0.49
sigma <- ar1_cor(n_comp, rho)*diag(nu_k)


X <- array(NA, dim = c(n, m, p))
for (i in 1:n) {
  # z_ij <- diag(nu_k) %*% rmvnorm(n_comp, sigma = sigma)
  z_ij <- rmvnorm(n_comp, sigma = sigma)
  eps <- rmvnorm(1, sigma = sigma_eps) %>% 
    as.numeric()
  X[i, , ] <- (basis_ftn %*% z_ij) + mu + eps
}
matplot(t(X[, , 1]), type = "l", lty = 1, col = "gray")



mu <- 4*gr
# mu <- 5*sin(2*pi*(gr-0.3))
# mu <- 0.5*sin(p/2*pi*gr) + 0.5*cos(p/2*pi*gr)
# mu <- 0.5*sin(p/2*pi*(gr-0.15)) + 0.5*cos(p/2*pi*(gr-0.15))
mu <- 4*gr + 0.3*sin(p/2*pi*(gr-0.15)) + 0.3*cos(p/2*pi*(gr-0.15))

n_comp <- 9
basis_obj <- create.fourier.basis(nbasis = n_comp)
basis_ftn <- eval.basis(gr, basis_obj) %>% 
  apply(2, function(x){ x / sqrt(sum(x^2)) })
X <- array(NA, dim = c(n, m, p))
# nu_k <- sqrt(exp(-(1:n_comp)/4))
for (i in 1:n) {
  # z_ij <- diag(nu_k) %*% rmvnorm(n_comp, sigma = sigma)
  z_ij <- rmvnorm(n_comp, sigma = sigma)
  # eps <- rmvnorm(m, sigma = sigma_eps)
  eps <- rmvnorm(1, sigma = sigma_eps) %>% 
    as.numeric()
  X[i, , ] <- (basis_ftn %*% z_ij) + mu + eps
}
matlines(t(X[, , 1]), type = "l", lty = 1, col = 2)






#-----------------------------------------------------
covfunexp <- function(gridpoints, alpha, beta, nu){
  d_matrix <- as.matrix(dist(gridpoints, upper = T, diag = T))
  return(alpha*exp(-beta*(d_matrix^nu)))
}


library(mvtnorm)

n <- 100
m <- 51
p <- 20

gr <- seq(0, 1, length.out = m)

n <- 100
m <- 51
p <- 20

rho <- 0.3
rho <- 0

# mu <- 5*sin(2*pi*gr)
# mu <- 5*(gr-1)^2
# mu <- 4*gr + 0.5*sin(p/4*pi*gr) + 0.5*cos(p/4*pi*gr)
# mu <- 0.5*sin(p/2*pi*gr) + 0.5*cos(p/2*pi*gr)
mu <- 0
n_comp <- 15
# basis_obj <- create.fourier.basis(nbasis = n_comp)
basis_obj <- create.bspline.basis(nbasis = n_comp)
basis_ftn <- eval.basis(gr, basis_obj) %>% 
  apply(2, function(x){ x / sqrt(sum(x^2)) })

# mu <- 4*gr
# nu_k <- sqrt(exp(-(1:n_comp)/4))
# nu <- (n_comp + 1 - (1:n_comp)) / n_comp * 5
nu <- (n_comp + 1 - (1:n_comp)) / n_comp
# sigma_eps <- diag(runif(p, 0.1, 0.3))
# sigma_eps <- ar1_cor(m, rho)*0.1
# sigma <- ar1_cor(p, rho)*0.49
sigma <- ar1_cor(p, rho)
sigma_eps <- covfunexp(gr, alpha = 1, beta = 1, nu = 1)
X <- array(NA, dim = c(n, m, p))
for (i in 1:n) {
  z_ij <- sapply(nu, function(nu_k){
    rmvnorm(1, sigma = sigma*nu_k)
  })
  eps <- rmvnorm(1, sigma = sigma_eps) %>%
    as.numeric()
  X[i, , ] <- (basis_ftn %*% t(z_ij)) + mu + eps
  # X[i, , ] <- (basis_ftn %*% t(z_ij)) + mu
}
matplot(t(X[, , 1]), type = "l", lty = 1, col = "gray")


# mu <- 4*gr + 0.1*sin(p*pi*(gr-0.15)) + 0.1*cos(p*pi*(gr-0.15))
# mu <- 4*gr + 0.5*sin(2*p*pi*gr) + 0.5*cos(2*p*pi*gr)
# mu <- 4*gr
# mu <- 0.5*sin(2*p*pi*gr) + 0.5*cos(2*p*pi*gr)
# sigma_eps <- diag(runif(p, 0.4, 0.6))
X <- array(NA, dim = c(n, m, p))
# nu <- (n_comp + 1 - (1:n_comp)) / n_comp * 0.1
# nu <- (n_comp + 1 - (1:n_comp)) / n_comp
# nu <- nu*5
sigma_eps <- covfunexp(gr, alpha = 0.5, beta = 2, nu = 0.2)
# sigma_eps <- covfunexp(gr, alpha = 2, beta = 2, nu = 0.5)
for (i in 1:10) {
  z_ij <- sapply(nu, function(nu_k){
    rmvnorm(1, sigma = sigma*nu_k)
    # rmvt(1, sigma = sigma*nu_k, df = 2)
  })
  eps <- rmvnorm(1, sigma = sigma_eps) %>%
    as.numeric()
  # eps <- rmvnorm(p, sigma = sigma_eps) %>% t()
  X[i, , ] <- (basis_ftn %*% t(z_ij)) + mu + eps
  # X[i, , ] <- (basis_ftn %*% t(z_ij)) + mu
  
  # # t_i <- runif(1, 0.1, 0.9)
  # for (j in 1:p) {
  #   t_i <- runif(1, 0.1, 0.9)
  #   # X[i, which(gr >= t_i), j] <- X[i, which(gr >= t_i), j] + sample(c(-1, 1), 1) * 2
  #   X[i, which(gr >= t_i & gr <= t_i + 0.05), j] <- X[i, which(gr >= t_i & gr <= t_i + 0.05), j] + sample(c(-1, 1), 1) * 1.5
  # }
  # # t_i <- runif(p, 0.1, 0.9)
  # # for (j in 1:p) {
  # #   X[i, which(gr >= t_i[j]), j] <- X[i, which(gr >= t_i[j]), j] + sample(c(-1, 1), 1) * 2
  # #   # X[i, which(gr >= t_i[j] & gr <= t_i[j] + 0.05), j] <- X[i, which(gr >= t_i[j] & gr <= t_i[j] + 0.05), j] + sample(c(-1, 1), 1) * 3
  # # }
}
matlines(t(X[, , 1]), type = "l", lty = 1, col = 2)




# matplot(X[1, , ], type = "l", lty = 1, col = "gray")






#-----------------------------------------------------


library(mvtnorm)

n <- 100
m <- 51
p <- 20

gr <- seq(0, 1, length.out = m)

rho <- 0.3
rho <- 0

mu <- 0.5*sin(p/2*pi*gr) + 0.5*cos(p/2*pi*gr)
sigma <- ar1_cor(p, rho)
sigma_eps <- covfunexp(gr, alpha = 1, beta = 1, nu = 1)
X <- array(NA, dim = c(n, m, p))
for (i in 1:n) {
  eps1 <- sapply(1:m, function(j){
    rmvnorm(1, sigma = sigma)
  })
  eps2 <- rmvnorm(1, sigma = sigma_eps) %>%
    as.numeric()
  X[i, , ] <- t(eps1) + mu + eps2
}
matplot(t(X[, , 1]), type = "l", lty = 1, col = "gray")


mu <- 0.5*sin(p*pi*(gr-0.3)) + 0.5*cos(p*pi*(gr-0.3))
X <- array(NA, dim = c(n, m, p))
sigma_eps <- covfunexp(gr, alpha = 1, beta = 2, nu = 0.2)
# sigma_eps <- covfunexp(gr, alpha = 2, beta = 2, nu = 0.5)
for (i in 1:10) {
  eps1 <- sapply(1:m, function(j){
    rmvnorm(1, sigma = sigma)
  })
  eps2 <- rmvnorm(1, sigma = sigma_eps) %>%
    as.numeric()
  X[i, , ] <- t(eps1) + mu + eps2
  
  # # t_i <- runif(1, 0.1, 0.9)
  # for (j in 1:p) {
  #   t_i <- runif(1, 0.1, 0.9)
  #   # X[i, which(gr >= t_i), j] <- X[i, which(gr >= t_i), j] + sample(c(-1, 1), 1) * 2
  #   X[i, which(gr >= t_i & gr <= t_i + 0.05), j] <- X[i, which(gr >= t_i & gr <= t_i + 0.05), j] + sample(c(-1, 1), 1) * 1.5
  # }
  # # t_i <- runif(p, 0.1, 0.9)
  # # for (j in 1:p) {
  # #   X[i, which(gr >= t_i[j]), j] <- X[i, which(gr >= t_i[j]), j] + sample(c(-1, 1), 1) * 2
  # #   # X[i, which(gr >= t_i[j] & gr <= t_i[j] + 0.05), j] <- X[i, which(gr >= t_i[j] & gr <= t_i[j] + 0.05), j] + sample(c(-1, 1), 1) * 3
  # # }
}
matlines(t(X[, , 1]), type = "l", lty = 1, col = 2)



# matplot(X[1, , ], type = "l", lty = 1, col = "gray")








data <- data.table::fread("../../../adhd200_preprocessed_phenotypics.tsv")
data <- as.data.frame(data)
head(data)

table(data$QC_Athena)
table(data$QC_NIAK)

df <- data %>% 
  filter(Site == 1) %>% 
  # filter(QC_Athena == 1) %>% 
  dplyr::select(`ScanDir ID`, Gender, Age, DX) %>% 
  mutate(DX = ifelse(DX == 0, 0, 1))

table(df$DX)


df$`ScanDir ID` %>% sort()

4334113 %in% df$`ScanDir ID`



data <- rbind(
  read_csv("/Users/hyunsung/Desktop/ADHD200_training_pheno/Peking_1/Peking_1_phenotypic.csv"),
  read_csv("/Users/hyunsung/Desktop/ADHD200_training_pheno/Peking_2/Peking_2_phenotypic.csv"),
  read_csv("/Users/hyunsung/Desktop/ADHD200_training_pheno/Peking_3/Peking_3_phenotypic.csv")
)
head(data)
dim(data)
data$DX %>% table()
data$QC_Anatomical_1 %>% table()
data$QC_Rest_4 %>% table()

idx <- data$`ScanDir ID`
data <- read_tsv("../../../adhd200_preprocessed_phenotypics.tsv")
head(data)
data[which(data$`ScanDir ID` %in% idx), ]$QC_Athena %>% table()
data[which(data$`ScanDir ID` %in% idx), ]$QC_NIAK %>% table()






summary(nonconform_score_calib)
summary(nonconform_score_test)

mean_calib <- colMeans(nonconform_score_calib)
sd_calib <- apply(nonconform_score_calib, 2, sd)

nonconform_score_calib_scaled <- t( (t(nonconform_score_calib) - mean_calib) / sd_calib )
nonconform_score_test_scaled <- t( (t(nonconform_score_test) - mean_calib) / sd_calib )
summary(nonconform_score_calib_scaled)
summary(nonconform_score_test_scaled)

par(mfrow = c(1, 2))
hist(nonconform_score_calib_scaled[, 1])
hist(nonconform_score_test_scaled[, 1])
hist(nonconform_score_calib_scaled[, 2])
hist(nonconform_score_test_scaled[, 2])
hist(nonconform_score_calib_scaled[, 3])
hist(nonconform_score_test_scaled[, 3])

nonconform_score_calib2 <- apply(nonconform_score_calib_scaled, 1, mean)
nonconform_score_test2 <- apply(nonconform_score_test_scaled, 1, mean)
plot(nonconform_score_test2)
points(nonconform_score_calib2, col = 2)

nonconform_score_calib2 <- apply(nonconform_score_calib, 1, mean)
nonconform_score_test2 <- apply(nonconform_score_test, 1, mean)
plot(nonconform_score_test2)
points(nonconform_score_calib2, col = 2)

apply(nonconform_score_calib, 2, sd)
apply(nonconform_score_calib, 2, range)

apply(nonconform_score_test, 2, sd)
apply(nonconform_score_test, 2, range)

apply(apply(arr_calib_trans, c(2,3), sd), 1, mean) %>% range
apply(apply(arr_test_trans, c(2,3), sd), 1, mean) %>% range


# Conformal p-value (marginal)
conf_pvalue_marg <- sapply(nonconform_score_test2, function(s){
  (1 + sum(nonconform_score_calib2 >= s)) / (n_calib + 1)
})

which(p.adjust(conf_pvalue_marg, method = "BH") < 0.1)


apply(conf_pvalue, 2, rank)[idx_bh$marginal, ]
conf_pvalue[idx_bh$marginal, ]




### fMRI
# 1st, 2nd derivatives
X_deriv_1 <- array(NA, c(m-1, n, p))
X_deriv_2 <- array(NA, c(m-2, n, p))
for (i in 1:p) {
  X_deriv_1[, , i] <- apply(X[, , i], 1, diff)
  X_deriv_2[, , i] <- apply(X[, , i], 1, function(x){ diff(diff(x)) })
}



par(mfrow = c(3, 1))
matplot(t(X[which(y == 0)[1:5], , 1]), type = "l")
matlines(t(X[c(33, 123), , 1]), type = "l", lwd = 2)
matplot(X_deriv_1[, which(y == 0)[1:5], 1], type = "l")
matlines(X_deriv_1[, c(33, 123), 1], type = "l", lwd = 2)
matplot(X_deriv_2[, which(y == 0)[1:5], 1], type = "l")
matlines(X_deriv_2[, c(33, 123), 1], type = "l", lwd = 2)




type_depth <- "projdepth"
type_depth <- "hdepth"

depth_values_1 <- mfd(aperm(X, c(2,1,3)), type = type_depth,
                      depthOptions = list(type = "Rotation"))
depth_values_1$MFDdepthX[c(33, 123)]
obj <- boxplot(depth_values_1$MFDdepthX, plot = FALSE)
idx_1 <- which(depth_values_1$MFDdepthX < obj$stats[1, 1])

depth_values_2 <- mfd(X_deriv_1, type = type_depth,
                      depthOptions = list(type = "Rotation"))
depth_values_2$MFDdepthX[c(33, 123)]
obj <- boxplot(depth_values_2$MFDdepthX, plot = FALSE)
idx_2 <- which(depth_values_2$MFDdepthX < obj$stats[1, 1])

depth_values_3 <- mfd(X_deriv_2, type = type_depth,
                      depthOptions = list(type = "Rotation"))
depth_values_3$MFDdepthX[c(33, 123)]
obj <- boxplot(depth_values_3$MFDdepthX, plot = FALSE)
idx_3 <- which(depth_values_3$MFDdepthX < obj$stats[1, 1])

par(mfrow = c(1, 3))
plot(depth_values_1$MFDdepthX, col = y+1, ylab = "")
points(idx_adhd_cand, depth_values_1$MFDdepthX[idx_adhd_cand], col = "blue", pch = 10)
points(c(33, 123), depth_values_1$MFDdepthX[c(33, 123)], col = "green", pch = 10)

plot(depth_values_2$MFDdepthX, col = y+1, ylab = "")
points(idx_adhd_cand, depth_values_2$MFDdepthX[idx_adhd_cand], col = "blue", pch = 10)
points(c(33, 123), depth_values_2$MFDdepthX[c(33, 123)], col = "green", pch = 10)

plot(depth_values_3$MFDdepthX, col = y+1, ylab = "")
points(idx_adhd_cand, depth_values_3$MFDdepthX[idx_adhd_cand], col = "blue", pch = 10)
points(c(33, 123), depth_values_3$MFDdepthX[c(33, 123)], col = "green", pch = 10)


cor(cbind(depth_values_1$MFDdepthX,
          depth_values_2$MFDdepthX,
          depth_values_3$MFDdepthX))

cbind(
  order(depth_values_1$MFDdepthX),
  order(depth_values_2$MFDdepthX),
  order(depth_values_3$MFDdepthX),
  order((depth_values_1$MFDdepthX + depth_values_2$MFDdepthX + depth_values_3$MFDdepthX)/3)
) %>%
  head(20) %>%
  apply(2, sort)

cbind(
  order(depth_values_1$MFDdepthX),
  order(depth_values_2$MFDdepthX),
  order(depth_values_3$MFDdepthX),
  order((depth_values_1$MFDdepthX + depth_values_2$MFDdepthX + depth_values_3$MFDdepthX)/3)
) %>%
  head(20) %>%
  as.integer() %>%
  unique() %>%
  length()


cbind(
  summary(depth_values_1$MFDdepthX),
  summary(depth_values_2$MFDdepthX),
  summary(depth_values_3$MFDdepthX)
)
cbind(
  sd(depth_values_1$MFDdepthX),
  sd(depth_values_2$MFDdepthX),
  sd(depth_values_3$MFDdepthX)
)


idx_1
idx_2
idx_3

obj <- boxplot(depth_values$MFDdepthX, plot = FALSE)
cand <- which(depth_values$MFDdepthX > obj$stats[5, 1])
cand[cand %in% idx_adhd]

c(idx_1, idx_2, idx_3) %>% unique()





depth_values <- c(
  mfd(aperm(X[y == 0, , ], c(2,1,3)), type = type_depth,
      depthOptions = list(type = "Rotation"))$MFDdepthX,
  mfd(aperm(X[y == 1, , ], c(2,1,3)), type = type_depth,
      depthOptions = list(type = "Rotation"))$MFDdepthX
)
plot(depth_values, col = c(rep(1, sum(y == 0)), rep(2, sum(y == 1))))
# abline(h = mean(depth_values[y == 0]), col = 1)
# abline(h = mean(depth_values[y == 1]), col = 2)
abline(h = median(depth_values[y == 0]), col = 1)
abline(h = median(depth_values[y == 1]), col = 2)

depth_values <- c(
  mfd(X_deriv_1[, y == 0, ], type = type_depth,
      depthOptions = list(type = "Rotation"))$MFDdepthX,
  mfd(X_deriv_1[, y == 1, ], type = type_depth,
      depthOptions = list(type = "Rotation"))$MFDdepthX
)
plot(depth_values, col = c(rep(1, sum(y == 0)), rep(2, sum(y == 1))))
# abline(h = mean(depth_values[y == 0]), col = 1)
# abline(h = mean(depth_values[y == 1]), col = 2)
abline(h = median(depth_values[y == 0]), col = 1)
abline(h = median(depth_values[y == 1]), col = 2)

depth_values <- c(
  mfd(X_deriv_2[, y == 0, ], type = type_depth,
      depthOptions = list(type = "Rotation"))$MFDdepthX,
  mfd(X_deriv_2[, y == 1, ], type = type_depth,
      depthOptions = list(type = "Rotation"))$MFDdepthX
)
plot(depth_values, col = c(rep(1, sum(y == 0)), rep(2, sum(y == 1))))
# abline(h = mean(depth_values[y == 0]), col = 1)
# abline(h = mean(depth_values[y == 1]), col = 2)
abline(h = median(depth_values[y == 0]), col = 1)
abline(h = median(depth_values[y == 1]), col = 2)


x <- X[, 2, ]
dim(x)
obj <- hdepth(x[1:70, ], options = list(type = "Rotation"))
obj <- hdepth(x)
obj <- hdepth(x, options = list(type = "Rotation"))
summary(obj$depthX)
plot(obj$depthX)







### Other real dataset

df <- read.table("/Users/hyunsung/Desktop/secom/secom.data", header = F)
dim(df)
y <- read.table("/Users/hyunsung/Desktop/secom/secom_labels.data", header = F)
dim(y)



df <- data.table::fread("/Users/hyunsung/Desktop/household_power_consumption.txt") %>% 
  as_tibble()

df <- data.table::fread("/Users/hyunsung/Desktop/real data/LD2011_2014.txt") 
dim(df)
colnames(df)
head(df[, 1:5])
tail(df[, 1:5])
rm(df)

df <- data.table::fread("/Users/hyunsung/Desktop/real data/Wafer/Wafer_TRAIN.txt") 
dim(df)
colnames(df)
head(df)
df$V1 %>% table
plot(as.numeric(df[1, -1]), type = "l")


df <- read.table("/Users/hyunsung/Desktop/human+activity+recognition+using+smartphones/UCI HAR Dataset/train/X_train.txt", header = F)
dim(df)
y <- read.table("/Users/hyunsung/Desktop/human+activity+recognition+using+smartphones/UCI HAR Dataset/train/y_train.txt", header = F)
y <- unlist(y)
table(y)
id <- read.table("/Users/hyunsung/Desktop/human+activity+recognition+using+smartphones/UCI HAR Dataset/train/subject_train.txt", header = F)
table(id)

idx <- which(id == 1 & y == 3)
df[idx, 1] %>% plot(type = "l")



"https://anson.ucdavis.edu/~mueller/data/medfly1000.txt"


path <- "/Users/hyunsung/Desktop/259.전력신사업을 위한 전라남도 지역 전력소비패턴 데이터 구축/01-1.정식개방데이터/Training/01.원천데이터/TS_1.suncheon/"
f_names <- list.files(path)

X <- matrix(NA, 100, 365)
for (i in 1:100) {
  df <- read_csv(paste0(path, f_names[i]), show_col_types = F) %>% 
    mutate(date = format(mrdDt, "%Y%m%d")) %>% 
    group_by(date) %>% 
    summarise(pwr = mean(pwrQrt)) %>% 
    mutate(time = row_number())
  # X[i, ] <- smooth.spline(df$time, df$pwr)$y
  X[i, ] <- df$pwr
}
X <- X[1:40, ]
matplot(t(X), type = "l")


df <- read_csv(paste0(path, f_names[1]))
plot(df$mrdDt, df$pwrQrt, type = "l")

df3 <- df %>% 
  mutate(date = format(mrdDt, "%Y%m%d")) %>% 
  group_by(date) %>% 
  summarise(pwr = mean(pwrQrt)) %>% 
  mutate(time = row_number())
plot(df3$time, df3$pwr, type = "l")
fit <- smooth.spline(df3$time, df3$pwr)
fit$x %>% length
lines(fit, col = 2, lwd = 2)

df2 <- df %>% 
  mutate(date = format(mrdDt, "%Y%m%d"))
plot(df2$pwrQrt[df2$date == "20210902"], type = "l")
plot(df2$pwrQrt[df2$date == "20210901"], type = "l")



df <- readxl::read_excel("/Users/hyunsung/Desktop/Dataset Population-Level Wearable Device Information/Dataset_Health.xlsx")
dim(df)
head(df)

column_names <- colnames(df)
var_names <- sapply(strsplit(column_names, "_"), function(x){ x[1] })
data <- list(
  df[, which(var_names == "Steps")],
  df[, which(var_names == "HeartRate")],
  df[, which(var_names == "SleepDuration")]
)
data <- lapply(data, as.matrix)
data[[2]] %>% dim

matplot(t(data[[1]][1:30, ]), type = "l", lty = 1)
matplot(t(data[[2]][1:30, ]), type = "l", lty = 1)
matplot(t(data[[3]][1:30, ]), type = "l", lty = 1)

n <- nrow(data[[1]])
m <- ncol(data[[1]])
p <- 3


obj <- get_clean_null(X = data, type = "depth_transform", type_depth = "projdepth")
setdiff(1:n, obj$idx_clean_null)

obj <- get_clean_null(X = data, type = "esssup")
setdiff(1:n, obj$idx_clean_null)

arr_df <- abind::abind(data, along = 3)
# Outlier detection for mixed training set
idx_ms <- msplot(dts = arr_df, plot = F, seed = b)$outliers
idx_seq <- seq_transform(arr_df, sequence = "O", depth_method = "erld",
                         erld_type = "one_sided_right", seed = b)$outliers$O



matplot(t(data[[1]][setdiff(1:n, obj$idx_clean_null), ]), type = "l", lty = 1, col ="red")
matlines(t(data[[1]][1:30, ]), type = "l", lty = 1, col = "gray")



set.seed(10)
### Split data into training and test set
alpha <- 0.2  # coverage level
prop_train <- 0.5  # proportion of training set

n_train <- round(n * prop_train)   # size of training set
n_test <- n - n_train   # size of test set

# Split training and test data
idx_train <- sample(1:n, n_train)
idx_test <- setdiff(1:n, idx_train)

data_train <- lapply(data, function(x){ x[idx_train, ] })
data_test <- lapply(data, function(x){ x[idx_test, ] })


b <- 1
cp_obj <- split_conformal_fd(X = data_train, X_test = data_test,
                             type = "depth_transform", type_depth = "projdepth",
                             train_type = "mixed",
                             seed = b)
conf_pvalue <- cp_obj$conf_pvalue

# BH procedure
idx_bh <- apply(conf_pvalue, 2, function(x){ 
  if (sum(sort(x) < (1:n_test)/n_test * alpha) == 0) {
    return(NA)
  } else {
    order(x)[1:max(which(sort(x) < (1:n_test)/n_test * alpha))]
  }
}, simplify = F)
idx_bh


matplot(t(data_test[[1]][idx_bh$marginal, ]), type = "l", lty = 1, col ="red")
matlines(t(data_test[[1]][1:50, ]), type = "l", lty = 1, col = "gray")











N <- 2e2
P <- 1e3
grid <- seq(0, 1, length.out = P)
Cov <- exp_cov_function(grid, alpha = 0.3, beta = 0.4)
Data <- list()
Data[[1]] <- generate_gauss_fdata(
  N,
  centerline = sin(2 * pi * grid),
  Cov = Cov
)
Data[[2]] <- generate_gauss_fdata(
  N,
  centerline = sin(2 * pi * grid),
  Cov = Cov
)
names <- paste0("id_", 1:nrow(Data[[1]]))
DG1 <- depthgram(Data, marginal_outliers = T)
fD <- fData(grid, Data[[1]])
DG2 <- depthgram(fD, marginal_outliers = TRUE, ids = names)
mfD <- mfData(grid, Data)
DG3 <- depthgram(mfD, marginal_outliers = TRUE, ids = names)
DG1$mag.out.det
DG1$shp.out.det


depth_values <- sapply(1:p, function(i){
  MBD_relative(data_test[[i]], data_train[[i]])
})
apply(depth_values, 1, mean)


MBD_relative(data_test[[1]], data_train[[1]])



library(roahd)
obj <- depthgram(data_test, marginal_outliers = T)
# # obj <- depthgram(mfData(gr, data_test), ids = names, marginal_outliers = F)
# obj$shp.out.det
# obj$mag.out.det
# 
# idx_test[obj$shp.out.det %>% 
#   unlist() %>% 
#   unique() %>% 
#   sort()]
# idx_test[obj$mag.out.det %>% 
#   unlist() %>% 
#   unique() %>% 
#   sort()]
get_fdr(idx_test[unique(c(unlist(obj$shp.out.det), unlist(obj$mag.out.det)))],
        idx_outliers)
get_tpr(idx_test[unique(c(unlist(obj$shp.out.det), unlist(obj$mag.out.det)))],
        idx_outliers)


obj <- depthgram(data, marginal_outliers = T)
# plot(obj)
# obj$shp.out.det %>% 
#   unlist() %>% 
#   unique() %>% 
#   sort()
# obj$mag.out.det %>% 
#   unlist() %>% 
#   unique() %>% 
#   sort()
get_fdr(unique(c(unlist(obj$shp.out.det), unlist(obj$mag.out.det))),
        idx_outliers)
get_tpr(unique(c(unlist(obj$shp.out.det), unlist(obj$mag.out.det))),
        idx_outliers)


seqobj <- seq_transform(df, sequence = "O", depth_method = "erld",
                        erld_type = "one_sided_right", save_data = TRUE)

dts <- fdaoutlier:::outlyingness_curves(df,  n_projections = 200L, seed = NULL)
t0_outliers <- functional_boxplot(dts, depth_method = "erld",
                                  central_region = 0.5,
                                  emp_factor = 1.5,
                                  erld_type = NULL,
                                  dq_quantiles = NULL)
t0_outliers$depth_values

# multivaraite으로 말고 각 variable마다 depth 계산하고, 얘들을 aggregate하는 방법을 생각해보자
# - depth의 scale이 비슷한지를 먼저 확인해야함
#   - depth가 작을수록 outlier로 볼 수 있을듯???
# - 각 calib data를 하나씩 추가해가면서 그 값의 depth 계산




plot(extreme_rank_length(dts))


par(mfrow = c(2, 2))
for (i in 1:4) {
  # depth_values <- modified_band_depth(df[, , i])
  # depth_values <- total_variation_depth(df[, , i])$tvd
  # depth_values <- extremal_depth(df[, , i])
  # depth_values <- -dir_out(df[, , i], return_distance = T)$distance
  # depth_values <- linfinity_depth(df[, , i])
  # depth_values <- band_depth(df[, , i])
  depth_values <- -directional_quantile(df[, , i])
  
  plot(depth_values, col = ifelse(1:n %in% idx_outliers, 2, 1))
}

# Compute univariate functional depth for calibration set
nonconform_score_calib <- matrix(NA, n_calib, p)
for (i in 1:p) {
  nonconform_score_calib[, i] <- apply(X_calib[[i]], 1, function(x){
    # directional_quantile(rbind(X_train[[i]], x))[n_train+1]
    -linfinity_depth(rbind(X_train[[i]], x))[n_train+1]
  })
}

# Compute univariate functional depth for test set
nonconform_score_test <- matrix(NA, n_test, p)
for (i in 1:p) {
  nonconform_score_test[, i] <- apply(X_test[[i]], 1, function(x){
    # directional_quantile(rbind(X_train[[i]], x))[n_train+1]
    -linfinity_depth(rbind(X_train[[i]], x))[n_train+1]
  })
}

agg_score_calib <- apply(nonconform_score_calib, 1, max)
agg_score_test <- apply(nonconform_score_test, 1, max)
conf_pvalue_marg <- sapply(agg_score_test, function(s){
  (1 + sum(agg_score_calib >= s)) / (n_calib + 1)
})
conf_pvalue_marg[which(idx_test %in% idx_outliers)]


depth_values[1, ] %>% plot
depth_values[2, ] %>% plot
depth_values[3, ] %>% plot

par(mfrow = c(2, 2))
depth_values[, 1] %>% plot(col = ifelse(idx_test %in% idx_outliers, 2, 1))
depth_values[, 2] %>% plot(col = ifelse(idx_test %in% idx_outliers, 2, 1))
depth_values[, 3] %>% plot(col = ifelse(idx_test %in% idx_outliers, 2, 1))

par(mfrow = c(2, 2))
apply(depth_values, 1, mean) %>% plot(col = ifelse(idx_test %in% idx_outliers, 2, 1))
apply(depth_values, 1, min) %>% plot(col = ifelse(idx_test %in% idx_outliers, 2, 1))





apply(X_test[[i]], 1, function(x){
  directional_quantile(rbind(X_train[[i]], x))[n_train+1]
})



# library(ddalpha)

library(mrfDepth)
depth_values <- mfd(aperm(df, c(2, 1, 3)), type = "hdepth")
depth_values <- mfd(aperm(df, c(2, 1, 3)), type = "projdepth")
depth_values <- mfd(aperm(df, c(2, 1, 3)), type = "dprojdepth")
plot(depth_values$MFDdepthX, col = ifelse(1:n %in% idx_outliers, 2, 1))


depth_values <- mfd(aperm(df[1:100, , ], c(2, 1, 3)), 
                    aperm(df[101:200, , ], c(2, 1, 3)),
                    type = "projdepth")$MFDdepthZ
depth_values[1:5, ]
mfd(aperm(df[1:100, , ], c(2, 1, 3)), 
    aperm(df[101:105, , ], c(2, 1, 3)),
    type = "projdepth")$MFDdepthZ




library(mclust)
fit_xy <- densityMclust(cbind(y, x), G = 2, plot = FALSE)  # f(x,y) 
fit_x <- densityMclust(x, G = 2, plot = FALSE)  # f(x)
f_hat <- fit_xy$density / ifelse(fit_x$density < 1e-8, 1e-8, fit_x$density)   # f(y|x) = f(x,y) / f(x)
range(f_hat)






load("RData/sim_2.RData")

lapply(res2, function(sim){
  sub <- paste0(
    rbind(fdr = colMeans(sim$fdr),
          tpr = colMeans(sim$tpr)) %>% 
      round(3) %>% 
      format(nsmall = 3),
    " (",
    rbind(fdr = apply(sim$fdr, 2, sd),
          tpr = apply(sim$tpr, 2, sd)) %>% 
      round(3) %>% 
      format(nsmall = 3),
    ")"
  )
  dim(sub) <- c(2, ncol(sim$fdr))
  rownames(sub) <- c("FDR","TPR")
  colnames(sub) <- colnames(sim$fdr)
  sub <- data.frame(sub)
  sub %>% xtable::xtable()
}) 






# Risk-controliling Prediction Set (RCPS)
rcps <- list()
for (j in 1:length(which_partition_test)) {
  # Calibration set belongs to A_i
  idx_calib_in_A_i <- which(which_partition_train == which_partition_test[j])
  y_calib_A_i <- y_calib[idx_calib_in_A_i]
  
  # CDE for i-th test data
  f_hat <- data.frame(x = pred_test$z,
                      y = pred_test$CDE[j, ])
  
  # # Remove vanishing part of the CDE
  # f_hat <- f_hat[which(f_hat$y > 0), ]
  
  # Nested set-valued predictors using Greedy algorithm
  T_lambda <- list()
  d <- 0.1   # step size rate
  zeta <- 3  # bound of random variables in Hoeffding's inequality
  lambda_list <- -seq(1e-8, max(f_hat$y)-(1e-8), length.out = 20)  # candidates of lambda
  for (i in 1:length(lambda_list)) {
    lambda <- lambda_list[i]
    T_lambda_i <- c()
    
    # iter <- 0
    zeta_update <- zeta
    while (zeta_update > -lambda) {
      # iter <- iter + 1
      zeta_update <- zeta_update - d*zeta_update
      T_lambda_i <- c(T_lambda_i, 
                      f_hat$x[(f_hat$y > zeta_update) & !(f_hat$x %in% T_lambda_i)])
    }
    T_lambda[[i]] <- sort(T_lambda_i)
    # print(iter)
  }
  # sapply(T_lambda, length)
  
  
  # Split the dis-connected region
  T_lambda_split <- lapply(T_lambda, function(tol_band){
    # Find connected components of bands
    dif <- round(diff(tol_band), 5)
    split_idx <- which(dif > min(dif))
    if (length(split_idx) > 0) {
      tol_band_split <- list()
      for (i in 1:length(split_idx)) {
        if (i == 1) {
          tol_band_split[[1]] <- tol_band[1:split_idx[1]]
        } else {
          tol_band_split[[i]] <- tol_band[(split_idx[i-1]+1):split_idx[i]]
        }
      }
      tol_band_split[[length(split_idx)+1]] <- tol_band[(split_idx[length(split_idx)]+1):length(tol_band)]
    } else {
      tol_band_split <- list(tol_band)
    }
    
    # Obtain the range of each band
    tol_band_split <- lapply(tol_band_split, range)
    
    return(tol_band_split)
  })
  
  # # Visualization
  # par(mfrow = c(1, 2))
  # plot(f_hat, type = "l")
  # for (i in 1:length(lambda_list)) {
  #   tol_band_split <- T_lambda_split[[i]]
  #   for (j in 1:length(tol_band_split)) {
  #     lines(tol_band_split[[j]], rep(-lambda_list[[i]], 2), col = i, lwd = 3)
  #   }
  # }
  # plot(f_hat[which(f_hat$y > 0), ], type = "l")
  
  
  ## Find UCB using the concentration inequality
  # Empirical risk
  R_hat_lambda <- sapply(1:length(T_lambda_split), function(i){
    number <- sapply(T_lambda_split[[i]], function(interval){
      sum(y_calib_A_i >= interval[1] & y_calib_A_i <= interval[2])
    })
    1 - sum(number)/length(y_calib_A_i)
  })
  # Bound of concentration inequality
  bound <- sqrt(log(1/delta)/(2*length(y_calib_A_i)))  # Hoeffding's inequality
  
  # Find lambda_hat (We use max since lambda has decreasing order!)
  ucb_rule <- which(R_hat_lambda < alpha - bound)
  if (length(ucb_rule) == 0) {
    lambda_hat_order <- NA
  } else {
    lambda_hat_order <- max(ucb_rule)
  }
  # lambda_hat_order
  
  # RCPS with lambda_hat
  tol_band <- T_lambda_split[[lambda_hat_order]]
  rcps[[j]] <- tol_band
}


par(mfrow = c(1, 2))
plot(x_test, y_test, col = which_partition_test)
idx <- which(sapply(rcps, is.null))
plot(x_test[idx], y_test[idx])










