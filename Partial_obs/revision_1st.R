### Boxplot of eigenfunction angle
par(mfrow = c(2, 2))
load("RData/Delaigle-normal-prop0.RData")
boxplot(mse_eigen2[, c(1,2,3,4,7)],
        ylab = "Eigenfunction angle",
        main = "Gaussian")

load("RData/Delaigle-tdist.RData")
boxplot(mse_eigen2[, c(1,2,3,4,7)],
        ylab = "Eigenfunction angle",
        main = "t(3)")

load("RData/Delaigle-normal-prop1.RData")
boxplot(mse_eigen2[, c(1,2,3,4,7)],
        ylab = "Eigenfunction angle",
        main = "10% contaminated")

load("RData/Delaigle-normal-prop2.RData")
boxplot(mse_eigen2[, c(1,2,3,4,7)],
        ylab = "Eigenfunction angle",
        main = "20% contaminated")


subspace(c(1,1), c(-1,-1))
work.grid <- seq(0, 1, length.out = 51)
eig.true <- get_delaigle_eigen(work.grid, model = 2) 
eig.est <- pca.est[[1]]$pca.obj$pca.ogk.sm.MM.obj$eig.fun
subspace(eig.true, eig.est)
subspace(eig.true, check_eigen_sign(eig.est, eig.true))

i <- 1
subspace(eig.true[, i], eig.est[, i])


### 각 eigenfunction마다 계산하면 작게 나오는데, 전체 subspace로 하면 나쁘게 나오네...
load("RData/Delaigle-normal-prop0.RData")
k <- 4
sapply(pca.est[[1]]$pca.obj, function(obj) {
  subspace(eig.true[, k], obj$eig.fun[, k])
})



# black : True
# red : yao
# green : kraus
# blue : proposed
par(mfrow = c(2, 2))
for (i in 1:4) {
  matplot(cbind(eig.true[, i],
                check_eigen_sign(pca.est[[1]]$pca.obj$pca.yao.obj$eig.fun[, i], eig.true[, i]),
                check_eigen_sign(pca.est[[1]]$pca.obj$pca.kraus.obj$eig.fun[, i], eig.true[, i]),
                check_eigen_sign(pca.est[[1]]$pca.obj$pca.ogk.sm.MM.obj$eig.fun[, i], eig.true[, i])),
          type = "l", 
          ylab = "", main = paste("Eigen", i))
}



