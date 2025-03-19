library(tidyverse)

#######################################
### Triceps data
#######################################
triceps <- MultiKink::triceps
dim(triceps)
head(triceps)

# with(triceps, plot(age, lntriceps, ylim = c(0, 5)))
# lines(smooth.spline(triceps$age, triceps$lntriceps, cv = T), col = 2, lwd = 3)

n <- nrow(triceps)

# Training + validation set
x <- matrix(triceps$age, ncol = 1)
y <- triceps$lntriceps

#######################################
### Gallaxy data
#######################################
data <- read.table("data/photoz_catalogues/Happy/happy_A", skip = 1, header = F)
dim(data)
colnames(data) <- c("id", "mag_r", "u-g", "g-r", "r-i", "i-z", "z_spec",
                    "feat1", "feat2", "feat3", "feat4", "feat5")
head(data)
y <- data$z_spec   # spectroscopic redshift
x <- data[, c("mag_r", "feat1", "feat2", "feat3", "feat4", "feat5")] %>%
  as.matrix()

# Split data
set.seed(1000)
idx_test <- sample(1:nrow(x), 5000)
idx_calib <- sample(setdiff(1:nrow(x), idx_test), 5000)
idx_train <- setdiff(1:nrow(x), c(idx_test, idx_calib))

x_train <- x[idx_train, ]
y_train <- y[idx_train]

x_calib <- x[idx_calib, ]
y_calib <- y[idx_calib]

x_test <- x[idx_test, ]
y_test <- y[idx_test]

x <- rbind(x_train, x_calib)
y <- c(y_train, y_calib)

n <- nrow(x)

idx <- sample(1:nrow(x), 5000)
x <- rbind(x_train, x_calib)[idx, ]
y <- c(y_train, y_calib)[idx]
n <- nrow(x)


####################################################
### Tolerance interval for nonparametric regression
####################################################

# Tolerance level
alpha <- 0.1    # 1 - proportion of the sampled population
delta <- 0.05   # 1 - confidence level (coverage level)

### Algorithm 1 (Bootstrap Algorithm for Simultaneous k-Factor)
# (1)
if (ncol(x) == 1) {
  # Smoothing spline for univariate covariate
  fit <- smooth.spline(x, y, cv = T, keep.stuff = T)
  lambda <- fit$lambda
  pred <- fitted(fit)
  eps <- y - pred
  
  # Smoother matrix of spline (hat matrix)
  L_lambda <- apply(diag(n), 2, function(y_i) {
    fitted(smooth.spline(x, y_i, lambda = lambda))
  })
  # L_lambda <- sfsmisc::hatMat(x, lambda = fit$lambda)
  # sum((L_lambda %*% y - fitted(fit))^2)
  # par(mfrow = c(1, 2))
  # plot(eps)
  # plot( (diag(n) - L_lambda) %*% y )
  # sum((diag(L_lambda) - hatvalues(fit))^2)
} else {
  # GAM for smoothing spline for multivariate covariate
  library(mgcv)
  df <- data.frame(y = y,
                   x = x)
  formula_fit <- as.formula( paste0("y ~ ", paste0("s(", colnames(df)[-1], ", bs='cr')", collapse = " + ")) )
  fit <- gam(formula_fit, data = df)
  pred <- fitted(fit)
  eps <- y - pred
  
  # Smoother matrix of GAM (hat matrix)
  # # sparse matrix로 바꾸기
  # z <- Matrix::sparseMatrix(3, 3)
  # z[1, ]
  system.time({
    L_lambda <- apply(diag(n), 2, function(y_i) {
      fitted(gam(formula_fit, data = data.frame(y = y_i,
                                                x = x),
                 sp = fit$sp))
    })
  })
  # n = 1000
  #   user  system elapsed 
  # 25.221   4.243  29.602 
  # n = 5000
  #    user  system elapsed 
  # 296.039  82.355 382.064 
  # n = 10000
  #     user   system  elapsed 
  # 1032.973  364.942 1400.869 
  # sum((L_lambda %*% y - fitted(fit))^2)
  # sum((diag(L_lambda) - fit$hat)^2)
  
  
  
  # # Find bandwidth of NW estimator (Do not use product kernel)
  # system.time({
  #   fit <- np::npregbw(xdat = x, ydat = y, regtype = "lc", ckertype = "epanechnikov")
  # })
  # # n = 1000
  # #    user  system elapsed 
  # # 112.719   0.680 114.243 
  # # pred <- np::npreg(bws = fit, exdat = x)
  # # plot(x[, 1], pred$mean)
  # # points(x[, 1], y, col = 2)
  # 
  # # Epanichinikov kernel
  # # x: matrix, x_eval: vector, bw: vector
  # epan_prod <- function(x, x_eval, bw) {
  #   kern_value <- apply(x, 1, function(x_i){
  #     u <- (x_i - x_eval) / bw
  #     u[abs(u) > 1] <- 1  # make kernel estimate = 0
  #     prod(3/4 * (1 - u^2))
  #   })
  #   return(kern_value)
  # }
  # 
  # # Nadaraya-Watson estimator using product kernel
  # # x: matrix, x_eval: matrix, bw: vector
  # nw_est <- function(x, x_eval, bw) {
  #   # Kernel weight using product kernel
  #   k <- apply(x, 1, function(x_i){ epan_prod(x = x, x_eval = x_i, bw = fit$bw) })   # n_eval x n
  #   w <- apply(k, 1, function(k_i){ k_i / sum(k_i) })   # n x n_eval
  #   
  #   # NW estimator
  #   est <- t(w) %*% y
  #   
  #   return(list(y_hat = est,
  #               hat_matrix = t(w)))
  # }

}
# L2 norm of each row of hat matrix
l_x_norm <- sqrt(rowSums(L_lambda^2))


# formula_fit <- as.formula( paste0("y ~ ", paste0(colnames(df)[-1], collapse = " + ")) )
# fit <- loess(formula_fit, data = df)
# z <- x^2+rnorm(n)
# df <- data.frame(y = y,
#                  x = x,
#                  z = z)
# fit <- loess(y ~ x + z, data = df)
# L_lambda <- apply(diag(n), 2, function(y_i) {
#   # fitted(loess(y ~ x, data = data.frame(y = y_i,
#   #                                       x = x)))
#   fitted(loess(y ~ x + z, data = data.frame(y = y_i,
#                                             x = x,
#                                             z = z)))
# })
# sum(abs(L_lambda %*% y - fitted(fit)))
# 
# par(mfrow = c(2, 3))
# for (i in 1:6) {
#   plot(x[, i], y)
# }


# sigma^2 estimate
sub <- t(diag(n) - L_lambda) %*% (diag(n) - L_lambda)
tr_sub <- sum(diag(sub))
sigma2 <- sum(eps^2) / tr_sub

# Pointwise k factor
if ("gam" %in% class(fit)) {
  # GAM
  nu <- fit$df.residual
} else {
  nu <- tr_sub^2 / sum(diag(sub %*% sub))
  # nu <- n - fit$df
  # nu <- tr_sub
}
k_pt <- sqrt( nu * qchisq(1-alpha, df = 1, ncp = l_x_norm^2) / qchisq(delta, df = nu) )

# (2)
set.seed(100)
B <- 1000
k_cand <- seq(1, 5, length.out = 30)
min_content <- matrix(NA, B, length(k_cand))
exp_f_hat <- L_lambda %*% pred  # E(\hat{f}^{(b)}(x)) using \hat{f} (assume true f is \hat{f})
for (b in 1:B) {
  # Bootstrap sample
  idx_boot <- sample(1:n, n, replace = T)
  y_boot <- pred + eps[idx_boot]
  if (ncol(x) == 1) {
    # Smoothing spline for univariate covariate 
    fit_boot <- smooth.spline(x, y_boot, lambda = lambda)
  } else {
    # GAM for smoothing spline for multivariate covariate
    fit_boot <- gam(formula_fit, 
                    data = data.frame(y = y_boot,
                                      x = x),
                    sp = fit$sp)
  }
  pred_boot <- fitted(fit_boot)
  sigma2_boot <- sum((y_boot - pred_boot)^2) / tr_sub
  
  
  # (3)
  # Minimum content using CLT
  min_content[b, ] <- sapply(k_cand, function(k){
    # z_value_lb <- (pred_boot - k*sqrt(sigma2_boot) - pred) / (sqrt(sigma2) * l_x_norm)
    # z_value_ub <- (pred_boot + k*sqrt(sigma2_boot) - pred) / (sqrt(sigma2) * l_x_norm)
    # z_value_lb <- (pred_boot - k*sqrt(sigma2_boot) - pred) / sqrt(sigma2)
    # z_value_ub <- (pred_boot + k*sqrt(sigma2_boot) - pred) / sqrt(sigma2)
    z_value_lb <- (pred_boot - exp_f_hat - k*sqrt(sigma2_boot)) / sqrt(sigma2)
    z_value_ub <- (pred_boot - exp_f_hat + k*sqrt(sigma2_boot)) / sqrt(sigma2)
    # z_value_lb <- (pred_boot - exp_f_hat - k*sigma2_boot) / sigma2
    # z_value_ub <- (pred_boot - exp_f_hat + k*sigma2_boot) / sigma2
    min( pnorm(z_value_ub) - pnorm(z_value_lb) )
  })
}
# (4)
gamma <- apply(min_content, 2, function(col){
  mean(col >= 1 - alpha)
})
# (5)
k_sim <- k_cand[ which.min(abs(gamma - (1 - delta))) ]

# Width of tolerance band
2*k_sim*sqrt(sigma2)
mean(2*k_pt*sqrt(sigma2))

# Simultaneous tolerance band
tol_band <- data.frame(
  lb = pred - k_sim*sqrt(sigma2),
  ub = pred + k_sim*sqrt(sigma2)
)
print(paste("Coverage (simultaneous):", 
            round(mean(y >= tol_band$lb & y <= tol_band$ub), 3)))
print(paste("Size (simultaneous):", 
            round(mean(tol_band$ub - tol_band$lb), 3)))

# Pointwise tolerance band
tol_band_pt <- data.frame(
  lb = pred - k_pt*sqrt(sigma2),
  ub = pred + k_pt*sqrt(sigma2)
)
print(paste("Coverage (pointwise):", 
            round(mean(y >= tol_band_pt$lb & y <= tol_band_pt$ub), 3)))
print(paste("Size (pointwise):", 
            round(mean(tol_band_pt$ub - tol_band_pt$lb), 3)))


# Coverage for test set - ??? 이거 맞나???
pred_test <- predict(fit, newdata = data.frame(x = x_test,
                                               y = y_test))
print(paste("Coverage (simultaneous) - test:", 
            round(mean(y_test >= pred_test - k_sim*sqrt(sigma2) 
                       & y_test <= pred_test + k_sim*sqrt(sigma2)), 3)))
# print(paste("Coverage (pointwise) - test:", 
#             round(mean(y_test >= tol_band_pt$lb & y_test <= tol_band_pt$ub), 3)))



df <- data.frame(
  x = triceps$age,
  y = triceps$lntriceps,
  pred = pred,
  lb = tol_band$lb,
  ub = tol_band$ub,
  lb_pt = tol_band_pt$lb,
  ub_pt = tol_band_pt$ub
)

ggplot() +
  geom_ribbon(data = df, aes(x = x, ymin = lb, ymax = ub), 
              color = "gray", fill = "gray", alpha = 0.6, linewidth = 0.7) +
  geom_ribbon(data = df, aes(x = x, ymin = lb_pt, ymax = ub_pt),
              color = "darkgray", fill = "darkgray", alpha = 1, linewidth = 0.7) +
  geom_point(data = df, aes(x = x, y = y)) +
  geom_line(data = df, aes(x = x, y = pred), color = "red", linewidth = 1) +
  scale_y_continuous(limits = c(0, 4.5)) +
  labs(x = "Age", y = "Triceps", title = "Nonparametric Tolerance region") +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        legend.position = "bottom")
# ggsave("./Rmd/2024_1022/nonpara_tol.pdf", width = 7, height = 7)






# source("R/npreg_tol.R")
# 
# fit_npreg_tol <- npreg_tol(x, y, B = 1000, alpha = 0.1, delta = 0.05)
# 
# df <- data.frame(
#   x = triceps$age,
#   y = triceps$lntriceps,
#   pred = predict(fit_npreg_tol$fit, x)$y,
#   lb = fit_npreg_tol$tol_band_sim$lb,
#   ub = fit_npreg_tol$tol_band_sim$ub,
#   lb_pt = fit_npreg_tol$tol_band_pt$lb,
#   ub_pt = fit_npreg_tol$tol_band_pt$ub
# )
# 
# ggplot() +
#   geom_ribbon(data = df, aes(x = x, ymin = lb, ymax = ub), 
#               color = "gray", fill = "gray", alpha = 0.6, linewidth = 0.7) +
#   geom_ribbon(data = df, aes(x = x, ymin = lb_pt, ymax = ub_pt),
#               color = "darkgray", fill = "darkgray", alpha = 1, linewidth = 0.7) +
#   geom_point(data = df, aes(x = x, y = y)) +
#   geom_line(data = df, aes(x = x, y = pred), color = "red", linewidth = 1) +
#   scale_y_continuous(limits = c(0, 4.5)) +
#   labs(x = "Age", y = "Triceps", title = "Nonparametric Tolerance region") +
#   theme_light() +
#   theme(plot.title = element_text(hjust = 0.5),
#         legend.title = element_blank(),
#         legend.position = "bottom")



source("R/npreg_tol.R")
p_list <- list()

fit_npreg_tol <- npreg_tol(x, y, B = 1000, alpha = 0.1, delta = 0.05)
df <- data.frame(
  x = triceps$age,
  y = triceps$lntriceps,
  pred = fitted(fit_npreg_tol$fit)*fit_npreg_tol$scale_est,
  lb = fit_npreg_tol$tol_band_sim$lb,
  ub = fit_npreg_tol$tol_band_sim$ub,
  lb_pt = fit_npreg_tol$tol_band_pt$lb,
  ub_pt = fit_npreg_tol$tol_band_pt$ub
)
p_list[[1]] <- ggplot() +
  geom_ribbon(data = df, aes(x = x, ymin = lb, ymax = ub), 
              color = "gray", fill = "gray", alpha = 0.6, linewidth = 0.7) +
  geom_ribbon(data = df, aes(x = x, ymin = lb_pt, ymax = ub_pt),
              color = "darkgray", fill = "darkgray", alpha = 1, linewidth = 0.7) +
  geom_point(data = df, aes(x = x, y = y)) +
  geom_line(data = df, aes(x = x, y = pred), color = "red", linewidth = 1) +
  scale_y_continuous(limits = c(0.5, 4)) +
  labs(x = "Age", y = "Triceps", title = "") +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        legend.position = "bottom")

fit_npreg_tol <- npreg_tol(x, y, B = 1000, alpha = 0.1, delta = 0.05, scale = T)
df <- data.frame(
  x = triceps$age,
  y = triceps$lntriceps,
  pred = fitted(fit_npreg_tol$fit)*fit_npreg_tol$scale_est,
  lb = fit_npreg_tol$tol_band_sim$lb,
  ub = fit_npreg_tol$tol_band_sim$ub,
  lb_pt = fit_npreg_tol$tol_band_pt$lb,
  ub_pt = fit_npreg_tol$tol_band_pt$ub
)
p_list[[2]] <- ggplot() +
  geom_ribbon(data = df, aes(x = x, ymin = lb, ymax = ub), 
              color = "gray", fill = "gray", alpha = 0.6, linewidth = 0.7) +
  geom_ribbon(data = df, aes(x = x, ymin = lb_pt, ymax = ub_pt),
              color = "darkgray", fill = "darkgray", alpha = 1, linewidth = 0.7) +
  geom_point(data = df, aes(x = x, y = y)) +
  geom_line(data = df, aes(x = x, y = pred), color = "red", linewidth = 1) +
  scale_y_continuous(limits = c(0.5, 4)) +
  labs(x = "Age", y = "Triceps", title = "") +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        legend.position = "bottom")
p <- gridExtra::grid.arrange(grobs = p_list, nrow = 1)

# ggsave("./Rmd/2024_1029/npreg_tol.pdf", p, width = 10, height = 6)



sign(df$lb - df$lb_pt) %>% table()
sign(df$ub - df$ub_pt) %>% table()
