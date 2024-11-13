
conformal <- function(type = "split") {
  
  if (type == "split") {
    
  } else if (type == "full") {
    conformal.split()
  }
}


conformal.split <- function(x, y, alpha = 0.1, rho = 0.5, s.type = "std", train.idx = NULL) {
  
  n <- length(y)   # number of data
  m <- round(n * rho)   # size of training set
  l <- n - m   # size of calibration set
  
  # Split data
  idx_train <- sample(1:n, m)
  idx_calib <- setdiff(1:n, idx_train)
  
  data <- data.frame(y = y,
                     x = x)
  data_train <- data[idx_train, ]
  data_calib <- data[idx_calib, ]
  
  # Point predictor
  fit <- lm(y ~ x, data = data_train)
  pred <- predict(fit, data_calib)
  
  # Scale estimator function
  if (s.type == "std") {
    # Scaling using standard deviation given X = x
    df <- data.frame(resid = abs(fit$residuals),
                     x = data_train$x)
    fit_local <- loess(resid ~ x, df)
    s_ftn <- predict(fit_local, data_calib$x)
  } else if (s.type == "identity") {
    # Not scaled
    s_ftn <- 1
  }
  
  # Non-conformity score
  nonconform_score <- abs(data_calib$y - pred) / s_ftn
  idx_cutoff <- ceiling((1 - alpha) * (l + 1))
  k_s <- sort(nonconform_score)[idx_cutoff]
  
  # Empirical coverage
  sum(nonconform_score <= k_s) / l
  
  # Conformal prediction band
  x_new <- data.frame(x = seq(min(data_train$x), max(data_train$x), length.out = 100))
  pred_x_new <- predict(fit, x_new)
  lb <- pred_x_new - k_s*s_ftn
  ub <- pred_x_new + k_s*s_ftn
  
  # Output
  out <- list(
    pred = pred,
    lb = lb,
    ub = ub,
    split_obj = list(
      idx_train = idx_train,
      idx_calib = idx_calib
    )
  )
  
  return(out)
}



### Split conformal prediction
n <- 1000

set.seed(1000)
x <- rnorm(n)
# y <- x + rnorm(n)
y <- x + rnorm(n, sd = exp(0.5*x))
data <- data.frame(y = y,
                   x = x)

alpha <- 0.1  # coverage level
rho <- 0.5  # split proportion

m <- round(n * rho)   # size of training set
l <- n - m   # size of calibration set

# Split data
idx_train <- sample(1:n, m)
idx_calib <- setdiff(1:n, idx_train)

data_train <- data[idx_train, ]
data_calib <- data[idx_calib, ]

# Point predictor
fit <- lm(y ~ x, data = data_train)
pred <- predict(fit, data_calib)

# Non-conformity score using calibration set
nonconform_score <- abs(data_calib$y - pred)
idx_cutoff <- ceiling((1 - alpha) * (l + 1))
k <- sort(nonconform_score)[idx_cutoff]

# Coverage
sum(nonconform_score <= k) / l

# Visualize conformal prediction band
par(mfrow = c(1, 2))
plot(data_calib$x, data_calib$y)
x_new <- data.frame(x = seq(min(data_train$x), max(data_train$x), length.out = 100))
pred_x_new <- predict(fit, x_new)
polygon(c(x_new$x, rev(x_new$x)),
        c(pred_x_new + k, rev(pred_x_new - k)),
        border = F, col = "gray")
points(data_calib$x, data_calib$y)
abline(fit, col = 2, lwd = 2)


### Locally adaptive prdiction band
# Scale estimator function
df <- data.frame(resid = abs(fit$residuals),
                 x = data_train$x)
fit_local <- loess(resid ~ x, df)
sigma <- predict(fit_local, data_calib$x)

# Non-conformity score
nonconform_score <- abs(data_calib$y - pred) / sigma
idx_cutoff <- ceiling((1 - alpha) * (l + 1))
k_s <- sort(nonconform_score)[idx_cutoff]

# Coverage
sum(nonconform_score <= k_s) / l

# Visualize conformal prediction band
plot(data_calib$x, data_calib$y)
pred_sigma <- predict(fit_local, x_new)
polygon(c(x_new$x, rev(x_new$x)),
        c(pred_x_new + k_s*pred_sigma, rev(pred_x_new - k_s*pred_sigma)),
        border = F, col = "lightblue")
points(data_calib$x, data_calib$y)
abline(fit, col = 2, lwd = 2)
