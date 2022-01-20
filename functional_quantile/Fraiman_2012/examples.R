# Example of empirical quantile contours

library(MASS)
source("Fraiman_2012/code.R")

# Spherical distribution in R2
# ================================
d <- 2               # dimension
n <- 2000            # sample size
mu <- c(0, 0)         # Mean vector
Sigma <- diag(1, 2)
set.seed(1000)
x <- mvrnorm(
  n = n,
  mu = mu,
  Sigma = Sigma,
  tol = 1e-6,
  empirical = FALSE
)

par(mfrow = c(1, 2))
matplot(t(x), type = "l")
dire <- 100 # Number of directions to be considered
theta <- seq(0, 2 * pi, le = dire)
u <- cbind(cos(theta), sin(theta)) # Matrix with unit vectors (bidimensional case)
out <- mvquantile(x, alpha = 0.9, u)
plot(x)
points(out[, 1, 1], out[, 2, 1], col = 2)



pqc.fit <- qpc(x,
               alpha = 0.5,
               np = 2,
               ct = "median",
               depth = "FM")


# Elliptical distribution in R2
# ================================
d <- 2               # dimension
n <- 2000            # sample size
mu <- c(0, 0)         # Mean vector
Sigma <- matrix(c(1, 0.75, 0.75, 1),
                nr = 2,
                nc = 2,
                byr = T)
set.seed(1000)
x <- mvrnorm(
  n = n,
  mu = mu,
  Sigma = Sigma,
  tol = 1e-6,
  empirical = FALSE
)

par(mfrow = c(1, 2))
matplot(t(x), type = "l")
dire <- 100 # Number of directions to be considered
theta <- seq(0, 2 * pi, le = dire)
u <- cbind(cos(theta), sin(theta)) # Matrix with unit vectors (bidimensional case)
out <- mvquantile(x, alpha = 0.9, u)
plot(x)
points(out[, 1, 1], out[, 2, 1], col = 2)
