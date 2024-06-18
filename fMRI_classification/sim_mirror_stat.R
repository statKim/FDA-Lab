# Simulation: mirror statistic centering해도 될까??
# X_0i ~ sqrt(1-rho) * y_i + sqrt(rho) * z,  y_i ~ iid N(0,1), z ~ N(0,1)
# y_i and z are independent
# X_1i ~ N(3, 1)

library(tidyverse)
source("R/homo_test.R")

n <- 10000
p <- 0.9
n1 <- n*p
n2 <- n*(1-p)

rho <- 0.9
z <- rnorm(1)
y <- rnorm(n1)
x0 <- sqrt(1-rho)*y + sqrt(rho)*z
x1 <- rnorm(n2, 3, 1)
x <- c(x0, x1)
z

par(mfrow = c(2, 2))
hist(x0, probability = T, breaks = 30,
     xlab = "", main = latex2exp::TeX("$X_0$"))
lines(density(x0), col = 2, lwd = 2)
hist(x1, probability = T, breaks = 30,
     xlab = "", main = latex2exp::TeX("$X_1$"))
lines(density(x1), col = 2, lwd = 2)
hist(x, probability = T, breaks = 30,
     xlab = "", main = latex2exp::TeX("$X$"))
lines(density(x), col = 2, lwd = 2)
# lines(density(x0), col = 4, lwd = 2)
# lines(density(x1), col = 5, lwd = 2)
# lines(seq(-4, 4, by = 0.1), dnorm(seq(-4, 4, by = 0.1)), col = 3, lwd = 2)
# legend("topright",
#        c(latex2exp::TeX("$X$"), latex2exp::TeX("$X_0$"), latex2exp::TeX("$X_1$"), "N(0,1)"),
#        lty = c(1, 1, 1, 1), lwd = c(2, 2, 2, 2), col = c(2, 4, 5, 3))


# obj <- mirror_stat_cutoff(x, q = 0.1, est = TRUE)
# fdp <- obj$fdp
# cutoff <- obj$cutoff
# 
# fdp <- sum(x0 > cutoff) / sum(x > cutoff)
# tpr <- sum(x1 > cutoff) / n2
# fdp
# tpr

sort(x, decreasing = T) %>% head()
sort(x0, decreasing = T) %>% head()
sort(x1, decreasing = T) %>% head()

# Centering between mode
dens <- density(x)
center <- dens$x[which.max(dens$y)]

hist(x - center, probability = T, breaks = 30,
     xlab = "", main = latex2exp::TeX("$X_{centered}$"))
lines(density(x - center), col = 2, lwd = 2)

cutoff <- mirror_stat_cutoff(x - center, q = 0.1)
cutoff

sum(x0 - center > cutoff) / sum(x - center > cutoff)
sum(x1 - center > cutoff) / n2



# Simulate 1000 times
B <- 1000
fdp <- rep(NA, B)
tpr <- rep(NA, B)
for (b in 1:B) {
  print(b)
  set.seed(b)
  
  n <- 10000
  p <- 0.9
  n1 <- n*p
  n2 <- n*(1-p)
  
  rho <- 0.9
  z <- rnorm(1)
  y <- rnorm(n1)
  x0 <- sqrt(1-rho)*y + sqrt(rho)*z
  x1 <- rnorm(n2, 3, 1)
  x <- c(x0, x1)
  
  # Centering between mode
  dens <- density(x)
  center <- dens$x[which.max(dens$y)]
  
  cutoff <- mirror_stat_cutoff(x - center, q = 0.1)
  # cutoff
  
  fdp[b] <- sum(x0 - center > cutoff) / sum(x - center > cutoff)
  tpr[b] <- sum(x1 - center > cutoff) / n2
}
mean(fdp)
mean(tpr)





