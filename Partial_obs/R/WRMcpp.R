### Weighted repeated median smoothing
# It is the same algorithm of "wrm.smooth()" in robfilter pacakges.
# For fast computation, transform the R code to C++.
# The parameters of a function are the same with "wrm.smooth()" except on "kernel".
# "kernel" option is provided "gauss" and "epanechnikov" only.
wrm_smooth <- function(x, y, h, xgrid = NULL, kernel = "epanechnikov") {
  
  if (is.null(xgrid)) {
    xgrid <- sort(x)
    if (length(xgrid) > 100) {
      xgrid <- seq(min(x), max(x), length.out = 100)
    }
  }
  
  res <- wrm_smooth_cpp(x, y, h, xgrid, kernel)
  
  return(res)
}

library(rbenchmark)
data(faithful) # Old Faithful Geyser data
fit1 <- wrm.smooth(faithful$w, faithful$e, h=4)
fit2 <- wrm_smooth(faithful$w, faithful$e, h=4)
all.equal(fit1$level, fit2$mu)
all.equal(fit1$slope, fit2$beta)

benchmark(robfilter = wrm.smooth(faithful$w, faithful$e, h=4), 
          Cpp = wrm_smooth(faithful$w, faithful$e, h=4),
          replications = 100,
          columns = c("test", "replications", "elapsed", "relative"))

gr <- seq(min(faithful$w), max(faithful$w), length.out = 100)
wrm_smooth_cpp(faithful$w, faithful$e, 4, gr, "epanechnikov")

x <- faithful$w
y <- faithful$e
h <- 4
xgrid <- seq(min(faithful$w), max(faithful$w), length.out = 100)
kernel <- "epanechnikov"

i <- 1
weights = get_kernel(x, xgrid[i], h, kernel)
fitted = wrm_fit(x, y, xgrid[i], weights)
