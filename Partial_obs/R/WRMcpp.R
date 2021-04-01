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


# ### benchmark between "robfilter" and C++
# library(rbenchmark)
# data(faithful) # Old Faithful Geyser data
# fit1 <- wrm.smooth(faithful$w, faithful$e, h=4)
# fit2 <- wrm_smooth(faithful$w, faithful$e, h=4)
# all.equal(fit1$level, fit2$mu)
# all.equal(fit1$slope, fit2$beta)
# 
# benchmark(robfilter = wrm.smooth(faithful$w, faithful$e, h=4), 
#           Cpp = wrm_smooth(faithful$w, faithful$e, h=4),
#           replications = 100,
#           columns = c("test", "replications", "elapsed", "relative"))
# 
# gr <- seq(min(faithful$w), max(faithful$w), length.out = 100)
# wrm_smooth_cpp(faithful$w, faithful$e, 4, gr, "epanechnikov")
# 
# x <- faithful$w
# y <- faithful$e
# h <- 4
# xgrid <- seq(min(faithful$w), max(faithful$w), length.out = 100)
# kernel <- "epanechnikov"
# 
# i <- 1
# weights = get_kernel(x, xgrid[i], h, kernel)
# fitted = wrm_fit(x, y, xgrid[i], weights)
# 
# 
# ### Comparison with "robfilter" and C++
# load("RData/sim_20210323.RData")
# 
# err <- numeric(100)
# par(mfrow = c(4, 4))
# for (i in 1:100) {
#   print(i)
#   
#   x <- sim.obj$Out_3$data.list[[i]]$x
#   gr <- sim.obj$Out_3$data.list[[i]]$gr
#   
#   t <- unlist(x$Lt)
#   y <- unlist(x$Ly)
#   
#   # R function in "robfilter"
#   start_time <- Sys.time()
#   fit.R <- wrm.smooth(t, y, h = 0.1, weight = 2, xgrid = gr)
#   # fit.R <- wrm.smooth(t, y, h = 0.1, weight = 3, xgrid = gr)
#   end_time <- Sys.time()
#   print(paste0("R: ", round(end_time - start_time, 3)))
#   
#   # C++ function
#   start_time <- Sys.time()
#   fit.Cpp <- wrm_smooth(t, y, h = 0.1, kernel = "epanechnikov", xgrid = gr)
#   # fit.Cpp <- wrm_smooth(t, y, h = 0.1, kernel = "gauss", xgrid = gr)  
#   end_time <- Sys.time()
#   print(paste0("C++: ", round(end_time - start_time, 3)))
#   
#   err[i] <- sum((fit.R$level - fit.Cpp$mu)^2)
#   
#   matplot(gr, cbind(fit.R$level,
#                     fit.Cpp$mu),
#           type = "l",
#           xlab = "", ylab = "",
#           main = i)
# }
# 
# benchmark(robfilter = wrm.smooth(t, y, h = 0.1, weight = 2, xgrid = gr), 
#           Cpp = wrm_smooth(t, y, h = 0.1, kernel = "epanechnikov", xgrid = gr)  ,
#           replications = 100,
#           columns = c("test", "replications", "elapsed", "relative"))



