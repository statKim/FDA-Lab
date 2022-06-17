
#####################################
### Find optimal delta for Huber loss
### - k <- 1 / max(abs(y))
#####################################
# load("RData/sim3-2_20210204.RData")
load("RData/sim_20210317.RData")
data.list <- sim.obj[["Out_2"]]$data.list
cov.est <- sim.obj[["Out_2"]]$cov.est

# generated된 데이터의 range가 가장 큰 3개의 데이터 select
y_range <- sapply(data.list, function(x){
  return( diff(range(x$x$y)) )
})
ind_top3 <- sort(y_range, index.return = T, decreasing = T)$ix[1:3]


# sim <- 20
sim <- ind_top3[1]
model.cov <- 2   # covariance function setting of the paper (1, 2)
kernel <- "gauss"
bw <- 0.2

fig_delta <- list()
fig_bw <- list()

for (j in 1:length(ind_top3)) {
  sim <- ind_top3[j]
  
  # Get simulation data
  x <- data.list[[sim]]$x
  gr <- data.list[[sim]]$gr
  x.2 <- list(Ly = x$y,
              Lt = x$t)
  
  ### Covariance estimation
  work.grid <- cov.est[[sim]]$work.grid
  cov.true <- cov.est[[sim]]$cov$true
  cov.yao <- cov.est[[sim]]$cov$yao
  cov.lin <- cov.est[[sim]]$cov$lin
  cov.huber <- cov.est[[sim]]$cov$huber
  cov.wrm <- cov.est[[sim]]$cov$wrm
  
  
  # fixed tuning parameters
  bw_fixed_mu <- cov.est[[sim]]$cov.obj$huber$mu$bw
  bw_fixed_var <- cov.est[[sim]]$cov.obj$huber$sig2x$obj$bw
  
  kernel <- "gauss"
  # k_cand <- c(1 / max(abs(unlist(x.2$Ly))), 1e-4, 0.001, 0.01, 0.1, 0.5, 0.7, 1, 1.345)
  k_cand <- union(1.345,
                  10^seq(-3, 2, length.out = 10) * (1 - 0)/3)
  var_est <- matrix(0, length(work.grid), length(k_cand)+1)
  var_est[, 1] <- diag(cov.true)
  system.time({
    for (i in 1:length(k_cand)) {
      print(k_cand[i])
      
      mu.huber.obj <- meanfunc.rob(x.2$Lt, x.2$Ly, method = "huber", kernel = kernel, 
                                   bw = bw_fixed_mu, k2 = k_cand[i])
      var.huber.obj <- varfunc.rob(x.2$Lt, x.2$Ly, mu = mu.huber.obj,  
                                   method = "huber", kernel = kernel, bw = bw_fixed_var, k2 = k_cand[i])
      var_est[, i+1] <- predict(var.huber.obj, work.grid)
    }
  }) 
  df_delta <- data.frame(t = rep(work.grid, length(k_cand)+3),
                         y = c(as.numeric(var_est),
                               diag(cov.yao),
                               diag(cov.lin)),
                         k_huber = rep(factor(c("True", round(k_cand, 3), "Yao", "Lin"), 
                                              levels = c("True", round(k_cand, 3), "Yao", "Lin")), 
                                       each = length(work.grid)))
  p1 <- ggplot(df_delta, 
               aes(t, y, color = k_huber, linetype = k_huber)) +
    geom_line(size = 1) +
    labs(x = TeX("$t$"), y = TeX("$\\sigma_X^2$")) +
    scale_linetype_manual(breaks = c("True", round(k_cand, 3), "Yao", "Lin"),
                          values = c("solid", rep("dashed", 11), rep("solid", 2))) +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.title = element_blank())
  p2 <- p1 + 
    ylim(0, 5) 
  fig_delta[[2*(j-1)+1]] <- p1
  fig_delta[[2*j]] <- p2
  
  # # Compare CV results and true optimal delta
  # var.huber.obj <- varfunc.rob(x.2$Lt, x.2$Ly, mu = mu.huber.obj,
  #                              method = "huber", kernel = kernel, bw = bw_fixed_var, k2 = NULL,
  #                              cv_delta_loss = "L1", cv_K = 10)
  # var.huber.obj$obj$k2
  
  ise_var <- df_delta %>% 
    group_by(k_huber) %>% 
    summarise(ise = get_ise(y, diag(cov.true), gr))
  
  # cov.est[[sim]]$cov.obj$huber$sig2x$obj$k2
  # k_cand[which.min(ise_var$ise[-1])]
  
  
  #####################################
  ### bandwidth test for Huber loss
  #####################################
  bw_cand <- 10^seq(-2, 0, length.out = 10) * (1 - 0)/3
  var_est <- matrix(0, length(gr), length(bw_cand)+1)
  var_est[, 1] <- diag(cov.true)
  delta_selected <- cov.est[[sim]]$cov.obj$huber$sig2x$obj$k2
  delta_selected <- k_cand[which.min(ise_var$ise[-1])]
  
  fig <- list()
  delta_selected <- c(cov.est[[sim]]$cov.obj$huber$sig2x$obj$k2,
                      k_cand[which.min(ise_var$ise[-1])])
  for (k in 1:length(delta_selected)) {
    for (i in 1:length(bw_cand)) {
      print(bw_cand[i])
      
      # k <- 1 / max(abs(unlist(x.2$Ly)))
      # k <- k_cand[which.min(ise_var$ise[-1])]   # optimal delta for Huber loss from above procedure
      kernel <- "gauss"
      tryCatch({
        mu.huber.obj <- meanfunc.rob(x.2$Lt, x.2$Ly, method = "Huber", kernel = kernel, 
                                     bw = bw_cand[i], k2 = delta_selected[k])
        var.huber.obj <- varfunc.rob(x.2$Lt, x.2$Ly, mu = mu.huber.obj,  
                                     method = "Huber", kernel = kernel, 
                                     bw = bw_cand[i], k2 = delta_selected[k])
        var_est[, i+1] <- predict(var.huber.obj, gr)
      }, error = function(e){
        print(e)
      })
    }
    
    df_bw <- data.frame(t = rep(gr, length(bw_cand)+1),
                        y = as.numeric(var_est),
                        bw = rep(factor(c("True",round(bw_cand, 3)),
                                        levels = c("True",round(bw_cand, 3))), 
                                 each = length(gr)))
    fig_bw[[2*(j-1)+k]] <- ggplot(df_bw, aes(t, y, color = bw, linetype = bw)) +
      geom_line(size = 1) +
      labs(x = TeX("$t$"), y = TeX("$\\sigma_X^2$"), 
           title = round(delta_selected[k], 4)) +
      # scale_linetype_manual(breaks = c("True", round(bw_cand, 3)),
      #                       values = c("solid", rep("dashed", 10))) +
      theme_bw() +
      theme(legend.position = "bottom",
            legend.title = element_blank())
  }
}
save(list = c("fig_delta","fig_bw"),
     file = "RData/20210317_fig_huber.RData")

gridExtra::grid.arrange(grobs = fig_delta, 
                        nrow = 3)
gridExtra::grid.arrange(grobs = fig_bw, 
                        nrow = 3)


#####################################
### bandwidth test for WRM
#####################################
bw_cand <- 10^seq(-2, 0, length.out = 10) * (1 - 0)/3
df_wrm <- list()
for (sim in 1:length(ind_top3)) {
  # Get simulation data
  x <- data.list[[sim]]$x
  gr <- data.list[[sim]]$gr
  x.2 <- list(Ly = x$y,
              Lt = x$t)
  
  # estimated covariances
  work.grid <- cov.est[[sim]]$work.grid
  cov.true <- cov.est[[sim]]$cov$true
  # cov.yao <- cov.est[[sim]]$cov$yao
  # cov.lin <- cov.est[[sim]]$cov$lin
  # cov.huber <- cov.est[[sim]]$cov$huber
  # cov.wrm <- cov.est[[sim]]$cov$wrm
  
  var_est <- matrix(0, length(gr), length(bw_cand)+1)
  var_est[, 1] <- diag(cov.true)
  system.time({
    for (i in 1:length(bw_cand)) {
      print(bw_cand[i])
      
      kernel <- "gauss"
      
      # mu.wrm.obj <- meanfunc.rob(x.2$Lt, x.2$Ly, method = "WRM", kernel = kernel, bw = bw_cand[i])
      mu.wrm.obj <- cov.est[[sim]]$mu.obj$wrm
      var.wrm.obj <- varfunc.rob(x.2$Lt, x.2$Ly, mu = mu.wrm.obj,  
                                 method = "WRM", kernel = kernel, bw = bw_cand[i])
      var_est[, i+1] <- predict(var.wrm.obj, gr)
    }
  }) 
  df_bw_wrm <- data.frame(t = rep(gr, length(bw_cand)+1),
                          y = as.numeric(var_est),
                          bw = rep(factor(c("True", round(bw_cand, 3)),
                                          levels = c("True", round(bw_cand, 3))),
                                   each = length(gr)))
  
  df_wrm[[sim]] <- df_bw_wrm
}
save(list = c("bw_cand","df_wrm"),
     file = "RData/20210317_fig_wrm.RData")

fig_list <- list()
for (i in 1:3) {
  p <- ggplot(df_wrm[[i]], aes(t, y, color = bw, linetype = bw)) +
    geom_line(size = 1) +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.title = element_blank())
  
  fig_list[[i]] <- p
}
gridExtra::grid.arrange(grobs = fig_list,   # mget(): get multiple objects
                        nrow = 2)


# save(list = c("k_cand","df_delta","bw_cand","df_bw","df_bw_wrm"),
#      file = "RData/20210310_fig_1st.RData")
# save(list = c("k_cand","df_delta","bw_cand","df_bw","df_bw_wrm"),
#      file = "RData/20210310_fig_20th.RData")


# selected bandwidth
for (sim in 1:length(ind_top3)) {
  # print(cov.est[[sim]]$cov.obj$wrm$sig2x$obj$bw)
  print(cov.est[[sim]]$cov.obj$huber$sig2x$obj$bw)
  
}


