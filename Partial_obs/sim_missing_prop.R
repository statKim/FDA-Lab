###################################################
### Sensitivity analysis for Simulations
###################################################
### Download required packages
# devtools::install_github("statKim/robfpca")
# devtools::install_github("statKim/mcfda.rob")

### Load packages
library(robfpca)   # proposed methods and data generating
library(mcfda.rob)   # R-Kraus
library(tidyverse)
library(fdapace)
library(doParallel)   # parallel computing
library(doRNG)   # set.seed for foreach
library(pracma)   # subspace

# Codes can be obtained from Kraus(2015), JRSS-B.
source("sim_utills/pred.missfd.R")
source("sim_utills/simul.missfd.R")

# R-Kraus
source("sim_utills/robust_Kraus.R")

# For Boente et al.(2021), you may install the package from the follow code.
# devtools::install_github('msalibian/sparseFPCA', ref = "master")
library(sparseFPCA)   # Boente (2021) 
source("sim_utills/Boente_cov.R")


#####################################
### Simulation Model setting
### - Model 1 : "Delaigle"
### - Model 2 : "Kraus"
#####################################

### Model 1
setting <- "Delaigle"
K <- 4   # fixed number of PCs (If NULL, it is selected by PVE)
pve <- 0.95   # Not used if K is given
bw_cand <- seq(0.3, 0.4, length.out = 10)

### Model 2
setting <- "Kraus"
K <- 3   # fixed number of PCs (If NULL, it is selected by PVE)
pve <- 0.95   # Not used if K is given
bw_cand <- seq(0.01, 0.1, length.out = 10)



#####################################
### Outlier setting
### - Case 1 : Not-contaminated
### - Case 2 : t-distribution
### - Case 3 : 10% contamination
### - Case 4 : 20% contamination
#####################################

### Case 1
dist_type <- "normal"
out_type <- 1   # type of outliers (fixed; Do not change)
out_prop <- 0   # proportion of outliers

### Case 2
dist_type <- "tdist"
out_prop <- 0   # proportion of outliers

### Case 3
dist_type <- "normal"
out_type <- 1   # type of outliers (fixed; Do not change)
out_prop <- 0.1   # proportion of outliers

### Case 4
dist_type <- "normal"
out_type <- 1   # type of outliers (fixed; Do not change)
out_prop <- 0.2   # proportion of outliers


if (dist_type == "tdist") {
  file_name <- paste0("RData/", setting, "-missing-", dist_type, ".RData")
} else {
  file_name <- paste0("RData/", setting, "-missing-", dist_type,
                      "-prop", out_prop*10, ".RData")
}
file_name


#####################################
### Simulation Parameters
#####################################
num_sim <- 100    # number of simulations
data_type <- "partial"   # type of functional data
kernel <- "epanechnikov"   # kernel function for local smoothing
n_cores <- 12   # number of threads for parallel computing


#####################################
### Simulation
#####################################
f <- seq(0.1, 0.4, by = 0.05)   # candidate of missingness parameter
# d=1.4, f=0.0 => 0.002745098
# d=1.4, f=0.05 => 0.02788235
# d=1.4, f=0.1 => 0.05305686
# d=1.4, f=0.15 => 0.07828431
# d=1.4, f=0.2 => 0.1035333
# d=1.4, f=0.25 => 0.1289373
# d=1.4, f=0.3 => 0.1545
# d=1.4, f=0.35 => 0.1805157
# d=1.4, f=0.4 => 0.2064902
# d=1.4, f=0.5 => 0.2588078
mse_eigen <- matrix(NA, num_sim, length(f))
mse_eigen2 <- matrix(NA, num_sim, length(f))
mse_reconstr <- matrix(NA, num_sim, length(f))
mse_completion <- matrix(NA, num_sim, length(f))
pve_res <- matrix(NA, num_sim, length(f))
K_res <- matrix(NA, num_sim, length(f))
time_d <- matrix(NA, num_sim, length(f)) 
p_missing <- matrix(NA, num_sim, length(f)) 

colnames(mse_eigen) <- paste0("f=", f)
colnames(mse_eigen2) <- colnames(mse_eigen)
colnames(mse_reconstr) <- colnames(mse_eigen)
colnames(mse_completion) <- colnames(mse_eigen)
colnames(pve_res) <- colnames(mse_eigen)
colnames(time_d) <- colnames(mse_eigen)



### Simulation
pca.est <- list()   # pca objects
num.sim <- 0   # number of simulations
seed <- 0   # current seed
sim.seed <- rep(NA, num_sim)   # collection of seed with no error occurs
while (num.sim < num_sim) {
  seed <- seed + 1
  print(paste0("Seed: ", seed))
  
  pca.obj <- list()
  for (j in 1:length(f)) {
    set.seed(seed)
    
    #############################
    ### Data generation
    #############################
    n <- 100
    n.grid <- 51
    if (setting == 'Kraus') {
      x.2 <- sim_kraus(n = n, 
                       type = data_type,  
                       out.prop = out_prop, 
                       out.type = out_type, 
                       dist = dist_type,
                       f = f[j])
    } else if (setting == 'Delaigle') {
      x.2 <- sim_delaigle(n = n,  
                          type = data_type, 
                          out.prop = out_prop, 
                          out.type = out_type, 
                          dist = dist_type,
                          f = f[j]) 
    }
    
    x <- list2matrix(x.2)
    # matplot(t(x), type = "l")
    
    
    #############################
    ### Covariance estimation
    #############################
    skip_sim <- FALSE   # if skip_sim == TRUE, pass this seed
    work.grid <- seq(0, 1, length.out = n.grid)
    
    ### Proposed method
    start_time <- Sys.time()
    registerDoRNG(seed)
    tryCatch({
      cov.sm.obj.cv <- cv.cov_ogk(x,  
                                  K = 5, 
                                  bw_cand = bw_cand,
                                  MM = TRUE,
                                  type = 'huber')
      print(cov.sm.obj.cv$selected_bw)
      cov.obj <- cov_ogk(x,   
                         type = "huber",
                         MM = TRUE,
                         smooth = T, 
                         bw = cov.sm.obj.cv$selected_bw)
      mu.ogk.sm <- cov.obj$mean
      cov.ogk.sm <- cov.obj$cov
    }, error = function(e) { 
      print("Huber cov error")
      print(e)
      skip_sim <<- TRUE
    })
    if (skip_sim == TRUE) {
      next
    }
    end_time <- Sys.time()
    time_d[num.sim + 1, j] <- round(difftime(end_time, 
                                             start_time, 
                                             units = "secs"), 3)
    print(paste0("f = ", f[j], " : ", 
                 time_d[num.sim + 1, j],
                 " secs"))
    
    
    ### Principal component analysis
    pca.ogk.sm.obj <- funPCA(x.2$Lt, x.2$Ly,
                             mu.ogk.sm, cov.ogk.sm, sig2 = 0,
                             work.grid, PVE = pve, K = K)
    
    ### Eigen function - Compute for fixed K
    if (is.null(K)) {
      mse_eigen[num.sim + 1, j] <- rep(NA, 4)
    } else {
      if (setting == 'Delaigle') {
        eig.true <- get_delaigle_eigen(work.grid, model = 2) 
      } else if (setting == 'Kraus') {
        eig.true <- get_kraus_eigen(work.grid) 
      }
      
      # Eigen MISE
      mse_eigen[num.sim + 1, j] <- c(
        mean((check_eigen_sign(pca.ogk.sm.obj$eig.fun, eig.true) - eig.true)^2)
      )
      
      
      # Eigne angle
      mse_eigen2[num.sim + 1, j] <- c(
        mean(
          sapply(1:K, function(i){
            subspace(pca.ogk.sm.obj$eig.fun[, i], eig.true[, i])
          })
        )
      )
      
    }
    
    
    ### Curve reconstruction via PCA
    # reconstructed curves
    pred_reconstr <- list(
      predict(pca.ogk.sm.obj, K = K)
    )
    
    # MISE of reconstruction
    Not_out_ind <- which(x.2$out.ind == 0)
    sse_reconstr <- sapply(pred_reconstr, function(method){
      if (is.matrix(method)) {
        return( mean((method[Not_out_ind, ] - x.2$x.full[Not_out_ind, ])^2) )
      } else {
        return(NA)
      }
    })
    
    
    # index of non-outlying curves having missing values (Only non-outlier index)
    cand <- which(
      (apply(x, 1, function(x){ sum(is.na(x)) }) > 0) & (x.2$out.ind == 0)
    )
    
    sse_completion <- matrix(NA, length(cand), 1)
    
    for (i in 1:length(cand)) {
      ind <- cand[i]
      
      # prediction for missing parts
      pred_comp <- list(
        pred_reconstr[[1]][ind, ]
      )
      
      
      # # ISE for reconstruction of overall interval
      # sse_reconstr[i, ] <- sapply(pred_reconstr, function(method){
      #   if (is.matrix(method)) {
      #     return( mean((method[ind, ] - x.2$x.full[ind, ])^2) )
      #   } else {
      #     return(NA)
      #   }
      # })
      
      # ISE for completion
      NA_ind <- which(is.na(x[ind, ]))   # index of missing periods
      sse_completion[i, ] <- sapply(pred_comp, function(method){
        mean((method[NA_ind] - x.2$x.full[ind, NA_ind])^2)
      })
    }
    
    
    
    mse_reconstr[num.sim + 1, j] <- sse_reconstr
    mse_completion[num.sim + 1, j] <- colMeans(sse_completion)
    
    pve_res[num.sim + 1, j] <- c(
      pca.ogk.sm.obj$PVE
    )
    
    K_res[num.sim + 1, j] <- c(
      pca.ogk.sm.obj$K
    )
    
    p_missing[num.sim + 1, j] <- sum(is.na(x)) / length(x)
    
    pca.obj[[j]] <- pca.ogk.sm.obj
  }
  
  # Update number of simulations and save seed which does not occur errors
  num.sim <- num.sim + 1
  sim.seed[num.sim] <- seed
  print(paste0("Total # of simulations: ", num.sim))
  
  print(colMeans(mse_eigen, na.rm = T))
  # print(colMeans(mse_eigen2, na.rm = T))
  print(colMeans(mse_completion, na.rm = T))
  
  
  ### Save the objects
  pca.est[[num.sim]] <- pca.obj
}
save(pca.est, mse_eigen, mse_eigen2,
     mse_reconstr, mse_completion,
     K_res, pve_res, time_d, 
     p_missing,
     file = file_name)


# load("RData/Delaigle-missing-tdist.RData")
# load("RData/Kraus-missing-tdist.RData")

### Summary results
if (is.null(K)) {
  PVE_K <- K_res
} else {
  PVE_K <- pve_res
}

res <- data.frame(Method = colnames(mse_eigen)) %>% 
  # Proportion of the missingness
  left_join(data.frame(
    Method = colnames(PVE_K),
    "P(missing)" = format(round(colMeans(p_missing), 3), 3)
  ), by = "Method") %>% 
  # Eigen MISE
  left_join(data.frame(
    Method = colnames(mse_eigen),
    "Eigen MISE" = paste0(
      format(round(colMeans(mse_eigen), 3), 3),
      " (",
      format(round(apply(mse_eigen, 2, sd), 3), 3),
      ")"
    )
  ), by = "Method") %>% 
  # Eigen Angle
  left_join(data.frame(
    Method = colnames(mse_eigen2),
    "Eigen angle" = paste0(
      format(round(colMeans(mse_eigen2), 3), 3),
      " (",
      format(round(apply(mse_eigen2, 2, sd), 3), 3),
      ")"
    )
  ), by = "Method") %>% 
  # Reconstruction MISE
  left_join(data.frame(
    Method = colnames(mse_reconstr),
    "Recon MISE" = paste0(
      format(round(colMeans(mse_reconstr), 3), 3),
      " (",
      format(round(apply(mse_reconstr, 2, sd), 3), 3),
      ")"
    )
  ), by = "Method") %>% 
  # Completion MISE
  left_join(data.frame(
    Method = colnames(mse_completion),
    "Comp MISE" = paste0(
      format(round(colMeans(mse_completion), 3), 3),
      " (",
      format(round(apply(mse_completion, 2, sd), 3), 3),
      ")"
    )
  ), by = "Method")
print(res)

# Make results to LaTeX code
library(xtable)
xtable(res[, -(1:2)])




### Plot missingness (MISE vs f)
for (setting in c("Delaigle","Kraus")) {
  file_name <- paste0("RData/", setting,
                      c("-missing-normal-prop0.RData",
                        "-missing-tdist.RData",
                        "-missing-normal-prop1.RData",
                        "-missing-normal-prop2.RData"))
  model_name <- c("(i) Gaussian process",
                  "(ii) t(3) process",
                  "(iii) Gaussian process with 10% contamination",
                  "(iii) Gaussian process with 20% contamination")

  load(file_name[1])
  df <- data.frame(
    f = f,
    model = model_name[1],
    Eigenfunction = colMeans(mse_eigen),
    Completion = colMeans(mse_completion)
  )
  for (i in 2:length(model_name)) {
    load(file_name[i])
    df <- rbind(
      df,
      data.frame(
        f = f,
        model = model_name[i],
        Eigenfunction = colMeans(mse_eigen),
        Completion = colMeans(mse_completion)
      )
    )
  }
  df <- df %>%
    gather(case, mise, -f, -model)

  fig_list <- list()
  case_name <- c("Eigenfunction","Completion")
  for (i in 1:2) {
    fig_list[[i]] <- df %>%
      filter(case == case_name[i]) %>%
      ggplot(aes(x = f,
                 y = mise,
                 color = model)) +
      geom_line(size = 1.4) +
      geom_point(size = 3) +
      ggtitle( paste(case_name[i], "MISE") ) +
      labs(y = "MISE") +
      theme_bw() +
      theme(
        legend.position = "none",
        legend.title = element_blank(),
        text = element_text(size = 15),   # text size
        plot.title = element_text(hjust = 0.5),   # title center
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        # panel.border = element_rect(colour = "black", fill = NA, size = 2),   # panel border
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        # panel.grid.major = element_blank(), # get rid of major grid
        # panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        # legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
      )
  }

  library(ggpubr)
  p <- ggarrange(plotlist = fig_list,
                 ncol = 2, nrow = 1,
                 common.legend = TRUE, legend = "bottom")
  p
  ggsave(p, width = 20, height = 10,
         filename = paste0("figure/MISE_missing_", setting, ".eps"),
         bg = "transparent")
}






# load("RData/Delaigle-missing-normal-prop0.RData")
# df <- data.frame(
#   f = f,
#   model = "(i) Gaussian process",
#   Eigenfunction = colMeans(mse_eigen),
#   Completion = colMeans(mse_completion)
# )
# load("RData/Delaigle-missing-tdist.RData")
# df <- rbind(
#   df,
#   data.frame(
#     f = f,
#     model = "(ii) t(3) process",
#     Eigenfunction = colMeans(mse_eigen),
#     Completion = colMeans(mse_completion)
#   )
# )
# df <- df %>%
#   gather(case, mise, -f, -model)
# 
# p1 <- df %>%
#   filter(case == "Eigenfunction") %>%
#   ggplot(aes(x = f,
#              y = mise,
#              color = model)) +
#   geom_line(size = 1.4) +
#   geom_point(size = 3) +
#   ggtitle("Eigenfunction MISE") +
#   labs(y = "MISE") +
#   theme_bw() +
#   theme(
#     legend.position = "none",
#     legend.title = element_blank(),
#     text = element_text(size = 15),   # text size
#     plot.title = element_text(hjust = 0.5),   # title center
#     panel.background = element_rect(fill = "transparent"), # bg of the panel
#     # panel.border = element_rect(colour = "black", fill = NA, size = 2),   # panel border
#     plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
#     # panel.grid.major = element_blank(), # get rid of major grid
#     # panel.grid.minor = element_blank(), # get rid of minor grid
#     legend.background = element_rect(fill = "transparent"), # get rid of legend bg
#     # legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
#   )
# p2 <- df %>%
#   filter(case == "Completion") %>%
#   ggplot(aes(x = f,
#              y = mise,
#              color = model)) +
#   geom_line(size = 1.4) +
#   geom_point(size = 3) +
#   ggtitle("Completion MISE") +
#   labs(y = "MISE") +
#   theme_bw() +
#   theme(
#     legend.position = "none",
#     legend.title = element_blank(),
#     text = element_text(size = 15),   # text size
#     plot.title = element_text(hjust = 0.5),   # title center
#     panel.background = element_rect(fill = "transparent"), # bg of the panel
#     # panel.border = element_rect(colour = "black", fill = NA, size = 2),   # panel border
#     plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
#     # panel.grid.major = element_blank(), # get rid of major grid
#     # panel.grid.minor = element_blank(), # get rid of minor grid
#     legend.background = element_rect(fill = "transparent"), # get rid of legend bg
#     # legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
#   )
# 
# library(ggpubr)
# p <- ggarrange(p1, p2,
#                ncol = 2, nrow = 1,
#                common.legend = TRUE, legend = "bottom")
# p
# ggsave(p, width = 20, height = 10,
#        filename = "/Users/hyunsung/Desktop/Rproject/FDA-Lab/Partial_obs/figure/MISE_missing_Delaigle.eps",
#        bg = "transparent")











# load("RData/Delaigle-missing-normal-prop0.RData")
# df <- data.frame(
#   f = f,
#   Eigenfunction = colMeans(mse_eigen),
#   Completion = colMeans(mse_completion)
# ) %>% 
#   mutate(Completion = Completion/1000 + 0.035) %>% 
#   gather(case, mise, -f)
# p1 <- df %>%
#   ggplot(aes(x = f,
#              y = mise,
#              color = case)) +
#   geom_line(size = 1.4) +
#   geom_point(size = 3) +
#   scale_y_continuous(
#     name = "Eigen MISE",
#     sec.axis = sec_axis(~(.-0.035) * 1000, name = "Comp. MISE")
#   ) +
#   ggtitle("t(3)") +
#   theme_bw() +
#   # theme(legend.position = "bottom",
#   #       legend.title = element_blank(),
#   #       # text = element_text(size = 15),
#   #       plot.title = element_text(hjust = 0.5))
#   theme(
#     legend.position = "none",
#     legend.title = element_blank(),
#     # text = element_text(size = 15),   # text size
#     plot.title = element_text(hjust = 0.5),   # title center
#     panel.background = element_rect(fill = "transparent"), # bg of the panel
#     # panel.border = element_rect(colour = "black", fill = NA, size = 2),   # panel border
#     plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
#     # panel.grid.major = element_blank(), # get rid of major grid
#     # panel.grid.minor = element_blank(), # get rid of minor grid
#     legend.background = element_rect(fill = "transparent"), # get rid of legend bg
#     # legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
#   )
# 
# load("RData/Delaigle-missing-tdist.RData")
# df <- data.frame(
#   f = f,
#   Eigenfunction = colMeans(mse_eigen),
#   Completion = colMeans(mse_completion)
# ) %>% 
#   mutate(Completion = Completion/1000 + 0.035) %>% 
#   gather(case, mise, -f)
# p2 <- df %>%
#   ggplot(aes(x = f,
#              y = mise,
#              color = case)) +
#   geom_line(size = 1.4) +
#   geom_point(size = 3) +
#   scale_y_continuous(
#     name = "Eigen MISE",
#     sec.axis = sec_axis(~(.-0.035) * 1000, name = "Comp. MISE")
#   ) +
#   ggtitle("t(3)") +
#   theme_bw() +
#   # theme(legend.position = "bottom",
#   #       legend.title = element_blank(),
#   #       # text = element_text(size = 15),
#   #       plot.title = element_text(hjust = 0.5))
#   theme(
#     legend.position = "none",
#     legend.title = element_blank(),
#     # text = element_text(size = 15),   # text size
#     plot.title = element_text(hjust = 0.5),   # title center
#     panel.background = element_rect(fill = "transparent"), # bg of the panel
#     # panel.border = element_rect(colour = "black", fill = NA, size = 2),   # panel border
#     plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
#     # panel.grid.major = element_blank(), # get rid of major grid
#     # panel.grid.minor = element_blank(), # get rid of minor grid
#     legend.background = element_rect(fill = "transparent"), # get rid of legend bg
#     # legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
#   )
# 
# library(ggpubr)
# p <- ggarrange(p1, p2,
#                ncol = 2, nrow = 1,
#                common.legend = TRUE, legend = "bottom")
# p

