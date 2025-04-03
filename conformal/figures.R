########################################
### Figure for Simulation Examples
########################################
library(fdaoutlier)

n <- 100
m <- 51
p <- 20
outlier_rate <- 0.1

# # Function for simulated data from `fdaoutlier` package
# sim_ftn_list <- list(
#   function(...){ simulation_model1(q = 2, ...) },
#   function(...){ simulation_model2(q = 2, ...) },
#   function(...){ simulation_model3(q = 1.5, ...) },
#   function(...){ simulation_model5(cov_alpha2 = 0.5, ...) }
# )
# 
# pdf(file = "figures/sim_setting.pdf", width = 10, height = 3)
# par(mfrow = c(1, 4), mar = c(3, 2, 3, 1))
# set.seed(100)
# x <- sim_ftn_list[[1]](n, m, outlier_rate = outlier_rate, 
#                        plot = T, plot_title = "Model 1", xlabel = "")
# x <- sim_ftn_list[[2]](n, m, outlier_rate = outlier_rate, 
#                        plot = T, plot_title = "Model 2", xlabel = "")
# x <- sim_ftn_list[[3]](n, m, outlier_rate = outlier_rate, 
#                        plot = T, plot_title = "Model 3", xlabel = "")
# x <- sim_ftn_list[[4]](n, m, outlier_rate = outlier_rate, 
#                        plot = T, plot_title = "Model 4", xlabel = "")
# dev.off()

pdf(file = "figures/sim_setting.pdf", width = 10, height = 3)
par(mfrow = c(1, 4), mar = c(3, 2, 3, 1))
set.seed(100)
for (i in 1:4) {
  plot.foutlier_sim_mfd(foutlier_sim_mfd(n = n, m = m, p = p, 
                                         outlier_rate = outlier_rate, 
                                         model = i),
                        p = 1,
                        xlabel = "", ylabel = "", plot_title = paste("Model", i))
}
dev.off()



########################################
### Boxplots for Simulation Results
########################################
library(tidyverse)
i <- 2
for (k in 1:length(res[[i]]$bh$fdr)) {
  sub_fdr <- res[[i]]$bh$fdr[[k]] %>% 
    gather(key = "Type") %>% 
    mutate(Type = factor(Type, levels = c("marg","simes","asymp")),
           Method = names(res[[i]]$bh$fdr)[k])
  sub_tpr <- res[[i]]$bh$tpr[[k]] %>% 
    gather(key = "Type") %>% 
    mutate(Type = factor(Type, levels = c("marg","simes","asymp")),
           Method = names(res[[i]]$bh$tpr)[k])
  
  if (k == 1) {
    df_fdr <- sub_fdr
    df_tpr <- sub_tpr
  } else {
    df_fdr <- rbind(df_fdr, sub_fdr)
    df_tpr <- rbind(df_tpr, sub_tpr)
  }
}
fig_fdr <- ggplot(df_fdr, aes(x = Method, y = value, fill = Type)) +
  geom_boxplot(alpha = 0.3) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4)) +
  lims(y = c(0, 0.4)) +
  labs(title = "FDR") +
  theme_minimal()
fig_tpr <- ggplot(df_tpr, aes(x = Method, y = value, fill = Type)) +
  geom_boxplot(alpha = 0.3) +
  labs(title = "Power") +
  theme_minimal()
fig_list <- list(fig_fdr, fig_tpr)
gridExtra::grid.arrange(fig_fdr, fig_tpr, nrow = 1)


df_fdr$Measure <- "FDR"
df_tpr$Measure <- "Power"
df <- rbind(df_fdr, df_tpr) 
ggplot(df, aes(x = Method, y = value, fill = Type)) +
  geom_boxplot(alpha = 0.3) +
  facet_grid(.~Measure) +
  # theme_minimal() +
  theme(legend.position = "bottom",
        # panel.background=element_rect(fill="white")
        )



########################################
### Figures for preprocessing of ADHD-200
########################################
load("RData/ADHD-200_Schaefer17_400regions.RData")
dim(X)
n <- dim(X)[1]   # number of curves
m <- dim(X)[2]   # number of timepoints
p <- dim(X)[3]  # number of functional variables


### Make additional outliers from ADHD group
# 1st, 2nd derivatives
X_deriv_1 <- array(NA, c(m-1, n, p))
X_deriv_2 <- array(NA, c(m-2, n, p))
for (i in 1:p) {
  X_deriv_1[, , i] <- apply(X[, , i], 1, diff)
  X_deriv_2[, , i] <- apply(X[, , i], 1, function(x){ diff(diff(x)) })
}


# Candidates of outliers from ADHD group
# - Choose depth based outliers
idx_adhd <- which(y == 1)
type_depth <- "projdepth"
depth_values <- list(
  mfd(aperm(X[idx_adhd, , ], c(2,1,3)), type = type_depth, 
      depthOptions = list(type = "Rotation"))$MFDdepthX,
  mfd(X_deriv_1[, idx_adhd, ], type = type_depth, 
      depthOptions = list(type = "Rotation"))$MFDdepthX,
  mfd(X_deriv_2[, idx_adhd, ], type = type_depth, 
      depthOptions = list(type = "Rotation"))$MFDdepthX
)
idx_adhd_idx <- lapply(depth_values, function(x){ order(x)[1:20] }) %>% 
  unlist() %>% 
  unique()
idx_adhd_cand <- idx_adhd[idx_adhd_idx]


region_idx <- 1  # region index

# Fast Fourier Transform with smoothing splines
X_fft <- apply(X[, , region_idx], 1, function(x) {
  Mod(fft(x)) * (2/m)
}) %>% 
  t()

X_fft_sm <- apply(X[, , region_idx], 1, function(x) {
  smooth.spline(Mod(fft(x)) * (2/m))$y
}) %>% 
  t()


plot_curves <- function(data, y, title = "title", legend_pos = "bottomright",
                        x_axis_labels = seq(0, 1, length.out = 5), legend = TRUE, ylim = NULL) {
  x_col <- ifelse(y == 0, "grey61", "#D55E00")
  m <- ncol(data)
  
  if (is.null(ylim)) {
    ylim <- range(data) + c(-.5*sd(data), .5*sd(data))
  }
  
  matplot(t(data), type = "l", lty = 1, lwd = 1.5, col = x_col, xlab = "",
          ylim = ylim,
          col.lab = "gray20", axes = F)
  grid()
  matlines(t(data), lty = 1, lwd = 1.5, col = x_col)
  axis(1, at = seq(1, m, length.out = length(x_axis_labels)), labels = x_axis_labels, 
       col = "white", col.ticks = "grey61",
       lwd.ticks = .5, tck = -0.025,
       cex.axis = 0.9, col.axis = "gray30")
  axis(2, col = "white", col.ticks = "grey61",
       lwd.ticks = .5, tck = -0.025,
       cex.axis = 0.9, col.axis = "gray30")
  box(col = "grey51")
  mtext(title, 3, adj = 0.5, line = 1, cex = 1.5, col = "gray20")
  
  if (isTRUE(legend)) {
    legend(legend_pos, legend = c("normal", "outlier"),
           lty = c("solid", "solid"),
           lwd = c(.4, 1.3),
           col = c("grey61", "#D55E00"),
           text.col = "gray40", bty = "n",
           box.lwd = .1, xjust = 0, inset = .01)
  }
}

set.seed(123)
set.seed(1234)
set.seed(12345)
set.seed(1000)
set.seed(500)
set.seed(555)
idx <- c(
  sample(which(y == 0), 5),
  # sample(which(y == 1), 5)
  sample(idx_adhd_cand, 5)
)
x_sample <- X[idx, , region_idx]
x_fft_sample <- X_fft[idx, ]
x_fft_sm_sample <- X_fft_sm[idx, ]
y_sample <- y[idx]

# pdf(file = "figures/adhd200_fft.pdf", width = 10, height = 5)
par(mfrow = c(2, 1), mar = c(3, 2, 3, 1))
plot_curves(x_sample, y_sample, title = "Raw")
plot_curves(x_fft_sm_sample, y_sample, title = "FFT-smooth",
            x_axis_labels = round(seq(1, m, length.out = 5)))
# dev.off()


# set.seed(123)
# idx <- c(
#   sample(which(y == 0), 3),
#   sample(idx_adhd_cand, 1)
# )
# x_sample <- X[idx, , region_idx]
# y_sample <- y[idx]
# 
# # pdf("figures/test.pdf", width = 10, height = 5)
# # par(mfrow = c(1, 1), mar = c(3, 2, 3, 1))
# plot_curves(x_sample, y_sample, title = "")
# # dev.off()
# par(mfrow = c(2, 2))
# for (i in 1:4) {
#   plot_curves(matrix(x_sample[i, ], nrow = 1), y_sample[i], title = "", 
#               legend = F, ylim = range(x_sample) + c(-.5*sd(x_sample), .5*sd(x_sample)))  
# }



########################################
### Boxplots for ADHD-200 Results
########################################
library(tidyverse)
load("RData/adhd_200_res.RData")

preproc_list <- c("Raw","FFT","FFT-smooth")
for (i in c(1,3)) {
  for (k in 1:length(res[[i]]$bh$fdr)) {
    sub_fdr <- res[[i]]$bh$fdr[[k]] %>% 
      gather(key = "Type") %>% 
      filter(Type == "marg") %>% 
      mutate(Method = names(res[[i]]$bh$fdr)[k],
             Type = preproc_list[i])
    sub_tpr <- res[[i]]$bh$tpr[[k]] %>% 
      gather(key = "Type") %>% 
      filter(Type == "marg") %>% 
      mutate(Method = names(res[[i]]$bh$tpr)[k],
             Type = preproc_list[i])
    
    if (i == 1 & k == 1) {
      df_fdr <- sub_fdr
      df_tpr <- sub_tpr
    } else {
      df_fdr <- rbind(df_fdr, sub_fdr)
      df_tpr <- rbind(df_tpr, sub_tpr)
    }
  }
  
  for (k in 1:length(res[[i]]$comparison$fdr)) {
    sub_fdr <- res[[i]]$comparison$fdr %>% 
      gather(key = "Method") %>% 
      mutate(Type = preproc_list[i])
    sub_tpr <- res[[i]]$comparison$tpr %>% 
      gather(key = "Method") %>% 
      mutate(Type = preproc_list[i])
    
    df_fdr <- rbind(df_fdr, sub_fdr)
    df_tpr <- rbind(df_tpr, sub_tpr)
  }
}



df_fdr <- df_fdr %>% 
  filter(Method %in% c("T_projdepth","projdepth",
                       # "T_hdepth","hdepth",
                       "esssup","seq")) %>% 
  mutate(Method = ifelse(Method == "T_projdepth", "MFD-agg",
                         ifelse(Method == "projdepth", "MFD",
                                ifelse(Method == "esssup", "Esssup", "Seq-trans")))) %>% 
  mutate(Method = factor(Method, levels = c("MFD-agg","MFD",
                                            "Esssup","Seq-trans")),
         Type = factor(Type, levels = c("Raw","FFT-smooth")))
  # mutate(Method = ifelse(Method == "T_projdepth", "MFD-agg (projdepth)",
  #                        ifelse(Method == "projdepth", "MFD (projdepth)",
  #                               ifelse(Method == "T_hdepth", "MFD-agg (hdepth)",
  #                                      ifelse(Method == "hdepth", "MFD (hdepth)",
  #                                             ifelse(Method == "esssup", "Esssup", "Seq-trans")))))) %>% 
  # mutate(Method = factor(Method, levels = c("MFD-agg (projdepth)","MFD (projdepth)",
  #                                           "MFD-agg (hdepth)","MFD (hdepth)",
  #                                           "Esssup","Seq-trans")),
  #        Type = factor(Type, levels = c("Raw","FFT-smooth")))
df_tpr <- df_tpr %>% 
  filter(Method %in% c("T_projdepth","projdepth",
                       # "T_hdepth","hdepth",
                       "esssup","seq")) %>% 
  mutate(Method = ifelse(Method == "T_projdepth", "MFD-agg",
                         ifelse(Method == "projdepth", "MFD",
                                ifelse(Method == "esssup", "Esssup", "Seq-trans")))) %>% 
  mutate(Method = factor(Method, levels = c("MFD-agg","MFD",
                                            "Esssup","Seq-trans")),
         Type = factor(Type, levels = c("Raw","FFT-smooth")))
  # mutate(Method = ifelse(Method == "T_projdepth", "MFD-agg (projdepth)",
  #                        ifelse(Method == "projdepth", "MFD (projdepth)",
  #                               ifelse(Method == "T_hdepth", "MFD-agg (hdepth)",
  #                                      ifelse(Method == "hdepth", "MFD (hdepth)",
  #                                             ifelse(Method == "esssup", "Esssup", "Seq-trans")))))) %>% 
  # mutate(Method = factor(Method, levels = c("MFD-agg (projdepth)","MFD (projdepth)",
  #                                           "MFD-agg (hdepth)","MFD (hdepth)",
  #                                           "Esssup","Seq-trans")),
  #        Type = factor(Type, levels = c("Raw","FFT-smooth")))


# fig_fdr <- ggplot(df_fdr, aes(x = Method, y = value, fill = Type)) +
#   geom_boxplot(alpha = 0.3) +
#   # scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4)) +
#   # lims(y = c(0, 0.4)) +
#   labs(title = "FDR", x = "", y = "") +
#   theme_bw() +
#   theme(plot.title = element_text(hjust = 0.5, size = 17),
#         axis.text.x = element_text(size = 13, face = "plain"),
#         axis.text.y = element_text(size = 13, face = "plain"),
#         legend.text = element_text(size = 13, face = "plain"),
#         legend.title = element_blank()) +
#   guides(fill = guide_legend(nrow = 1))
# fig_tpr <- ggplot(df_tpr, aes(x = Method, y = value, fill = Type)) +
#   geom_boxplot(alpha = 0.3) +
#   labs(title = "Power", x = "", y = "") +
#   theme_bw() +
#   theme(plot.title = element_text(hjust = 0.5, size = 17),
#         axis.text.x = element_text(size = 13, face = "plain"),
#         axis.text.y = element_text(size = 13, face = "plain"),
#         legend.text = element_text(size = 13, face = "plain"),
#         legend.title = element_blank()) +
#   guides(fill = guide_legend(nrow = 1))
# fig_list <- list(fig_fdr, fig_tpr)
# ggpubr::ggarrange(plotlist = fig_list, 
#                   nrow = 1,
#                   common.legend = TRUE, legend = "bottom")


df <- df_fdr %>% 
  mutate(Measure = "FDR") %>% 
  rbind(df_tpr %>% 
          mutate(Measure = "Power")) %>% 
  filter(Type == "FFT-smooth")
fig <- ggplot(df, aes(x = Method, y = value, fill = Type, color = Type)) +
  geom_boxplot(alpha = 0.3, size = 0.7) +
  facet_wrap(. ~ Measure, scales = "free") +
  # scale_fill_manual(values = c("Raw" = "#E69F00", "FFT-smooth" = "#56B4E9")) +
  # scale_color_manual(values = c("Raw" = "#E69F00", "FFT-smooth" = "#56B4E9")) +
  # scale_fill_manual(values = c("Raw" = "#F8766D", "FFT-smooth" = "#00BFC4")) +
  # scale_color_manual(values = c("Raw" = "#F8766D", "FFT-smooth" = "#00BFC4")) +
  scale_fill_manual(values = c("Raw" = "#D55E00", "FFT-smooth" = "#0072B2")) +
  scale_color_manual(values = c("Raw" = "#D55E00", "FFT-smooth" = "#0072B2")) +
  labs(title = "", x = "", y = "") +
  theme_bw() +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 17),
        axis.text.x = element_text(size = 14, face = "plain"),
        axis.text.y = element_text(size = 13, face = "plain"),
        legend.text = element_text(size = 14, face = "plain"),
        strip.text =  element_text(size = 15, face = "plain"),
        legend.title = element_blank()) +
  guides(fill = guide_legend(nrow = 1))
fig
ggsave("figures/res_adhd200.pdf", width = 10, height = 6)
