# Summary the results
num_decimal <- 3
res2 <- list()
for (i in 1:length(res)) {
  res2[[i]] <- list(
    fdr = cbind(res[[i]]$bh$fdr,
                res[[i]]$comparison$fdr),
    tpr = cbind(res[[i]]$bh$tpr,
                res[[i]]$comparison$tpr)
  )
}

res3 <- lapply(res2, function(sim){
  sub <- paste0(
    rbind(fdr = colMeans(sim$fdr),
          tpr = colMeans(sim$tpr)) %>% 
      round(num_decimal) %>% 
      format(nsmall = num_decimal),
    " (",
    rbind(fdr = apply(sim$fdr, 2, sd),
          tpr = apply(sim$tpr, 2, sd)) %>% 
      round(num_decimal) %>% 
      format(nsmall = num_decimal),
    ")"
  )
  dim(sub) <- c(2, ncol(sim$fdr))
  rownames(sub) <- c("FDR","TPR")
  colnames(sub) <- colnames(sim$fdr)
  sub <- data.frame(sub)
  # sub[, c("esssup.marg","hdepth.marg","projdepth.marg","seq_trans.marg","ms","seq")]
  # sub[, c("T_projdepth.marg","T_hdepth.marg","T_mbd.marg",
  #         "esssup.marg","hdepth.marg","projdepth.marg",
  #         "ms","seq","ms_all","seq_all")]
  # sub[, c("T_projdepth.marg","T_hdepth.marg",
  #         "esssup.marg","projdepth.marg","hdepth.marg","ms","seq")]
  sub[, c("T_projdepth.marg","T_projdepth.simes","T_projdepth.asymp",
          "esssup.marg","esssup.simes","esssup.asymp",
          "ms","seq")] %>% 
    t()
})

cbind(res3[[1]], res3[[2]], res3[[3]], res3[[4]]) %>% 
  xtable::xtable()


par(mfrow = c(1, 2))
matplot(apply(conf_pvalue, 2, sort), pch = 1, ylab = "")
lines((1:n_test)/n_test * alpha, col = "blue", lwd = 2)
matplot(apply(conf_pvalue, 2, sort)[1:100, ], pch = 1, ylab = "")
lines((1:n_test)/n_test * alpha, col = "blue", lwd = 2)
legend("topleft", c("Marginal","Simes","Asymptotic","BH cutoff"),
       col = c(1,2,3,"blue"), pch = 1)





### Clean
# [[1]]
#       esssup.marg  esssup.simes  esssup.asymp projdepth.marg projdepth.simes projdepth.asymp
# FDR 0.077 (0.048) 0.000 (0.000) 0.001 (0.007)  0.090 (0.057)   0.017 (0.026)   0.016 (0.025)
# TPR 0.757 (0.178) 0.000 (0.000) 0.016 (0.110)  1.000 (0.000)   1.000 (0.000)   1.000 (0.000)
#     T_projdepth.marg T_projdepth.simes T_projdepth.asymp            ms           seq
# FDR    0.079 (0.049)     0.019 (0.026)     0.019 (0.026) 0.168 (0.047) 0.000 (0.000)
# TPR    1.000 (0.000)     1.000 (0.002)     1.000 (0.002) 1.000 (0.000) 0.761 (0.063)
# 
# [[2]]
#       esssup.marg  esssup.simes  esssup.asymp projdepth.marg projdepth.simes projdepth.asymp
# FDR 0.038 (0.091) 0.000 (0.000) 0.000 (0.000)  0.000 (0.000)   0.000 (0.000)   0.000 (0.000)
# TPR 0.049 (0.111) 0.000 (0.000) 0.000 (0.000)  0.000 (0.000)   0.000 (0.000)   0.000 (0.000)
#     T_projdepth.marg T_projdepth.simes T_projdepth.asymp            ms           seq
# FDR    0.076 (0.055)     0.018 (0.031)     0.017 (0.024) 0.636 (0.119) 0.010 (0.100)
# TPR    1.000 (0.000)     1.000 (0.000)     1.000 (0.000) 0.113 (0.042) 0.001 (0.004)
# 
# [[3]]
#       esssup.marg  esssup.simes  esssup.asymp projdepth.marg projdepth.simes projdepth.asymp
# FDR 0.034 (0.082) 0.000 (0.000) 0.000 (0.000)  0.072 (0.061)   0.000 (0.000)   0.003 (0.017)
# TPR 0.053 (0.111) 0.000 (0.000) 0.000 (0.000)  0.575 (0.169)   0.000 (0.000)   0.022 (0.128)
#     T_projdepth.marg T_projdepth.simes T_projdepth.asymp            ms           seq
# FDR    0.076 (0.055)     0.018 (0.031)     0.017 (0.024) 0.174 (0.047) 0.001 (0.008)
# TPR    1.000 (0.000)     0.999 (0.003)     0.999 (0.003) 0.931 (0.041) 0.110 (0.046)
# 
# [[4]]
#       esssup.marg  esssup.simes  esssup.asymp projdepth.marg projdepth.simes projdepth.asymp
# FDR 0.000 (0.000) 0.000 (0.000) 0.000 (0.000)  0.000 (0.000)   0.000 (0.000)   0.000 (0.000)
# TPR 0.000 (0.000) 0.000 (0.000) 0.000 (0.000)  0.000 (0.000)   0.000 (0.000)   0.000 (0.000)
#     T_projdepth.marg T_projdepth.simes T_projdepth.asymp            ms           seq
# FDR    0.082 (0.054)     0.019 (0.025)     0.019 (0.025) 1.000 (0.000) 0.010 (0.100)
# TPR    1.000 (0.000)     1.000 (0.000)     1.000 (0.000) 0.000 (0.000) 0.000 (0.000)



### Mixed
# [[1]]
#       esssup.marg  esssup.simes  esssup.asymp projdepth.marg projdepth.simes projdepth.asymp
# FDR 0.367 (0.082) 0.209 (0.137) 0.225 (0.109)  0.085 (0.057)   0.022 (0.032)   0.022 (0.032)
# TPR 0.989 (0.016) 0.737 (0.417) 0.890 (0.248)  1.000 (0.000)   0.520 (0.502)   1.000 (0.000)
#     T_projdepth.marg T_projdepth.simes T_projdepth.asymp            ms           seq
# FDR    0.076 (0.054)     0.019 (0.026)     0.019 (0.026) 0.169 (0.052) 0.000 (0.000)
# TPR    1.000 (0.000)     0.480 (0.502)     0.999 (0.005) 1.000 (0.000) 0.754 (0.075)
# 
# [[2]]
#       esssup.marg  esssup.simes  esssup.asymp projdepth.marg projdepth.simes projdepth.asymp
# FDR 0.434 (0.088) 0.038 (0.132) 0.108 (0.190)  0.000 (0.000)   0.000 (0.000)   0.000 (0.000)
# TPR 0.659 (0.135) 0.060 (0.206) 0.163 (0.288)  0.000 (0.000)   0.000 (0.000)   0.000 (0.000)
#     T_projdepth.marg T_projdepth.simes T_projdepth.asymp            ms           seq
# FDR    0.082 (0.057)     0.026 (0.030)     0.026 (0.030) 0.689 (0.144) 0.010 (0.100)
# TPR    1.000 (0.000)     0.660 (0.476)     1.000 (0.000) 0.068 (0.036) 0.000 (0.002)
# 
# [[3]]
#       esssup.marg  esssup.simes  esssup.asymp projdepth.marg projdepth.simes projdepth.asymp
# FDR 0.431 (0.078) 0.026 (0.106) 0.113 (0.188)  0.001 (0.007)   0.000 (0.000)   0.000 (0.000)
# TPR 0.638 (0.119) 0.040 (0.160) 0.165 (0.275)  0.023 (0.073)   0.000 (0.000)   0.000 (0.000)
#     T_projdepth.marg T_projdepth.simes T_projdepth.asymp            ms           seq
# FDR    0.081 (0.058)     0.024 (0.030)     0.024 (0.030) 0.172 (0.046) 0.003 (0.033)
# TPR    1.000 (0.000)     0.600 (0.492)     0.999 (0.005) 0.931 (0.039) 0.031 (0.025)
# 
# [[4]]
#       esssup.marg  esssup.simes  esssup.asymp projdepth.marg projdepth.simes projdepth.asymp
# FDR 0.969 (0.171) 0.010 (0.100) 0.040 (0.197)  0.000 (0.000)   0.000 (0.000)   0.000 (0.000)
# TPR 0.001 (0.006) 0.000 (0.000) 0.000 (0.000)  0.000 (0.000)   0.000 (0.000)   0.000 (0.000)
#     T_projdepth.marg T_projdepth.simes T_projdepth.asymp            ms           seq
# FDR    0.082 (0.057)     0.021 (0.025)     0.021 (0.025) 1.000 (0.000) 0.000 (0.000)
# TPR    1.000 (0.000)     0.560 (0.499)     1.000 (0.000) 0.000 (0.000) 0.000 (0.000)



### fMRI (30번)

### Raw X
# T_projdepth.marg T_hdepth.marg   esssup.marg projdepth.marg   hdepth.marg            ms           seq
# FDR    0.116 (0.123) 0.125 (0.109) 0.116 (0.103)  0.099 (0.110) 0.113 (0.108) 0.028 (0.108) 0.059 (0.080)
# TPR    0.597 (0.364) 0.569 (0.313) 0.683 (0.287)  0.719 (0.268) 0.647 (0.293) 0.100 (0.063) 0.450 (0.099)

#     T_projdepth.marg T_hdepth.marg
# FDR    0.112 (0.120) 0.123 (0.108)
# TPR    0.619 (0.345) 0.569 (0.315)


### FFT
# T_projdepth.marg T_hdepth.marg   esssup.marg projdepth.marg   hdepth.marg            ms           seq
# FDR    0.107 (0.109) 0.108 (0.096) 0.102 (0.074)  0.107 (0.100) 0.078 (0.111) 0.000 (0.000) 0.052 (0.066)
# TPR    0.819 (0.103) 0.783 (0.186) 0.742 (0.149)  0.836 (0.080) 0.431 (0.333) 0.117 (0.068) 0.569 (0.122)

#     T_projdepth.marg T_hdepth.marg
# FDR    0.111 (0.124) 0.119 (0.123)
# TPR    0.647 (0.287) 0.542 (0.296)

### FFT-sm
# T_projdepth.marg T_hdepth.marg   esssup.marg projdepth.marg   hdepth.marg            ms           seq
# FDR    0.118 (0.107) 0.095 (0.094) 0.094 (0.102)  0.088 (0.102) 0.064 (0.127) 0.000 (0.000) 0.100 (0.090)
# TPR    0.636 (0.170) 0.672 (0.199) 0.614 (0.206)  0.728 (0.120) 0.192 (0.338) 0.136 (0.089) 0.639 (0.103)

#     T_projdepth.marg T_hdepth.marg
# FDR    0.102 (0.102) 0.066 (0.130)
# TPR    0.731 (0.140) 0.183 (0.330)

