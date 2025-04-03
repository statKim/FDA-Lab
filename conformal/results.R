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



############################
### Simulation Results
############################

### Clean (Not scaling for T_; average)
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



### Mixed (Not scaling for T_; average)
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


### Clean (Scaling for T_; weighted)
# [[1]]
#         marg         simes         asymp
# FDR 0.068 (0.066) 0.000 (0.000) 0.000 (0.000)
# TPR 0.399 (0.224) 0.000 (0.000) 0.000 (0.000)
# 
# [[2]]
#         marg         simes         asymp
# FDR 0.071 (0.044) 0.013 (0.017) 0.013 (0.017)
# TPR 1.000 (0.000) 1.000 (0.000) 1.000 (0.000)
# 
# [[3]]
#         marg         simes         asymp
# FDR 0.071 (0.044) 0.013 (0.017) 0.013 (0.017)
# TPR 1.000 (0.000) 1.000 (0.000) 1.000 (0.000)
# 
# [[4]]
#         marg         simes         asymp
# FDR 0.078 (0.049) 0.016 (0.029) 0.014 (0.023)
# TPR 1.000 (0.000) 1.000 (0.000) 1.000 (0.000)






############################
### Real Data Results
############################

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



### fMRI (30번) - random하게 20개
# T_projdepth.marg T_hdepth.marg   esssup.marg projdepth.marg   hdepth.marg            ms           seq
# FDR    0.075 (0.144) 0.027 (0.087) 0.064 (0.130)  0.083 (0.146) 0.043 (0.122) 0.075 (0.219) 0.047 (0.111)
# TPR    0.083 (0.152) 0.047 (0.109) 0.158 (0.185)  0.117 (0.185) 0.058 (0.120) 0.067 (0.063) 0.247 (0.086)

### fMRI (100번) - raw, 1st, 2nd derivative 20개
# T_projdepth.marg T_hdepth.marg   esssup.marg projdepth.marg   hdepth.marg            ms           seq
# FDR    0.102 (0.137) 0.092 (0.141) 0.105 (0.151)  0.093 (0.130) 0.109 (0.161) 0.032 (0.143) 0.065 (0.116)
# TPR    0.236 (0.238) 0.196 (0.203) 0.321 (0.223)  0.256 (0.245) 0.229 (0.219) 0.068 (0.070) 0.338 (0.102)

### fMRI-FFT (30번) - raw, 1st, 2nd derivative 20개
# T_projdepth.marg T_hdepth.marg   esssup.marg projdepth.marg   hdepth.marg            ms           seq
# FDR    0.103 (0.125) 0.057 (0.128) 0.100 (0.120)  0.117 (0.134) 0.152 (0.174) 0.000 (0.000) 0.080 (0.100)
# TPR    0.444 (0.167) 0.114 (0.229) 0.392 (0.167)  0.486 (0.161) 0.286 (0.250) 0.111 (0.080) 0.403 (0.103)

### fMRI-FFT-sm (30번) - raw, 1st, 2nd derivative 20개
# T_projdepth.marg T_hdepth.marg   esssup.marg projdepth.marg   hdepth.marg            ms           seq
# FDR    0.094 (0.128) 0.049 (0.133) 0.093 (0.127)  0.070 (0.114) 0.072 (0.146) 0.058 (0.204) 0.123 (0.121)
# TPR    0.442 (0.145) 0.094 (0.216) 0.356 (0.162)  0.447 (0.144) 0.097 (0.190) 0.122 (0.061) 0.450 (0.104)


### fMRI (20번) - raw, 1st, 2nd derivative를 depth 기준으로 20개
# T_projdepth.marg T_hdepth.marg   esssup.marg projdepth.marg   hdepth.marg            ms           seq
# FDR    0.095 (0.096) 0.110 (0.088) 0.112 (0.092)  0.101 (0.094) 0.093 (0.091) 0.000 (0.000) 0.063 (0.105)
# TPR    0.688 (0.294) 0.637 (0.320) 0.700 (0.198)  0.767 (0.305) 0.696 (0.340) 0.062 (0.066) 0.462 (0.116)




##############################
### 2025.02.25 ~
##############################

##### derivative가 도움이 되는지 확인

### fMRI (30번) - 기존 outlier 방법 (seq는 c("O","D1","D2"))
#     T_projdepth.marg projdepth.marg projdepth_1d.marg projdepth_2d.marg T_hdepth.marg   hdepth.marg
# FDR    0.121 (0.124)  0.107 (0.110)     0.125 (0.117)     0.111 (0.127) 0.138 (0.106) 0.121 (0.114)
# TPR    0.606 (0.318)  0.753 (0.219)     0.619 (0.296)     0.472 (0.325) 0.586 (0.254) 0.631 (0.307)
#     hdepth_1d.marg hdepth_2d.marg   esssup.marg            ms           seq
# FDR  0.124 (0.117)  0.144 (0.133) 0.127 (0.116) 0.024 (0.093) 0.072 (0.082)
# TPR  0.478 (0.340)  0.394 (0.254) 0.686 (0.268) 0.100 (0.091) 0.486 (0.124)


### fMRI (30번) - raw, 1st, 2nd derivative를 depth 기준으로 20개 (seq는 c("O","D1","D2"))
#     T_projdepth.marg projdepth.marg projdepth_1d.marg projdepth_2d.marg T_hdepth.marg   hdepth.marg
# FDR    0.094 (0.094)  0.093 (0.092)     0.100 (0.093)     0.102 (0.109) 0.105 (0.096) 0.091 (0.095)
# TPR    0.628 (0.353)  0.731 (0.321)     0.581 (0.367)     0.475 (0.339) 0.575 (0.355) 0.686 (0.321)
#     hdepth_1d.marg hdepth_2d.marg   esssup.marg            ms           seq
# FDR  0.108 (0.112)  0.111 (0.119) 0.108 (0.092) 0.033 (0.183) 0.078 (0.099)
# TPR  0.506 (0.344)  0.394 (0.307) 0.692 (0.236) 0.072 (0.072) 0.478 (0.118)

### fMRI (100번) - raw, 1st, 2nd derivative를 depth 기준으로 20개 (seq는 c("O","D1","D2"))
#     T_projdepth.marg projdepth.marg projdepth_1d.marg projdepth_2d.marg T_hdepth.marg   hdepth.marg
# FDR    0.105 (0.091)  0.104 (0.087)     0.108 (0.096)     0.114 (0.111) 0.111 (0.098) 0.101 (0.097)
# TPR    0.684 (0.319)  0.807 (0.237)     0.637 (0.335)     0.507 (0.321) 0.618 (0.301) 0.708 (0.289)
#     hdepth_1d.marg hdepth_2d.marg   esssup.marg            ms           seq
# FDR  0.103 (0.108)  0.116 (0.120) 0.111 (0.101) 0.033 (0.159) 0.070 (0.095)
# TPR  0.492 (0.325)  0.428 (0.280) 0.672 (0.254) 0.086 (0.076) 0.504 (0.114)


### fMRI (30번) - align 새로 한 것(scale = F) - raw, 1st, 2nd derivative를 depth 기준으로 20개
# T_projdepth.marg projdepth.marg projdepth_1d.marg projdepth_2d.marg T_hdepth.marg   hdepth.marg hdepth_1d.marg
# FDR    0.094 (0.127)  0.101 (0.123)     0.107 (0.138)     0.094 (0.135) 0.102 (0.131) 0.103 (0.137)  0.105 (0.134)
# TPR    0.347 (0.353)  0.463 (0.379)     0.327 (0.345)     0.277 (0.294) 0.387 (0.367) 0.427 (0.384)  0.367 (0.345)
# hdepth_2d.marg   esssup.marg            ms           seq
# FDR  0.085 (0.133) 0.065 (0.111) 0.172 (0.351) 0.182 (0.207)
# TPR  0.237 (0.305) 0.257 (0.385) 0.073 (0.069) 0.300 (0.136)

### fMRI (30번) - align 새로 한 것(scale = T) - raw, 1st, 2nd derivative를 depth 기준으로 20개
# T_projdepth.marg projdepth.marg projdepth_1d.marg projdepth_2d.marg T_hdepth.marg   hdepth.marg
# FDR    0.080 (0.123)  0.093 (0.132)     0.090 (0.125)     0.080 (0.134) 0.079 (0.133) 0.107 (0.129)
# TPR    0.350 (0.436)  0.377 (0.434)     0.350 (0.414)     0.210 (0.309) 0.273 (0.365) 0.363 (0.373)
# hdepth_1d.marg hdepth_2d.marg   esssup.marg            ms           seq
# FDR  0.104 (0.149)  0.090 (0.149) 0.080 (0.159) 0.113 (0.251) 0.165 (0.133)
# TPR  0.237 (0.345)  0.167 (0.266) 0.123 (0.242) 0.153 (0.086) 0.483 (0.178)


### 시뮬레이션 (30번) - max로 aggregate
# [[1]]
#     T_projdepth.marg T_hdepth.marg   esssup.marg projdepth.marg projdepth_1d.marg projdepth_2d.marg
# FDR    0.097 (0.068) 0.000 (0.000) 0.082 (0.045)  0.097 (0.068)     0.000 (0.000)     0.027 (0.146)
# TPR    1.000 (0.000) 0.000 (0.000) 0.729 (0.216)  1.000 (0.000)     0.000 (0.000)     0.001 (0.007)
#       hdepth.marg hdepth_1d.marg hdepth_2d.marg
# FDR 0.000 (0.000)  0.000 (0.000)  0.000 (0.000)
# TPR 0.000 (0.000)  0.000 (0.000)  0.000 (0.000)
# 
# [[2]]
#     T_projdepth.marg T_hdepth.marg   esssup.marg projdepth.marg projdepth_1d.marg projdepth_2d.marg
# FDR    0.079 (0.045) 0.000 (0.000) 0.045 (0.105)  0.000 (0.000)     0.080 (0.060)     0.081 (0.057)
# TPR    1.000 (0.000) 0.000 (0.000) 0.063 (0.125)  0.000 (0.000)     1.000 (0.000)     1.000 (0.000)
#       hdepth.marg hdepth_1d.marg hdepth_2d.marg
# FDR 0.000 (0.000)  0.034 (0.088)  0.055 (0.081)
# TPR 0.000 (0.000)  0.036 (0.083)  0.212 (0.210)
# 
# [[3]]
#     T_projdepth.marg T_hdepth.marg   esssup.marg projdepth.marg projdepth_1d.marg projdepth_2d.marg
# FDR    0.084 (0.052) 0.000 (0.000) 0.039 (0.090)  0.086 (0.052)     0.080 (0.060)     0.081 (0.057)
# TPR    0.717 (0.168) 0.000 (0.000) 0.057 (0.120)  0.618 (0.114)     1.000 (0.000)     1.000 (0.000)
#       hdepth.marg hdepth_1d.marg hdepth_2d.marg
# FDR 0.000 (0.000)  0.013 (0.073)  0.023 (0.090)
# TPR 0.000 (0.000)  0.004 (0.022)  0.015 (0.060)
# 
# [[4]]
#     T_projdepth.marg T_hdepth.marg   esssup.marg projdepth.marg projdepth_1d.marg projdepth_2d.marg
# FDR    0.081 (0.054) 0.000 (0.000) 0.000 (0.000)  0.000 (0.000)     0.078 (0.049)     0.083 (0.058)
# TPR    1.000 (0.000) 0.000 (0.000) 0.000 (0.000)  0.000 (0.000)     1.000 (0.000)     1.000 (0.000)
#       hdepth.marg hdepth_1d.marg hdepth_2d.marg
# FDR 0.000 (0.000)  0.062 (0.046)  0.062 (0.056)
# TPR 0.000 (0.000)  1.000 (0.000)  1.000 (0.000)


### 시뮬레이션 (30번) - average
# [[1]]
#     T_projdepth.marg T_hdepth.marg   esssup.marg projdepth.marg projdepth_1d.marg projdepth_2d.marg
# FDR    0.081 (0.053) 0.000 (0.000) 0.082 (0.045)  0.097 (0.068)     0.000 (0.000)     0.027 (0.146)
# TPR    1.000 (0.000) 0.000 (0.000) 0.729 (0.216)  1.000 (0.000)     0.000 (0.000)     0.001 (0.007)
#       hdepth.marg hdepth_1d.marg hdepth_2d.marg
# FDR 0.000 (0.000)  0.000 (0.000)  0.000 (0.000)
# TPR 0.000 (0.000)  0.000 (0.000)  0.000 (0.000)
# 
# [[2]]
#     T_projdepth.marg T_hdepth.marg   esssup.marg projdepth.marg projdepth_1d.marg projdepth_2d.marg
# FDR    0.069 (0.047) 0.053 (0.087) 0.045 (0.105)  0.000 (0.000)     0.080 (0.060)     0.081 (0.057)
# TPR    1.000 (0.000) 0.112 (0.148) 0.063 (0.125)  0.000 (0.000)     1.000 (0.000)     1.000 (0.000)
#       hdepth.marg hdepth_1d.marg hdepth_2d.marg
# FDR 0.000 (0.000)  0.034 (0.088)  0.055 (0.081)
# TPR 0.000 (0.000)  0.036 (0.083)  0.212 (0.210)
# 
# [[3]]
#     T_projdepth.marg T_hdepth.marg   esssup.marg projdepth.marg projdepth_1d.marg projdepth_2d.marg
# FDR    0.069 (0.047) 0.037 (0.103) 0.039 (0.090)  0.086 (0.052)     0.080 (0.060)     0.081 (0.057)
# TPR    1.000 (0.000) 0.022 (0.057) 0.057 (0.120)  0.618 (0.114)     1.000 (0.000)     1.000 (0.000)
#       hdepth.marg hdepth_1d.marg hdepth_2d.marg
# FDR 0.000 (0.000)  0.013 (0.073)  0.023 (0.090)
# TPR 0.000 (0.000)  0.004 (0.022)  0.015 (0.060)
# 
# [[4]]
#     T_projdepth.marg T_hdepth.marg   esssup.marg projdepth.marg projdepth_1d.marg projdepth_2d.marg
# FDR    0.076 (0.050) 0.000 (0.000) 0.000 (0.000)  0.000 (0.000)     0.078 (0.049)     0.083 (0.058)
# TPR    1.000 (0.000) 0.000 (0.000) 0.000 (0.000)  0.000 (0.000)     1.000 (0.000)     1.000 (0.000)
#       hdepth.marg hdepth_1d.marg hdepth_2d.marg
# FDR 0.000 (0.000)  0.062 (0.046)  0.062 (0.056)
# TPR 0.000 (0.000)  1.000 (0.000)  1.000 (0.000)


### 시뮬레이션 (10번) - average
# [[1]]
#     T_projdepth.marg T_projdepth.simes T_projdepth.asymp T_hdepth.marg T_hdepth.simes T_hdepth.asymp
# FDR    0.068 (0.050)     0.027 (0.022)     0.027 (0.022) 0.000 (0.000)  0.000 (0.000)  0.000 (0.000)
# TPR    1.000 (0.000)     1.000 (0.000)     1.000 (0.000) 0.000 (0.000)  0.000 (0.000)  0.000 (0.000)
#       esssup.marg  esssup.simes  esssup.asymp projdepth.marg projdepth.simes projdepth.asymp
# FDR 0.098 (0.051) 0.000 (0.000) 0.005 (0.015)  0.094 (0.068)   0.012 (0.016)   0.012 (0.016)
# TPR 0.744 (0.181) 0.000 (0.000) 0.078 (0.247)  1.000 (0.000)   1.000 (0.000)   1.000 (0.000)
#     projdepth_1d.marg projdepth_1d.simes projdepth_1d.asymp projdepth_2d.marg projdepth_2d.simes
# FDR     0.000 (0.000)      0.000 (0.000)      0.000 (0.000)     0.080 (0.253)      0.000 (0.000)
# TPR     0.000 (0.000)      0.000 (0.000)      0.000 (0.000)     0.004 (0.013)      0.000 (0.000)
#     projdepth_2d.asymp   hdepth.marg  hdepth.simes  hdepth.asymp hdepth_1d.marg hdepth_1d.simes
# FDR      0.000 (0.000) 0.000 (0.000) 0.000 (0.000) 0.000 (0.000)  0.000 (0.000)   0.000 (0.000)
# TPR      0.000 (0.000) 0.000 (0.000) 0.000 (0.000) 0.000 (0.000)  0.000 (0.000)   0.000 (0.000)
#     hdepth_1d.asymp hdepth_2d.marg hdepth_2d.simes hdepth_2d.asymp            ms           seq
# FDR   0.000 (0.000)  0.000 (0.000)   0.000 (0.000)   0.000 (0.000) 1.000 (0.000) 1.000 (0.000)
# TPR   0.000 (0.000)  0.000 (0.000)   0.000 (0.000)   0.000 (0.000) 0.000 (0.000) 0.000 (0.000)
# 
# [[2]]
#     T_projdepth.marg T_projdepth.simes T_projdepth.asymp T_hdepth.marg T_hdepth.simes T_hdepth.asymp
# FDR    0.078 (0.035)     0.023 (0.019)     0.023 (0.019) 0.034 (0.074)  0.000 (0.000)  0.000 (0.000)
# TPR    1.000 (0.000)     1.000 (0.000)     1.000 (0.000) 0.056 (0.120)  0.000 (0.000)  0.000 (0.000)
#       esssup.marg  esssup.simes  esssup.asymp projdepth.marg projdepth.simes projdepth.asymp
# FDR 0.070 (0.120) 0.000 (0.000) 0.000 (0.000)  0.000 (0.000)   0.000 (0.000)   0.000 (0.000)
# TPR 0.114 (0.157) 0.000 (0.000) 0.000 (0.000)  0.000 (0.000)   0.000 (0.000)   0.000 (0.000)
#     projdepth_1d.marg projdepth_1d.simes projdepth_1d.asymp projdepth_2d.marg projdepth_2d.simes
# FDR     0.085 (0.074)      0.020 (0.033)      0.020 (0.033)     0.070 (0.054)      0.021 (0.030)
# TPR     1.000 (0.000)      1.000 (0.000)      1.000 (0.000)     1.000 (0.000)      1.000 (0.000)
#     projdepth_2d.asymp   hdepth.marg  hdepth.simes  hdepth.asymp hdepth_1d.marg hdepth_1d.simes
# FDR      0.021 (0.030) 0.000 (0.000) 0.000 (0.000) 0.000 (0.000)  0.101 (0.132)   0.000 (0.000)
# TPR      1.000 (0.000) 0.000 (0.000) 0.000 (0.000) 0.000 (0.000)  0.086 (0.113)   0.000 (0.000)
#     hdepth_1d.asymp hdepth_2d.marg hdepth_2d.simes hdepth_2d.asymp            ms           seq
# FDR   0.000 (0.000)  0.067 (0.094)   0.000 (0.000)   0.000 (0.000) 1.000 (0.000) 0.200 (0.422)
# TPR   0.000 (0.000)  0.192 (0.194)   0.000 (0.000)   0.000 (0.000) 0.000 (0.000) 0.000 (0.000)
# 
# [[3]]
#     T_projdepth.marg T_projdepth.simes T_projdepth.asymp T_hdepth.marg T_hdepth.simes T_hdepth.asymp
# FDR    0.078 (0.035)     0.023 (0.019)     0.023 (0.019) 0.033 (0.105)  0.000 (0.000)  0.000 (0.000)
# TPR    1.000 (0.000)     1.000 (0.000)     1.000 (0.000) 0.016 (0.051)  0.000 (0.000)  0.000 (0.000)
#       esssup.marg  esssup.simes  esssup.asymp projdepth.marg projdepth.simes projdepth.asymp
# FDR 0.075 (0.111) 0.000 (0.000) 0.000 (0.000)  0.107 (0.061)   0.000 (0.000)   0.000 (0.000)
# TPR 0.120 (0.165) 0.000 (0.000) 0.000 (0.000)  0.614 (0.106)   0.000 (0.000)   0.000 (0.000)
#     projdepth_1d.marg projdepth_1d.simes projdepth_1d.asymp projdepth_2d.marg projdepth_2d.simes
# FDR     0.085 (0.074)      0.020 (0.033)      0.020 (0.033)     0.070 (0.054)      0.021 (0.030)
# TPR     1.000 (0.000)      1.000 (0.000)      1.000 (0.000)     1.000 (0.000)      1.000 (0.000)
#     projdepth_2d.asymp   hdepth.marg  hdepth.simes  hdepth.asymp hdepth_1d.marg hdepth_1d.simes
# FDR      0.021 (0.030) 0.000 (0.000) 0.000 (0.000) 0.000 (0.000)  0.040 (0.126)   0.000 (0.000)
# TPR      1.000 (0.000) 0.000 (0.000) 0.000 (0.000) 0.000 (0.000)  0.012 (0.038)   0.000 (0.000)
#     hdepth_1d.asymp hdepth_2d.marg hdepth_2d.simes hdepth_2d.asymp            ms           seq
# FDR   0.000 (0.000)  0.030 (0.095)   0.000 (0.000)   0.000 (0.000) 1.000 (0.000) 1.000 (0.000)
# TPR   0.000 (0.000)  0.014 (0.044)   0.000 (0.000)   0.000 (0.000) 0.000 (0.000) 0.000 (0.000)
# 
# [[4]]
#     T_projdepth.marg T_projdepth.simes T_projdepth.asymp T_hdepth.marg T_hdepth.simes T_hdepth.asymp
# FDR    0.052 (0.030)     0.012 (0.016)     0.012 (0.016) 0.000 (0.000)  0.000 (0.000)  0.000 (0.000)
# TPR    1.000 (0.000)     1.000 (0.000)     1.000 (0.000) 0.000 (0.000)  0.000 (0.000)  0.000 (0.000)
#       esssup.marg  esssup.simes  esssup.asymp projdepth.marg projdepth.simes projdepth.asymp
# FDR 0.000 (0.000) 0.000 (0.000) 0.000 (0.000)  0.000 (0.000)   0.000 (0.000)   0.000 (0.000)
# TPR 0.000 (0.000) 0.000 (0.000) 0.000 (0.000)  0.000 (0.000)   0.000 (0.000)   0.000 (0.000)
#     projdepth_1d.marg projdepth_1d.simes projdepth_1d.asymp projdepth_2d.marg projdepth_2d.simes
# FDR     0.089 (0.055)      0.010 (0.016)      0.010 (0.016)     0.068 (0.054)      0.027 (0.041)
# TPR     1.000 (0.000)      1.000 (0.000)      1.000 (0.000)     1.000 (0.000)      1.000 (0.000)
#     projdepth_2d.asymp   hdepth.marg  hdepth.simes  hdepth.asymp hdepth_1d.marg hdepth_1d.simes
# FDR      0.027 (0.041) 0.000 (0.000) 0.000 (0.000) 0.000 (0.000)  0.079 (0.073)   0.027 (0.041)
# TPR      1.000 (0.000) 0.000 (0.000) 0.000 (0.000) 0.000 (0.000)  1.000 (0.000)   1.000 (0.000)
#     hdepth_1d.asymp hdepth_2d.marg hdepth_2d.simes hdepth_2d.asymp            ms           seq
# FDR   0.027 (0.041)  0.039 (0.034)   0.006 (0.009)   0.006 (0.009) 1.000 (0.000) 0.000 (0.000)
# TPR   1.000 (0.000)  1.000 (0.000)   1.000 (0.000)   1.000 (0.000) 0.000 (0.000) 0.000 (0.000)



##############################
### 2025.03.04 ~
##############################

### Sequential transformation 논문의 세팅에 적용한 결과 (10번 반복; 단, 여기서는 multivariate)
# [[1]]
# T_projdepth.marg T_projdepth.simes T_projdepth.asymp T_hdepth.marg T_hdepth.simes T_hdepth.asymp
# FDR    0.068 (0.050)     0.027 (0.022)     0.027 (0.022) 0.000 (0.000)  0.000 (0.000)  0.000 (0.000)
# TPR    1.000 (0.000)     1.000 (0.000)     1.000 (0.000) 0.032 (0.101)  0.000 (0.000)  0.000 (0.000)
# esssup.marg  esssup.simes  esssup.asymp projdepth.marg projdepth.simes projdepth.asymp
# FDR 0.094 (0.047) 0.014 (0.016) 0.014 (0.016)  0.094 (0.068)   0.012 (0.016)   0.012 (0.016)
# TPR 0.996 (0.008) 0.880 (0.310) 0.880 (0.310)  1.000 (0.000)   0.998 (0.006)   0.998 (0.006)
# projdepth_1d.marg projdepth_1d.simes projdepth_1d.asymp projdepth_2d.marg projdepth_2d.simes
# FDR     0.097 (0.059)      0.006 (0.013)      0.006 (0.013)     0.073 (0.047)      0.024 (0.041)
# TPR     1.000 (0.000)      1.000 (0.000)      1.000 (0.000)     1.000 (0.000)      1.000 (0.000)
# projdepth_2d.asymp   hdepth.marg  hdepth.simes  hdepth.asymp hdepth_1d.marg hdepth_1d.simes
# FDR      0.024 (0.041) 0.000 (0.000) 0.000 (0.000) 0.000 (0.000)  0.098 (0.209)   0.000 (0.000)
# TPR      1.000 (0.000) 0.000 (0.000) 0.000 (0.000) 0.000 (0.000)  0.028 (0.062)   0.000 (0.000)
# hdepth_1d.asymp hdepth_2d.marg hdepth_2d.simes hdepth_2d.asymp            ms           seq
# FDR   0.000 (0.000)  0.000 (0.000)   0.000 (0.000)   0.000 (0.000) 1.000 (0.000) 1.000 (0.000)
# TPR   0.000 (0.000)  0.000 (0.000)   0.000 (0.000)   0.000 (0.000) 0.000 (0.000) 0.000 (0.000)
# 
# [[2]]
# T_projdepth.marg T_projdepth.simes T_projdepth.asymp T_hdepth.marg T_hdepth.simes T_hdepth.asymp
# FDR    0.068 (0.050)     0.027 (0.022)     0.027 (0.022) 0.010 (0.032)  0.000 (0.000)  0.000 (0.000)
# TPR    1.000 (0.000)     1.000 (0.000)     1.000 (0.000) 0.156 (0.165)  0.000 (0.000)  0.000 (0.000)
# esssup.marg  esssup.simes  esssup.asymp projdepth.marg projdepth.simes projdepth.asymp
# FDR 0.090 (0.052) 0.004 (0.014) 0.009 (0.019)  0.000 (0.000)   0.000 (0.000)   0.000 (0.000)
# TPR 0.914 (0.128) 0.088 (0.278) 0.174 (0.367)  0.000 (0.000)   0.000 (0.000)   0.000 (0.000)
# projdepth_1d.marg projdepth_1d.simes projdepth_1d.asymp projdepth_2d.marg projdepth_2d.simes
# FDR     0.097 (0.059)      0.006 (0.013)      0.006 (0.013)     0.073 (0.047)      0.024 (0.041)
# TPR     1.000 (0.000)      1.000 (0.000)      1.000 (0.000)     1.000 (0.000)      1.000 (0.000)
# projdepth_2d.asymp   hdepth.marg  hdepth.simes  hdepth.asymp hdepth_1d.marg hdepth_1d.simes
# FDR      0.024 (0.041) 0.000 (0.000) 0.000 (0.000) 0.000 (0.000)  0.091 (0.123)   0.000 (0.000)
# TPR      1.000 (0.000) 0.000 (0.000) 0.000 (0.000) 0.000 (0.000)  0.134 (0.175)   0.000 (0.000)
# hdepth_1d.asymp hdepth_2d.marg hdepth_2d.simes hdepth_2d.asymp            ms           seq
# FDR   0.000 (0.000)  0.036 (0.052)   0.000 (0.000)   0.000 (0.000) 1.000 (0.000) 0.300 (0.483)
# TPR   0.000 (0.000)  0.370 (0.288)   0.000 (0.000)   0.000 (0.000) 0.000 (0.000) 0.000 (0.000)
# 
# [[3]]
# T_projdepth.marg T_projdepth.simes T_projdepth.asymp T_hdepth.marg T_hdepth.simes T_hdepth.asymp
# FDR    0.052 (0.030)     0.012 (0.016)     0.012 (0.016) 0.061 (0.093)  0.000 (0.000)  0.000 (0.000)
# TPR    1.000 (0.000)     1.000 (0.000)     1.000 (0.000) 0.414 (0.091)  0.000 (0.000)  0.000 (0.000)
# esssup.marg  esssup.simes  esssup.asymp projdepth.marg projdepth.simes projdepth.asymp
# FDR 0.000 (0.000) 0.000 (0.000) 0.000 (0.000)  0.000 (0.000)   0.000 (0.000)   0.000 (0.000)
# TPR 0.000 (0.000) 0.000 (0.000) 0.000 (0.000)  0.000 (0.000)   0.000 (0.000)   0.000 (0.000)
# projdepth_1d.marg projdepth_1d.simes projdepth_1d.asymp projdepth_2d.marg projdepth_2d.simes
# FDR     0.089 (0.055)      0.010 (0.016)      0.010 (0.016)     0.068 (0.054)      0.027 (0.041)
# TPR     1.000 (0.000)      1.000 (0.000)      1.000 (0.000)     1.000 (0.000)      1.000 (0.000)
# projdepth_2d.asymp   hdepth.marg  hdepth.simes  hdepth.asymp hdepth_1d.marg hdepth_1d.simes
# FDR      0.027 (0.041) 0.000 (0.000) 0.000 (0.000) 0.000 (0.000)  0.079 (0.073)   0.027 (0.041)
# TPR      1.000 (0.000) 0.000 (0.000) 0.000 (0.000) 0.000 (0.000)  1.000 (0.000)   1.000 (0.000)
# hdepth_1d.asymp hdepth_2d.marg hdepth_2d.simes hdepth_2d.asymp            ms           seq
# FDR   0.027 (0.041)  0.039 (0.034)   0.006 (0.009)   0.006 (0.009) 1.000 (0.000) 0.000 (0.000)
# TPR   1.000 (0.000)  1.000 (0.000)   1.000 (0.000)   1.000 (0.000) 0.000 (0.000) 0.000 (0.000)
# 
# [[4]]
# T_projdepth.marg T_projdepth.simes T_projdepth.asymp T_hdepth.marg T_hdepth.simes T_hdepth.asymp
# FDR    0.077 (0.046)     0.015 (0.022)     0.015 (0.022) 0.033 (0.081)  0.000 (0.000)  0.000 (0.000)
# TPR    1.000 (0.000)     1.000 (0.000)     1.000 (0.000) 0.040 (0.085)  0.000 (0.000)  0.000 (0.000)
# esssup.marg  esssup.simes  esssup.asymp projdepth.marg projdepth.simes projdepth.asymp
# FDR 0.102 (0.067) 0.034 (0.052) 0.021 (0.024)  0.048 (0.035)   0.014 (0.016)   0.014 (0.016)
# TPR 1.000 (0.000) 1.000 (0.000) 1.000 (0.000)  1.000 (0.000)   1.000 (0.000)   1.000 (0.000)
# projdepth_1d.marg projdepth_1d.simes projdepth_1d.asymp projdepth_2d.marg projdepth_2d.simes
# FDR     0.081 (0.058)      0.013 (0.030)      0.013 (0.030)     0.000 (0.000)      0.000 (0.000)
# TPR     1.000 (0.000)      1.000 (0.000)      1.000 (0.000)     0.000 (0.000)      0.000 (0.000)
# projdepth_2d.asymp   hdepth.marg  hdepth.simes  hdepth.asymp hdepth_1d.marg hdepth_1d.simes
# FDR      0.000 (0.000) 0.000 (0.000) 0.000 (0.000) 0.000 (0.000)  0.053 (0.069)   0.000 (0.000)
# TPR      0.000 (0.000) 0.000 (0.000) 0.000 (0.000) 0.000 (0.000)  0.442 (0.157)   0.000 (0.000)
# hdepth_1d.asymp hdepth_2d.marg hdepth_2d.simes hdepth_2d.asymp            ms           seq
# FDR   0.000 (0.000)  0.000 (0.000)   0.000 (0.000)   0.000 (0.000) 1.000 (0.000) 1.000 (0.000)
# TPR   0.000 (0.000)  0.000 (0.000)   0.000 (0.000)   0.000 (0.000) 0.000 (0.000) 0.000 (0.000)
# 
# [[5]]
# T_projdepth.marg T_projdepth.simes T_projdepth.asymp T_hdepth.marg T_hdepth.simes T_hdepth.asymp
# FDR    0.900 (0.000)     0.900 (0.000)     0.900 (0.000) 0.532 (0.369)  0.000 (0.000)  0.000 (0.000)
# TPR    1.000 (0.000)     1.000 (0.000)     1.000 (0.000) 0.140 (0.121)  0.000 (0.000)  0.000 (0.000)
# esssup.marg  esssup.simes  esssup.asymp projdepth.marg projdepth.simes projdepth.asymp
# FDR 0.900 (0.000) 0.900 (0.000) 0.900 (0.000)  0.900 (0.000)   0.900 (0.000)   0.900 (0.000)
# TPR 1.000 (0.000) 1.000 (0.000) 1.000 (0.000)  1.000 (0.000)   1.000 (0.000)   1.000 (0.000)
# projdepth_1d.marg projdepth_1d.simes projdepth_1d.asymp projdepth_2d.marg projdepth_2d.simes
# FDR     0.430 (0.182)      0.229 (0.220)      0.230 (0.202)     0.000 (0.000)      0.000 (0.000)
# TPR     0.896 (0.114)      0.500 (0.448)      0.578 (0.419)     0.000 (0.000)      0.000 (0.000)
# projdepth_2d.asymp   hdepth.marg  hdepth.simes  hdepth.asymp hdepth_1d.marg hdepth_1d.simes
# FDR      0.000 (0.000) 0.900 (0.000) 0.720 (0.379) 0.720 (0.379)  0.053 (0.142)   0.000 (0.000)
# TPR      0.000 (0.000) 1.000 (0.000) 0.800 (0.422) 0.800 (0.422)  0.066 (0.148)   0.000 (0.000)
# hdepth_1d.asymp hdepth_2d.marg hdepth_2d.simes hdepth_2d.asymp            ms           seq
# FDR   0.000 (0.000)  0.000 (0.000)   0.000 (0.000)   0.000 (0.000) 1.000 (0.000) 1.000 (0.000)
# TPR   0.000 (0.000)  0.000 (0.000)   0.000 (0.000)   0.000 (0.000) 0.000 (0.000) 0.000 (0.000)
# 
# [[6]]
# T_projdepth.marg T_projdepth.simes T_projdepth.asymp T_hdepth.marg T_hdepth.simes T_hdepth.asymp
# FDR    0.097 (0.078)     0.015 (0.020)     0.015 (0.020) 0.000 (0.000)  0.000 (0.000)  0.000 (0.000)
# TPR    1.000 (0.000)     1.000 (0.000)     1.000 (0.000) 0.000 (0.000)  0.000 (0.000)  0.000 (0.000)
# esssup.marg  esssup.simes  esssup.asymp projdepth.marg projdepth.simes projdepth.asymp
# FDR 0.072 (0.049) 0.010 (0.016) 0.010 (0.016)  0.100 (0.082)   0.017 (0.019)   0.017 (0.019)
# TPR 1.000 (0.000) 1.000 (0.000) 1.000 (0.000)  1.000 (0.000)   1.000 (0.000)   1.000 (0.000)
# projdepth_1d.marg projdepth_1d.simes projdepth_1d.asymp projdepth_2d.marg projdepth_2d.simes
# FDR     0.090 (0.072)      0.015 (0.020)      0.015 (0.020)     0.105 (0.102)      0.016 (0.015)
# TPR     1.000 (0.000)      1.000 (0.000)      1.000 (0.000)     1.000 (0.000)      1.000 (0.000)
# projdepth_2d.asymp   hdepth.marg  hdepth.simes  hdepth.asymp hdepth_1d.marg hdepth_1d.simes
# FDR      0.016 (0.015) 0.000 (0.000) 0.000 (0.000) 0.000 (0.000)  0.000 (0.000)   0.000 (0.000)
# TPR      1.000 (0.000) 0.000 (0.000) 0.000 (0.000) 0.000 (0.000)  0.000 (0.000)   0.000 (0.000)
# hdepth_1d.asymp hdepth_2d.marg hdepth_2d.simes hdepth_2d.asymp            ms           seq
# FDR   0.000 (0.000)  0.000 (0.000)   0.000 (0.000)   0.000 (0.000) 1.000 (0.000) 1.000 (0.000)
# TPR   0.000 (0.000)  0.000 (0.000)   0.000 (0.000)   0.000 (0.000) 0.000 (0.000) 0.000 (0.000)



##############################
### 2025.03.17 ~ 
##############################
### ADHD-200 (100 regions; 30번 반복; ms는 seq 결과임) raw; FFT
# T_projdepth.marg projdepth.marg projdepth_1d.marg projdepth_2d.marg T_hdepth.marg   hdepth.marg
# FDR    0.112 (0.153)  0.102 (0.141)     0.118 (0.145)     0.134 (0.166) 0.097 (0.149) 0.096 (0.147)
# TPR    0.367 (0.317)  0.397 (0.362)     0.353 (0.304)     0.290 (0.252) 0.313 (0.340) 0.353 (0.350)
# hdepth_1d.marg hdepth_2d.marg   esssup.marg            ms           seq
# FDR  0.109 (0.156)  0.098 (0.165) 0.094 (0.136) 0.153 (0.156) 0.153 (0.156)
# TPR  0.287 (0.315)  0.177 (0.242) 0.287 (0.345) 0.490 (0.158) 0.490 (0.158)
# 
# T_projdepth.marg projdepth.marg projdepth_1d.marg projdepth_2d.marg T_hdepth.marg   hdepth.marg
# FDR    0.118 (0.113)  0.125 (0.130)     0.108 (0.108)     0.113 (0.119) 0.134 (0.121) 0.074 (0.180)
# TPR    0.543 (0.331)  0.447 (0.310)     0.593 (0.339)     0.523 (0.353) 0.543 (0.352) 0.090 (0.254)
# hdepth_1d.marg hdepth_2d.marg   esssup.marg            ms           seq
# FDR  0.110 (0.120)  0.122 (0.122) 0.133 (0.143) 0.214 (0.113) 0.214 (0.113)
# TPR  0.480 (0.377)  0.473 (0.342) 0.307 (0.305) 0.590 (0.149) 0.590 (0.149)


### ADHD-200 (82 regions; 30번 반복) 
# T_projdepth.marg projdepth.marg projdepth_1d.marg projdepth_2d.marg T_hdepth.marg   hdepth.marg
# FDR    0.112 (0.156)  0.110 (0.144)     0.114 (0.149)     0.137 (0.170) 0.103 (0.151) 0.135 (0.150)
# TPR    0.280 (0.306)  0.387 (0.366)     0.263 (0.299)     0.240 (0.242) 0.380 (0.370) 0.387 (0.388)
# hdepth_1d.marg hdepth_2d.marg   esssup.marg            ms           seq
# FDR  0.092 (0.152)  0.101 (0.170) 0.104 (0.141) 0.083 (0.265) 0.142 (0.140)
# TPR  0.270 (0.301)  0.197 (0.243) 0.333 (0.348) 0.033 (0.055) 0.493 (0.134)


# ### ADHD-200 (200 regions; 30번 반복) 
# T_projdepth.marg projdepth.marg projdepth_1d.marg projdepth_2d.marg T_hdepth.marg   hdepth.marg
# FDR    0.100 (0.162)  0.095 (0.155)     0.084 (0.147)     0.086 (0.157) 0.059 (0.140) 0.098 (0.178)
# TPR    0.183 (0.235)  0.213 (0.271)     0.147 (0.218)     0.133 (0.197) 0.123 (0.237) 0.153 (0.264)
# hdepth_1d.marg hdepth_2d.marg   esssup.marg            ms           seq
# FDR  0.049 (0.123)  0.083 (0.179) 0.056 (0.117) 0.199 (0.228) 0.199 (0.228)
# TPR  0.093 (0.176)  0.073 (0.144) 0.127 (0.249) 0.340 (0.138) 0.340 (0.138)


# ### ADHD-200 (400 regions; 30번 반복) raw; FFT; FFT-sm
# T_projdepth.marg projdepth.marg projdepth_1d.marg projdepth_2d.marg T_hdepth.marg   hdepth.marg
# FDR    0.124 (0.142)  0.119 (0.136)     0.111 (0.138)     0.101 (0.148) 0.118 (0.146) 0.126 (0.144)
# TPR    0.343 (0.330)  0.440 (0.355)     0.313 (0.335)     0.217 (0.257) 0.273 (0.316) 0.427 (0.375)
# hdepth_1d.marg hdepth_2d.marg   esssup.marg            ms           seq
# FDR  0.085 (0.129)  0.069 (0.137) 0.091 (0.130) 0.157 (0.164) 0.157 (0.164)
# TPR  0.223 (0.275)  0.140 (0.227) 0.287 (0.385) 0.410 (0.142) 0.410 (0.142)
# focsvm.marg (n_basis=10) n_basis=5      n_basis=30    <- scaling 안한 경우
# FDR 0.115 (0.111)       0.117 (0.108)   0.114 (0.111)
# TPR 0.667 (0.361)       0.690 (0.338)   0.670 (0.362)
# focsvm_5.marg focsvm_10.marg focsvm_20.marg focsvm_30.marg    <- scaling한 경우
# FDR 0.106 (0.115)  0.120 (0.120)  0.116 (0.115)  0.113 (0.117)
# TPR 0.757 (0.280)  0.783 (0.328)  0.800 (0.325)  0.800 (0.325)

# T_projdepth.marg projdepth.marg projdepth_1d.marg projdepth_2d.marg T_hdepth.marg   hdepth.marg
# FDR    0.116 (0.124)  0.104 (0.132)     0.122 (0.136)     0.130 (0.139) 0.094 (0.110) 0.112 (0.157)
# TPR    0.570 (0.310)  0.547 (0.257)     0.613 (0.292)     0.567 (0.341) 0.473 (0.269) 0.230 (0.305)
# hdepth_1d.marg hdepth_2d.marg   esssup.marg            ms           seq
# FDR  0.098 (0.130)  0.100 (0.124) 0.104 (0.142) 0.159 (0.134) 0.159 (0.134)
# TPR  0.413 (0.331)  0.350 (0.325) 0.320 (0.339) 0.520 (0.145) 0.520 (0.145)

# T_projdepth.marg projdepth.marg projdepth_1d.marg projdepth_2d.marg T_hdepth.marg   hdepth.marg
# FDR    0.096 (0.128)  0.119 (0.130)     0.090 (0.130)     0.086 (0.119) 0.095 (0.122) 0.058 (0.129)
# TPR    0.623 (0.253)  0.613 (0.213)     0.583 (0.270)     0.410 (0.354) 0.590 (0.220) 0.163 (0.261)
# hdepth_1d.marg hdepth_2d.marg   esssup.marg            ms           seq
# FDR  0.093 (0.125)  0.109 (0.141) 0.133 (0.150) 0.207 (0.121) 0.207 (0.121)
# TPR  0.443 (0.346)  0.487 (0.336) 0.303 (0.347) 0.600 (0.134) 0.600 (0.134)
# focsvm.marg (n_basis=20)
# FDR 0.111 (0.176)
# TPR 0.190 (0.287)


### ADHD-200 (400 regions; 100번 반복) raw; FFT; FFT-sm
# T_projdepth.marg projdepth.marg projdepth_1d.marg projdepth_2d.marg T_hdepth.marg   hdepth.marg
# FDR    0.117 (0.141)  0.127 (0.141)     0.104 (0.141)     0.112 (0.149) 0.104 (0.138) 0.112 (0.129)
# TPR    0.323 (0.310)  0.451 (0.353)     0.276 (0.305)     0.259 (0.240) 0.270 (0.298) 0.379 (0.360)
# hdepth_1d.marg hdepth_2d.marg   esssup.marg   focsvm.marg         seq
# FDR  0.092 (0.132)  0.090 (0.150) 0.101 (0.141) 0.117 (0.108)  0.162 (0.156)
# TPR  0.235 (0.289)  0.196 (0.241) 0.279 (0.352) 0.833 (0.271)  0.428 (0.142)

# T_projdepth.marg projdepth.marg projdepth_1d.marg projdepth_2d.marg T_hdepth.marg   hdepth.marg
# FDR    0.101 (0.118)  0.094 (0.122)     0.104 (0.125)     0.109 (0.128) 0.087 (0.110) 0.085 (0.154)
# TPR    0.566 (0.275)  0.531 (0.246)     0.587 (0.293)     0.564 (0.319) 0.475 (0.295) 0.202 (0.302)
# hdepth_1d.marg hdepth_2d.marg   esssup.marg   focsvm.marg           seq
# FDR  0.088 (0.117)  0.089 (0.117) 0.075 (0.122) 0.107 (0.129) 0.191 (0.142)
# TPR  0.385 (0.335)  0.351 (0.331) 0.227 (0.305) 0.446 (0.381) 0.512 (0.149)

# T_projdepth.marg projdepth.marg projdepth_1d.marg projdepth_2d.marg T_hdepth.marg   hdepth.marg
# FDR    0.105 (0.136)  0.112 (0.137)     0.106 (0.125)     0.123 (0.136) 0.113 (0.132) 0.049 (0.109)
# TPR    0.613 (0.257)  0.597 (0.246)     0.585 (0.274)     0.429 (0.314) 0.567 (0.245) 0.148 (0.256)
# hdepth_1d.marg hdepth_2d.marg   esssup.marg   focsvm.marg           seq
# FDR  0.099 (0.127)  0.103 (0.127) 0.113 (0.149) 0.120 (0.135) 0.206 (0.138)
# TPR  0.471 (0.328)  0.480 (0.331) 0.250 (0.323) 0.482 (0.370) 0.573 (0.146)


### fMRI (100번 반복) raw; FFT; FFT-sm
# T_projdepth.marg projdepth.marg projdepth_1d.marg projdepth_2d.marg T_hdepth.marg   hdepth.marg
# FDR    0.105 (0.091)  0.104 (0.087)     0.108 (0.096)     0.114 (0.111) 0.111 (0.098) 0.101 (0.097)
# TPR    0.684 (0.319)  0.807 (0.237)     0.637 (0.335)     0.507 (0.321) 0.618 (0.301) 0.708 (0.289)
# hdepth_1d.marg hdepth_2d.marg   esssup.marg   focsvm.marg            ms           seq
# FDR  0.103 (0.108)  0.116 (0.120) 0.111 (0.101) 0.100 (0.095) 0.033 (0.159) 0.070 (0.095)
# TPR  0.492 (0.325)  0.428 (0.280) 0.672 (0.254) 0.869 (0.144) 0.086 (0.076) 0.504 (0.114)

# T_projdepth.marg projdepth.marg projdepth_1d.marg projdepth_2d.marg T_hdepth.marg   hdepth.marg
# FDR    0.100 (0.086)  0.103 (0.089)     0.115 (0.092)     0.105 (0.092) 0.107 (0.094) 0.077 (0.109)
# TPR    0.871 (0.090)  0.889 (0.110)     0.833 (0.153)     0.840 (0.137) 0.858 (0.182) 0.408 (0.357)
# hdepth_1d.marg hdepth_2d.marg   esssup.marg   focsvm.marg            ms           seq
# FDR  0.099 (0.096)  0.104 (0.105) 0.111 (0.099) 0.105 (0.095) 0.014 (0.080) 0.077 (0.093)
# TPR  0.646 (0.320)  0.688 (0.316) 0.717 (0.201) 0.847 (0.186) 0.135 (0.071) 0.620 (0.122)

# T_projdepth.marg projdepth.marg projdepth_1d.marg projdepth_2d.marg T_hdepth.marg   hdepth.marg
# FDR    0.104 (0.112)  0.099 (0.104)     0.107 (0.123)     0.116 (0.125) 0.092 (0.111) 0.048 (0.101)
# TPR    0.557 (0.170)  0.627 (0.160)     0.517 (0.171)     0.482 (0.181) 0.574 (0.178) 0.138 (0.259)
# hdepth_1d.marg hdepth_2d.marg   esssup.marg   focsvm.marg            ms           seq
# FDR  0.096 (0.116)  0.101 (0.107) 0.105 (0.112) 0.104 (0.109) 0.028 (0.134) 0.099 (0.111)
# TPR  0.536 (0.171)  0.515 (0.226) 0.545 (0.229) 0.597 (0.215) 0.122 (0.092) 0.598 (0.127)


### Mixed training에 대한 outlier detection (100번 반복; test 데이터 사용 안함)
### - FOCSVM은 clean training이 있어야 해서 적용 불가
# [[1]]
# T_projdepth      T_hdepth        esssup        focsvm     projdepth  projdepth_1d  projdepth_2d
# FDR 0.000 (0.000) 0.000 (0.000) 0.004 (0.010)    NA (   NA) 0.000 (0.000) 0.782 (0.309) 0.734 (0.352)
# TPR 0.859 (0.076) 0.000 (0.000) 0.188 (0.058)    NA (   NA) 0.984 (0.015) 0.003 (0.004) 0.003 (0.004)
# hdepth     hdepth_1d     hdepth_2d            ms           seq
# FDR 0.000 (0.000) 0.000 (0.000) 0.000 (0.000) 0.024 (0.011) 0.000 (0.000)
# TPR 0.000 (0.000) 0.000 (0.000) 0.000 (0.000) 1.000 (0.000) 0.821 (0.033)
# 
# [[2]]
# T_projdepth      T_hdepth        esssup        focsvm     projdepth  projdepth_1d  projdepth_2d
# FDR 0.000 (0.000) 0.000 (0.000) 0.107 (0.038)    NA (   NA) 0.175 (0.359) 0.000 (0.000) 0.000 (0.000)
# TPR 1.000 (0.001) 0.998 (0.003) 0.287 (0.031)    NA (   NA) 0.002 (0.003) 1.000 (0.000) 1.000 (0.000)
# hdepth     hdepth_1d     hdepth_2d            ms           seq
# FDR 0.000 (0.000) 0.000 (0.000) 0.000 (0.000) 0.427 (0.211) 0.000 (0.000)
# TPR 0.000 (0.000) 0.999 (0.002) 1.000 (0.000) 0.024 (0.015) 0.000 (0.000)
# 
# [[3]]
# T_projdepth      T_hdepth        esssup        focsvm     projdepth  projdepth_1d  projdepth_2d
# FDR 0.000 (0.000) 0.000 (0.000) 0.105 (0.055)    NA (   NA) 0.000 (0.000) 0.000 (0.000) 0.000 (0.000)
# TPR 0.954 (0.023) 0.000 (0.000) 0.152 (0.030)    NA (   NA) 0.007 (0.007) 1.000 (0.001) 1.000 (0.000)
# hdepth     hdepth_1d     hdepth_2d            ms           seq
# FDR 0.000 (0.000) 0.000 (0.000) 0.000 (0.000) 0.025 (0.012) 0.000 (0.000)
# TPR 0.000 (0.000) 0.016 (0.011) 0.018 (0.014) 0.984 (0.009) 0.001 (0.003)
# 
# [[4]]
# T_projdepth      T_hdepth        esssup        focsvm     projdepth  projdepth_1d  projdepth_2d
# FDR 0.000 (0.000) 0.000 (0.000) 0.999 (0.005)    NA (   NA) 0.010 (0.100) 0.000 (0.000) 0.000 (0.000)
# TPR 1.000 (0.000) 1.000 (0.002) 0.000 (0.001)    NA (   NA) 0.000 (0.000) 1.000 (0.000) 1.000 (0.000)
# hdepth     hdepth_1d     hdepth_2d            ms           seq
# FDR 0.000 (0.000) 0.000 (0.000) 0.000 (0.000) 1.000 (0.000) 0.000 (0.000)
# TPR 0.000 (0.000) 1.000 (0.000) 1.000 (0.000) 0.000 (0.000) 0.000 (0.000)


### Mixed setting에 대한 conformal outlier detection (100번 반복)
### - FOCSVM은 clean training이 있어야 해서 적용 불가
# [[1]]
# T_projdepth.marg T_projdepth.simes T_projdepth.asymp   esssup.marg  esssup.simes  esssup.asymp
# FDR    0.076 (0.054)     0.019 (0.026)     0.019 (0.026) 0.043 (0.042) 0.000 (0.000) 0.001 (0.011)
# TPR    1.000 (0.000)     0.480 (0.502)     0.999 (0.005) 0.507 (0.269) 0.000 (0.000) 0.008 (0.082)
# projdepth.marg projdepth.simes projdepth.asymp projdepth_1d.marg projdepth_1d.simes
# FDR  0.084 (0.058)   0.018 (0.030)   0.018 (0.030)     0.000 (0.000)      0.000 (0.000)
# TPR  1.000 (0.000)   0.450 (0.500)   0.997 (0.014)     0.000 (0.000)      0.000 (0.000)
# projdepth_1d.asymp projdepth_2d.marg projdepth_2d.simes projdepth_2d.asymp            ms
# FDR      0.000 (0.000)     0.000 (0.000)      0.000 (0.000)      0.000 (0.000) 1.000 (0.000)
# TPR      0.000 (0.000)     0.000 (0.000)      0.000 (0.000)      0.000 (0.000) 0.000 (0.000)
# seq
# FDR 1.000 (0.000)
# TPR 0.000 (0.000)
# 
# [[2]]
# T_projdepth.marg T_projdepth.simes T_projdepth.asymp   esssup.marg  esssup.simes  esssup.asymp
# FDR    0.082 (0.057)     0.026 (0.030)     0.026 (0.030) 0.156 (0.141) 0.000 (0.000) 0.000 (0.000)
# TPR    1.000 (0.000)     0.660 (0.476)     1.000 (0.000) 0.204 (0.173) 0.000 (0.000) 0.000 (0.000)
# projdepth.marg projdepth.simes projdepth.asymp projdepth_1d.marg projdepth_1d.simes
# FDR  0.000 (0.000)   0.000 (0.000)   0.000 (0.000)     0.081 (0.053)      0.017 (0.020)
# TPR  0.000 (0.000)   0.000 (0.000)   0.000 (0.000)     1.000 (0.000)      0.520 (0.502)
# projdepth_1d.asymp projdepth_2d.marg projdepth_2d.simes projdepth_2d.asymp            ms
# FDR      0.017 (0.020)     0.085 (0.057)      0.018 (0.027)      0.018 (0.027) 1.000 (0.000)
# TPR      1.000 (0.000)     1.000 (0.000)      0.410 (0.494)      1.000 (0.000) 0.000 (0.000)
# seq
# FDR 0.330 (0.473)
# TPR 0.000 (0.000)
# 
# [[3]]
# T_projdepth.marg T_projdepth.simes T_projdepth.asymp   esssup.marg  esssup.simes  esssup.asymp
# FDR    0.081 (0.058)     0.024 (0.030)     0.024 (0.030) 0.153 (0.157) 0.000 (0.000) 0.000 (0.000)
# TPR    1.000 (0.000)     0.600 (0.492)     0.999 (0.005) 0.157 (0.152) 0.000 (0.000) 0.000 (0.000)
# projdepth.marg projdepth.simes projdepth.asymp projdepth_1d.marg projdepth_1d.simes
# FDR  0.085 (0.077)   0.000 (0.000)   0.000 (0.000)     0.082 (0.053)      0.017 (0.021)
# TPR  0.585 (0.201)   0.000 (0.000)   0.000 (0.000)     1.000 (0.000)      0.480 (0.502)
# projdepth_1d.asymp projdepth_2d.marg projdepth_2d.simes projdepth_2d.asymp            ms
# FDR      0.017 (0.021)     0.085 (0.057)      0.016 (0.025)      0.016 (0.025) 1.000 (0.000)
# TPR      0.940 (0.239)     1.000 (0.000)      0.390 (0.490)      0.978 (0.141) 0.000 (0.000)
# seq
# FDR 0.780 (0.416)
# TPR 0.000 (0.000)
# 
# [[4]]
# T_projdepth.marg T_projdepth.simes T_projdepth.asymp   esssup.marg  esssup.simes  esssup.asymp
# FDR    0.082 (0.057)     0.021 (0.025)     0.021 (0.025) 0.370 (0.485) 0.000 (0.000) 0.000 (0.000)
# TPR    1.000 (0.000)     0.560 (0.499)     1.000 (0.000) 0.000 (0.000) 0.000 (0.000) 0.000 (0.000)
# projdepth.marg projdepth.simes projdepth.asymp projdepth_1d.marg projdepth_1d.simes
# FDR  0.000 (0.000)   0.000 (0.000)   0.000 (0.000)     0.082 (0.059)      0.016 (0.023)
# TPR  0.000 (0.000)   0.000 (0.000)   0.000 (0.000)     1.000 (0.000)      0.470 (0.502)
# projdepth_1d.asymp projdepth_2d.marg projdepth_2d.simes projdepth_2d.asymp            ms
# FDR      0.016 (0.023)     0.087 (0.055)      0.021 (0.028)      0.021 (0.028) 1.000 (0.000)
# TPR      1.000 (0.000)     1.000 (0.000)      0.520 (0.502)      1.000 (0.000) 0.000 (0.000)
# seq
# FDR 0.070 (0.256)
# TPR 0.000 (0.000)



##############################
### 2025.04.01 ~ 
##############################
### ADHD-200 mixed setting (10번 반복)
# Raw
# T_projdepth.marg projdepth.marg   esssup.marg           seq
# FDR    0.073 (0.125)  0.083 (0.136) 0.117 (0.193) 0.134 (0.163)
# TPR    0.136 (0.198)  0.110 (0.159) 0.057 (0.083) 0.219 (0.068)
# 
# FFT-sm
# T_projdepth.marg projdepth.marg   esssup.marg           seq
# FDR    0.031 (0.068)  0.042 (0.100) 0.072 (0.118) 0.158 (0.139)
# TPR    0.227 (0.339)  0.128 (0.179) 0.072 (0.122) 0.200 (0.065)

### Training outlier detection
# T_projdepth   projdepth      esssup         seq 
# FDR  0.3883170   0.3762021   0.9074074   0.6228215 
# TPR  1.00        1.00        1.00        0.36 
#
# FFT-sm
# T_projdepth   projdepth      esssup         seq 
# FDR  0.2243590   0.2243590   0.9099099   0.5284737 
# TPR  1.00        1.00        1.00        0.71 


### 시뮬레이션 mixed setting
# [[1]]
# T_projdepth.marg T_projdepth.simes T_projdepth.asymp   esssup.marg  esssup.simes  esssup.asymp
# FDR    0.076 (0.054)     0.019 (0.026)     0.019 (0.026) 0.043 (0.042) 0.000 (0.000) 0.001 (0.011)
# TPR    1.000 (0.000)     0.480 (0.502)     0.999 (0.005) 0.507 (0.269) 0.000 (0.000) 0.008 (0.082)
# projdepth.marg projdepth.simes projdepth.asymp            ms           seq
# FDR  0.085 (0.057)   0.022 (0.032)   0.022 (0.032) 0.168 (0.051) 0.000 (0.000)
# TPR  1.000 (0.000)   0.520 (0.502)   1.000 (0.000) 1.000 (0.000) 0.765 (0.071)
# 
# [[2]]
# T_projdepth.marg T_projdepth.simes T_projdepth.asymp   esssup.marg  esssup.simes  esssup.asymp
# FDR    0.082 (0.057)     0.026 (0.030)     0.026 (0.030) 0.156 (0.141) 0.000 (0.000) 0.000 (0.000)
# TPR    1.000 (0.000)     0.660 (0.476)     1.000 (0.000) 0.204 (0.173) 0.000 (0.000) 0.000 (0.000)
# projdepth.marg projdepth.simes projdepth.asymp            ms           seq
# FDR  0.000 (0.000)   0.000 (0.000)   0.000 (0.000) 0.687 (0.149) 0.005 (0.050)
# TPR  0.000 (0.000)   0.000 (0.000)   0.000 (0.000) 0.068 (0.035) 0.009 (0.014)
# 
# [[3]]
# T_projdepth.marg T_projdepth.simes T_projdepth.asymp   esssup.marg  esssup.simes  esssup.asymp
# FDR    0.081 (0.058)     0.024 (0.030)     0.024 (0.030) 0.153 (0.157) 0.000 (0.000) 0.000 (0.000)
# TPR    1.000 (0.000)     0.600 (0.492)     0.999 (0.005) 0.157 (0.152) 0.000 (0.000) 0.000 (0.000)
# projdepth.marg projdepth.simes projdepth.asymp            ms           seq
# FDR  0.001 (0.007)   0.000 (0.000)   0.000 (0.000) 0.175 (0.049) 0.003 (0.033)
# TPR  0.023 (0.073)   0.000 (0.000)   0.000 (0.000) 0.935 (0.038) 0.030 (0.025)
# 
# [[4]]
# T_projdepth.marg T_projdepth.simes T_projdepth.asymp   esssup.marg  esssup.simes  esssup.asymp
# FDR    0.082 (0.057)     0.021 (0.025)     0.021 (0.025) 0.370 (0.485) 0.000 (0.000) 0.000 (0.000)
# TPR    1.000 (0.000)     0.560 (0.499)     1.000 (0.000) 0.000 (0.000) 0.000 (0.000) 0.000 (0.000)
# projdepth.marg projdepth.simes projdepth.asymp            ms           seq
# FDR  0.000 (0.000)   0.000 (0.000)   0.000 (0.000) 1.000 (0.000) 0.000 (0.000)
# TPR  0.000 (0.000)   0.000 (0.000)   0.000 (0.000) 0.000 (0.000) 0.002 (0.006)