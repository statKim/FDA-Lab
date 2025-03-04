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

### fMRI (30번) - 기존 outlier 방법
#     T_projdepth.marg projdepth.marg projdepth_1d.marg projdepth_2d.marg T_hdepth.marg   hdepth.marg
# FDR    0.115 (0.120)  0.099 (0.110)     0.117 (0.120)     0.114 (0.130) 0.125 (0.109) 0.113 (0.108)
# TPR    0.597 (0.364)  0.719 (0.268)     0.583 (0.343)     0.472 (0.336) 0.569 (0.313) 0.647 (0.293)
#     hdepth_1d.marg hdepth_2d.marg   esssup.marg
# FDR  0.139 (0.122)  0.131 (0.128) 0.116 (0.103)
# TPR  0.553 (0.321)  0.408 (0.281) 0.683 (0.287)

### fMRI (30번) - raw, 1st, 2nd derivative를 depth 기준으로 20개
#     T_projdepth.marg projdepth.marg projdepth_1d.marg projdepth_2d.marg T_hdepth.marg   hdepth.marg
# FDR    0.094 (0.094)  0.093 (0.092)     0.100 (0.093)     0.102 (0.109) 0.105 (0.096) 0.091 (0.095)
# TPR    0.601 (0.349)  0.710 (0.319)     0.558 (0.362)     0.435 (0.332) 0.607 (0.383) 0.714 (0.338)
#     hdepth_1d.marg hdepth_2d.marg   esssup.marg
# FDR  0.108 (0.112)  0.111 (0.119) 0.108 (0.092)
# TPR  0.530 (0.369)  0.402 (0.326) 0.646 (0.260)
### - score 계산시 max로 함
#     T_projdepth.marg T_hdepth.marg
# FDR    0.091 (0.088) 0.115 (0.089)
# TPR    0.648 (0.312) 0.642 (0.348)
#         ms           seq-c("O","D1","D2")
# FDR  0.033 (0.183) 0.056 (0.094)
# TPR  0.066 (0.064) 0.372 (0.139)


### fMRI (100번) - raw, 1st, 2nd derivative를 depth 기준으로 20개
# T_projdepth.marg projdepth.marg projdepth_1d.marg projdepth_2d.marg T_hdepth.marg   hdepth.marg hdepth_1d.marg
# FDR    0.105 (0.091)  0.104 (0.087)     0.108 (0.096)     0.114 (0.111) 0.111 (0.098) 0.101 (0.097)  0.103 (0.108)
# TPR    0.654 (0.323)  0.785 (0.241)     0.608 (0.337)     0.462 (0.320) 0.652 (0.330) 0.744 (0.305)  0.512 (0.353)
# hdepth_2d.marg   esssup.marg            ms           seq
# FDR  0.116 (0.120) 0.111 (0.101) 0.033 (0.159) 0.055 (0.088)
# TPR  0.431 (0.299) 0.622 (0.275) 0.071 (0.072) 0.390 (0.135)


### fMRI (30번) - align 새로 한 것 - raw, 1st, 2nd derivative를 depth 기준으로 20개
# T_projdepth.marg projdepth.marg projdepth_1d.marg projdepth_2d.marg T_hdepth.marg   hdepth.marg hdepth_1d.marg
# FDR    0.094 (0.127)  0.101 (0.123)     0.107 (0.138)     0.094 (0.135) 0.102 (0.131) 0.103 (0.137)  0.105 (0.134)
# TPR    0.347 (0.353)  0.463 (0.379)     0.327 (0.345)     0.277 (0.294) 0.387 (0.367) 0.427 (0.384)  0.367 (0.345)
# hdepth_2d.marg   esssup.marg            ms           seq
# FDR  0.085 (0.133) 0.065 (0.111) 0.172 (0.351) 0.182 (0.207)
# TPR  0.237 (0.305) 0.257 (0.385) 0.073 (0.069) 0.300 (0.136)




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


### 시뮬레이션 (30번) - average (ccv도 포함한 결과; marginal은 위와 동일)
# [[1]]
#     T_projdepth.marg T_projdepth.simes T_projdepth.asymp T_hdepth.marg T_hdepth.simes T_hdepth.asymp
# FDR    0.081 (0.053)     0.033 (0.033)     0.033 (0.033) 0.000 (0.000)  0.000 (0.000)  0.000 (0.000)
# TPR    1.000 (0.000)     1.000 (0.000)     1.000 (0.000) 0.000 (0.000)  0.000 (0.000)  0.000 (0.000)
#       esssup.marg  esssup.simes  esssup.asymp projdepth.marg projdepth.simes projdepth.asymp projdepth_1d.marg
# FDR 0.082 (0.045) 0.000 (0.000) 0.002 (0.009)  0.097 (0.068)   0.018 (0.035)   0.017 (0.031)     0.000 (0.000)
# TPR 0.729 (0.216) 0.000 (0.000) 0.026 (0.142)  1.000 (0.000)   1.000 (0.000)   1.000 (0.000)     0.000 (0.000)
# projdepth_1d.simes projdepth_1d.asymp projdepth_2d.marg projdepth_2d.simes projdepth_2d.asymp   hdepth.marg
# FDR      0.000 (0.000)      0.000 (0.000)     0.027 (0.146)      0.000 (0.000)      0.000 (0.000) 0.000 (0.000)
# TPR      0.000 (0.000)      0.000 (0.000)     0.001 (0.007)      0.000 (0.000)      0.000 (0.000) 0.000 (0.000)
# hdepth.simes  hdepth.asymp hdepth_1d.marg hdepth_1d.simes hdepth_1d.asymp hdepth_2d.marg hdepth_2d.simes
# FDR 0.000 (0.000) 0.000 (0.000)  0.000 (0.000)   0.000 (0.000)   0.000 (0.000)  0.000 (0.000)   0.000 (0.000)
# TPR 0.000 (0.000) 0.000 (0.000)  0.000 (0.000)   0.000 (0.000)   0.000 (0.000)  0.000 (0.000)   0.000 (0.000)
# hdepth_2d.asymp            ms           seq
# FDR   0.000 (0.000)    NA (   NA)    NA (   NA)
# TPR   0.000 (0.000)    NA (   NA)    NA (   NA)
# 
# [[2]]
# T_projdepth.marg T_projdepth.simes T_projdepth.asymp T_hdepth.marg T_hdepth.simes T_hdepth.asymp
# FDR    0.069 (0.047)     0.025 (0.032)     0.023 (0.024) 0.053 (0.087)  0.000 (0.000)  0.000 (0.000)
# TPR    1.000 (0.000)     1.000 (0.000)     1.000 (0.000) 0.112 (0.148)  0.000 (0.000)  0.000 (0.000)
# esssup.marg  esssup.simes  esssup.asymp projdepth.marg projdepth.simes projdepth.asymp projdepth_1d.marg
# FDR 0.045 (0.105) 0.000 (0.000) 0.000 (0.000)  0.000 (0.000)   0.000 (0.000)   0.000 (0.000)     0.080 (0.060)
# TPR 0.063 (0.125) 0.000 (0.000) 0.000 (0.000)  0.000 (0.000)   0.000 (0.000)   0.000 (0.000)     1.000 (0.000)
# projdepth_1d.simes projdepth_1d.asymp projdepth_2d.marg projdepth_2d.simes projdepth_2d.asymp   hdepth.marg
# FDR      0.015 (0.025)      0.015 (0.025)     0.081 (0.057)      0.018 (0.025)      0.018 (0.025) 0.000 (0.000)
# TPR      1.000 (0.000)      1.000 (0.000)     1.000 (0.000)      1.000 (0.000)      1.000 (0.000) 0.000 (0.000)
# hdepth.simes  hdepth.asymp hdepth_1d.marg hdepth_1d.simes hdepth_1d.asymp hdepth_2d.marg hdepth_2d.simes
# FDR 0.000 (0.000) 0.000 (0.000)  0.034 (0.088)   0.000 (0.000)   0.000 (0.000)  0.055 (0.081)   0.000 (0.000)
# TPR 0.000 (0.000) 0.000 (0.000)  0.036 (0.083)   0.000 (0.000)   0.000 (0.000)  0.212 (0.210)   0.000 (0.000)
# hdepth_2d.asymp            ms           seq
# FDR   0.000 (0.000)    NA (   NA)    NA (   NA)
# TPR   0.000 (0.000)    NA (   NA)    NA (   NA)
# 
# [[3]]
# T_projdepth.marg T_projdepth.simes T_projdepth.asymp T_hdepth.marg T_hdepth.simes T_hdepth.asymp
# FDR    0.069 (0.047)     0.025 (0.032)     0.023 (0.024) 0.037 (0.103)  0.000 (0.000)  0.000 (0.000)
# TPR    1.000 (0.000)     0.999 (0.004)     0.999 (0.004) 0.022 (0.057)  0.000 (0.000)  0.000 (0.000)
# esssup.marg  esssup.simes  esssup.asymp projdepth.marg projdepth.simes projdepth.asymp projdepth_1d.marg
# FDR 0.039 (0.090) 0.000 (0.000) 0.000 (0.000)  0.086 (0.052)   0.000 (0.000)   0.004 (0.022)     0.080 (0.060)
# TPR 0.057 (0.120) 0.000 (0.000) 0.000 (0.000)  0.618 (0.114)   0.000 (0.000)   0.025 (0.135)     1.000 (0.000)
# projdepth_1d.simes projdepth_1d.asymp projdepth_2d.marg projdepth_2d.simes projdepth_2d.asymp   hdepth.marg
# FDR      0.015 (0.025)      0.015 (0.025)     0.081 (0.057)      0.018 (0.025)      0.018 (0.025) 0.000 (0.000)
# TPR      1.000 (0.000)      1.000 (0.000)     1.000 (0.000)      1.000 (0.000)      1.000 (0.000) 0.000 (0.000)
# hdepth.simes  hdepth.asymp hdepth_1d.marg hdepth_1d.simes hdepth_1d.asymp hdepth_2d.marg hdepth_2d.simes
# FDR 0.000 (0.000) 0.000 (0.000)  0.013 (0.073)   0.000 (0.000)   0.000 (0.000)  0.023 (0.090)   0.000 (0.000)
# TPR 0.000 (0.000) 0.000 (0.000)  0.004 (0.022)   0.000 (0.000)   0.000 (0.000)  0.015 (0.060)   0.000 (0.000)
# hdepth_2d.asymp            ms           seq
# FDR   0.000 (0.000)    NA (   NA)    NA (   NA)
# TPR   0.000 (0.000)    NA (   NA)    NA (   NA)
# 
# [[4]]
# T_projdepth.marg T_projdepth.simes T_projdepth.asymp T_hdepth.marg T_hdepth.simes T_hdepth.asymp
# FDR    0.076 (0.050)     0.023 (0.028)     0.023 (0.028) 0.000 (0.000)  0.000 (0.000)  0.000 (0.000)
# TPR    1.000 (0.000)     1.000 (0.000)     1.000 (0.000) 0.000 (0.000)  0.000 (0.000)  0.000 (0.000)
# esssup.marg  esssup.simes  esssup.asymp projdepth.marg projdepth.simes projdepth.asymp projdepth_1d.marg
# FDR 0.000 (0.000) 0.000 (0.000) 0.000 (0.000)  0.000 (0.000)   0.000 (0.000)   0.000 (0.000)     0.078 (0.049)
# TPR 0.000 (0.000) 0.000 (0.000) 0.000 (0.000)  0.000 (0.000)   0.000 (0.000)   0.000 (0.000)     1.000 (0.000)
# projdepth_1d.simes projdepth_1d.asymp projdepth_2d.marg projdepth_2d.simes projdepth_2d.asymp   hdepth.marg
# FDR      0.008 (0.017)      0.008 (0.017)     0.083 (0.058)      0.015 (0.027)      0.015 (0.027) 0.000 (0.000)
# TPR      1.000 (0.000)      1.000 (0.000)     1.000 (0.000)      1.000 (0.000)      1.000 (0.000) 0.000 (0.000)
# hdepth.simes  hdepth.asymp hdepth_1d.marg hdepth_1d.simes hdepth_1d.asymp hdepth_2d.marg hdepth_2d.simes
# FDR 0.000 (0.000) 0.000 (0.000)  0.062 (0.046)   0.016 (0.026)   0.016 (0.026)  0.062 (0.056)   0.012 (0.018)
# TPR 0.000 (0.000) 0.000 (0.000)  1.000 (0.000)   0.967 (0.183)   0.967 (0.183)  1.000 (0.000)   0.900 (0.305)
# hdepth_2d.asymp            ms           seq
# FDR   0.012 (0.018)    NA (   NA)    NA (   NA)
# TPR   0.900 (0.305)    NA (   NA)    NA (   NA)



##############################
### 2025.03.04 ~
##############################

