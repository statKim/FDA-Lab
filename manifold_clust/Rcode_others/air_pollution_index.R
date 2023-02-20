library(tidyverse)

# ### 측정소 위치 데이터 - 어디서 받은지 모르겠음...
# # 여기서 위경도 정보 필요 없음 (air pollutant variable이 spherical data가 될 것임)
# full_data <- readRDS("/Users/hyunsung/GoogleDrive/Lab/KHS/partiall_obs/real_data/KoreaPM10.RDS")
# class(full_data)
# full_data$place
# dim(full_data)
# data_loc <- tibble(full_data$place[, c(3, 6, 7)])
# colnames(data_loc) <- c("ID","LON","LAT")
# head(data_loc)

######################################################
### Daily data for 2021.01 ~ 2021.10
######################################################

### Air pollution data
# 최종확정 측정자료 조회 -> 확정자료 다운로드
# https://www.airkorea.or.kr/web/last_amb_hour_data?pMENU_NO=123
path <- "/Users/hyunsung/Desktop/air_pollution_data/2021/"
file_list <- list.files(path)
file_list
paste0(path, file_list[4])
data <- readxl::read_excel("/Users/hyunsung/Desktop/2021/2021년 4월.xlsx")

data <- rbind(
    readxl::read_excel("/Users/hyunsung/Desktop/2021/2021년 3월.xlsx"),
    readxl::read_excel("/Users/hyunsung/Desktop/2021/2021년 4월.xlsx"),
    readxl::read_excel("/Users/hyunsung/Desktop/2021/2021년 5월.xlsx")
)
head(data)
tail(data)
dim(data)

# Remove redundant features
data <- data[, -c(2, 4, 12)]
data
# Rename features with Korean to English
colnames(data)[1:3] <- c("Region","ID","Time")
data

data %>% 
    distinct(Region, ID) %>% 
    mutate(Si = substr(Region, 1, 2)) %>% 
    select(Si) %>% 
    table() %>% 
    t() %>% t()
substr(data$Region, 1, 2)

# Combine location data
data2 <- data %>% 
    # filter(substr(Region, 1, 2) == "서울") %>%
    mutate(
        ID = as.integer(ID),
        Day = strptime(Time, "%Y %m %d")
    ) %>% 
    group_by(ID, Day) %>% 
    summarise(
        SO2 = mean(SO2, na.rm = T),
        CO = mean(CO, na.rm = T),
        O3 = mean(O3, na.rm = T),
        NO2 = mean(NO2, na.rm = T),
        PM10 = mean(PM10, na.rm = T),
        PM25 = mean(PM25, na.rm = T)
    ) %>% 
    mutate(Time = (0:(n()-1)) / (n()-1)) %>% 
    ungroup()
data2
tail(data2)

# Interpolation to replace NA values
# - There exists All NA values per days
# - `rule = 2` means that boundary extrapolation is performed
api_data <- data2 %>% 
    group_by(ID) %>% 
    filter(sum(is.na(SO2)) != n() &
               sum(is.na(CO)) != n() &
               sum(is.na(O3)) != n() &
               sum(is.na(NO2)) != n() &
               sum(is.na(PM10)) != n() &
               sum(is.na(PM25)) != n()) %>% 
    mutate(
        SO2 = approx(Time, SO2, xout = Time, rule = 2)$y,
        CO = approx(Time, CO, xout = Time, rule = 2)$y,
        O3 = approx(Time, O3, xout = Time, rule = 2)$y,
        NO2 = approx(Time, NO2, xout = Time, rule = 2)$y,
        PM10 = approx(Time, PM10, xout = Time, rule = 2)$y,
        PM25 = approx(Time, PM25, xout = Time, rule = 2)$y
    )
# which(is.na(data3$CO))
# data3[c(1097:1107), ]
# sum(is.na(data3))
length(unique(api_data$ID))



# # Make API (Air Pollution Index) and transform to compositional data
# # Using 5 variables:
# # - SO2, CO, O3, NO2, PM10
# api_data <- data3 %>% 
#     select(-PM25) %>% 
#     mutate(
#         sum_pollution = SO2 + CO + O3 + NO2 + PM10,
#         SO2 = SO2 / sum_pollution,
#         CO = CO / sum_pollution,
#         O3 = O3 / sum_pollution,
#         NO2 = NO2 / sum_pollution,
#         PM10 = PM10 / sum_pollution
#     ) %>% 
#     select(-sum_pollution)
# api_data
# # Transform compositional data to spherical data
# api_data[, 3:7] <- sqrt(api_data[, 3:7])
# rowSums(api_data[3:7]^2) %>% table()


# Make API (Air Pollution Index) and transform to compositional data
# Using 5 variables:
# - SO2, CO, O3, NO2, PM10
### Make data to input format
id <- unique(api_data$ID)
Ly <- lapply(id, function(i) {
    as.matrix( api_data[api_data$ID == i, 3:7] )
})
Lt <- lapply(id, function(i) {
    as.numeric( unlist( api_data[api_data$ID == i, "Time"] ) )
})


# ### Pre-smoothing for regular grids using local linear smoother
# n <- length(id)
# Ly <- lapply(1:n, function(i) {
#     y <- Ly[[i]]
#     t <- Lt[[i]]
#     bw <- min(diff(t))   # very small bandwidth
#     
#     # kernel smoothing with 51 regular grids
#     y <- apply(y, 2, function(col) {
#         stats::ksmooth(x = t, 
#                        y = col, 
#                        kernel = "normal", 
#                        bandwidth = bw,
#                        n.points = 51)$y
#     })
#     
#     # make spherical data
#     apply(y, 1, function(row){ sqrt(row / sum(row)) })
# })
# Lt <- rep(list(seq(0, 1, length.out = 51)), n)
# Ly[[1]]
# Lt[[1]]
# apply(Ly[[1]], 2, function(x){ sum(x^2) })   # check that it is the unit sphere

### Make compositional data and transform to spherical data
n <- length(id)
Ly <- lapply(1:n, function(i) {
    y <- t( Ly[[i]] )
    # make spherical data
    apply(y, 2, function(row){ sqrt(row / sum(row)) })
})
apply(Ly[[1]], 2, function(x){ sum(x^2) })   # check that it is the unit sphere


### Make "Others" group which combines SO2 + CO + NO2
Ly <- lapply(Ly, function(y) {
    rbind(
        y[c(5, 3), ],
        Others = sqrt(colSums(y[c(1,2,4), ]^2))
    )
})
Ly[[1]]
apply(Ly[[1]], 2, function(x){ sum(x^2) })   # check that it is the unit sphere




######################################################
### Monthly data for 2017.01 ~ 2020.12
######################################################
### Air pollution data
# 최종확정 측정자료 조회 -> 확정자료 다운로드
# https://www.airkorea.or.kr/web/last_amb_hour_data?pMENU_NO=123
file_list <- lapply(paste0("/Users/hyunsung/Desktop/air_pollution_data/", 2017:2020, "/"), function(dir) {
    paste0(dir, list.files(dir))
}) %>% 
    unlist()
file_list
# path <- "/Users/hyunsung/Desktop/air_pollution_data/"
# file_list <- list.files(path)
# file_list
# paste0(path, file_list)

for (file_name in file_list) {
    # file_name <- paste0(path, fname)
    data_month <- readxl::read_excel(file_name)
    
    # Remove redundant features
    cols <- c("지역","측정소코드","측정일시","SO2","CO","O3","NO2","PM10","PM25")
    # data_month <- data_month[, -c(2, 4, 12)]
    data_month <- data_month[, cols]
    data_month
    
    # Rename features with Korean to English
    colnames(data_month)[1:3] <- c("Region","ID","Time")
    data_month
    
    # Combine location data
    data_month <- data_month %>% 
        # filter(substr(Region, 1, 2) == "서울") %>%
        mutate(
            ID = as.integer(ID),
            Month = as.integer(substr(Time, 1, 6))
        ) %>% 
        group_by(ID, Month) %>% 
        summarise(
            SO2 = mean(SO2, na.rm = T),
            CO = mean(CO, na.rm = T),
            O3 = mean(O3, na.rm = T),
            NO2 = mean(NO2, na.rm = T),
            PM10 = mean(PM10, na.rm = T),
            PM25 = mean(PM25, na.rm = T)
        ) %>% 
        ungroup()
    data_month
    tail(data_month)
    
    print(file_name)
    print(data_month)
    
    # data_month[which(is.na(data_month), arr.ind = T)[, 1], ]
    
    # if (fname == file_list[1]) {
    if (file_name == file_list[1]) {
        data <- data_month
    } else {
        data <- rbind(data, data_month)
    }
}
head(data)
tail(data)
dim(data)

data[which(is.na(data), arr.ind = T)[, 1], ] %>% 
    mutate(sum_NA = is.na(SO2) + is.na(CO) + is.na(O3) + is.na(NO2) + is.na(PM10)) %>% 
    filter(sum_NA == 5)

table(data$Month)   # 202110까지만 데이터 존재...
length(unique(data$Month))   # 58 months (2017.01 ~ 2021.10)

data %>% 
    group_by(ID) %>% 
    summarise(n = n()) %>% 
    filter(n >= 58)


### Data Preprocessing
n_month <- length(unique(data$Month))   # 58
# gr <- seq(0, 1, length.out = n_month)
api_data <- data %>% 
    arrange(ID, Month) %>% 
    group_by(ID) %>% 
    filter(
        # At least 58 months are observed
        n() == n_month
    ) %>% 
    mutate(
        # Timepoints between 0 and 1
        Time = (order(Month) - 1) / (n_month-1)
    ) %>% 
    filter(
        # # At most 2 observations are NA
        # sum(is.na(SO2)) <= 2 &
        #     sum(is.na(CO)) <= 2 &
        #     sum(is.na(O3)) <= 2 &
        #     sum(is.na(NO2)) <= 2 &
        #     sum(is.na(PM10)) <= 2 &
        # Contains boundary timepoints
        first(Time) == 0 & 
            last(Time) == 1 &
            # Contains boundary values
            !is.na(first(SO2)) & !is.na(last(SO2)) &
            !is.na(first(CO)) & !is.na(last(CO)) &
            !is.na(first(O3)) & !is.na(last(O3)) &
            !is.na(first(NO2)) & !is.null(last(NO2)) &
            !is.na(first(PM10)) & !is.na(last(PM10))
    ) %>% 
    ungroup()
api_data
length(unique(api_data$ID))



# Make API (Air Pollution Index) and transform to compositional data
# Using 5 variables:
# - SO2, CO, O3, NO2, PM10
### Make data to input format
id <- unique(api_data$ID)
Ly <- lapply(id, function(i) {
    as.matrix( api_data[api_data$ID == i, 3:7] )
})
Lt <- lapply(id, function(i) {
    as.numeric( unlist( api_data[api_data$ID == i, "Time"] ) )
})

Ly[sapply(Ly, function(x){sum(is.na(x))}) > 0]


### Interpolation and make compositional data, transform to spherical data
### Pre-smoothing for regular grids using local linear smoother
n <- length(id)
Ly <- lapply(1:n, function(i) {
    y <- Ly[[i]]
    t <- Lt[[i]]
    bw <- max(diff(t))   # very small bandwidth
    
    # kernel smoothing with 51 regular grids
    y <- apply(y, 2, function(col) {
        col <- approx(t, col, xout = t)$y
        # stats::ksmooth(x = t,
        #                y = col,
        #                kernel = "normal",
        #                bandwidth = bw,
        #                n.points = 51)$y
    })
    
    # Transform compositional data into spherical data
    apply(y, 1, function(row){ sqrt(row / sum(row)) })
})
Lt <- rep(list(seq(0, 1, length.out = n_month)), n)
# Lt <- rep(list(seq(0, 1, length.out = 51)), n)
Ly[[1]]
Lt[[1]]
apply(Ly[[1]], 2, function(x){ sum(x^2) })   # check that it is the unit sphere

# ### Make compositional data and transform to spherical data
# n <- length(id)
# Ly <- lapply(1:n, function(i) {
#     y <- t( Ly[[i]] )
#     # make spherical data
#     apply(y, 2, function(row){ sqrt(row / sum(row)) })
# })
# apply(Ly[[1]], 2, function(x){ sum(x^2) })   # check that it is the unit sphere


### Make "Others" group which combines SO2 + CO + NO2
Ly <- lapply(Ly, function(y) {
    rbind(
        y[c(5, 3), ],
        Others = sqrt(colSums(y[c(1,2,4), ]^2))
    )
})
Ly[[1]]
apply(Ly[[1]], 2, function(x){ sum(x^2) })   # check that it is the unit sphere


sapply(Ly, ncol)
sapply(Ly, nrow)

sapply(Ly, function(x){ sum(is.na(x)) })
sapply(Lt, function(x){ sum(is.na(x)) })

Ly <- lapply(Ly, function(y){
    y[,  1:30]
})
Lt <- rep(list(seq(0, 1, length.out = 30)), n)



### Ternary plot
library(Ternary)
TernaryPlot("PM10", "O3", "Others")
AddToTernary(points, lapply(1:length(Lt[[1]]), function(i){ Ly[[1]][, i]^2 }))
AddToTernary(points, lapply(1:length(Lt[[1]]), function(i){ Ly[[2]][, i]^2 }), col = 2)

as.data.frame( t( Ly[[1]] ) ) %>% 
    mutate(Time = Lt[[1]]) %>% 
    gather("Source","Value", -Time) %>% 
    ggplot(aes(Time, Value, group = Source)) +
    geom_line() +
    theme_bw() +
    facet_grid(Source ~ .)

par(mfrow = c(2, 2))
for (i in 1:3) {
    plot(Lt[[1]], Ly[[1]][i, ], type = "l")
    for (j in 2:50) {
        lines(Lt[[j]], Ly[[j]][i, ], col = j)
    }
}

# compositional barplot
barplot(Ly[[1]]^2, col = c("gray","red","blue"), ylim = c(0.9, 1))



######################################################
### Clustering
######################################################
# devtools::install_github('CrossD/RFPCA')
# devtools::install_url("https://cran.r-project.org/src/contrib/Archive/Funclustering/Funclustering_1.0.2.tar.gz")
library(RFPCA)    # RFPCA and MFPCA
library(mclust)   # clustering measure
library(Funclustering)   # funclust (Currently, it is not supported by cran.)
library(funHDDC)   # funHDDC
library(gmfd)   # gmfd
source("functions.R")

### Model parameters
seed <- 1000
k <- 2    # number of clusters
num.pc.method <- "FVE"   # using FVE thresholds
# num.pc.method <- 2     # fixed number
if (num.pc.method == "FVE") {
    FVEthresholdSW <- 0.90
    FVEthresholdCS <- 0.70
    maxK <- Inf
} else if (as.integer(num.pc.method)) {
    FVEthresholdSW <- 1
    FVEthresholdCS <- 1
    maxK <- num.pc.method
}
init_method <- "kmeans"
# init_method <- "cmeans"


### kCFC with Riemannian metric
system.time({
    fit.kCFC.Riemann <- kCRFC(y = Ly, 
                              t = Lt, 
                              k = k,
                              initMethod = init_method,
                              kSeed = seed, 
                              maxIter = 125, 
                              optnsSW = list(mfdName = "Sphere",
                                             FVEthreshold = FVEthresholdSW,
                                             maxK = maxK),
                              optnsCS = list(mfdName = "Sphere",
                                             FVEthreshold = FVEthresholdCS,
                                             maxK = maxK))
})
clust.kCFC.Riemann <- fit.kCFC.Riemann$cluster   # clustering index
clust.kmeans.Riemann <- fit.kCFC.Riemann$clustConf0   # initial cluster
# fit.kCFC.Riemann$clustConf0   # initial clustering index from k-means


### kCFC with Euclidean metric (multivariate FPCA)
system.time({
    fit.kCFC.L2 <- kCRFC(y = Ly, 
                         t = Lt, 
                         k = k,
                         initMethod = init_method,
                         kSeed = seed, 
                         maxIter = 125, 
                         optnsSW = list(mfdName = "Euclidean",
                                        FVEthreshold = FVEthresholdSW,
                                        maxK = maxK),
                         optnsCS = list(mfdName = "Euclidean",
                                        FVEthreshold = FVEthresholdCS,
                                        maxK = maxK))
})
clust.kCFC.L2 <- fit.kCFC.L2$cluster   # clustering index
clust.kmeans.L2 <- fit.kCFC.L2$clustConf0   # initial cluster




table(clust.kCFC.Riemann, clust.kCFC.L2)
table(clust.kCFC.Riemann, clust.kmeans.Riemann)
table(clust.kCFC.Riemann, clust.kmeans.L2)


par(mfrow = c(4, 3))
col <- as.integer(clust.kCFC.Riemann) + 1
for (i in 1:3) {
    ylim <- range(sapply(Ly, function(y){ y[i, ]^2 }))
    plot(Lt[[1]], Ly[[1]][i, ]^2, type = "l", col = col[1], ylim = ylim)
    for (j in 2:length(Ly)) {
        lines(Lt[[j]], Ly[[j]][i, ]^2, col = col[j])
    }
}
col <- as.integer(clust.kCFC.L2) + 1
for (i in 1:3) {
    ylim <- range(sapply(Ly, function(y){ y[i, ]^2 }))
    plot(Lt[[1]], Ly[[1]][i, ]^2, type = "l", col = col[1], ylim = ylim)
    for (j in 2:length(Ly)) {
        lines(Lt[[j]], Ly[[j]][i, ]^2, col = col[j])
    }
}
col <- as.integer(clust.kmeans.Riemann) + 1
for (i in 1:3) {
    ylim <- range(sapply(Ly, function(y){ y[i, ]^2 }))
    plot(Lt[[1]], Ly[[1]][i, ]^2, type = "l", col = col[1], ylim = ylim)
    for (j in 2:length(Ly)) {
        lines(Lt[[j]], Ly[[j]][i, ]^2, col = col[j])
    }
}
col <- as.integer(clust.kmeans.L2) + 1
for (i in 1:3) {
    ylim <- range(sapply(Ly, function(y){ y[i, ]^2 }))
    plot(Lt[[1]], Ly[[1]][i, ]^2, type = "l", col = col[1], ylim = ylim)
    for (j in 2:length(Ly)) {
        lines(Lt[[j]], Ly[[j]][i, ]^2, col = col[j])
    }
}

par(mfrow = c(4, 3))
col <- as.integer(clust.kCFC.Riemann) + 1
for (i in 1:3) {
    matplot(cbind(rowMeans(sapply(Ly[col == 2], function(y){ y[i, ]^2 })),
                  rowMeans(sapply(Ly[col == 3], function(y){ y[i, ]^2 }))),
            type = "l", lty = 1)
}
col <- as.integer(clust.kCFC.L2) + 1
for (i in 1:3) {
    matplot(cbind(rowMeans(sapply(Ly[col == 2], function(y){ y[i, ]^2 })),
                  rowMeans(sapply(Ly[col == 3], function(y){ y[i, ]^2 }))),
            type = "l", lty = 1)
}
col <- as.integer(clust.kmeans.Riemann) + 1
for (i in 1:3) {
    matplot(cbind(rowMeans(sapply(Ly[col == 2], function(y){ y[i, ]^2 })),
                  rowMeans(sapply(Ly[col == 3], function(y){ y[i, ]^2 }))),
            type = "l", lty = 1)
}
col <- as.integer(clust.kmeans.L2) + 1
for (i in 1:3) {
    matplot(cbind(rowMeans(sapply(Ly[col == 2], function(y){ y[i, ]^2 })),
                  rowMeans(sapply(Ly[col == 3], function(y){ y[i, ]^2 }))),
            type = "l", lty = 1)
}






