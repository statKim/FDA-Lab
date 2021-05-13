#########################
### Combine AMI data
#########################
library(tidyverse)
library(data.table)

path <- "C:/Users/user/Google 드라이브/Lab/KHS/partiall_obs/real_data/AMI자료/house_2012_2014/"
flist <- list.files(path)
flist <- setdiff(flist, "desktop.ini")

for (name in flist) {
  print(name)
  
  tmp <- fread(paste0(path, name))
  
  # data has different variable names...
  if (name %in% c("y2014_01.csv","y2014_10.csv")) {
    tmp <- tmp[, -c("순번","고객명")]
  } else {
    tmp <- tmp[, -c("순번","계약전력")]
  }
  
  num_day_hour <- ncol(tmp) - 2   # 31 * 24
  num_day <- num_day_hour / 24   # number of day
  yyyymm <- as.numeric(tmp[1, 1])   # YYYYMM
  
  # transform time variable to "YYYYMMDDHH" format
  df <- data.frame(yymm = yyyymm,
                   day = rep(1:num_day, each = 24),
                   hour = rep(1:24, num_day)) %>%  
    mutate(
      time1 = paste0(yymm, "-", day, "-", hour),
      time2 = format(as.POSIXct(time1, format = "%Y%m-%d-%H"), 
                     "%Y%m%d%H"),
      time = as.numeric(time2)
    ) %>% 
    dplyr::select("time")
  
  colnames(tmp) <- c("yyyymm","ID",df$time)
  tmp <- tmp %>% 
    dplyr::select(-c("yyyymm")) %>% 
    gather(key = Time, value = Consumption, -ID) %>% 
    tibble
  
  if (name == flist[1]) {
    AMI <- tmp
  } else {
    AMI <- rbind(AMI, tmp)
  }
  
  rm(tmp)
}
save(list = c("AMI"),
     file = "AMI_data/AMI.RData")
