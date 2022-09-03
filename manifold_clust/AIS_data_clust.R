### Ship Trajectories - AIS data
### - https://marinecadastre.gov/ais/
### Variable Description: 
### - https://coast.noaa.gov/data/marinecadastre/ais/2018DataDictionary.png

data <- read_csv("~/Desktop/2019/AIS_2019_01_01.csv")
dim(data)
data
tail(data)

unique(data$MMSI) %>% length

data %>% 
    filter(MMSI == 367528080)

data$TransceiverClass %>% table

data %>% 
    filter(!is.na(CallSign))

data$BaseDateTime[1:10]
class(data$BaseDateTime[1])

as.POSIXct(data$BaseDateTime[1:10], "GMT")
as.character(data$BaseDateTime[1:10])

unique(data$BaseDateTime)

# 시간별 time 뽑기 => 이거로 시간별 데이터로 뽑아보자!!
as.POSIXct(strptime(as.character(data$BaseDateTime[c(1, 10000, 1000000, 100000)]), "%Y-%m-%d %H"))

# Data per hours
data2 <- data %>% 
    # filter(MMSI == 367528080) %>% 
    arrange(BaseDateTime) %>% 
    mutate(
        Hour = as.POSIXct(strptime(as.character(BaseDateTime), "%Y-%m-%d %H"))
    ) %>% 
    # dplyr::select(MMSI, Hour, LAT, LON) %>% 
    group_by(MMSI, Hour) %>% 
    filter(row_number() == 1) %>% 
    ungroup()


range(data2$LAT)

world <- ne_countries(scale = "medium", returnclass = "sf")
map_bg <- ggplot(data = world) +
    geom_sf() +
    coord_sf(xlim = range(data2$LON) + c(-10, 10),
             ylim = range(data2$LAT) + c(-10, 10),
             expand = FALSE)
map_bg +
    geom_path(
        data = data2,
        aes(
            x = LON,
            y = LAT,
            # group = MMSI,
            color = factor(MMSI)
        ),
        size = 0.3
    ) +
    theme_bw() +
    theme(legend.position = "none")


df <- data2 %>% 
    group_by(MMSI) %>% 
    arrange(MMSI) %>% 
    mutate(dif = abs(LON - lag(LON, 1, LON[1])) + abs(LAT - lag(LAT, 1, LAT[1]))) %>% 
    ungroup() %>% 
    filter(dif > 0.5) %>% 
    filter(VesselType == 70 & Length > 260 & Width > 32)
df

summary(df$Length, na.rm = T)
summary(df$Width, na.rm = T)


world <- ne_countries(scale = "medium", returnclass = "sf")
map_bg <- ggplot(data = world) +
    geom_sf() +
    coord_sf(xlim = range(df$LON) + c(-10, 10),
             ylim = range(df$LAT) + c(-10, 10),
             expand = FALSE)
map_bg +
    geom_path(
        data = df,
        aes(
            x = LON,
            y = LAT,
            # group = MMSI,
            color = factor(MMSI)
        ),
        size = 0.3
    ) +
    theme_bw() +
    theme(legend.position = "none")


data2 %>% 
    group_by(MMSI) %>% 
    filter(row_number() == 1) %>% 
    dplyr::select(VesselType) %>% 
    unlist() %>% 
    table






#####################################################################
### Monthly data
### - 이 데이터의 경우, 미국 근처에서의 항로만 나옴...jgu
#####################################################################
path <- "~/Desktop/2019/"
file_list <- list.files(path)

##### 월 바꿔가면서 해보자...
file_list <- file_list[which(substr(file_list, 10, 11) == "05")]
file_list <- paste0(path, file_list)

for (file_name in file_list) {
    print(file_name)
    9
    data <- read_csv(file_name)
    
    # Data per hours
    data2 <- data %>% 
        # filter(MMSI == 367528080) %>% 
        arrange(BaseDateTime) %>% 
        mutate(
            Hour = as.POSIXct(strptime(as.character(BaseDateTime), "%Y-%m-%d %H"))
        ) %>% 
        # dplyr::select(MMSI, Hour, LAT, LON) %>% 
        group_by(MMSI, Hour) %>% 
        filter(row_number() == 1) %>% 
        ungroup()
    
    # Filter Cargo ships
    df <- data2 %>% 
        # group_by(MMSI) %>% 
        # arrange(MMSI) %>% 
        # mutate(dif = abs(LON - lag(LON, 1, LON[1])) + abs(LAT - lag(LAT, 1, LAT[1]))) %>% 
        # ungroup() %>% 
        # # ship is moved
        # filter(dif > 0.5) %>% 
        # dplyr::select(-dif) %>% 
        # select huge ships
        filter(VesselType == 70 & Length > 260 & Width > 32)
    
    if (file_name == file_list[1]) {
        ais_data <- df
    } else {
        ais_data <- rbind(ais_data, df)
    }
}
ais_data
tail(ais_data)
length(unique(ais_data$MMSI))
# save(ais_data, file = "~/GoogleDrive/Lab/KHS/manifold_clust/real_data/ais.RData")


df <- ais_data %>% 
    filter(SOG > 0)

world <- ne_countries(scale = "medium", returnclass = "sf")
map_bg <- ggplot(data = world) +
    geom_sf() +
    coord_sf(xlim = range(ais_data$LON) + c(-10, 10),
             ylim = range(ais_data$LAT) + c(-10, 10),
             expand = FALSE)
map_bg +
    geom_path(
        data = ais_data,
        aes(
            x = LON,
            y = LAT,
            # group = MMSI,
            color = factor(MMSI)
        ),
        size = 0.3
    ) +
    theme_bw() +
    theme(legend.position = "none")


ais_data %>% 
    filter(MMSI == 218614000) %>% 
    # dplyr::select(BaseDateTime, LAT, LON) %>% 
    as.data.frame()





