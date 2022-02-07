#####################################################################
### Aircraft Localization Competition
### - https://www.aicrowd.com/challenges/cyd-campus-aircraft-localization-competition
### - https://github.com/radoslawkrolikowski/adsb-flight-localization
#####################################################################
library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

### Load data
# ### 1.3GB - 56 secs
# system.time({
#     data <- read.csv("/Users/hyunsung/Desktop/round2/round2_training1.csv", header = T)
# })

### 1.3GB - 4 secs
flight <- data.table::fread("/Users/hyunsung/Desktop/round2/round2_training1.csv")
flight <- as.data.frame(flight)
dim(flight)   # 6535444 9
head(flight)

flight <- flight[, -9]   # remove 9th columns
flight$aircraft <- factor(flight$aircraft)

flight_id <- unique(flight$aircraft)   # unique id
length(flight_id)   # 2888

range(flight$timeAtServer)


### Trajectories
df <- flight %>% 
    filter(aircraft %in% flight_id[1:200])

world <- ne_countries(scale = "medium", returnclass = "sf")
map_bg <- ggplot(data = world) +
    geom_sf() +
    coord_sf(xlim = range(df$longitude) + c(-5, 5), 
             ylim = range(df$latitude) + c(-5, 5), 
             expand = FALSE)

map_bg +
    geom_path(data = df, aes(x = longitude, 
                             y = latitude, 
                             col = aircraft)) +
    theme_bw() +
    theme(legend.position = "none")




### sensors data
sensors <- read.csv("/Users/hyunsung/Desktop/round2/round2_sensors.csv", header = T)
dim(sensors)
head(sensors)















library(arrow)
df <- read_parquet("/Users/hyunsung/Desktop/flights_data5_develop_sample_2019-09-01.snappy.parquet", as_tibble = TRUE)
df <- as.data.frame(df)
head(df)
dim(df)
length(unique(df$icao24))
df[df$icao24 == "44d2b8", ]

as.Date("1567296000", format='%m%d%Y')

date("1567296000")


