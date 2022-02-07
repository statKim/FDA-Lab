#####################################################################
### Typhoon Trajectory data
### - International Best Track Archive for Climate Stewardship (IBTrACS)
### - https://www.ncei.noaa.gov/products/international-best-track-archive?name=gisSLD
### - Column documentation
###   https://www.ncei.noaa.gov/sites/default/files/2021-07/IBTrACS_v04_column_documentation.pdf
### - 참고 : https://wsyang.com/2014/02/typhoon-trajectories/
#####################################################################
library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

### Load data while skipping 2nd row
all_content <- readLines("/Users/hyunsung/Desktop/ibtracs.WP.list.v04r00.csv")
skip_second <- all_content[-2]
typhoon <- read.csv(textConnection(skip_second), header = TRUE)
dim(typhoon)   # 240690 163
head(typhoon)

typhoon <- typhoon[, -c(11:161)]   # remove unnecessary variables
typhoon$STORM_SPEED <- typhoon$STORM_SPEED * 0.5144   # knot -> m/s
head(typhoon)

### Typhoon trajectories
df <- typhoon[typhoon$SEASON %in% 2003:2004, ]   # select years
df$SID <- factor(df$SID)

world <- ne_countries(scale = "medium", returnclass = "sf")
map_bg <- ggplot(data = world) +
    geom_sf() +
    coord_sf(xlim = range(df$LON) + c(-10, 10),
             ylim = range(df$LAT) + c(-10, 10), 
             expand = FALSE)
# df <- df %>% 
#     arrange(order(as.POSIXlt(df$ISO_TIME)))
map_bg + 
    geom_path(data = df, 
              aes(x = LON, y = LAT, color = STORM_SPEED)) +
    theme_bw()


