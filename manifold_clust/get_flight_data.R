#####################################################################
### Aircraft trajectory data from OpenSky Network
### - We use the historical database using impala shell.
### - To access the database, you need to log-in and 
###   submit the "Historical Data Access Application".
### - https://opensky-network.org/data/apply
### - The useful frameworks exists to access the database :
###     - https://opensky-network.org/data/data-tools
### - We use the SQL query to get data.
#####################################################################
library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
source("OpenSkyNetwork/connect_impala_shell.R")

### Access to Impala shell
session <- osn_connect("statKim")

### 2019-03-01 ~ 2019-06-30
start_period <- date2unixtime("2019-07-01 00:00:00")
end_period <- date2unixtime("2019-12-31 23:59:59")
arrival <- "EGLL"   # arrival airport icao24

# ### Throw query to Impalla DB
# query <- paste0(
#     "SELECT
# icao24, firstseen, estdepartureairport, lastseen, estarrivalairport, callsign, day
# FROM flights_data4
# WHERE day >= ", start_period, " and day <= ", end_period, "
# and nonnullvalue(callsign)",
# "and nonnullvalue(estarrivalairport)",
# "and nonnullvalue(estdepartureairport)",
# "and estarrivalairport = '", arrival, "'",
# # "and estarrivalairport = 'EGLL'",
# # "and estarrivalairport = 'RKSI'",
# # "and estdepartureairport = 'EGLL'",
# # "and estdepartureairport = 'WSSS'",
# ";"
# )
# query <- gsub("\n", " ", query)   # remove "\n"
# query
# 
# ### Get data from database
# system.time({
#     flight <- impala_query(session, query)    
# })
# flight
# 
# ### Check the number of flight per departure airport
# flight %>% 
#     filter(!(substr(estdepartureairport, 1, 1) %in% c("E","L","K"))) %>%
#     dplyr::select(estdepartureairport) %>% 
#     # dplyr::select(estarrivalairport) %>% 
#     table() %>% 
#     sort()


### Select specefic departure airport
# LLBG: 이스라엘, WSSS: 싱가포르, LD, LG
departure <- "KJFK"   # departure airport icao24
# ind <- which(flight$estdepartureairport == departure)
flight %>%
    mutate(id = row_number()) %>%
    filter(estdepartureairport == departure) %>%
    mutate(takeoff = unixtime2date(firstseen),
           landing = unixtime2date(lastseen),
           airline = substr(callsign, 1, 3)) %>%
    arrange(takeoff) %>%
    filter(takeoff >= as.POSIXct("2019-10-01") & takeoff <= as.POSIXct("2019-10-31")) %>%
    dplyr::select(airline) %>%
    table
ind <- flight %>%
    mutate(id = row_number()) %>%
    filter(estdepartureairport == departure) %>%
    mutate(takeoff = unixtime2date(firstseen),
           landing = unixtime2date(lastseen),
           airline = substr(callsign, 1, 3)) %>%
    arrange(takeoff) %>%
    filter(airline %in% c("AAL","DAL","JBU")) %>%
    filter(takeoff >= as.POSIXct("2019-10-01") & takeoff <= as.POSIXct("2019-10-05")) %>%
    dplyr::select(id) %>%
    unlist() %>%
    as.numeric()


### Total 667 flights
for (j in 1:length(ind)) {
    print(j)
    i <- ind[j]
    
    query <- paste0(
        "SELECT *
FROM state_vectors_data4
WHERE icao24 = '", flight$icao24[i], "'",
" AND hour >= ", flight$firstseen[i] - (flight$firstseen[i] %% 3600), 
" AND hour <= ", flight$lastseen[i] - (flight$lastseen[i] %% 3600),
" AND time >= ", flight$firstseen[i], " AND time <= ", flight$lastseen[i], 
# " AND callsign = '", flight$callsign[i], "'",
" AND NONNULLVALUE(lat) AND NONNULLVALUE(lon)",
";"
    )
    query <- gsub("\n", " ", query)   # remove "\n"
    query
    
    ### Get data from database
    t1 <- Sys.time()
    data <- impala_query(session, query)    
    t2 <- Sys.time()
    print(paste(flight$icao24[i], ":", round(t2 - t1, 2)))
    
    # data %>% 
    #     filter(callsign == flight$callsign[i])
    
    
    # flight[i, ] %>% 
    #     mutate(firstseen = unixtime2date(firstseen),
    #            lastseen = unixtime2date(lastseen),
    #            daty = unixtime2date(day))
    
    df <- data %>% 
        filter(callsign == flight$callsign[i]) %>%
        dplyr::select(time, icao24, callsign, lat, lon, hour) %>% 
        filter(time >= flight$firstseen[i] & time <= flight$lastseen[i]) %>% 
        arrange(time) %>% 
        mutate(time = unixtime2date(time),
               hour = unixtime2date(hour),
               minute = substr(time, 1, 16),
               day = substr(time, 1, 10),
               flight_id = i)
    df <- df %>%
        group_by(minute) %>%
        filter(row_number() == 1) %>% 
        ungroup()
    df
    
    if (j == 1) {
        flight_df <- df
        # plot(df$lon, df$lat, col = i)
    } else {
        flight_df <- rbind(flight_df,
                           df)
        # points(df$lon, df$lat, col = i)
    }
    
    
    # world <- ne_countries(scale = "medium", returnclass = "sf")
    # map_bg <- ggplot(data = world) +
    #     geom_sf() +
    #     coord_sf(xlim = range(df$lon) + c(-10, 10),
    #              ylim = range(df$lat) + c(-10, 10),
    #              expand = FALSE)
    # map_bg +
    #     geom_point(
    #         data = df,
    #         aes(
    #             x = lon,
    #             y = lat,
    #             # group = day,
    #         ),
    #         color = "red"
    #     ) +
    #     theme_bw()
}

flight_df <- flight_df %>%
    mutate(airline = substr(callsign, 1, 3))

flight_df %>%
    group_by(flight_id) %>%
    filter(row_number() == 1) %>%
    ungroup() %>%
    group_by(airline) %>%
    summarise(n = n())

save(flight_df, file = paste0("RData/", departure, "_to_", arrival, ".RData"))


# world <- ne_countries(scale = "medium", returnclass = "sf")
# map_bg <- ggplot(data = world) +
#     geom_sf() +
#     coord_sf(xlim = range(flight_df$lon) + c(-10, 10),
#              ylim = range(flight_df$lat) + c(-10, 10), 
#              expand = FALSE)
# map_bg + 
#     geom_point(
#         data = flight_df, 
#         aes(
#             x = lon, 
#             y = lat, 
#             group = flight_id,
#             color = airline,
#         ),
#         size = 0.3
#     ) +
#     theme_bw()
