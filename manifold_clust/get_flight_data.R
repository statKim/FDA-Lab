#############################################
### Test!!!
#############################################
library(tidyverse)
# library(osn)
source("OpenSkyNetwork/connect_impala_shell.R")

### Access to Impala shell
session <- osn_connect("statKim")

### 2019-03-01 ~ 2019-06-30
start_period <- date2unixtime("2019-01-01 00:00:00")
end_period <- date2unixtime("2019-12-31 00:00:00")
query <- paste0(
    "SELECT
icao24, firstseen, estdepartureairport, lastseen, estarrivalairport, callsign, day
FROM flights_data4
WHERE day >= ", start_period, " and day <= ", end_period, "
and nonnullvalue(callsign)",
"and nonnullvalue(estarrivalairport)",
"and nonnullvalue(estdepartureairport)",
"and estarrivalairport = 'EGLL'",
# "and estarrivalairport = 'RKSI'",
# "and estdepartureairport = 'EGLL'",
# "and estdepartureairport = 'WSSS'",
";"
)
query <- gsub("\n", " ", query)   # remove "\n"
query

### Get data from database
system.time({
    flight <- impala_query(session, query)    
})
flight

# query <- "SELECT * FROM state_vectors_data4 LIMIT 5;"
# df <- impala_query(session, query)    
# df$time %>% unixtime2date()

flight %>% 
    filter(!(substr(estdepartureairport, 1, 1) %in% c("E","L"))) %>%
    dplyr::select(estdepartureairport) %>% 
    # dplyr::select(estarrivalairport) %>% 
    table() %>% 
    sort()

# LLBG: 이스라엘, WSSS: 싱가포르
ind <- which(flight$estdepartureairport == "OERK")
ind <- sample(ind, 30)
i <- ind[1]
flight[i, ]


# plot(1, type = "n",
#      xlim = c(-10, 110), ylim = c(-10, 60))
for (i in ind[1:3]) {
    print(i)
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
    print(paste(flight_id$icao24[i], ":", round(t2 - t1, 2), "mins"))
    
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
    
    if (i == ind[1]) {
        flight_df <- df
        plot(df$lon, df$lat, col = i)
    } else {
        flight_df <- rbind(flight_df,
                           df)
        points(df$lon, df$lat, col = i)
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
    mutate(airline = substr(callsign, 1, 3)) %>% 
    filter(airline %in% c("BAW","SIA"))

flight_df %>% 
    group_by(flight_id) %>% 
    filter(row_number() == 1) %>% 
    ungroup() %>% 
    group_by(airline) %>% 
    summarise(n = n())

world <- ne_countries(scale = "medium", returnclass = "sf")
map_bg <- ggplot(data = world) +
    geom_sf() +
    coord_sf(xlim = range(flight_df$lon) + c(-10, 10),
             ylim = range(flight_df$lat) + c(-10, 10), 
             expand = FALSE)
map_bg + 
    geom_point(
        data = flight_df, 
        aes(
            x = lon, 
            y = lat, 
            group = flight_id,
            color = airline,
        ),
        size = 0.3
    ) +
    theme_bw()
