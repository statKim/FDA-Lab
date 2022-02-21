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
source("OpenSkyNetwork/connect_impala_shell.R")

### Access to Impala shell
session <- osn_connect("statKim")

query <- "show tables;"
query <- "describe state_vectors_data4;"
query <- "describe flights_data4;"
query <- "select * FROM allcall_replies_data4 LIMIT 10"
query <- "SELECT * FROM state_vectors_data4 WHERE icao24='a0d724' AND time>=1480760100 AND time<=1480764600 AND hour>=1480759200 AND hour<=1480762800;"
query <- "SELECT
icao24, firstseen, estdepartureairport, lastseen, estarrivalairport, callsign 
FROM flights_data4 
WHERE icao24='a0d724' AND hour>=1480759200 AND hour<=1480762800;"
query <- "SELECT
icao24, firstseen, estdepartureairport, lastseen, estarrivalairport, callsign, day
FROM flights_data4
WHERE day = 1612483200;"
query <- "SELECT
icao24, firstseen, estdepartureairport, lastseen, estarrivalairport, callsign, day
FROM flights_data4
WHERE day = 1612483200 and estdepartureairport = 'KRIV';"
query <- "SELECT
icao24, firstseen, estdepartureairport, lastseen, estarrivalairport, callsign, day
FROM flights_data4
WHERE day = 1612483200 and nonnullvalue(estdepartureairport);"


#####################################################################
### Get icao24 and airline variables from "flights_data4"
### - Reference for SQL function for impala shell
###     - https://impala.apache.org/docs/build/html/topics/impala_conditional_functions.html#conditional_functions__nonnullvalue
#####################################################################

### 2019-03-01 ~ 2019-06-30
start_period <- date2unixtime("2019-03-01 00:00:00")
end_period <- date2unixtime("2019-07-01 00:00:00")
query <- paste0(
    "SELECT
icao24, firstseen, estdepartureairport, lastseen, estarrivalairport, callsign, day
FROM flights_data4
WHERE day >= ", start_period, " and day <= ", end_period, "
and nonnullvalue(estarrivalairport) and nonnullvalue(callsign)
and estdepartureairport = 'RKSI';"
)
query <- gsub("\n", " ", query)   # remove "\n"
query

### Get data from database
system.time({
    flight <- impala_query(session, query)    
})
flight

unixtime2date(flight$day %>% unique %>% sort)
apply(flight, 2, function(col){ sum(is.na(col)) })


### Remove NA
df2 <- flight %>% 
    drop_na() %>% 
    # filter(!is.na(callsign) & 
    #            !is.na(estarrivalairport) & 
    #            (estdepartureairport != estarrivalairport)) %>% 
    group_by(icao24, callsign, estdepartureairport, estarrivalairport) %>% 
    select(icao24, estdepartureairport, estarrivalairport, callsign) %>% 
    ungroup() %>% 
    distinct()
df2
length(unique(df2$icao24))   # 2042


### Count per arrival
### - KSFO(San Francisco), KLAX(LA), KJFX(New York), KSEA(Seatle), KDFW(Dallas), KORD(Chicago)
df2 %>% 
    filter(substr(estarrivalairport, 1, 1) == "K") %>% 
    group_by(estarrivalairport) %>% 
    summarise(n = n()) %>% 
    arrange(desc(n))
#   estarrivalairport     n
# 1 KHHR                 96
# 2 KSFO                 83
# 3 KLAX                 70
# 4 KJFK                 38
# 5 KDFW                 36
# 6 KSEA                 34
# 7 KORD                 24


### Count per airline
df2 %>% 
    mutate(airline = substr(callsign, 1, 3)) %>% 
    group_by(airline) %>% 
    summarise(n = n()) %>% 
    arrange(desc(n))
#   airline     n
# 1 KAL      1954
# 2 AAR      1005
# 3 JJA       717
# 4 JNA       494

### Unique flights of Incheon to LA
flight_to_LA <- flight %>% 
    drop_na() %>% 
    group_by(icao24, callsign, estdepartureairport, estarrivalairport) %>% 
    select(icao24, estdepartureairport, estarrivalairport, callsign, day) %>% 
    ungroup() %>% 
    filter(substr(estarrivalairport, 1, 4) == "KLAX") %>% 
    mutate(day = substr(unixtime2date(day), 1, 10),
           airline = substr(callsign, 1, 3)) %>% 
    distinct()


### Unique icao24, airline
flight_id <- flight_to_LA %>% 
    select(icao24, airline) %>% 
    distinct()
flight_id
table(flight_id$airline)
# AAR GTI KAL 
# 23   4  23 
flight_id$icao24

flight_id <- flight_id %>% 
    filter(airline %in% c("KAL","AAR"))


#####################################################################
### (Test) Get trajectories from "state_vectors_data4"
#####################################################################

### 2019-03-01 ~ 2019-06-30
### => Too much time...Do not use it!!!
start_period <- date2unixtime("2019-03-01 00:00:00")
end_period <- date2unixtime("2019-07-01 00:00:00")
query <- paste0(
    "SELECT *
FROM state_vectors_data4
WHERE icao24 IN (", 
paste0("'", flight_id$icao24, "'", collapse = ", "),
")",
" AND hour >= ", start_period, " AND hour <= ", end_period, 
" AND time >= ", start_period, " AND time <= ", end_period, 
" AND SUBSTR(callsign, 1, 3) IN (", 
paste0("'", unique(flight_id$airline), "'", collapse = ", "),
")",
" AND NONNULLVALUE(lat) AND NONNULLVALUE(lon)
AND NONNULLVALUE(callsign);"
)


### Get 1st icao24 for test
start_period <- date2unixtime("2019-03-01 00:00:00")
end_period <- date2unixtime("2019-07-01 00:00:00")
query <- paste0(
    "SELECT *
FROM state_vectors_data4
WHERE icao24 = '", 
flight_id$icao24[1],
"' AND hour >= ", start_period, " AND hour <= ", end_period, 
" AND time >= ", start_period, " AND time <= ", end_period, 
" AND SUBSTR(callsign, 1, 3) IN (", 
paste0("'", unique(flight_id$airline), "'", collapse = ", "),
")",
" AND NONNULLVALUE(lat) AND NONNULLVALUE(lon)
AND NONNULLVALUE(callsign);"
)
query <- gsub("\n", " ", query)   # remove "\n"
query


### Get data from database
system.time({
    flight_state <- impala_query(session, query)    
})
# user   system  elapsed 
# 35.223    8.214 1433.202 
flight_state
# save(flight, flight_state, file = "RData/flight.RData")

unique(flight_state$icao24)
unixtime2date(flight_state$hour %>% unique %>% sort)
apply(flight_state, 2, function(col){ sum(is.na(col)) })

flight_state$callsign %>% unique
table(flight_state$onground)


flight_to_LA %>% 
    filter(icao24 == flight_id$icao24[1])

sub <- flight_state %>% 
    # filter(vertrate == 0) %>% 
    arrange(time) %>% 
    select(time, icao24, lat, lon, vertrate, callsign, onground, hour) %>% 
    mutate(time = unixtime2date(time),
           hour = unixtime2date(hour))

sub2 <- sub %>% 
    mutate(minute = substr(time, 1, 16),
           day = substr(time, 1, 10)) %>% 
    group_by(icao24, callsign, minute) %>% 
    filter(row_number() == n()) %>%    # select last row
    ungroup()

sub3 <- sub2 %>% 
    filter(substr(time, 1, 10) == "2019-06-09") %>% 
    # filter(callsign == "KAL595" & substr(time, 1, 10) == "2019-03-01") %>% 
    as.data.frame



sub3 <- flight_to_LA %>% 
    filter(icao24 == flight_id$icao24[1]) %>% 
    select(icao24, callsign, day, airline) %>% 
    left_join(sub2, by = c("icao24", "day", "callsign"))
    

plot(sub3$lon, sub3$lat, type = "l")

sub3 %>% 
    # filter(onground == TRUE) %>% 
    select(icao24, callsign, day, minute, lat, lon) %>% 
    group_by(icao24, callsign, day) %>% 
    filter(row_number() == n())

### Flight trajectories example
library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

world <- ne_countries(scale = "medium", returnclass = "sf")
map_bg <- ggplot(data = world) +
    geom_sf() +
    coord_sf(xlim = range(sub3$lon) + c(-10, 10),
             ylim = range(sub3$lat) + c(-10, 10), 
             expand = FALSE)
map_bg + 
    geom_path(
        data = sub3, 
        aes(
            x = lon, 
            y = lat, 
            group = day,
            color = airline
        )
    ) +
    theme_bw()



#####################################################################
### Get trajectories from "state_vectors_data4"
#####################################################################

### Access to Impala shell
session <- osn_connect("statKim")

### Get all data - 50 icao24
start_period <- date2unixtime("2019-03-01 00:00:00")
end_period <- date2unixtime("2019-07-01 00:00:00")
for (i in 1:length(flight_id$icao24)) {
    print(i)
    
    query <- paste0(
        "SELECT *
FROM state_vectors_data4
WHERE icao24 = '", 
flight_id$icao24[i],
"' AND hour >= ", start_period, " AND hour <= ", end_period, 
" AND time >= ", start_period, " AND time <= ", end_period, 
" AND SUBSTR(callsign, 1, 3) IN (", 
paste0("'", unique(flight_id$airline), "'", collapse = ", "),
")",
" AND NONNULLVALUE(lat) AND NONNULLVALUE(lon)
AND NONNULLVALUE(callsign);"
    )
    query <- gsub("\n", " ", query)   # remove "\n"
    query
    
    ### Get data from database
    t1 <- Sys.time()
    data <- impala_query(session, query)    
    t2 <- Sys.time()
    print(paste(flight_id$icao24[i], ":", round(t2 - t1, 2), "mins"))
    # user   system  elapsed 
    # 35.223    8.214 1433.202 
    
    ### Save data
    if (i == 1) {
        flight_state <- data     
    } else {
        flight_state <- rbind(flight_state,
                              data)
    }
    save(flight, flight_state, file = "RData/flight.RData")
}
dim(flight_state)   # 69756853       17

### Terminate the connection
osn_disconnect(session)



#####################################################################
### Data preprocessing
#####################################################################

### Check the number of missings
### Caution!! Because of too large data, it takes a lot of times. (8 mins)
apply(flight_state, 2, function(col){ sum(is.na(col)) })
# time        icao24           lat           lon      velocity       heading      vertrate      callsign      onground 
# 0             0             0             0        720501        720501        720373             0             0 
# alert           spi        squawk  baroaltitude   geoaltitude lastposupdate   lastcontact          hour 
# 0             0       4525967       1079393       2664876             0             0             0 

### Transform unix time to date
### Create minute and day variables
### 7 mins
flight_state_to_LA <- flight_state %>% 
    arrange(time) %>% 
    select(time, icao24, lat, lon, vertrate, callsign, onground, hour) %>% 
    mutate(time = unixtime2date(time),
           hour = unixtime2date(hour),
           minute = substr(time, 1, 16),
           day = substr(time, 1, 10)) %>% 
    group_by(icao24, callsign, minute) %>% 
    filter(row_number() == n()) %>%    # select last row
    ungroup()
dim(flight_state_to_LA)   # 1188867      10
flight_state_to_LA
save(flight_to_LA, flight_state_to_LA, file = "RData/flight_reduced.RData")


### Test trajectories
sub3 <- flight_to_LA %>% 
    # filter(icao24 == flight_id$icao24[1]) %>% 
    select(icao24, callsign, day, airline) %>% 
    left_join(flight_state_to_LA, by = c("icao24", "day", "callsign")) %>% 
    drop_na()   # remove to trajectories of icao24 for such day
dim(sub3)   # 11610    11
# apply(sub3, 2, function(col){ sum(is.na(col)) })
# sub3[which(is.na(sub3$time)), ]

world <- ne_countries(scale = "medium", returnclass = "sf")
map_bg <- ggplot(data = world) +
    geom_sf() +
    coord_sf(xlim = range(sub3$lon) + c(-10, 10),
             ylim = range(sub3$lat) + c(-10, 10), 
             expand = FALSE)
map_bg + 
    geom_path(
        data = sub3, 
        aes(
            x = lon, 
            y = lat, 
            group = day,
            color = airline
        )
    ) +
    theme_bw()

sub3[which(sub3$lat > 50), ]
sub3 %>% 
    filter(icao24 == "71be27" & day == "2019-06-23") %>% 
    as.data.frame()

### Remove un-usual trajectory (it seems to be malfunction of the ADS-B reciever)
sub3 <- sub3 %>% 
    filter( !(icao24 == "71be27" & day == "2019-06-23") ) 
dim(sub3)   # 11536    11
sub3


### RKSI - lat: 37.46, lon: 126.44
### KLAX - lat: 33.94, lon: -118.41
sub3 %>% 
    group_by(icao24, callsign, day) %>% 
    arrange(time) %>% 
    filter(row_number() == 1) %>%     # select 1st row
    ungroup() %>% 
    select(icao24, callsign, time, lat, lon, vertrate) %>% 
    mutate(lat = round(lat, 2),
           lon = round(lon, 2)) %>% 
    as.data.frame()
sub3 %>% 
    group_by(icao24, callsign, day) %>% 
    arrange(time) %>% 
    filter(row_number() == n()) %>%     # select last row
    ungroup() %>% 
    select(icao24, callsign, time, lat, lon, vertrate) %>% 
    mutate(lat = round(lat, 2),
           lon = round(lon, 2)) %>% 
    as.data.frame()


sub3 %>% 
    filter(round(lat, 1) == 37.5 & round(lon, 1) == 126.4) %>% 
    select(icao24, callsign, day, time, lat, lon) %>% 
    group_by(icao24, callsign, day) %>% 
    filter(row_number() == 1) %>% 
    as.data.frame()


sub3 %>% 
    filter(round(lat, 1) == 33.9 & round(lon, 1) == -118.4) %>% 
    select(icao24, callsign, day, time, lat, lon) %>% 
    group_by(icao24, callsign, day) %>% 
    filter(row_number() == 1) %>% 
    as.data.frame()
# count > 2이상인거??



### Obtain Take-off and Landing time
### - Incheon(RKSI) -> lat: 37.46, lon: 126.44
### - LA(KLAX) -> lat: 33.94, lon: -118.41

# Take-off time for RKSI location
takeoff <- flight_state_to_LA %>% 
    filter( (lat > 36.46 & lat < 38.46) & (lon > 125.44 & lon < 127.44) ) %>% 
    arrange(time) %>%
    select(time, icao24, lat, lon, vertrate, callsign, onground, hour) %>% 
    mutate(time = unixtime2date(time),
           hour = unixtime2date(hour),
           minute = substr(time, 1, 16),
           day = substr(time, 1, 10)) %>% 
    group_by(icao24, callsign, day) %>% 
    filter(row_number() == n()) %>% 
    ungroup() %>% 
    select(time, icao24, callsign, minute, day) %>% 
    mutate(day = as.Date(day),
           takeoff = time)
takeoff

# Landing time for KLAX location
landing <- flight_state_to_LA %>% 
    filter( (lat > 32.94 & lat < 34.94) & (lon > -119.41 & lon < -117.41) ) %>% 
    arrange(time) %>%
    select(time, icao24, lat, lon, vertrate, callsign, onground, hour) %>% 
    mutate(time = unixtime2date(time),
           hour = unixtime2date(hour),
           minute = substr(time, 1, 16),
           day = substr(time, 1, 10)) %>% 
    group_by(icao24, callsign, day) %>% 
    filter(row_number() == 1) %>% 
    ungroup() %>% 
    select(time, icao24, callsign, minute, day) %>% 
    mutate(day = as.Date(day),
           landing = time)
landing


# Each flight flies on once per 1 day
takeoff %>% 
    group_by(icao24, callsign, day) %>% 
    summarise(n = n()) %>% 
    filter(n > 1)
landing %>% 
    group_by(icao24, callsign, day) %>% 
    summarise(n = n()) %>% 
    filter(n > 1)

# Take-off before 15:00
takeoff_landing_1 <- takeoff %>% 
    left_join(landing, by = c("icao24", "callsign", "day")) %>% 
    drop_na() %>% 
    select(icao24, callsign, takeoff, landing) %>% 
    mutate(duration = landing - takeoff) %>% 
    filter(duration > 0)   # remove negative flying time
takeoff_landing_1
takeoff_landing_1$duration

# Take-off after 15:00 (Day is changed!!)
takeoff_landing_2 <- takeoff %>% 
    mutate(day = as.Date(day) + 1) %>% 
    left_join(landing, by = c("icao24", "callsign", "day")) %>% 
    drop_na() %>% 
    select(icao24, callsign, takeoff, landing) %>% 
    mutate(duration = landing - takeoff) %>% 
    filter(duration < 15)   # remove too large flying time
takeoff_landing_2
takeoff_landing_2$duration

# Combine
takeoff_landing <- rbind(takeoff_landing_1,
                         takeoff_landing_2) %>% 
    arrange(takeoff) %>% 
    mutate(flight_id = row_number(),
           airline = substr(callsign, 1, 3))
takeoff_landing   # 590  5
takeoff_landing$icao24 %>% 
    unique()
table(takeoff_landing$airline)
# AAR GTI KAL 
# 231  53 306 
# Asiana Airlines, Atlas Air, Korean Air


### Get data for each takeoff, landing
flight_df <- takeoff_landing %>% 
    left_join(flight_state_to_LA, 
              by = c("icao24","callsign")) %>% 
    filter(time >= takeoff & time <= landing) %>% 
    select(flight_id, time, icao24, airline, callsign, lat, lon, takeoff, landing, duration)
flight_df   # 43326  10


### Draw trajectories again
world <- ne_countries(scale = "medium", returnclass = "sf")
map_bg <- ggplot(data = world) +
    geom_sf() +
    coord_sf(xlim = range(flight_df$lon) + c(-10, 10),
             ylim = range(flight_df$lat) + c(-5, 5), 
             expand = FALSE)
map_bg + 
    geom_path(
        data = flight_df, 
        aes(
            x = lon, 
            y = lat, 
            group = flight_id,
            color = airline
        ),
        size = 0.3
    ) +
    theme_bw()








