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
icao24, firstseen, estdepartureairport, lastseen, estarrivalairport, callsign, 
day
FROM flights_data4
WHERE DAY = 1612483200;"

### Get data from database
query <- gsub("\n", " ", query)   # remove "\n"
system.time({
    df <- impala_query(session, query)    
})
df

unixtime2date(df$day %>% unique)
dim(df)
length(unique(df$icao24))
df[df$icao24 == "a04c67", ]

sub <- df %>% 
    filter(!is.na(callsign) & !is.na(estdepartureairport) & !is.na(estarrivalairport)) %>% 
    group_by(icao24, callsign, estdepartureairport, estarrivalairport) %>% 
    select(icao24, estdepartureairport, estarrivalairport, callsign)
dim(sub)
length(unique(sub$icao24))
dim(unique(sub))

sub %>% 
    group_by(icao24) %>% 
    summarise(n = n()) %>% 
    filter(n > 1)
df %>% 
    filter(icao24 == "008de1") %>% 
    mutate(firstseen = unixtime2date(firstseen),
           lastseen = unixtime2date(lastseen))


date2unixtime("2021-02-05 09:00:00")


### Terminate the connection
osn_disconnect(session)











# devtools::install_github("espinielli/osn")
library(logger)
log_threshold(TRACE, namespace = 'osn')

library(osn)
session <- osn_connect("statkim")

# EDDF
df <- state_vector(
    session,
    icao24 = c("3c6589", "3c6757"),
    wef = "2019-01-01 09:00:00",
    til = "2019-01-01 10:00:00",
    # bbox = c(xmin = 7.553013, ymin = 49.378819,  xmax = 9.585482, ymax = 50.688044)
)
df <- state_vector(
    session,
    icao24 = c("718935","718940","718933","718a24","718a49","718936","718a22","718934","718a23"),
    wef = "2019-01-01 00:00:00",
    til = "2019-01-30 00:00:00",
    # bbox = c(xmin = 7.553013, ymin = 49.378819,  xmax = 9.585482, ymax = 50.688044)
)
df <- as.data.frame(df)
dim(df)
tail(df)
df$lat[!is.na(df$lat)] %>% order %>% tail
df[168:200, ]

# Convert unix timestamp to YYMMDD hhmmss
as.POSIXct(df$time, origin = "1970-01-01")[1:5]

as.POSIXct(df$hour %>% unique, origin = "1970-01-01")



plot(df$lon, df$lat, type = "l")

# 30분 정도 걸림...
system.time({
    # arriv <- arrivals(session, "EDDF", "2019-04-22 00:00:00", til=NULL)
    arriv <- arrivals_state_vector(session, "EDDF", "2019-04-22 00:00:00", til=NULL)
})
arriv <- as.data.frame(arriv)
dim(arriv)
head(arriv)
unique(arriv$icao24)
as.POSIXct(max(arriv$lastseen), origin = "1970-01-01")


sub <- arriv %>% 
    filter(onground == TRUE & !is.na(callsign) & !is.na(estdepartureairport)) %>% 
    group_by(icao24, callsign, estdepartureairport, estarrivalairport) %>% 
    slice(1)
unique(sub$icao24) %>% length
as.POSIXct(sub$hour %>% unique %>% sort, origin = "1970-01-01")

arriv %>% 
    filter(onground == TRUE & !is.na(callsign) & !is.na(estdepartureairport)) %>% 
    group_by(icao24, callsign, estdepartureairport, estarrivalairport) %>% 
    summarise(n = n())

# callsign 앞의 3개 문자로 항공사 알 수 있음
# https://en.wikipedia.org/wiki/List_of_airline_codes



aircraft <- data.table::fread("/Users/hyunsung/Desktop/aircraftDatabase-2020-11.csv")
aircraft <- as.data.frame(aircraft)
dim(aircraft)
head(aircraft)

aircraft[substr(aircraft$icao24, 1, 3) == "718", ]$icao24


as.Date(as.POSIXct(arriv$time, origin = "1970-01-01")) %>% unique()



