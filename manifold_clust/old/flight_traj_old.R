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












library(arrow)
df <- read_parquet("/Users/hyunsung/Desktop/flights_data5_develop_sample_2019-09-01.snappy.parquet", as_tibble = TRUE)
df <- as.data.frame(df)
head(df)
dim(df)
length(unique(df$icao24))
df[df$icao24 == "44d2b8", ]

as.Date("1567296000", format='%m%d%Y')

date("1567296000")








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

source("connect_impala_shell.R")

### Access to Impala shell
session <- osn_connect("statKim")

query <- "select * FROM allcall_replies_data4 LIMIT 10"
query <- "SELECT * FROM state_vectors_data4 WHERE icao24='a0d724' AND time>=1480760100 AND time<=1480764600 AND hour>=1480759200 AND hour<=1480762800;"

### Get data from database
df <- impala_query(session, query)
df

### Terminate the connection
osn_disconnect(session)







# # Install the development version from GitHub:
# # install.packages("devtools")
# # devtools::install_github("luisgasco/openskyr")
# 
library(openskyr)
state_vectors_df <- get_state_vectors()
dim(state_vectors_df)

username <- "statKim"
password <- ""
data_flight_track <- get_track_data(username = username, 
                                    password = password,
                                    icao24 = "494103",
                                    time = 1587126600)
dim(data_flight_track)





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




paste0("script log.txt ssh -o StrictHostKeyChecking=no -p 2230 -l ", username, " data.opensky-network.org")
password

"cat log.txt | grep '^|.*' | sed -e 's/\s*|\s*/,/g' -e 's/^,\|,$//g' -e 's/NULL//g' | awk '!seen[$0]++' >> log.csv"

script <- "ls -al"
script <- paste0("script log.txt ssh -o StrictHostKeyChecking=no -p 2230 -l ", username, " data.opensky-network.org")
system(script, intern = T)


query <- "show tables;"

query <- "SELECT * FROM state_vectors_data4 WHERE icao24='a0d724' AND time>=1480760100 AND time<=1480764600 AND hour>=1480759200 AND hour<=1480762800;"
username <- "statKim"
password <- ""

df <- impala_query(query, username, password)
head(df)


impala_query <- function(query, username, password) {
    session <- ssh::ssh_connect(paste0(username,"@data.opensky-network.org:2230"), passwd=password)
    lines <- base::rawToChar(ssh::ssh_exec_internal(session, paste("-q ", query))$stdout)
    if(lines == "") {
        return(NULL)
    }
    lines <- unlist(strsplit(lines, '\n'))
    colnames_line <- lines[2]
    colnames_line <- gsub("^.|.$", "", colnames_line)
    columnNames <- unlist(strsplit(colnames_line, '|', fixed=TRUE))
    columnNames <- trimws(columnNames)
    data_lines <- lines[4:(length(lines)-1)]
    data_lines <- gsub("^.|.$", "", data_lines)
    split_data_lines <- strsplit(data_lines, '|', fixed=TRUE)
    trimmed_data_lines <- lapply(split_data_lines, trimws)
    if(length(trimmed_data_lines) > 1024) {
        extra_pages_number <- length(trimmed_data_lines) %/% 1024
        repeated_header_lines <- sequence(nvec=rep(4, extra_pages_number), 
                                          from=(1:extra_pages_number)*1024 + 1
                                          + seq(from=0, by=4, length=extra_pages_number))
        trimmed_data_lines <- trimmed_data_lines[-repeated_header_lines]
    }
    data_matrix <- matrix(unlist(trimmed_data_lines), ncol=length(columnNames), byrow=TRUE)
    colnames(data_matrix) <- columnNames
    ssh::ssh_disconnect(session)
    
    # data_matrix <- as.data.frame(data_matrix)
    
    return(data_matrix)
}





library(openSkies)
## In the following example, we retrieve information about all flights
## that departed from Frankfurt Airport the 29th of January 2018 after 12 pm.
## It should be noted that such large requests can take a long time unless
## performed with the Impala shell
## It should also be noted that, in this and the following examples, the
## values for the username and password arguments should be substituted with
## personal credentials registered at the OpenSky Network
flights <- getAirportArrivals(airport="EDDF", startTime="2018-01-29 12:00:00",
                              endTime="2018-01-29 24:00:00", timeZone="Europe/Berlin",
                              includeStateVectors = TRUE, timeResolution = 60,
                              useImpalaShell = TRUE, username = username,
                              password = password)
## We can then easily check the amount of flights
length(flights) # 316 flights

## Trajectories of the 316 flights can be obtained by retrieving
## the set of state vectors of each flight
trajectories_frankfurt <- lapply(flights, function(f) f$state_vectors)

## It is also possible to retrieve all state vectors received from a
## given aircraft in a given period of time. In the following example, we
## obtain all the state vectors received for aircraft with ICAO24 code
## 403003 between 12 PM of the 8th of October, 2020 and 6 PM of the
## 9th of October, 2020
stateVectors <- getAircraftStateVectorsSeries("403003", startTime = "2020-10-08 12:00:00",
                                              endTime = "2020-10-09 18:00:00",
                                              timeZone="Europe/London",
                                              timeResolution = 60,
                                              useImpalaShell = TRUE,
                                              username = username,
                                              password = password)

## The ensemble of state vectors can then be split into individual
## flight instances, revealing that the aircraft performed 6 flights
## in the analyzed period
flights <- stateVectors$split_into_flights()
length(flights) # 6 flights

## Let us get some additional information about the flights performed by
## this aircraft. For example, the maximum speed that it reached in km/h
maxSpeed <- max(stateVectors$get_values("velocity", removeNAs = TRUE))/1000*3600

## The maximum speed that it reached is just above 210 km/h. This is well below
## the speed typically used by commercial passenger jets, usually around 800 km/h
## Investigation of the aircraft model confirms that it is indeed a small CESSNA 152
aircraftData <- getAircraftMetadata("403003")
aircraftData$model

## First, let us obtain the trajectories of the flights performed
## by the CESSNA 152
trajectories <- lapply(flights, function(f) f$state_vectors)

## Then, we create a color palette
library(viridis)
colors <- magma(length(trajectories))

## Then, the trajectories are plotted
plotRoutes(trajectories, pathColors = colors, lineSize = 1.2,
           lineAlpha = 0.8, includeArrows = TRUE,
           paddingFactor = 0.05)

## Firstly we retrieve the state vectors of all aircraft flying over Switzerland
## the 2nd of March, 2018 at 14.30 (London time)
vectors_switzerland <- getSingleTimeStateVectors(time="2018-03-02 14:30:00",
                                                 timeZone="Europe/London",
                                                 minLatitude=45.8389,
                                                 maxLatitude=47.8229,
                                                 minLongitude=5.9962,
                                                 maxLongitude=10.5226,
                                                 username=username,
                                                 password=password,
                                                 useImpalaShell = TRUE)
## Then, the aircraft are plotted
plotPlanes(vectors_switzerland)


