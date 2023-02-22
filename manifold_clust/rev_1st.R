map_bg +
    geom_path(
        data = bird[bird$id %in% c(16,18,21,23,25,27,29,31,46), ],
        aes(
            x = lon,
            y = lat,
            group = id,
            color = season
        ),
        size = 0.3
    ) +
    labs(x = "Longitude", y = "Latitude", color = "Season") +
    theme_bw()

bird$id %>% unique()
bird %>% 
    filter(id %in% 46:50) %>% 
    dplyr::select(id, season) %>% 
    distinct()

# spring 16,18,21,23,25,27,29,31,46 (총 9개)


### Start and end date for 
bird_period <- data.frame(
    id = c(
        rep("Agri", 2), "Aras", rep("Ardahan", 3), "Armenia2", "Armenia3",
        rep("Cabuk", 3), rep("Haydi", 3), rep("Igdir", 5), rep("Iste", 12),
        rep("Logiya", 16), rep("Orada", 3), "Serhat", rep("Tuzluca", 6)
    ),
    start_date = c(
        "2013-09-11","2014-03-21",
        "2012-09-04",
        "2013-08-22","2014-04-03","2014-09-05",
        "2015-08-31",
        "2015-09-20",
        "2014-09-08","2015-05-28","2015-10-06",
        "2014-09-19","2015-04-11","2015-09-20",
        "2012-09-19","2013-03-16","2013-09-25","2014-03-03","2014-10-06",
        "2014-09-19","2015-04-14","2015-09-17","2016-03-22","2016-09-25",
        "2017-03-13","2017-09-24","2018-03-16","2018-09-16","2019-03-20",
        "2019-09-20","2020-03-21",
        "2014-05-19","2014-08-21","2015-04-17","2015-09-07","2016-04-16",
        "2016-09-14","2017-04-06","2017-09-15","2018-03-22","2018-09-06",
        "2019-03-11","2019-09-11","2020-03-24","2020-09-16","2021-03-22",
        "2021-09-13",
        "2015-08-28","2016-05-12","2016-09-15",
        "2014-09-21",
        "2013-08-23","2014-04-07","2014-09-06","2015-03-15","2015-09-08",
        "2016-03-14"
    ),
    end_date = c(
        "2013-09-28","2014-04-08",
        "2012-09-24",
        "2013-09-16","2014-04-30","2014-09-18",
        "2015-09-16",
        "2015-10-09",
        "2014-09-22","2015-07-09","2015-10-29",
        "2014-10-25","2015-05-13","2015-10-16",
        "2012-10-12","2013-04-04","2013-10-14","2014-03-23","2014-11-04",
        "2014-10-04","2015-04-30","2015-09-30","2016-04-10","2016-10-09",
        "2017-03-31","2017-10-09","2018-04-02","2018-09-30","2019-04-06",
        "2019-10-05","2020-04-07",
        "2014-06-16","2014-09-10","2015-05-26","2015-09-24","2016-05-09",
        "2016-09-28","2017-04-29","2017-09-28","2018-04-10","2018-09-21",
        "2019-03-30","2019-09-26","2020-04-13","2020-09-27","2021-04-06",
        "2021-09-25",
        "2015-09-23","2016-06-04","2016-10-02",
        "2014-10-08",
        "2013-10-11","2014-04-30","2014-10-03","2015-04-05","2015-10-02",
        "2016-04-10"
    )
)
bird_period <- bird_period %>% 
    mutate(curve_id = 1:n(),
           start_date = as.POSIXct(paste(start_date, "00:00:00"), tz = "UTC"),
           end_date = as.POSIXct(paste(end_date, "23:59:59"), tz = "UTC"),
           season = ifelse(as.integer(substr(start_date, 6, 7)) <= 6, 
                           "Spring", "Fall"))
bird_period$id
bird_period$curve_id
bird$id %>% unique()


bird_period %>% 
    filter(curve_id %in% c(16,18,21,23,25,27,29,31,46))
bird_period %>% 
    filter(id %in% c("Igdir","Iste","Logiya"))


bird_period %>% 
    mutate(diff = end_date - start_date) %>% 
    filter(season == "Spring") %>% 
    filter(curve_id %in% c(16,18,21,23,25,27,29,31,46))
bird_period %>% 
    mutate(diff = end_date - start_date) %>% 
    filter(season == "Spring") %>% 
    filter(!(curve_id %in% c(16,18,21,23,25,27,29,31,46)))

bird_period %>% 
    mutate(diff = end_date - start_date) %>% 
    filter(season == "Fall")


rfpca.obj <- RFPCA(Lt = Lt,
                   Ly = Ly,
                   optns = list(mfdName = "Sphere",
                                FVEthreshold = 1,
                                userBwMu = "GCV", 
                                userBwCov = "GCV"))
par(mfrow = c(1, 2))
plot(rfpca.obj$xi[, 1:2], col = ifelse(cluster == "Fall", 1, 2))
plot(rfpca.obj$xi[, 1:2], col = ifelse(cluster == "Fall", 1, 2))
points(rfpca.obj$xi[c(16,18,21,23,25,27,29,31,46), 1:2], col = 3)



## change direction
Ly <- lapply(1:n, function(i){
    if (cluster[i] == "Fall") {
        Ly[[i]][, ncol(Ly[[i]]):1]
    } else {
        Ly[[i]]
    }
})
Ly[[1]][, 51:1]

