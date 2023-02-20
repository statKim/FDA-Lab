#####################################################################
### Access the Impala shell OpenSky Network
### - The functions are refered from R package "osn".
#####################################################################

### Converte "YYYY-MM-DD hh:mm:ss" to UNIX timestamp
date2unixtime <- function(YYYYMMDD_hhmmss) {
    # wef <- "2021-02-05 09:00:00"
    wef <- lubridate::as_datetime(YYYYMMDD_hhmmss)
    wef <- wef %>% as.integer()
    # floor to POSIX hour
    wefh <- wef - (wef %% 3600)
    # floor to POSIX day
    wefd <- wefh - (wefh %% 86400)
    
    return(wefd)
}

### Converte UNIX timestamp to "YYYY-MM-DD hh:mm:ss"
unixtime2date <- function(unixtime) {
    return( as.POSIXct(unixtime, origin = "1970-01-01") )
}


### Run the SQL query and get data from Impala shell
impala_query <- function(session, query) {
    # impala_query <- function(session, query, cols) {
    stopifnot(class(session) == "ssh_session")
    # stopifnot(!is.null(cols))
    lines <- ssh::ssh_exec_internal(
        session,
        stringr::str_glue("-q {query}", query = query)) %>%
        { rawToChar(.$stdout) }
    if (logger::log_threshold() == logger::TRACE) {
        lines %>%
            readr::write_lines("query_output.txt")
    }
    lines <- lines %>%
        parse_impala_query_output()
    
    # make a 1 line data so to return an empty tibble in case of empty Impala result
    if (length(lines) == 0) {
        stop("There is no available data!")
        # lines <- paste0(paste(names(cols$cols), collapse = "|"),
        #                                     "\n")
    }
    
    I(lines) %>%
        readr::read_delim(
            # col_types = cols,
            delim = "|",
            na = c("", "NULL"),
            trim_ws = TRUE
        )
}


#' Create an ssh session to OpenSky Network’s Impala shell.
#'
#' @param usr     user account
#' @param port    port to connect to
#' @inheritParams ssh::ssh_connect
#'
#' @return an SSH session
#' @export
#'
#' @examples
#' \dontrun{
#' # connect directly to OSN
#' session <- osn_connect("cucu", verbose = 2)
#'
#' # connect via SSH port forwarding
#' session <- osn_connect_osn(
#'   usr = Sys.getenv("OSN_USER"),
#'   passwd = Sys.getenv("OSN_PASSWORD"),
#'   port = 6666,
#'   host = "localhost"
#' )
#' }
osn_connect <- function(usr, passwd = askpass::askpass,
                        host = "data.opensky-network.org", port = 2230,
                        verbose = FALSE) {
    fullhost <- stringr::str_glue("{usr}@{host}:{port}")
    ssh::ssh_connect(fullhost, passwd = passwd, verbose = verbose)
}

#' Disconnect from OpenSky Network’s Impala shell.
#'
#' @inheritParams ssh::ssh_disconnect
#'
#' @return an SSH session
#' @export
#'
#' @examples
#' \dontrun{
#' session <- osn_connect("cucu", verbose = 2)
#' osn_disconnect(session)
#' }
osn_disconnect <- function(session) {
    ssh::ssh_disconnect(session)
}


### Utill function to make data from database to tibble
parse_impala_query_output <- function(lines) {
    lines %>%
        stringi::stri_split_lines() %>%
        purrr::flatten_chr() %>%
        # remove empty lines
        stringr::str_subset(pattern = "^$", negate = TRUE) %>%
        # remove delimiting lines
        stringr::str_subset(pattern = "^\\+-", negate = TRUE) %>%
        # remove blanks
        stringr::str_replace_all(pattern = "[ ][ ]*", "") %>%
        # remove leading/last '|'
        stringr::str_replace_all("^[|](.+)[|]$", "\\1") %>%
        # remove duplicated lines, i.e. repeated column names header
        unique()
}



