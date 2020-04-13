## Functions for working with remote sensing data
##
## Author: Thomas Bryce Kelly (tbk14 at fsu.edu)
## http://about.tkelly.org/
##
## Dept of Earth, Ocean & Atmospherical Sciences
## Florida State University
##
## Center for Ocean & Atmospheric Prediction Studies
## Florida State University
##
## National High Magnetic Field Laboratory
## Florida State University


###############################
#### Conversion functions #####
###############################

#' @title Convert Data to SST (Kahru)
#' @author Thomas Bryce Kelly
#' @description Convert
#' @param x The pixel values (or a vector or matrix)
#' @export
conv.kahru.sst = function(x) {
    ## Apply fix for signed vs unsigned int values
    if (any(x < 0)) {
        x[which(x < 0)] = x[which(x < 0)] + 256
    }
    ## Remove invalid data
    l = which(x <= 0 | x >= 255)
    x[l] = NA

    0.15 * x - 3
}

#' @title Convert Data to Chl (Kahru)
#' @author Thomas Bryce Kelly
#' @description Convert Kahru's data to mg Chl m-2
#' @param x Raw value
#' @export
conv.kahru.chl = function(x) {
    ## Apply fix for signed vs unsigned int values
    if (any(x < 0)) {
        x[which(x < 0)] = x[which(x < 0)] + 256
    }
    ## Remove invalid data
    l = which(x <= 0 | x >= 255)
    x[l] = NA

    10^(0.015 * x - 2)
}

#' @title Convert Data to NPP (Kahru)
#' @author Thomas Bryce Kelly
#' @description Converts kahru's data to npp
#' @param x raw value
#' @export
conv.kahru.npp = function(x) {
    x[x < 0] = NA
    x
}

#' @title Convert Data to Chl (NASA)
#' @author Thomas Bryce Kelly
#' @description convert
#' @param x raw value
#' @export
conv.nasa.chl = function(x) { 10 ^ x }


###########################
###### File SECTION #######
###########################

#' @title Read Satellite Data
#' @author Thomas Bryce Kelly
#' @description read
#' @param input.dir the input directory
#' @export
## Read nc file and do preliminary parsing/conversion
read.satellite = function(input.dir, file, conv = function(x){x}, entry = 1) {
    nc.file = nc_open(paste0(input.dir, file))

    x = ncvar_get(nc.file, varid = names(nc.file$var)[entry])
    attr = ncatt_get(nc = nc.file, varid = names(nc.file$var)[entry])

    ## Lat and Lon
    l.lat = which(names(nc.file$dim) == 'lat')
    l.lon = which(names(nc.file$dim) == 'lon')

    ## Determine Lat and Lon
    if (length(l.lat) > 0 & length(l.lon) > 0) {
        lat = ncvar_get(nc.file, varid = 'lat')
        lon = ncvar_get(nc.file, varid = 'lon')
    } else {

      ## Determine lat and lon based on presets
      if (dim(x)[1] == 4320) { ## Standard NASA 4km grid
        lat = seq(90, -90, length.out = dim(x)[2])
        lon = seq(-180, 180, length.out = dim(x)[1])
      }

      else if (dim(x)[1] == 540) {  # Standard 4k
        lat = seq(45, 30.03597, length.out = dim(x)[2])
        lon = seq(-140, -115.5454, length.out = dim(x)[1])
      }

      else if (dim(x)[1] == 588) {  ## Calcofi Fit
        lat = seq(37, 29.51349, length.out = dim(x)[2])
        lon = seq(-126.125, -116.6828, length.out = dim(x)[1])
      }

      else {  ## No match.
        warning('Remote sensing field size is not recognized.')
        lat = NULL
        lon = NULL
      }
    }
    nc_close(nc.file)
    x = conv(x)
    list(field = x, file = file, dir = input.dir, grid = function() {expand.grid(lon = lon, lat = lat)}, lon = lon, lat = lat,
         times = get.satellite.times(file), conv = conv, attr = function() {attr})
}


#' @title Tirm Satellite Field
#' @export
#' @param satellite a satellite object as returned from read.satellite()
#' @param lon a vector of length two indicating the limits desired. By default it is set to the range in longitude values in the satellite object.
#' @param lat same as above but for latitude
#'
trim.satellite = function(satellite, lon = NULL, lat = NULL) {
  grid = expand.grid(lon = satellite$lon, lat = satellite$lat)

  if (is.null(lon)) {lon = range(satellite$lon)}
  if (is.null(lat)) {lat = range(satellite$lat)}

  ## Replace lat and lon
  satellite$lat = satellite$lat[satellite$lat >= lat[1] & satellite$lat <= lat[2]]
  satellite$lon = satellite$lon[satellite$lon >= lon[1] & satellite$lon <= lon[2]]

  ## Replace field
  l = which(grid$lon >= lon[1] & grid$lon <= lon[2] &
              grid$lat >= lat[1] & grid$lat <= lat[2])

  satellite$field = matrix(satellite$field[l], ncol = length(satellite$lat), nrow = length(satellite$lon))
  satellite$grid = function() {expand.grid(lon = satellite$lon, lat = satellite$lat)}

  ## Return
  satellite
}

###########################
###### TIME SECTION #######
###########################


#' @title Get Satellite Times
#' @author Thomas Bryce Kelly
#' @description time
#' @param x raw
#' @export
get.satellite.times = function(x) {
    start = rep(conv.time.unix(1), length(x))
    mid = rep(conv.time.unix(1), length(x))
    end = rep(conv.time.unix(1), length(x))

    for (i in 1:length(x)) {
        start[i] = as.POSIXct(as.numeric(substr(x[i], 6, 8))*86400,
                                   origin = paste0(substr(x[i], 2, 5), '-01-01'), tz='GMT')

        timeframe = strsplit(x, '_')[[1]][2]
        end[i] = start[i]
        if (timeframe == 'DAY') {end[i] = start[i] + 1 * 86400}
        if (timeframe == '8D') {end[i] = start[i] + 8 * 86400}
        if (timeframe == 'MO') {end[i] = start[i] + 365/12 * 86400}
        if (timeframe == 'YR') {end[i] = start[i] + 365 * 86400}

        mid[i] = mean(c(start[i], end[i]))
    }

    list(start = start, mid = mid, end = end)
}


#' @title Print NetCDF
#' @import ncdf4
#' @export
print.nc = function(file) {
  file = nc_open(file)
  print(file)
  nc_close(file)
}


#' @title Load NetCDF
#' @import ncdf4
#' @export
load.nc = function(file) {
  file = nc_open(file)
  data = file
  nc_close(file)

  ##return
  data
}



###########################
###### Index SECTION #######
###########################

#' @title Rebuild Satellite Index
#' @author Thomas Bryce Kelly
#' @description rebuild
#' @param verbose Flag to print diagnostics
#' @export
rebuild.satellite.index = function(dir = 'D:/Data/Data_Satellite/', verbose = FALSE) {
    files = list.files(pattern = '.nc', path = paste0(dir, 'Raw/'))

    summary = data.frame(File = files, Sensor = NA, Datetime = NA, Year = NA,
                         Julian.Day = NA, Month = NA, Day = NA, Level = NA,
                         Timeframe = NA, Param = NA, Resolution = NA, stringsAsFactors = FALSE)

    for (i in 1:length(files)) {
        if (verbose) {cat(i)}

          prime = strsplit(files[i], split = '\\.')[[1]]
          second = strsplit(prime[2], split = '_')[[1]]
          datetime = substr(prime[1], 2, 8)

          summary$Sensor[i] = substr(prime[1], 1, 1)

          summary$Year[i] = as.numeric(substr(datetime, 1, 4))
          summary$Julian.Day[i] = as.numeric(substr(datetime, 5, 7))

          summary$Datetime[i] = paste0(as.POSIXct(summary$Julian.Day[i]*86400,
                                                    origin = paste0(summary$Year[i], '-01-01'), tz = 'GMT'))
          summary$Month[i] = as.numeric(substr(summary$Datetime[i], 6, 7))
          summary$Day[i] = as.numeric(substr(summary$Datetime[i], 9, 10))

          summary$Level[i] = second[1]
          summary$Timeframe[i] = second[2]
          summary$Resolution[i] = second[length(second)]
          if (length(second) > 4) {
            summary$Param[i] = paste(second[3], second[4], sep = '.')
          } else {
            summary$Param[i] = second[3]
          }

    }
    summary$Datetime = as.POSIXct(summary$Datetime, tz = 'GMT')

    if (verbose) {
        print(unique(summary$Year))
        print(unique(summary$Sensor))
        print(unique(summary$Level))
        print(unique(summary$Timeframe))
        print(unique(summary$Param))
        print(unique(summary$Resolution))
    }

    summary
}
