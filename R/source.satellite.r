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

###########################
###### File SECTION #######
###########################

#' @title Read Satellite Data
#' @author Thomas Bryce Kelly
#' @description read
#' @export
## Read nc file and do preliminary parsing/conversion
read.satellite = function(file = NULL, entry = 1,
                          verbose = F, lon = NULL, lat = NULL) {

  sate = list()

  ## Load and trim satellite data products
  for (i in 1:length(file)) {
    if (verbose) { message(Sys.time(), ': Loading satellite file: ', file[i]) }
    sate[[i]] = load.satellite(file = file[i], entry = entry, verbose = verbose)
    if (!is.null(lon) | !is.null(lat)) {
      sate[[i]] = trim.satellite(sate[[i]], lon = lon, lat = lat)
    }
  }

  if (verbose) { message(' Running checks on satellite files...') }

  if (length(file) > 1) {
    ## Run checks
    for (i in 2:length(sate)) {

      ## Lon check
      if (any(sate[[i]]$lon != sate[[1]]$lon)) { message('  Longitudes do not match in files: ', file[1], ' and ', file[i], '. Code will likely fail now.') }

      ## lat check
      if (any(sate[[i]]$lat != sate[[1]]$lat)) { message('  Latitudes do not match in files: ', file[1], ' and ', file[i], '. Code will likely fail now.') }

      ## lat check
      if (length(sate[[i]]$field) != length(sate[[1]]$field)) { message('  Number of entries do not match in files: ', file[1], ' and ', file[i], '. Code will likely fail now.') }
    }
  }
  if (verbose) { message(' Calculating field statistics...') }

  ## Calculate average field
  if (length(file) > 1) {
    n = matrix(0, nrow = nrow(sate[[1]]$field), ncol = ncol(sate[[1]]$field))
    res = matrix(0, nrow = nrow(sate[[1]]$field), ncol = ncol(sate[[1]]$field))
    times = list(start = 0, mid = 0, end = 0)

    ## Add together all fields
    for (i in 1:length(sate)) {
      temp = sate[[i]]$field
      n = n + !is.na(temp)
      temp[is.na(temp)] = 0
      res = res + temp

      times$start = times$start + as.numeric(sate[[i]]$times$start)
      times$mid = times$mid + as.numeric(sate[[i]]$times$mid)
      times$end = times$end + as.numeric(sate[[i]]$times$end)
      sate[[1]]$file = c(sate[[1]]$file, sate[[i]]$file)
    }

    ## Calculate mean field value
    res = res / n
    res[!is.finite(res)] = NA
    sate[[1]]$field = res

    ## Update times and filenames
    sate[[1]]$file = 'composite file'
    sate[[1]]$times$start = conv.time.unix(times$start / length(sate))
    sate[[1]]$times$mid = conv.time.unix(times$mid / length(sate))
    sate[[1]]$times$end = conv.time.unix(times$end / length(sate))
  }

  ## Update metadata
  meta = list(min = min(as.numeric(sate[[1]]$field), na.rm = T),
              max = max(as.numeric(sate[[1]]$field), na.rm = T),
              mean = mean(as.numeric(sate[[1]]$field), na.rm = T),
              median = median(as.numeric(sate[[1]]$field), na.rm = T),
              n.na = sum(is.na(as.numeric(sate[[1]]$field))),
              n.fields = length(sate))
  sate[[1]]$meta = meta

  ## Return
  sate[[1]]
}


#' @title Read Satellite Data
#' @author Thomas Bryce Kelly
#' @description read
## Read nc file and do preliminary parsing/conversion
load.satellite = function(file, entry = 1, verbose = T) {

  nc.file = ncdf4::nc_open(file)

  x = ncdf4::ncvar_get(nc.file, varid = names(nc.file$var)[entry])

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

    ### Kahru data aproducts:
    else if (dim(x)[1] == 540) {  # Standard 4k
      lat = seq(45, 30.03597, length.out = dim(x)[2])
      lon = seq(-140, -115.5454, length.out = dim(x)[1])
    }

    else if (dim(x)[1] == 588) {  ## Calcofi Fit
      lat = seq(37, 29.51349, length.out = dim(x)[2])
      lon = seq(-126.125, -116.6828, length.out = dim(x)[1])
    }

    else {  ## No match.
      message(' Remote sensing field size is not recognized! Lat and Lon will not be available. Downstream code will likely fail.')
      lat = NULL
      lon = NULL
    }
  }
  ncdf4::nc_close(nc.file)

  meta = list(min = min(as.numeric(x), na.rm = T),
              max = max(as.numeric(x), na.rm = T),
              mean = mean(as.numeric(x), na.rm = T),
              median = median(as.numeric(x), na.rm = T),
              n.na = sum(is.na(as.numeric(x))),
              n.fields = 1)

  ## Return
  list(field = x,
       file = file,
       dir = file,
       grid = function() {expand.grid(lon = lon, lat = lat)},
       lon = lon,
       lat = lat,
       times = get.satellite.times(file, verbose = verbose),
       meta = meta)
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
  satellite$grid = expand.grid(lon = satellite$lon, lat = satellite$lat)

  ## Return
  satellite
}


load.sss = function(file, lon = NULL, lat = NULL) {
  ss = load.nc(file)
  ss$sss = ss$sss * ss$sss_qc ## Apply quality control

  if (!is.null(lon)) {
    l = which(ss$lon >= lon[1] & ss$lon <= lon[2])
    ss$sss = ss$sss[l,]
  }
  if (!is.null(lat)) {
    l = which(ss$lat >= lat[1] & ss$lat <= lat[2])
    ss$sss = ss$sss[,l]
  }
  ss$grid = expand.grid(lon = ss$lon, lat = ss$lat)

  ## Return
  list(lon = ss$lon, lat = ss$lat, time = ss$time, field = ss$sss)
}

###########################
###### TIME SECTION #######
###########################


#' @title Get Satellite Times
#' @author Thomas Bryce Kelly
#' @description time
#' @param x raw
#' @export
get.satellite.times = function(x, verbose = T) {
  start = rep(make.time(), length(x))
  mid = rep(make.time(1), length(x))
  end = rep(make.time(1), length(x))

  for (i in 1:length(x)) {
    temp = strsplit(x[i], split = '/')[[1]]
    x[i] = temp[length(temp)]
    day = as.numeric(substr(x[i], 6, 8))
    year = as.numeric(substr(x[i], 2, 5))

    if (verbose) { message(' Processing satellite times for file ', i, ' of ', length(x), ':\tday = ', day, '\tyear = ', year)}
    if (!is.na(day) & !is.na(year)) {
      start[i] = as.POSIXct(day * 86400,
                            origin = paste0(year, '-01-01'),
                            tz = 'GMT')

      timeframe = strsplit(x, '_')[[1]][2]
      end[i] = start[i]
      if (timeframe == 'DAY') {
        end[i] = start[i] + 1 * 86400
      }
      if (timeframe == '8D') {
        end[i] = start[i] + 8 * 86400
      }
      if (timeframe == 'MO') {
        end[i] = start[i] + 365 / 12 * 86400
      }
      if (timeframe == 'YR') {
        end[i] = start[i] + 365 * 86400
      }
      mid[i] = mean(c(start[i], end[i]))
    } else {
      start[i] = NA
      mid[i] = NA
      end[i] = NA
    }
  }

  ## Return
  list(start = start, mid = mid, end = end)
}


#' @title Load NetCDF
#' @import ncdf4
#' @param file The location of a .nc file.
#' @export
load.nc = function(file, var = NULL, test = F, verbose = T) {

  if (verbose) { message('Attempting to load data from ', file)}
  file = ncdf4::nc_open(file)
  if (verbose) { message(' File openned.')}

  all.var = c(names(file$var), names(file$dim))
  all.var = all.var[all.var != 'nv']

  if (is.null(var)) {
    var = all.var
    if (verbose) { message(' No variables specified, loading all ', length(var), ' entires.') }
  }

  ## Allow number of variable to be used as well
  if (is.numeric(var)) {
    var = names(file$var)[var]
  }

  ## Load data
  data = list()

  for (v in var) {

    if (v %in% all.var) {
      a = Sys.time()
      data[[v]] = ncdf4::ncvar_get(nc = file, varid = v)
      if (test) { data[[v]] = data[[v]][1] }
      if (verbose) {
        d = dim(data[[v]])
        if (is.null(d)) { d = length(data[[v]])}
        message(' Variable ', v, ' (', paste(d, collapse = 'x'),') loaded from file. \t', round(difftime(Sys.time(), a, units = 'secs')), 's')
      }
    } else {
      a = Sys.time()
      data[[v]] = NA
      if (verbose) {
        message(' Variable ', v, ' (X) does not exist in file. \t', round(difftime(Sys.time(), a, units = 'secs')), 's')
      }
    }
  }

  if (verbose) {message(' Closing file. Finished.')}
  ncdf4::nc_close(file)

  ##return
  data
}


#' @title Read NetCDF Variables
#' @import ncdf4
#' @param file The location of a .nc file.
#' @export
load.nc.vars = function(file, verbose = F) {

  if (verbose) { message('Attempting to load variables from ', file)}
  file = ncdf4::nc_open(file)
  if (verbose) { message(' File openned.')}

  var = c(names(file$var), names(file$dim))

  if (verbose) {message(' Closing file. Finished.')}
  ncdf4::nc_close(file)

  ##return
  var
}



###########################
###### Index SECTION #######
###########################

#' @title Rebuild Satellite Index
#' @author Thomas Bryce Kelly
#' @description rebuild
#' @param verbose Flag to print diagnostics
#' @export
rebuild.satellite.index = function(dir = 'Z:/Data/Data_Satellite/Raw', verbose = FALSE) {

  #### Find files
  if (verbose) { message(Sys.time(), ': Looking for available satellite files...')}

  # List files
  a = Sys.time()
  files = list.files(pattern = '.nc', path = dir, recursive = F)
  b = Sys.time()

  if (verbose) { message(Sys.time(), ': Found ', length(files), ' files! Search took ', round(as.numeric(b) - as.numeric(a)), ' seconds.')}


  #### Populate data fields
  summary = data.frame(File = c(files, NA), Sensor = NA, Datetime = NA, Year = NA,
                       Julian.Day = NA, Month = NA, Day = NA, Level = NA,
                       Timeframe = NA, Param = NA, Resolution = NA, stringsAsFactors = FALSE)

  if (length(files) > 0) {
    for (i in 1:length(files)) {
      if (file.size(paste0(dir, files[i])) > 0) {
        if (verbose) { message(Sys.time(), ': Attempting to scrub metadata from ', files[i], ' (', i, '/', length(files), ')... ', appendLF = F) }

        summary$File[i] = paste0(dir, files[i])

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
          summary$Param[i] = paste(second[-c(1:2)], collapse = '.')
        } else {
          summary$Param[i] = second[3]
        }

        if (verbose) { message('Success.') }
      }
    }
  }

  ## Remove bad files or folders
  summary = summary[!is.na(summary$Param),]

  c = Sys.time()

  #### Finish up
  summary$Datetime = as.POSIXct(summary$Datetime, tz = 'GMT')

  ## Recursively delve into subfolders:
  dirs = list.dirs(path = dir, recursive = F, full.names = T)
  if (length(dirs) > 0) {
    for (i in 1:length(dirs)) {
      if (verbose) { message(' Delving into ', dirs[i])}
      summary = rbind(summary, rebuild.satellite.index(dir = paste0(dirs[i], '/'), verbose = verbose))
    }
  }
  ## Return
  summary
}


###############################
#### Conversion functions #####
###############################

## Depreciated

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

