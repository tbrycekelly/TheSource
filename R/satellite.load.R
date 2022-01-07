#' @title Read Satellite Data
#' @author Thomas Bryce Kelly
#' @description read
#' @export
## Read nc file and do preliminary parsing/conversion
read.satellite = function(file = NULL, verbose = F, lon = NULL, lat = NULL) {

  sate = list()

  ## Load and trim satellite data products
  for (i in 1:length(file)) {
    if (verbose) { message(Sys.time(), ': Loading satellite file: ', file[i]) }
    sate[[i]] = load.satellite(file = file[i], verbose = verbose)

    if (!is.null(lon) | !is.null(lat)) {
      sate[[i]] = trim.satellite(sate[[i]], lon = lon, lat = lat, verbose = verbose)
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
      if (length(sate[[i]]$field[[1]]) != length(sate[[1]]$field[[1]])) { message('  Number of entries do not match in files: ', file[1], ' and ', file[i], '. Code will likely fail now.') }
    }
  }
  if (verbose) { message(' Calculating field statistics...') }

  ## Calculate average field
  if (length(file) > 1) {
    for (k in 1:sate[[1]]$meta$n.fields) {
      n = matrix(0, nrow = nrow(sate[[1]]$field[[k]]), ncol = ncol(sate[[1]]$field[[k]]))
      res = matrix(0, nrow = nrow(sate[[1]]$field[[k]]), ncol = ncol(sate[[1]]$field[[k]]))
      times = list(start = 0, mid = 0, end = 0)

      ## Add together all fields
      for (i in 1:length(sate)) {
        temp = sate[[i]]$field[[k]]
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
      sate[[1]]$field[[k]] = res
    }
    ## Update times and filenames
    sate[[1]]$file = 'composite file'
    sate[[1]]$times$start = conv.time.unix(times$start / length(sate))
    sate[[1]]$times$mid = conv.time.unix(times$mid / length(sate))
    sate[[1]]$times$end = conv.time.unix(times$end / length(sate))

  }

  ## Update metadata
  meta = list(min = min(as.numeric(sate[[1]]$field[[1]]), na.rm = T),
              max = max(as.numeric(sate[[1]]$field[[1]]), na.rm = T),
              mean = mean(as.numeric(sate[[1]]$field[[1]]), na.rm = T),
              median = median(as.numeric(sate[[1]]$field[[1]]), na.rm = T),
              n.na = length(as.numeric(sate[[1]]$field[[1]])),
              n.na = sum(is.na(as.numeric(sate[[1]]$field[[1]]))),
              n.fields = sate[[1]]$meta$n.fields)

  sate[[1]]$meta = meta

  ## Return
  sate[[1]]
}


#' @title Read Satellite Data
#' @author Thomas Bryce Kelly
#' @description read
## Read nc file and do preliminary parsing/conversion
load.satellite = function(file, verbose = T) {

  ## Load data
  data = load.nc(file = file, verbose = verbose)
  dimnames = names(data)
  sizes = sapply(data, function(x) {max(cumprod(c(1,dim(x))), na.rm = T)})

  ## Determine primary field by size
  l = which.max(sizes)
  field.size = dim(data[[l]])
  l = which(sizes == sizes[l])
  x = data[l] ## list

  ## Get Lat/lon if possible
  lat = NULL; lon = NULL
  if ('lat' %in% dimnames) { lat = data[['lat']] }
  if ('latitude' %in% dimnames) { lat = data[['latitude']] }
  if ('LAT' %in% dimnames) { lat = data[['LAT']] }
  if ('LATITUDE' %in% dimnames) { lat = data[['LATITUDE']] }
  if ('Lat' %in% dimnames) { lat = data[['Lat']] }
  if ('Latitude' %in% dimnames) { lat = data[['Latitude']] }

  if ('lon' %in% dimnames) { lon = data[['lon']] }
  if ('LON' %in% dimnames) { lon = data[['LON']] }
  if ('LONGITUDE' %in% dimnames) { lon = data[['LONGITUDE']] }
  if ('longitude' %in% dimnames) { lon = data[['longitude']] }
  if ('Lon' %in% dimnames) { lon = data[['Lon']] }
  if ('Longitude' %in% dimnames) { lon = data[['Longitude']] }

  ## time to guess lat and lon
  if (is.null(lat) & is.null(lon)) {

    ## Standard NASA 4km grid
    if (4320 %in% field.size) {
      lat = seq(90, -90, length.out = dim(x)[2])
      lon = seq(-180, 180, length.out = dim(x)[1])
    }
    if (540 %in% field.size) {  # Standard 4k
      lat = seq(45, 30.03597, length.out = field.size[2])
      lon = seq(-140, -115.5454, length.out = field.size[1])
    }
  }

  meta = list(min = min(as.numeric(x[[1]]), na.rm = T),
              max = max(as.numeric(x[[1]]), na.rm = T),
              mean = mean(as.numeric(x[[1]]), na.rm = T),
              median = median(as.numeric(x[[1]]), na.rm = T),
              n = length(as.numeric(x[[1]])),
              n.na = sum(is.na(as.numeric(x[[1]]))),
              n.fields = length(l),
              trim = list(
                x = -1,
                y = -1,
                z = -1,
                t = -1
              ),
              meta = list(
                time = Sys.time(),
                Source.version = packageVersion('TheSource'),
                R.version = R.version.string))

  times = NULL
  if (is.null(times)) { times = get.satellite.times(file, verbose = verbose) }

  ## Return
  list(field = x,
       file = file,
       grid = function() {expand.grid(lon = lon, lat = lat)},
       lon = lon,
       lat = lat,
       times = times,
       meta = meta)
}



#' @title Read Satellite Data
#' @author Thomas Bryce Kelly
#' @description read
## Read nc file and do preliminary parsing/conversion
load.satellite.old = function(file, entry = 1, verbose = T) {

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
              n.fields = 1,
              meta = list(
                time = Sys.time(),
                Source.version = packageVersion('TheSource'),
                R.version = R.version.string))

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
trim.satellite = function(satellite, lon = NULL, lat = NULL, verbose = T) {
  prior = object.size(satellite)

  satellite$lon = satellite$lon %% 360
  if (is.null(lon)) {lon = range(satellite$lon)}
  if (is.null(lat)) {lat = range(satellite$lat)}

  anti = lon[2] < lon[1]
  lon = lon %% 360
  grid = expand.grid(lon = satellite$lon, lat = satellite$lat)

  if (verbose) { message('Trimming satellite domain... \n\tAntimeridian:\t', anti,
                         '\n\tLon:\t\t', paste(round(lon, 2), collapse = '\t'),
                         '\n\tLat:\t\t', paste(round(lat, 2), collapse = '\t')) }


  ## Replace lat and lonD
  satellite$lat = satellite$lat[satellite$lat >= lat[1] & satellite$lat <= lat[2]]

  if (anti) {
    satellite$lon = satellite$lon[satellite$lon >= lon[1] | satellite$lon <= lon[2]]
    l = which(grid$lon >= lon[1] | grid$lon <= lon[2]  &
                grid$lat >= lat[1] & grid$lat <= lat[2])
  } else {
    satellite$lon = satellite$lon[satellite$lon >= lon[1] & satellite$lon <= lon[2]]
    l = which(grid$lon >= lon[1] & grid$lon <= lon[2] &
                grid$lat >= lat[1] & grid$lat <= lat[2])
  }


  ## Replace field

  if (class(satellite$field) == 'list') {
    for (i in 1:length(satellite$field)) {
      satellite$field[[i]] = matrix(as.numeric(satellite$field[[i]])[l], ncol = length(satellite$lat), nrow = length(satellite$lon))
    }
  } else {
    satellite$field = matrix(satellite$field[l], ncol = length(satellite$lat), nrow = length(satellite$lon))
  }
  satellite$grid = expand.grid(lon = satellite$lon, lat = satellite$lat)

  if (verbose) { message(' Saved ', round(100 * (1 - object.size(satellite) / prior)), '%.') }

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

  all.var = names(file$var)
  for (name in names(file$dim)) {
    if (file$dim[[name]]$dimvarid$id > 0) {
      all.var = c(all.var, name)
    }
  }
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
