#' @title Read Satellite Data
#' @author Thomas Bryce Kelly
#' @description read
## Read nc file and do preliminary parsing/conversion
read.satellite = function(file, lon = NULL, lat = NULL, verbose = T) {

  if (!is.null(lon) | !is.null(lat)) {
    pre.trimmed = T
    if (is.null(lon)) { lon = c(-180, 180)}
    if (is.null(lat)) { lat = c(-90,90)}
    
    # info = load.nc.dims(file)
    # 
    # ## See if lat/lon exists in the file
    # lat.check = c('lat', 'latitude', 'Lat', 'Latitude', 'LAT', 'LATITUDE', 'y', 'Y', 'fakeDim1') %in% names(info$dim)
    # lon.check = c('lon', 'longitude', 'Lon', 'Longitude', 'LON', 'LONGITUDE', 'x', 'X', 'fakeDim0') %in% names(file.con$dim)
    # 
    # if (any(lat.check) & any(lon.check)) {
    #   
    #   ## Valid lon/lat
    #   file.lons = info$dim[[which(lon.check)[1]]] %% 360
    #   file.lats = info$dim[[which(lat.check)[1]]]
    #   
    #   k = which(file.lons <= lon[2] & file.lons >= lon[1])
    #   
    #   ## Determine subset of lon
    #   if (max(diff(k.lon)) > 1) {
    #     lon.start = 1
    #     lon.count = length(file.lons)
    #   } else {
    #     file.lons = file.lons[k]
    #     lon.start = k[1]
    #     lon.count = length(k)
    #   }
    #   
    #   k = which(file.lats <= lat[2] & file.lats >= lat[1])
    #   
    #   ## Determine subset of lat
    #   if (max(diff(k.lat)) > 1) {
    #     lat.start = 1
    #     lat.count = length(file.lats)
    #   } else {
    #     file.lats = file.lats[k]
    #     lat.start = k[1]
    #     lat.count = length(k)
    #   }
    #   
    # } else {
    #   ## Lat/lon not found reliably enough...
    #   lon.count = -1
    #   lon.start = 1
    #   lat.start = 1
    #   lat.count = -1
    #   message('No coordiantes found.')
    # }
    # 
    # 
    
    ## Can we load only subsets of the data?
    file.con = ncdf4::nc_open(file)
    
    lat.check = c('lat', 'latitude', 'Lat', 'Latitude', 'LAT', 'LATITUDE') %in% names(file.con$dim)
    
    if (any(lat.check)) {
      ## Get available latitudes
      file.lats = file.con$dim[[which(lat.check)[1]]]$vals
      k = which(file.lats <= lat[2] & file.lats >= lat[1])
      
      if (max(diff(k)) == 1) {
        file.lats = file.lats[k]
        lat.count = k[length(k)] - k[1] + 1
        lat.start = k[1]
      } else {
        lat.count = -1
        lat.start = 1
      }
    } else {
      lat.count = -1
      lat.start = 1
      file.lats = NULL
    }
    
    lon.check = c('lon', 'longitude', 'Lon', 'Longitude', 'LON', 'LONGITUDE') %in% names(file.con$dim)
    
    if (any(lon.check)) {
      file.lons = file.con$dim$lon$vals
      k = which(file.lons <= lon[2] & file.lons >= lon[1])
      
      if (max(diff(k)) == 1) {
        file.lons = file.lons[k]
        lon.count = k[length(k)] - k[1] + 1
        lon.start = k[1]
      } else {
        lon.count = -1
        lon.start = 1
      }
    } else {
      lon.count = -1
      lon.start = 1
      file.lons = NULL
    }
    
    data = list()
    data[[names(file.con$var)[1]]] = ncdf4::ncvar_get(file.con, names(file.con$var)[1], start = c(lon.start, lat.start), count = c(lon.count, lat.count))
    data[['lon']] = file.lons
    data[['lat']] = file.lats
    
    ncdf4::nc_close(file.con)
      
  } else{
  
    ## Load data
    data = load.nc(file = file, verbose = verbose)
    pre.trimmed = F
  }
  
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
      lat = seq(90, -90, length.out = field.size[2]+1)[-1]
      lat = lat - diff(lat)[1]
      lon = seq(-180, 180, length.out = field.size[1]+1)[-1]
      lon = lon - diff(lon)[1]
    } else if (540 %in% field.size) {  # Mati Kahru's output grid
      lat = seq(45, 30.03597, length.out = field.size[2])
      lon = seq(-140, -115.5454, length.out = field.size[1])
    } else if (588 %in% field.size) {  ## mati Kahru's Calcofi Fit
        lat = seq(37, 29.51349, length.out = field.size[2])
        lon = seq(-126.125, -116.6828, length.out = field.size[1])
    } else {
      if (verbose) { message(' No known grids found, assuming global extent.')}
      lat = seq(90, -90, length.out = field.size[2]+1)[-1]
      lat = lat - diff(lat)[1]
      lon = seq(-180, 180, length.out = field.size[1]+1)[-1]
      lon = lon - diff(lon)[1]
    }
  }
  
  meta = list(n = length(x[[1]]),
              n.na = sum(is.na(x[[1]])),
              n.fields = length(l),
              trim = list(
                lon = NA,
                lat = NA
              ),
              meta = list(
                time = Sys.time(),
                Source.version = packageVersion('TheSource'),
                R.version = R.version.string))

  times = NULL
  if (is.null(times)) { times = get.satellite.times(file, verbose = verbose) }

  sate = list(field = x,
              file = file,
              grid = function() {expand.grid(lon = lon, lat = lat)},
              lon = lon,
              lat = lat,
              times = times,
              meta = meta)
  
  ## Trim
  if (!pre.trimmed & (!is.null(lon) | !is.null(lat))) {
    sate = trim.satellite(sate, lon = lon, lat = lat, verbose = verbose)
  }
  
  ## Return
  sate
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
  satellite$grid = function() {expand.grid(lon = satellite$lon, lat = satellite$lat)}
  satellite$meta$trim$lon = lon
  satellite$meta$trim$lat = lat

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
  mid = rep(make.time(), length(x))
  end = rep(make.time(), length(x))

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


#' @title Load NetCDF
#' @import ncdf4
#' @param file The location of a .nc file.
#' @export
load.nc.subset = function(file, var = NULL, x = NULL, y = NULL, test = F, verbose = T) {
  
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
      data[[v]] = ncdf4::ncvar_get(nc = file, varid = v, )
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


#' @title Read NetCDF Variables
#' @import ncdf4
#' @param file The location of a .nc file.
#' @export
load.nc.dims = function(file, verbose = F) {
  
  if (verbose) { message('Attempting to load dims from ', file)}
  file.con = ncdf4::nc_open(file)
  if (verbose) { message(' File openned.')}
  
  dim = list()
  var = names(file.con$dim)
  for (v in var) {
    dim[[v]] = file.con$dim[[v]]$vals
  }
  
  vars = list()
  var = names(file.con$var)
  
  for (v in var) {
    vars[[v]] = file.con$var[[v]]$size
    names(vars[[v]]) = file.con$var[[v]]$dimids + 1
  }
  
  if (verbose) {message(' Closing file. Finished.')}
  ncdf4::nc_close(file.con)
  
  ##return
  list(dim = dim, var = vars)
}
