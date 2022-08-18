library(tidync)

lons = c(-20, 20)
lats = c(30, 50)
verbose = T
{
  a = Sys.time()
  file = tidync::tidync(x = 'Z:/Data/Satellite/Raw/2012/V/V2012025.L4_8D_CBPM_4km.nc')
  
  ## Get dimension information
  ndim = length(file$transforms)
  dims = list()
  
  for (i in 1:ndim) {
    d.name = names(file$transforms)[i]
    ## Set val
    dims[[d.name]] = file$transforms[[i]][[1]]
    
    ## Determine if we know lat/lon values?
    
    if (names(file$transforms)[i] %in% c('lat', 'LAT', 'Lat', 'latitude', 'LATITUDE', 'Latitude')) {
      if (verbose) { message(' Identified Latitude axis.')}
      d.name = 'lat'
      
    } else if (names(file$transforms)[i] %in% c('lon', 'LON', 'Lon', 'longitude', 'LONGITUDE', 'Longitude')) {
      if (verbose) { message(' Identified Longitude axis.')}
      d.name = 'lon'
      
    } else if (length(file$transforms[[i]][[1]]) == 4320) {
      if (verbose) { message(' Identified NASA Longitude axis.')}
      d.name = 'lon'
      dims[[d.name]] = seq(-180, 180, length.out = 4321)[-1]
      dims[[d.name]] = dims[[d.name]] - diff(dims[[d.name]])[1]
      
    } else if (length(file$transforms[[i]][[1]]) == 2160) {
      if (verbose) { message(' Identified NASA Latitude axis.')}
      d.name = 'lat'
      dims[[d.name]] = seq(-90, 90, length.out = 2161)[-1]
      dims[[d.name]] = dims[[d.name]] - diff(dims[[d.name]])[1]
      
    } 
    
    ## Set subsetting if able:
    if (d.name == 'lon' & !is.null(lons)) {
      if (verbose) { message(' Filtering by Longitude range.')}
      file$transforms[[i]]$selected = F
      file$transforms[[i]]$selected[dims[[d.name]] >= lons[1] & dims[[d.name]] <= lons[2]] = T
      dims[[d.name]] = dims[[d.name]][file$transforms[[i]]$selected]
    }
    
    if (d.name == 'lat' & !is.null(lats)) {
      if (verbose) { message(' Filtering by Latitude range.')}
      file$transforms[[i]]$selected = F
      file$transforms[[i]]$selected[dims[[d.name]] >= lats[1] & dims[[d.name]] <= lats[2]] = T
      dims[[d.name]] = dims[[d.name]][file$transforms[[i]]$selected]
    }
    
    
  }
  
  ## Load data
  temp = as.data.frame(tidync::hyper_tibble(file))
  
  message(Sys.time() - a)
}



#' @title Read Satellite Data
#' @author Thomas Bryce Kelly
#' @description read
## Read nc file and do preliminary parsing/conversion
read.satellite2 = function(file, lon = NULL, lat = NULL, verbose = T) {
  
  file = tidync::tidync(file)
  
  ## Get dimension information
  ndim = length(file$transforms)
  dims = list()
  
  for (i in 1:ndim) {
    d.name = names(file$transforms)[i]
    ## Set val
    dims[[d.name]] = file$transforms[[i]][[1]]
    
    ## Determine if we know lat/lon values?
    
    if (names(file$transforms)[i] %in% c('lat', 'LAT', 'Lat', 'latitude', 'LATITUDE', 'Latitude')) {
      if (verbose) { message(' Identified Latitude axis.')}
      d.name = 'lat'
      
    } else if (names(file$transforms)[i] %in% c('lon', 'LON', 'Lon', 'longitude', 'LONGITUDE', 'Longitude')) {
      if (verbose) { message(' Identified Longitude axis.')}
      d.name = 'lon'
      
    } else if (length(file$transforms[[i]][[1]]) == 4320) {
      if (verbose) { message(' Identified NASA Longitude axis.')}
      d.name = 'lon'
      dims[[d.name]] = seq(-180, 180, length.out = 4321)[-1]
      dims[[d.name]] = dims[[d.name]] - diff(dims[[d.name]])[1]
      
    } else if (length(file$transforms[[i]][[1]]) == 2160) {
      if (verbose) { message(' Identified NASA Latitude axis.')}
      d.name = 'lat'
      dims[[d.name]] = seq(-90, 90, length.out = 2161)[-1]
      dims[[d.name]] = dims[[d.name]] - diff(dims[[d.name]])[1]
      
    } else {
      if (verbose) { message(' No standard axis was found for ', d.name, '.')}
    }
    
    ## Set subsetting if able:
    if (d.name == 'lon' & !is.null(lon)) {
      if (verbose) { message(' Filtering by Longitude range.')}
      file$transforms[[i]]$selected = F
      file$transforms[[i]]$selected[dims[[d.name]] >= lon[1] & dims[[d.name]] <= lon[2]] = T
      dims[[d.name]] = dims[[d.name]][file$transforms[[i]]$selected]
    }
    
    if (d.name == 'lat' & !is.null(lat)) {
      if (verbose) { message(' Filtering by Latitude range.')}
      file$transforms[[i]]$selected = F
      file$transforms[[i]]$selected[dims[[d.name]] >= lat[1] & dims[[d.name]] <= lat[2]] = T
      dims[[d.name]] = dims[[d.name]][file$transforms[[i]]$selected]
    }
  }
  
  ## Load data
  for (n in names(file$transforms)) {
    if (any(file$variable$name == n)) {
      file$variable$active[file$variable$name == n] = F
    }
  }
  #vars = file$variable$name[file$variable$ndims == max(file$variable$ndims)]
  #nvars = length(vars)
  #var = list()
  
  var = as.data.frame(tidync::hyper_tibble(file))
  var = var[,!names(var) %in% names(file$transforms)]
  
  
  meta = list(n = nrow(var),
              n.fields = ncol(var),
              trim = list(
                lon = lon,
                lat = lat
              ),
              meta = list(
                time = Sys.time(),
                Source.version = packageVersion('TheSource'),
                R.version = R.version.string))
  
  times = NULL
  if (is.null(times)) { times = get.satellite.times(file$source$source, verbose = verbose) }
  
  ## Repackage vars
  var.names = colnames(var)
  var = lapply(1:ncol(var), function(x) {array(var[,x], dim = c(length(dims[[1]]), length(dims[[2]])))})
  names(var) = var.names
  
  sate = list(field = var,
              file = file$source$source,
              grid = function() {expand.grid(lon = dims[[1]], lat = dims[[2]])},
              lon = dims[[1]],
              lat = dims[[2]],
              times = times,
              meta = meta)
  
  ## Return
  sate
}

{
  a = Sys.time()
  temp = read.satellite2(path, lon = c(-120, -80), lat = c(50, 100))
  message(Sys.time() - a)
}

{
  a = Sys.time()
  temp = read.satellite(path)
  message(Sys.time() - a)
}


