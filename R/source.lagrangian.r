## A general Lagrangian package for odeling particles in 2 and 3D.
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
##
## College of Fisheries and Oceanographic Science
## University of Alaska Fairbanks

#' @title Initialize Lagrangian Model
#' @param t.start start time for the Lagrangian model (POSIX)
#' @param t.end end time for the Lagrangian model (POSIX)
#' @param t.step time between internal model steps (sec)
#' @param save.freq number of time steps between particle saves
#' @param nparticles maximum number of particles considered
#' @param verbose output flag for setup data
#' @author Thomas Bryce Kelly
#' @export
init.lagrangian.model = function(t.start,
                                 t.end,
                                 t.step,
                                 save.freq,
                                 nparticles,
                                 verbose = T) {
  if (verbose) {message(Sys.time(), ': Setting up Lagrangian model meta parameters (e.g. timestepping values, number of save points, etc)')}

  meta = list(t.start = t.start,
       t.end = t.end,
       t.step = t.step,
       nstep = abs(floor(as.numeric(difftime(t.end, t.start, units = 'sec')) / t.step)),
       save.freq = save.freq,
       nsave = abs(ceiling(as.numeric(ceiling(difftime(t.end, t.start, units = 'sec')) / t.step) / save.freq)) + 1,
       nparticles = nparticles,
       params = c('lon', 'lat', 'depth', 'u', 'v', 'w', 'T', 'S', 'Alive'))

  hist = init.lagrangian.history(meta, verbose = verbose)

  list(meta = meta, hist = hist)
}


#' @title Initialize Lagrangian History
#' @param meta list object containing lagrangian model setup (see init.lagrangian.model())
#' @author Thomas Bryce Kelly
init.lagrangian.history = function(meta,
                                   extra.tracers = NULL,
                                   verbose = T) {
  ## 1) lon, 2) lat, 3) depth, 4) u, 5) v, 6) w, 7) T, 8) S, 9) Alive
  hist = array(NA, dim = c(meta$nparticles, meta$nsave, 9), dimnames = list(x = c(1:meta$nparticles),
                                                                            y = c(1:meta$nsave),
                                                                            z = c('lon', 'lat', 'depth', 'u', 'v', 'w', 'T', 'S', 'Alive')))

  if (!is.null(extra.tracers)) {
    if (class(extra.tracers) == 'data.frame') { ## Incomplete
      for (i in names(extra.tracers)) {
        hist[[i]] = NA
      }
    } else {
      for (i in extra.tracers) {
        hist[[i]] = NA
      }
    }
  }
  if (verbose) {
    message(Sys.time(), ': Initializing Lagrangian history file with z array indexes:\n 1) lon,\n 2) lat,\n 3) depth,\n 4) u,\n 5) v,\n 6) w,\n 7) T,\n 8) S,\n 9) Alive')
    message('Current history file has x = ', meta$nparticles, ' particles and y = ', meta$nsave, ' snapshots.')
  }
  hist
}


#' @title Initiate Lagrangian Particles
#' @author Thomas Bryce Kelly
#' @param lon Longitude values of the lagrangian floats
#' @param lat Latitude values of the lagrangian floats
#' @param depth Depths of the lagrangian floats
#' @param time POSIX datetime of lagrangian floats
#' @param sink Sinking velocity (m s-1) of the lagrangian float
#' @param const.depth Depth (m) at which to maintain regardless of vertical velocity
#' @param density Density (kg m-3) of the float in order to maintain isopycnal position
#' @param verbose Verbose flag, default = T
#' @export
init.lagrangian.particles = function(lon = NA,
                                     lat = NA,
                                     depth = NA,
                                     time = NA,
                                     sink = NA,
                                     const.depth = NA,
                                     density = NA,
                                     verbose = T) {

  if (verbose) {
    message('Initializing Langraign Particle Data Structure...')
    message('This function is helpful in organizing the initial data structure in order to be compatible for the advection scripts in this library. Optionally you can provide vectors (one entry per particle).')
    message('Optional Entries:')
    message(' sink:\t\t the sinking speed of the particle (m/s).')
    message(' const.depth:\t a depth for the particle to remain at (m)')
    message(' density:\t an isopycnal density for the float to stay at (kg m-3) -- Not yet implemented!')
    message('Note:\t\t Longitudes will be coerced to 0-360 for compatibility.\n')
  }

  a = data.frame(ID = c(1:length(lon)), lon = lon, lat = lat, depth = depth, start.time = time, sink = sink, const.depth = const.depth, density = density,
                 current.lon = lon, current.lat = lat, current.depth = depth, current.u = 0, current.v = 0, current.w = 0, alive = F)

  if (verbose){message('\nReturning particle file with ', nrow(a), ' particles.')}
  ## Return
  a
}


##############################
#### Loading Functions #######
##############################

#' @title Load Advection (MITgcm)
#' @author Thomas Bryce Kelly
#' @description A function used to load an mitgcm file set (offline data fields) and structure in a manner suitable for Lagrangian particle release.
load.advection.mitgcm = function(u.file, v.file, w.file = NULL) {

  ## U
  finfo = file.info(u.file)
  uin = file(u.file, 'rb')
  U = readbin(uin, double(), size = 4, n = finfo$size, endian = 'big')
  close(uin)

  ## V
  finfo = file.info(v.file)
  vin = file(v.file, 'rb')
  V = readbin(vin, double(), size = 4, n = finfo$size, endian = 'big')
  close(vin)

  ## W
  if (!is.null(w.file)) {
    finfo = file.info(w.file)
    win = file(w.file, 'rb')
    W = readbin(uwin, double(), size = 4, n = finfo$size, endian = 'big')
    close(win)
  } else {
    W = 0
  }

  ## TODO Add Lon, LAt, TIME, grid, DEPTH, etc
  list(u = U, v = V, w = W)
}



#' @title Load Advection (OSCAR)
#' @author Thomas Bryce Kelly
#' @param file OSCAR file (includign path) to read in
#' @param lats The latitude limits of the resulting model domain
#' @param lons The longitude of the resulting model domain
#' @param verbose A verbose flag
#' @description A function used to load an OSCAR circulation netcdf file and structure in a manner suitable for Lagrangian particle release.
#' @export
load.advection.oscar = function(file, lats = NULL, lons = NULL, get.vel = get.vel.oscar, verbose = F) {

  if (verbose) {
    message('Loading OSCAR Advection Product\n Loading OSCAR input file ', file, ' (', format(file.size(file)/2^20, digits = 3), ' MB)... ', appendLF = F)
    a = Sys.time()
  }

  f = ncdf4::nc_open(file)
  oscar = f
  u = ncdf4::ncvar_get(f, 'u')
  v = ncdf4::ncvar_get(f, 'v')
  ncdf4::nc_close(f)
  if (verbose) {message('Done.')}

  ## Set names
  times = as.POSIXct(oscar$var$u$dim[[4]]$vals * 86400, origin = '1992-10-05 00:00:00', tz = 'GMT')
  dimnames(u) = list(longitude = oscar$var$u$dim[[1]]$vals, latitude = oscar$var$u$dim[[2]]$vals, time = times)
  dimnames(v) = list(longitude = oscar$var$u$dim[[1]]$vals, latitude = oscar$var$u$dim[[2]]$vals, time = times)

  oscar = list(u = u, v = v, lon = oscar$var$u$dim[[1]]$vals, lat = oscar$var$u$dim[[2]]$vals, time = times,
               grid = expand.grid(lon = oscar$var$u$dim[[1]]$vals, lat = oscar$var$u$dim[[2]]$vals))

  if (!is.null(lats)) {
    l = which(oscar$lat >= lats[1] & oscar$lat <= lats[2])
    oscar$u = oscar$u[,l,]
    oscar$v = oscar$v[,l,]
    oscar$lat = oscar$lat[l]
    oscar$grid = expand.grid(lon = oscar$lon, lat = oscar$lat)
    if (verbose) {message(' Trimming OSCAR latitudes (', length(l), ' entries remain).')}
  }
  if (!is.null(lons)) {
    lons[lons<0] = lons[lons<0] + 360
    ## Filter
    l = which(oscar$lon >= lons[1] & oscar$lon <= lons[2])
    oscar$u = oscar$u[l,,]
    oscar$v = oscar$v[l,,]
    oscar$lon = oscar$lon[l]
    oscar$grid = expand.grid(lon = oscar$lon, lat = oscar$lat)
    if (verbose) {message(' Trimming OSCAR longitudes (', length(l), ' entries remain).')}
  }
  if (verbose){message(' OSCAR succesfully loaded (', format(Sys.time() - a, digits = 3), ' sec).')}

  oscar$w = matrix(0, nrow = nrow(oscar$u), ncol = ncol(oscar$u))
  oscar$get.vel = get.vel
  ## Return
  oscar
}







#' @title Load in ROMS Advection
#' @param path ROMS file (includign path) to read in
#' @param verbose A verbose flag
#' @export
#' @author Thomas Bryce Kelly
#' @import ncdf4
load.advection.roms = function(path, lats = NULL, lons = NULL, get.vel = get.vel.roms, verbose = T) {

  if (verbose) {
    message('Loading ROMS Advection Product\n Loading ROMS input file ', path, ' (', format(file.size(path)/2^20, digits = 3), ' MB)... ', appendLF = F)
    a = Sys.time()
  }

  file = ncdf4::nc_open(path, write = FALSE)

  #### Load variables
  ## Vertical Grid
  h = ncdf4::ncvar_get(file, 'h')
  hc = ncdf4::ncvar_get(file, 'hc')
  theta = ncdf4::ncvar_get(file, 'theta_s')
  N = 42 # Hard coded

  h = (h[2:dim(h)[1],] + h[1:(dim(h)[1] - 1),]) / 2 ## interpolate to center of each grid cell
  h = (h[,2:dim(h)[2]] + h[,1:(dim(h)[2] - 1)]) / 2

  ## Horizontal
  lat = ncdf4::ncvar_get(file, 'lat_psi')
  lon = ncdf4::ncvar_get(file, 'lon_psi') %% 360

  ## Advect
  u = ncdf4::ncvar_get(file, 'u')
  u = (u[,2:dim(u)[2],,] + u[,1:(dim(u)[2] - 1),,]) / 2

  v = ncdf4::ncvar_get(file, 'v')
  v = (v[2:dim(v)[1],,,] + v[1:(dim(v)[1] - 1),,,]) / 2

  w = ncdf4::ncvar_get(file, 'w')
  w = (w[,,2:dim(w)[3],] + w[,,1:(dim(w)[3] - 1),]) / 2 # depth averaged w
  w = (w[2:dim(w)[1],,,] + w[1:(dim(w)[1] - 1),,,]) / 2
  w = (w[,2:dim(w)[2],,] + w[,1:(dim(w)[2] - 1),,]) / 2

  ncdf4::nc_close(file)

  dims = c(dim(u)[1], dim(u)[2], dim(u)[3], dim(u)[4])
  dim(u) = dims
  dim(v) = dims
  dim(w) = dims

  ## Filter domain
  if (!is.null(lons)) {
    lons = lons %% 360
    l = which(lon[,1] >= lons[1] & lon[,1] <= lons[2]) # TODO Improve lat lon filtering for curvilinear grids!
    u = u[l,,,]
    v = v[l,,,]
    w = w[l,,,]
    lon = lon[l,]
    h = h[l,]
  }

  if (!is.null(lats)) {
    l = which(lat[1,] >= lats[1] & lat[1,] <= lats[2])
    u = u[,l,,]
    v = v[,l,,]
    w = w[,l,,]
    lat = lat[l,]
    h = h[,l]
  }

  #### Calculated
  grid = expand.grid(lon = lon, lat = lat)
  z = get.zlevel.roms(h, hc, theta, N = N)
  time = get.roms.time(path, TRUE)

  if (verbose){message(' ROMS succesfully loaded (', format(Sys.time() - a, digits = 3), ' sec).')}

  #### Return object
  list(lat = lat, lon = lon, grid = grid, h = h, z = z, u = u, v = v, w = w, time = time, get.vel = get.vel)
}


#' @title Load in CMEMS Advection
#' @param path CMEMS NC file (includign path) to read in
#' @param verbose A verbose flag
#' @export
#' @author Thomas Bryce Kelly
#' @import ncdf4
load.advection.cmems = function(path, lats = NULL, lons = NULL, get.vel = get.vel.cmems, verbose = T) {
  if (verbose) {message('Loading CMEMS Advection Product.', appendLF = F)}
  a = Sys.time()
  for (i in 1:length(path)) {
    if (verbose) {
      message('\n Loading CMEMS input file ', path[i], ' (', format(file.size(path[i])/2^20, digits = 3), ' MB): ', appendLF = F)
    }

    file = ncdf4::nc_open(path[i], write = FALSE)

    #### Load variables
    ## Vertical Grid
    z = ncdf4::ncvar_get(file, 'depth')
    #ssh = ncdf4::ncvar_get(file, 'zos')
    ssh = 'not implemented yet'

    ## Horizontal
    lat = ncdf4::ncvar_get(file, 'latitude')
    lon = ncdf4::ncvar_get(file, 'longitude') %% 360
    if (verbose) {message(' Grid.', appendLF = F)}

    ## Temporal
    time = ncdf4::ncvar_get(file, 'time')
    time = make.time(1950, 1, 1) + time * 3600 ## hours since 1950-01-01

    ## Advect
    if (verbose) {message(' U.', appendLF = F)}
    u = ncdf4::ncvar_get(file, 'uo')
    u = array(u, dim = c(length(lon), length(lat), length(z), length(time)))
    if (verbose) {message(' V.', appendLF = F)}
    v = ncdf4::ncvar_get(file, 'vo')
    v = array(v, dim = c(length(lon), length(lat), length(z), length(time)))

    ncdf4::nc_close(file)

    ## Filter domain
    if (!is.null(lons)) {
      if (verbose) {message(' Filtering Lon.', appendLF = F)}
      lons = lons %% 360
      l = which(lon >= lons[1] & lon <= lons[2]) # TODO Improve lat lon filtering for curvilinear grids!
      u = array(u[l,,,], dim = c(length(l), length(lat), length(z), length(time)))
      v = array(v[l,,,], dim = c(length(l), length(lat), length(z), length(time)))
      lon = lon[l]
      #ssh = array(ssh[l,,,], dim = c(length(l), length(lat), length(time), 1))
    }

    if (!is.null(lats)) {
      if (verbose) {message(' Filtering Lat.', appendLF = F)}
      l = which(lat >= lats[1] & lat <= lats[2])
      u = array(u[,l,,], dim = c(length(lon), length(l), length(z), length(time)))
      v = array(v[,l,,], dim = c(length(lon), length(l), length(z), length(time)))
      lat = lat[l]
      #ssh = array(ssh[,l,,], dim = c(length(lon), length(l), length(time), 1))
    }

    if (i == 1) {
      cmems = list(lat = lat, lon = lon, z = -z, ssh = ssh, u = u, v = v, w = array(0, dim = c(length(lon), length(lat), length(z), length(time))),
                   time = time) ## negate depth to be consistent with other products
    } else {
      dims = dim(cmems$u)
      cmems$u = array(c(cmems$u, u), dim = c(dims[1:3], dims[4]+1))
      cmems$v = array(c(cmems$u, u), dim = c(dims[1:3], dims[4]+1))
      cmems$time = c(cmems$time, time)
    }
  }
  if (verbose){message('\nCMEMS succesfully loaded (', format(Sys.time() - a, digits = 3), ', ', format(object.size(cmems) / 2^20, digits = 3), ' MB).')}
  #### Return object
  cmems$get.vel = get.vel
  cmems
}



#' @title Load in CMEMS Advection
#' @param path CMEMS NC file (includign path) to read in
#' @param verbose A verbose flag
#' @export
#' @author Thomas Bryce Kelly
#' @import ncdf4
load.advection.ecco = function(u.path, v.path, w.path, lats = NULL, lons = NULL, get.vel = NULL, verbose = T) {
  if (verbose) {message('Loading ECCO Advection Product.', appendLF = F)}
  a = Sys.time()

  if (verbose) {
    message('\n Loading ECCO input file ', u.path, ' (', format(file.size(u.path)/2^20, digits = 3), ' MB): ', appendLF = F)
    message('\n Loading ECCO input file ', v.path, ' (', format(file.size(v.path)/2^20, digits = 3), ' MB): ', appendLF = F)
    message('\n Loading ECCO input file ', w.path, ' (', format(file.size(w.path)/2^20, digits = 3), ' MB): ', appendLF = F)
  }

  u.file = ncdf4::nc_open(u.path, write = FALSE)
  v.file = ncdf4::nc_open(v.path, write = FALSE)
  w.file = ncdf4::nc_open(w.path, write = FALSE)

  #### Load variables
  ## Vertical Grid
  u = ncdf4::ncvar_get(u.file, 'UVEL')
  v = ncdf4::ncvar_get(v.file, 'VVEL')
  w = ncdf4::ncvar_get(w.file, 'WVEL')

  ## Domain info
  time = ncdf4::ncvar_get(u.file, 'TIME')* 86400 + make.time(1992, 1, 1)
  lat = ncdf4::ncvar_get(u.file, 'LATITUDE_T')
  lon = ncdf4::ncvar_get(u.file, 'LONGITUDE_T')
  depth = ncdf4::ncvar_get(u.file, 'DEPTH_T')

  ## Close files
  ncdf4::nc_close(u.file)
  ncdf4::nc_close(v.file)
  ncdf4::nc_close(w.file)

  ## Filter domain
  if (!is.null(lons)) {
    if (verbose) {message(' Filtering Lon.', appendLF = F)}
    lons = lons %% 360
    l = which(lon >= lons[1] & lon <= lons[2]) # TODO Improve lat lon filtering for curvilinear grids!
    u = u[l,,]
    v = v[l,,]
    w = w[l,,]
  }

  if (!is.null(lats)) {
    if (verbose) {message(' Filtering Lat.', appendLF = F)}
      l = which(lat >= lats[1] & lat <= lats[2])
      u = u[,l,]
      v = v[,l,]
      w = w[,l,]
    }

    #### Calculated
    #grid = expand.grid(lon = lon, lat = lat)

  ecco = list(lat = lat, lon = lon, z = depth,
              u = u, v = v, w = w, time = time) ## negate depth to be consistent with other products

  if (verbose){message('\nECCO succesfully loaded (', format(Sys.time() - a, digits = 3), ', ', format(object.size(ecco) / 2^20, digits = 3), ' MB).')}
  #### Return object
  ecco$get.vel = get.vel
  ecco
}



#########################################
#### Functions to retreive velocity #####
#########################################

#' @title Get Velocity (OSCAR)
#' @author Thomas Bryce Kelly
#' @export
get.vel.oscar = function(lon, lat, depth, time, advection) {

  grid = data.frame(lon = lon, lat = lat)
  grid$lon[grid$lon < 0] = grid$lon[grid$lon < 0] + 360
  grid$u = 0
  grid$v = 0

  if (time < min(advection$time) | time > max(advection$time)) {
    message('Improper advection product loaded for time = ', time)
    nil = rep(0, length(lon))

    return(list(u = nil, v = nil, w = nil))
  }

  t2 = min(which(advection$time > time))
  t1 = t2-1
  t1w = as.numeric(advection$time[t2] - time) / as.numeric(advection$time[t2] - advection$time[t1])

  x2 = sapply(grid$lon, function(x) {min(which(advection$lon > x), length(advection$lon))})
  y2 = sapply(grid$lat, function(x) {max(which(advection$lat > x), 1)})
  x1 = x2 - 1
  y1 = y2 + 1

  ## Calculate weights
  x1w = (advection$lon[x2] - grid$lon) / (advection$lon[x2] - advection$lon[x1])
  y1w = (advection$lat[y2] - grid$lat) / (advection$lat[y2] - advection$lat[y1])

  ## calculate u
  grid$u = t1w * ((advection$u[cbind(x1,y1,t1)] * x1w + advection$u[cbind(x2,y1,t1)] * (1 - x1w)) * y1w + (advection$u[cbind(x1,y2,t1)] * x1w + advection$u[cbind(x2,y2,t1)] * (1 - x1w)) * (1 - y1w)) +
    (1 - t1w) * ((advection$u[cbind(x1,y1,t2)] * x1w + advection$u[cbind(x2,y1,t2)] * (1 - x1w)) * y1w + (advection$u[cbind(x1,y2,t2)] * x1w + advection$u[cbind(x2,y2,t2)] * (1 - x1w)) * (1 - y1w))

  ## calculate v
  grid$v = t1w * ((advection$v[cbind(x1,y1,t1)] * x1w + advection$v[cbind(x2,y1,t1)] * (1 - x1w)) * y1w + (advection$v[cbind(x1,y2,t1)] * x1w + advection$v[cbind(x2,y2,t1)] * (1 - x1w)) * (1 - y1w)) +
    (1 - t1w) * ((advection$v[cbind(x1,y1,t2)] * x1w + advection$v[cbind(x2,y1,t2)] * (1 - x1w)) * y1w + (advection$v[cbind(x1,y2,t2)] * x1w + advection$v[cbind(x2,y2,t2)] * (1 - x1w)) * (1 - y1w))

  grid[is.na(grid)] = 0

  ## Return
  list(u = grid$u,
       v = grid$v,
       w = rep(0, length(lon)))
}


#' @title Get Velocity (CMEMS)
#' @author Thomas Bryce Kelly
#' @export
get.vel.cmems = function(lon, lat, depth, time, advection) {

  grid = data.frame(lon = lon, lat = lat, depth = depth)
  grid$lon[grid$lon < 0] = grid$lon[grid$lon < 0] + 360
  grid$u = NA
  grid$v = NA
  grid$w = NA

  ## check for timing conflict
  if (time < min(advection$time) | time > max(advection$time)) {
    message('Improper advection product loaded for time = ', time)
    nil = rep(0, length(lon))

    return(list(u = nil, v = nil, w = nil))
  }

  t2 = min(which(advection$time > time))
  t1 = t2-1
  t1w = as.numeric(advection$time[t2] - time) / as.numeric(advection$time[t2] - advection$time[t1])

  for (i in 1:nrow(grid)) {
    ## Find Indicies
    x2 = min(which(advection$lon > grid$lon[i]))
    x1 = x2 - 1                         ## Makes assumptions about no boundaries
    y2 = max(which(advection$lat > grid$lat[i]))
    y1 = y2 + 1
    z2 = min(which(advection$z < grid$depth[i]))
    z1 = z2 - 1

    if (x1 == 0 | y2 == 0 | t1 == 0 | z1 == 0 | x2 > length(advection$lon) | y1 > length(advection$lat) | z2 > length(advection$z) | t2 > length(advection$time)){
      message('Point out of bounds! x1 = ', x1, ', y1 = ', y1, ', t1 = ', t1)
      grid$u[i] = 0
      grid$v[i] = 0
      grid$w[i] = 0
    } else {

      ## Calculate weights
      x1w = (advection$lon[x2] - grid$lon[i]) / (advection$lon[x2] - advection$lon[x1])
      y1w = (advection$lat[y2] - grid$lat[i]) / (advection$lat[y2] - advection$lat[y1])
      z1w = (advection$z[z2] - grid$depth[i]) / (advection$z[z2] - advection$z[z1])

      ## calculate u
      u1 = z1w * ((advection$u[x1,y1,z1,t1] * x1w + advection$u[x2,y1,z1,t1] * (1 - x1w)) * y1w + (advection$u[x1,y2,z1,t1] * x1w + advection$u[x2,y2,z1,t1] * (1 - x1w)) * (1 - y1w)) +
          (1 - z1w) * ((advection$u[x1,y1,z2,t1] * x1w + advection$u[x2,y1,z2,t1] * (1 - x1w)) * y1w + (advection$u[x1,y2,z2,t1] * x1w + advection$u[x2,y2,z2,t1] * (1 - x1w)) * (1 - y1w))
      u2 = z1w * ((advection$u[x1,y1,z1,t2] * x1w + advection$u[x2,y1,z1,t2] * (1 - x1w)) * y1w + (advection$u[x1,y2,z1,t2] * x1w + advection$u[x2,y2,z1,t2] * (1 - x1w)) * (1 - y1w)) +
        (1 - z1w) * ((advection$u[x1,y1,z1,t2] * x1w + advection$u[x2,y1,z1,t2] * (1 - x1w)) * y1w + (advection$u[x1,y2,z1,t2] * x1w + advection$u[x2,y2,z1,t2] * (1 - x1w)) * (1 - y1w))
      grid$u[i] = u1 * t1w + u2 * (1 - t1w)

      ## calculate v
      v1 = z1w * ((advection$v[x1,y1,z1,t1] * x1w + advection$v[x2,y1,z1,t1] * (1 - x1w)) * y1w + (advection$v[x1,y2,z1,t1] * x1w + advection$v[x2,y2,z1,t1] * (1 - x1w)) * (1 - y1w)) +
        (1 - z1w) * ((advection$v[x1,y1,z2,t1] * x1w + advection$v[x2,y1,z2,t1] * (1 - x1w)) * y1w + (advection$v[x1,y2,z2,t1] * x1w + advection$v[x2,y2,z2,t1] * (1 - x1w)) * (1 - y1w))
      v2 = z1w * ((advection$v[x1,y1,z1,t2] * x1w + advection$v[x2,y1,z1,t2] * (1 - x1w)) * y1w + (advection$v[x1,y2,z1,t2] * x1w + advection$v[x2,y2,z1,t2] * (1 - x1w)) * (1 - y1w)) +
        (1 - z1w) * ((advection$v[x1,y1,z1,t2] * x1w + advection$v[x2,y1,z1,t2] * (1 - x1w)) * y1w + (advection$v[x1,y2,z1,t2] * x1w + advection$v[x2,y2,z1,t2] * (1 - x1w)) * (1 - y1w))
      grid$v[i] = v1 * t1w + v2 * (1 - t1w)

      ## calculate w
      w1 = z1w * ((advection$w[x1,y1,z1,t1] * x1w + advection$w[x2,y1,z1,t1] * (1 - x1w)) * y1w + (advection$w[x1,y2,z1,t1] * x1w + advection$w[x2,y2,z1,t1] * (1 - x1w)) * (1 - y1w)) +
        (1 - z1w) * ((advection$w[x1,y1,z2,t1] * x1w + advection$w[x2,y1,z2,t1] * (1 - x1w)) * y1w + (advection$w[x1,y2,z2,t1] * x1w + advection$w[x2,y2,z2,t1] * (1 - x1w)) * (1 - y1w))
      w2 = z1w * ((advection$w[x1,y1,z1,t2] * x1w + advection$w[x2,y1,z1,t2] * (1 - x1w)) * y1w + (advection$w[x1,y2,z1,t2] * x1w + advection$w[x2,y2,z1,t2] * (1 - x1w)) * (1 - y1w)) +
        (1 - z1w) * ((advection$w[x1,y1,z1,t2] * x1w + advection$w[x2,y1,z1,t2] * (1 - x1w)) * y1w + (advection$w[x1,y2,z1,t2] * x1w + advection$w[x2,y2,z1,t2] * (1 - x1w)) * (1 - y1w))
      grid$w[i] = w1 * t1w + w2 * (1 - t1w)
    }
  }

  grid[is.na(grid)] = 0

  ## Return
  list(u = grid$u,
       v = grid$v,
       w = grid$w)
}


#' @title Get Velocity (ROMS)
#' @author Thomas Bryce Kelly
#' @export
get.vel.roms = function(lon, lat, depth, time, advection) {

  ## Build return object
  velocity = data.frame(lon = lon, lat = lat, depth = depth, time = time, u = NA, v = NA, w = NA)
  velocity$lon = velocity$lon %% 360

  if (time < min(advection$time)) {
    warning(paste0('Improper advection product loaded for time = ', time, ' (too early)'))
    t2 = 1
    t1 = 1
    t1w = 0.5
  } else if (time > max(advection$time)) {
    warning(paste0('Improper advection product loaded for time = ', time, ' (too early)'))
    t2 = length(roms$time)
    t1 = length(roms$time)
    t1w = 0.5
  } else {
    ## Time weights: linear in time
    t2 = min(which(advection$time > time))
    t1 = t2-1
    t1w = as.numeric(advection$time[t2] - time) / as.numeric(advection$time[t2] - advection$time[t1])
  }

  for (i in 1:nrow(velocity)) {
    index = which.min((velocity$lon[i] - advection$lon)^2 + (velocity$lat[i] - advection$lat)^2)
    xx = index %% dim(advection$lon)[1]
    if (xx == 0) { xx = dim(advection$lon)[1] }
    yy = floor((index - 1) / dim(advection$lon)[1]) + 1


    if (advection$lon[xx, yy] > velocity$lon[i]) { x1 = max(1, xx-1); x2 = xx } else { x1 = xx; x2 = min(xx + 1, dim(advection$lon)[1])}
    if (advection$lat[xx, yy] > velocity$lat[i]) { y1 = max(1, yy-1); y2 = yy } else { y1 = yy; y2 = min(yy + 1, dim(advection$lon)[2])}


    depths = apply(roms$z[c(x1, x2), c(y1,y2),], 3, function(x){ mean(x[!is.na(x)])})
    z2 = max(which(depths < velocity$depth[i]))
    z1 = min(z2 + 1, dim(roms$u)[3])

    if (x1 == 0 | y2 == 0 | x2 > dim(advection$lon)[1] | y1 > dim(advection$lat)[2]){
      message('Point out of bounds! x1 = ', x1, ', y1 = ', y1, ', z1 = ', z1, ', t1 = ', t1)
      grid$u[i] = 0
      grid$v[i] = 0
      grid$w[i] = 0
    } else {

      ## Calculate weights
      x1w = (advection$lon[x2,y1] - velocity$lon[i]) / (advection$lon[x2,y1] - advection$lon[x1,y1])
      y1w = (advection$lat[x1,y2] - velocity$lat[i]) / (advection$lat[x1,y2] - advection$lat[x1,y1])
      z1w = (depths[z2] - velocity$depth[i]) / (depths[z2] - depths[z1])

      ## calculate u
      u1 = ((advection$u[x1,y1,z1,t1] * x1w + advection$u[x2,y1,z1,t1] * (1 - x1w)) * y1w + (advection$u[x1,y2,z1,t1] * x1w + advection$u[x2,y2,z1,t1] * (1 - x1w)) * (1 - y1w)) * z1w +
        ((advection$u[x1,y1,z1,t1] * x1w + advection$u[x2,y1,z1,t1] * (1 - x1w)) * y1w + (advection$u[x1,y2,z1,t1] * x1w + advection$u[x2,y2,z1,t1] * (1 - x1w)) * (1 - y1w)) * (1 - z1w)
      u2 = ((advection$u[x1,y1,z2,t2] * x1w + advection$u[x2,y1,z2,t2] * (1 - x1w)) * y1w + (advection$u[x1,y2,z2,t2] * x1w + advection$u[x2,y2,z2,t2] * (1 - x1w)) * (1 - y1w)) * z1w +
        ((advection$u[x1,y1,z2,t2] * x1w + advection$u[x2,y1,z2,t2] * (1 - x1w)) * y1w + (advection$u[x1,y2,z2,t2] * x1w + advection$u[x2,y2,z2,t2] * (1 - x1w)) * (1 - y1w)) * (1 - z1w)
      velocity$u[i] = u1 * t1w + u2 * (1 - t1w)

      ## calculate v
      v1 = ((advection$v[x1,y1,z1,t1] * x1w + advection$v[x2,y1,z1,t1] * (1 - x1w)) * y1w + (advection$v[x1,y2,z1,t1] * x1w + advection$v[x2,y2,z1,t1] * (1 - x1w)) * (1 - y1w)) * z1w +
        ((advection$v[x1,y1,z1,t1] * x1w + advection$v[x2,y1,z1,t1] * (1 - x1w)) * y1w + (advection$v[x1,y2,z1,t1] * x1w + advection$v[x2,y2,z1,t1] * (1 - x1w)) * (1 - y1w)) * (1 - z1w)
      v2 = ((advection$v[x1,y1,z1,t2] * x1w + advection$v[x2,y1,z1,t2] * (1 - x1w)) * y1w + (advection$v[x1,y2,z1,t2] * x1w + advection$v[x2,y2,z1,t2] * (1 - x1w)) * (1 - y1w)) * z1w +
        ((advection$v[x1,y1,z1,t2] * x1w + advection$v[x2,y1,z1,t2] * (1 - x1w)) * y1w + (advection$v[x1,y2,z1,t2] * x1w + advection$v[x2,y2,z1,t2] * (1 - x1w)) * (1 - y1w)) * (1 - z1w)
      velocity$v[i] = v1 * t1w + v2 * (1 - t1w)

      ## calculate w
      w1 = ((advection$w[x1,y1,z1,t1] * x1w + advection$w[x2,y1,z1,t1] * (1 - x1w)) * y1w + (advection$w[x1,y2,z1,t1] * x1w + advection$w[x2,y2,z1,t1] * (1 - x1w)) * (1 - y1w)) * z1w +
         ((advection$w[x1,y1,z1,t1] * x1w + advection$w[x2,y1,z1,t1] * (1 - x1w)) * y1w + (advection$w[x1,y2,z1,t1] * x1w + advection$w[x2,y2,z1,t1] * (1 - x1w)) * (1 - y1w)) * (1 - z1w)
      w2 = ((advection$w[x1,y1,z1,t2] * x1w + advection$w[x2,y1,z1,t2] * (1 - x1w)) * y1w + (advection$w[x1,y2,z1,t2] * x1w + advection$w[x2,y2,z1,t2] * (1 - x1w)) * (1 - y1w)) * z1w +
        ((advection$w[x1,y1,z1,t2] * x1w + advection$w[x2,y1,z1,t2] * (1 - x1w)) * y1w + (advection$w[x1,y2,z1,t2] * x1w + advection$w[x2,y2,z1,t2] * (1 - x1w)) * (1 - y1w)) * (1 - z1w)
      velocity$w[i] = w1 * t1w + w2 * (1 - t1w)
    }
  }
  velocity[is.na(velocity)] = 0

  ## Return
  list(u = velocity$u, v = velocity$v, w = velocity$w)
}



#######################
##### Advection Fn ####
#######################


#' @title Advect Lagrangian Particle
#' @author Thomas Bryce Kelly
#' @param particles a particle release dataframe such as that generated by `init.lagrangian.particles()`
#' @param model a model list such as that generated by `init.lagrangian.model()`
#' @param advection An advection list, such as one generated by `load.oscar.advection()`, that contains U V and W
#' @param dt Time step for RK4 advection scheme
#' @export
run.advection = function(particles, model, advection, zlim = c(-6e3, 0), verbose = T) {

  if (verbose) { message('Starting Lagrangian model run:')}
  a = Sys.time()
  b = 0
  for (timestep in 1:model$meta$nstep) {

    if (verbose) { message('\nTimestep = ', timestep, ' (', round(timestep / model$meta$nstep * 100), '%)', appendLF = F) }

    ## Save if it is a saving timestep!
    if (timestep %% model$meta$save.freq == 0) {

      c = Sys.time()

      savestep = timestep / model$meta$save.freq

      ## 1) lon, 2) lat, 3) depth, 4) u, 5) v, 6) w, 7) T, 8) S, 9) alive
      model$hist[,savestep,1] = particles$current.lon
      model$hist[,savestep,2] = particles$current.lat
      model$hist[,savestep,3] = particles$current.depth
      model$hist[,savestep,4] = particles$current.u
      model$hist[,savestep,5] = particles$current.v
      model$hist[,savestep,6] = particles$current.w
      model$hist[,savestep,9] = particles$alive

      ## TODO Add u, v, w, T, and S info!
      if (verbose) { message(' Saved savestep = ', savestep, appendLF = F)}
      b = as.numeric(difftime(Sys.time(), c, units = 'secs')) + b
    }


    if (model$meta$t.step > 0) {
      timet = model$meta$t.start + timestep * model$meta$t.step
      time2 = model$meta$t.start + (timestep + 0.5) * model$meta$t.step
      time3 = model$meta$t.start + (timestep + 0.5) * model$meta$t.step
      time4 = model$meta$t.start + (timestep + 1) * model$meta$t.step
    } else {
      timet = model$meta$t.end + timestep * model$meta$t.step
      time2 = model$meta$t.end + (timestep + 0.5) * model$meta$t.step
      time3 = model$meta$t.end + (timestep + 0.5) * model$meta$t.step
      time4 = model$meta$t.end + (timestep + 1) * model$meta$t.step
    }

    ## make any new particles alive
    if (verbose) {message('  Awaking particles.', appendLF = F)}
    if (model$meta$t.step > 0) {
      particles$alive[particles$start.time < time4 & particles$start.time >= timet] = TRUE
    } else {
      particles$alive[particles$start.time > time4 & particles$start.time <= timet] = TRUE
    }

    ## Advect all living particles:
    l = which(particles$alive)

    if (verbose) {message('  Advecting particles.', appendLF = F)}

    if (length(l) > 0) {
      ## Initial Guess
      latt = particles$current.lat[l]
      lont = particles$current.lon[l]
      deptht = particles$current.depth[l]
      vel = advection$get.vel(lont, latt, deptht, timet, advection)

      ## p2 based on velocity at final pos at meta$t.step/2
      dlat = vel$v * 1e-3 * model$meta$t.step/2 / 6378.14 # distance / earth radius
      lat2 = 180 / pi * asin(sin(latt / 180 * pi) * cos(dlat) + cos(latt / 180 * pi) * sin(dlat))
      dlon = vel$u * 1e-3 * model$meta$t.step/2 / 6378.14 # distance / earth radius
      lon2 = lont + 180 / pi * atan2(sin(dlon) * cos(latt / 180 * pi), cos(dlon) - sin(latt / 180 * pi) * sin(lat2 / 180 * pi))
      ddepth = vel$w * model$meta$t.step/2
      depth2 = deptht + ddepth

      vel2 = advection$get.vel(lon2, lat2, depth2, time2, advection)

      ## p3 based on velocity at p2 at meta$t.step/2
      dlat = vel2$v * 1e-3 * model$meta$t.step/2 / 6378.14 # distance / earth radius
      lat3 = 180 / pi * asin(sin(latt / 180 * pi) * cos(dlat) + cos(latt / 180 * pi) * sin(dlat))
      dlon = vel2$u * 1e-3 * model$meta$t.step/2 / 6378.14 # distance / earth radius
      lon3 = lont + 180 / pi * atan2(sin(dlon) * cos(latt / 180 * pi), cos(dlon) - sin(latt / 180 * pi) * sin(lat3 / 180 * pi))
      ddepth = vel2$w * model$meta$t.step/2
      depth3 = deptht + ddepth

      vel3 = advection$get.vel(lon3, lat3, depth3, time3, advection)

      ## p4 based on velocity at p2 at meta$t.step/2
      dlat = vel$v *  model$meta$t.step * 1e-3 / 6378.14 # distance / earth radius
      lat4 = 180 / pi * asin(sin(latt / 180 * pi) * cos(dlat) + cos(latt / 180 * pi) * sin(dlat))
      dlon = vel$u * 1e-3 * model$meta$t.step / 6378.14 # distance / earth radius
      lon4 = lont + 180 / pi * atan2(sin(dlon) * cos(latt / 180 * pi), cos(dlon) - sin(latt / 180 * pi) * sin(lat4 / 180 * pi))
      ddepth = vel$w * model$meta$t.step
      depth4 = deptht + ddepth

      vel4 = advection$get.vel(lon4, lat4, depth4, time4, advection)

      ## RK4 velocity
      vel$u = (vel$u + 2 * vel2$u + 2 * vel3$u + vel4$u)/ 6
      vel$v = (vel$v + 2 * vel2$v + 2 * vel3$v + vel4$v)/ 6
      vel$w = (vel$w + 2 * vel2$w + 2 * vel3$w + vel4$w)/ 6

      dlat = vel$v * 1e-3 * model$meta$t.step / 6378.14 # distance / earth radius
      lat.final = 180 / pi * asin(sin(latt / 180 * pi) * cos(dlat) + cos(latt / 180 * pi) * sin(dlat))
      dlon = vel$u * 1e-3 * model$meta$t.step / 6378.14 # distance / earth radius
      lon.final = lont + 180 / pi * atan2(sin(dlon) * cos(latt / 180 * pi), cos(dlon) - sin(latt / 180 * pi) * sin(lat.final / 180 * pi))
      ddepth = vel$w * model$meta$t.step

      for (p in l){
        #### Sinking
        if(!is.na(particles$sink[p])) { ddepth[p] = ddepth[p] + particles$sink[p] * model$meta$t.step }
        if (!is.na(particles$const.depth[p])) {ddepth[p] = 0}
      }
      depth.final = deptht + ddepth

      ## Update fields
      particles$current.lon[l] = lon.final
      particles$current.lat[l] = lat.final
      particles$current.depth[l] = depth.final
      particles$current.u[l] = vel$u
      particles$current.v[l] = vel$v
      particles$current.w[l] = vel$w
    }

    if (verbose) {message('  Killing particles.', appendLF = F)}
    particles$alive[!is.na(particles$current.depth) & (particles$current.depth > zlim[2] | particles$current.depth < zlim[1])] = F

  }

  if (verbose) {
    message('\nModel finished in\t', format(difftime(Sys.time(), a, units = 'sec'), digits = 4))
    message('\tSaving used:\t', format(b, digits = 4), ' secs')
  }

  ## Return
  list(particles = particles, hist = model$hist, meta = model$meta)
}


#' @title Advect Lagrangian Particles Fast
#' @author Thomas Bryce Kelly
#' @param particles a particle release dataframe such as that generated by `init.lagrangian.particles()`
#' @param model a model list such as that generated by `init.lagrangian.model()`
#' @param advection An advection list, such as one generated by `load.oscar.advection()`, that contains U V and W
#' @param dt Time step for forward euler advection
#' @export
run.advection.fast = function(particles, model, advection, zlim = c(-6e3, 0), verbose = T) {

  if (verbose) { message('Starting Lagrangian model run:')}
  a = Sys.time()
  b = 0
  for (timestep in 1:model$meta$nstep) {

    if (verbose) { message('\nTimestep = ', timestep, ' (', round(timestep / model$meta$nstep * 100), '%)', appendLF = F) }

    ## Save if it is a saving timestep!
    if (timestep %% model$meta$save.freq == 0) {

      c = Sys.time()

      savestep = timestep / model$meta$save.freq

      ## 1) lon, 2) lat, 3) depth, 4) u, 5) v, 6) w, 7) T, 8) S, 9) alive
      model$hist[,savestep,1] = particles$current.lon
      model$hist[,savestep,2] = particles$current.lat
      model$hist[,savestep,3] = particles$current.depth
      model$hist[,savestep,4] = particles$current.u
      model$hist[,savestep,5] = particles$current.v
      model$hist[,savestep,6] = particles$current.w
      model$hist[,savestep,9] = particles$alive

      ## TODO Add u, v, w, T, and S info!
      if (verbose) { message(' Saved savestep = ', savestep, appendLF = F)}
      b = as.numeric(difftime(Sys.time(), c, units = 'secs')) + b
    }


    if (model$meta$t.step > 0) {
      timet = model$meta$t.start + timestep * model$meta$t.step
      time2 = model$meta$t.start + (timestep + 0.5) * model$meta$t.step
      time3 = model$meta$t.start + (timestep + 0.5) * model$meta$t.step
      time4 = model$meta$t.start + (timestep + 1) * model$meta$t.step
    } else {
      timet = model$meta$t.end + timestep * model$meta$t.step
      time2 = model$meta$t.end + (timestep + 0.5) * model$meta$t.step
      time3 = model$meta$t.end + (timestep + 0.5) * model$meta$t.step
      time4 = model$meta$t.end + (timestep + 1) * model$meta$t.step
    }

    ## make any new particles alive
    if (verbose) {message('  Awaking particles.', appendLF = F)}
    if (model$meta$t.step > 0) {
      particles$alive[particles$start.time < time4 & particles$start.time >= timet] = TRUE
    } else {
      particles$alive[particles$start.time > time4 & particles$start.time <= timet] = TRUE
    }

    ## Advect all living particles:
    l = which(particles$alive)

    if (verbose) {message('  Advecting particles.', appendLF = F)}

    if (length(l) > 0) {
      ## Initial Guess
      latt = particles$current.lat[l]
      lont = particles$current.lon[l]
      deptht = particles$current.depth[l]
      vel = advection$get.vel(lont, latt, deptht, timet, advection)

      dlat = vel$v * 1e-3 * model$meta$t.step / 6378.14 # distance / earth radius
      lat.final = 180 / pi * asin(sin(latt / 180 * pi) * cos(dlat) + cos(latt / 180 * pi) * sin(dlat))
      dlon = vel$u * 1e-3 * model$meta$t.step / 6378.14 # distance / earth radius
      lon.final = lont + 180 / pi * atan2(sin(dlon) * cos(latt / 180 * pi), cos(dlon) - sin(latt / 180 * pi) * sin(lat.final / 180 * pi))
      ddepth = vel$w * model$meta$t.step

      for (p in l){
        #### Sinking
        if(!is.na(particles$sink[p])) { ddepth[p] = ddepth[p] + particles$sink[p] * model$meta$t.step }
        if (!is.na(particles$const.depth[p])) {ddepth[p] = 0}
      }
      depth.final = deptht + ddepth

      ## Update fields
      particles$current.lon[l] = lon.final
      particles$current.lat[l] = lat.final
      particles$current.depth[l] = depth.final
      particles$current.u[l] = vel$u
      particles$current.v[l] = vel$v
      particles$current.w[l] = vel$w
    }

    if (verbose) {message('  Killing particles.', appendLF = F)}
    particles$alive[!is.na(particles$current.depth) & (particles$current.depth > zlim[2] | particles$current.depth < zlim[1])] = F

  }

  if (verbose) {
    message('\nModel finished in\t', format(difftime(Sys.time(), a, units = 'sec'), digits = 4))
    message('\tSaving used:\t', format(b, digits = 4), ' secs')
  }

  ## Return
  list(particles = particles, hist = model$hist, meta = model$meta)
}


#' @title Advect Lagrangian Particle Model in Parallel
#' @author Thomas Bryce Kelly
#' @param lat Initial latitude of particle
#' @param lon Initial longitude of particle
#' @param depth Intial depth of particle
#' @param time Start time of particle
#' @param advection An advection list, such as one generated by `load.oscar.advection()`, that contains U V and W
#' @param dt Time step for RK4 advection scheme
#' @export
run.advection.parallel = function(particles, meta, hist, advection, get.vel, zlim = c(-6e3, 0), verbose = T) {

  if (verbose) { message('Starting Lagrangian model run:')}
  a = Sys.time()
  b = 0

  parallelCluster <- parallel::makeCluster(4,
                                 type = "fork",
                                 methods = FALSE)
  parallel::setDefaultCluster(parallelCluster)
  doParallel::registerDoParallel(parallelCluster)

  for (timestep in 1:meta$nstep) {

    if (verbose) { message('\nTimestep = ', timestep, ' (', round(timestep / meta$nstep * 100), '%)', appendLF = F) }

    ## Save if it is a saving timestep!
    if (timestep %% meta$save.freq == 0) {

      c = Sys.time()

      savestep = timestep / meta$save.freq

      ## 1) lon, 2) lat, 3) depth, 4) u, 5) v, 6) w, 7) T, 8) S, 9) alive
      hist[,savestep,1] = particles$current.lon
      hist[,savestep,2] = particles$current.lat
      hist[,savestep,3] = particles$current.depth
      hist[,savestep,4] = particles$current.u
      hist[,savestep,5] = particles$current.v
      hist[,savestep,6] = particles$current.w
      hist[,savestep,9] = particles$alive

      ## TODO Add u, v, w, T, and S info!
      if (verbose) { message(' Saved savestep = ', savestep, appendLF = F)}
      b = as.numeric(difftime(Sys.time(), c, units = 'secs')) + b
    }


    if (meta$t.step > 0) {
      timet = meta$t.start + timestep * meta$t.step
      time2 = meta$t.start + (timestep + 0.5) * meta$t.step
      time3 = meta$t.start + (timestep + 0.5) * meta$t.step
      time4 = meta$t.start + (timestep + 1) * meta$t.step
    } else {
      timet = meta$t.end + timestep * meta$t.step
      time2 = meta$t.end + (timestep + 0.5) * meta$t.step
      time3 = meta$t.end + (timestep + 0.5) * meta$t.step
      time4 = meta$t.end + (timestep + 1) * meta$t.step
    }

    ## make any new particles alive
    if (verbose) {message('  Awaking particles.', appendLF = F)}
    if (meta$t.step > 0) {
      particles$alive[particles$start.time < time4 & particles$start.time >= timet] = TRUE
    } else {
      particles$alive[particles$start.time > time4 & particles$start.time <= timet] = TRUE
    }

    ## Advect all living particles:
    l = which(particles$alive)
    l.split <- split(1:length(l), l)
    results <- foreach(currentRow = l.split, .combine = c) %dopar% {
      data$c[currentRow] = data$a[currentRow] + data$b[currentRow]
    }

    if (verbose) {message('  Advecting particles.', appendLF = F)}

    if (length(l) > 0) {
      res = foreach(l = l.split, .combine = c) %dopar% {
        ## Initial Guess
        latt = particles$current.lat[l]
        lont = particles$current.lon[l]
        deptht = particles$current.depth[l]
        vel = get.vel(lont, latt, deptht, timet, advection)

        ## p2 based on velocity at final pos at meta$t.step/2
        dlat = vel$v * 1e-3 * meta$t.step/2 / 6378.14 # distance / earth radius
        lat2 = 180 / pi * asin(sin(latt / 180 * pi) * cos(dlat) + cos(latt / 180 * pi) * sin(dlat))
        dlon = vel$u * 1e-3 * meta$t.step/2 / 6378.14 # distance / earth radius
        lon2 = lont + 180 / pi * atan2(sin(dlon) * cos(latt / 180 * pi), cos(dlon) - sin(latt / 180 * pi) * sin(lat2 / 180 * pi))
        ddepth = vel$w * meta$t.step/2
        depth2 = deptht + ddepth

        vel2 = get.vel(lon2, lat2, depth2, time2, advection)

        ## p3 based on velocity at p2 at meta$t.step/2
        dlat = vel2$v * 1e-3 * meta$t.step/2 / 6378.14 # distance / earth radius
        lat3 = 180 / pi * asin(sin(latt / 180 * pi) * cos(dlat) + cos(latt / 180 * pi) * sin(dlat))
        dlon = vel2$u * 1e-3 * meta$t.step/2 / 6378.14 # distance / earth radius
        lon3 = lont + 180 / pi * atan2(sin(dlon) * cos(latt / 180 * pi), cos(dlon) - sin(latt / 180 * pi) * sin(lat3 / 180 * pi))
        ddepth = vel2$w * meta$t.step/2
        depth3 = deptht + ddepth

        vel3 = get.vel(lon3, lat3, depth3, time3, advection)

        ## p4 based on velocity at p2 at meta$t.step/2
        dlat = vel$v *  meta$t.step * 1e-3 / 6378.14 # distance / earth radius
        lat4 = 180 / pi * asin(sin(latt / 180 * pi) * cos(dlat) + cos(latt / 180 * pi) * sin(dlat))
        dlon = vel$u * 1e-3 * meta$t.step / 6378.14 # distance / earth radius
        lon4 = lont + 180 / pi * atan2(sin(dlon) * cos(latt / 180 * pi), cos(dlon) - sin(latt / 180 * pi) * sin(lat4 / 180 * pi))
        ddepth = vel$w * meta$t.step
        depth4 = deptht + ddepth

        vel4 = get.vel(lon4, lat4, depth4, time4, advection)

        ## RK4 velocity
        vel$u = (vel$u + 2 * vel2$u + 2 * vel3$u + vel4$u)/ 6
        vel$v = (vel$v + 2 * vel2$v + 2 * vel3$v + vel4$v)/ 6
        vel$w = (vel$w + 2 * vel2$w + 2 * vel3$w + vel4$w)/ 6

        dlat = vel$v * 1e-3 * meta$t.step / 6378.14 # distance / earth radius
        lat.final = 180 / pi * asin(sin(latt / 180 * pi) * cos(dlat) + cos(latt / 180 * pi) * sin(dlat))
        dlon = vel$u * 1e-3 * meta$t.step / 6378.14 # distance / earth radius
        lon.final = lont + 180 / pi * atan2(sin(dlon) * cos(latt / 180 * pi), cos(dlon) - sin(latt / 180 * pi) * sin(lat.final / 180 * pi))
        ddepth = vel$w * meta$t.step

        for (p in l){
          #### Sinking
          if(!is.na(particles$sink[p])) { ddepth[p] = ddepth[p] + particles$sink[p] * meta$t.step }
          if (!is.na(particles$const.depth[p])) {ddepth[p] = 0}
        }
        depth.final = deptht + ddepth

        ## Update fields
        particles$current.lon[l] = lon.final
        particles$current.lat[l] = lat.final
        particles$current.depth[l] = depth.final
        particles$current.u[l] = vel$u
        particles$current.v[l] = vel$v
        particles$current.w[l] = vel$w
      }
    }

    if (verbose) {message('  Killing particles.', appendLF = F)}
    particles$alive[!is.na(particles$current.depth) & (particles$current.depth > zlim[2] | particles$current.depth < zlim[1])] = F

  }

  if (verbose) {
    message('\nModel finished in\t', format(difftime(Sys.time(), a, units = 'sec'), digits = 4))
    message('\tSaving used:\t', format(b, digits = 4), ' secs')
  }

  ## Return
  list(particles = particles, hist = hist, meta = meta)
}


