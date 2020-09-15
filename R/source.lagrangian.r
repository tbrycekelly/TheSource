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
init.lagrangian.particles = function(lon = NA, lat = NA, depth = NA, time = NA, sink = NULL, const.depth = NULL, density = NULL, verbose = T) {

  if (verbose) {
    message('Initializing Langraign Particle Data Structure...')
    message('\t This function is helpful in organizing the initial data structure in order to be compatible for the advection scripts in this library. Optionally you can provide vectors (one entry per particle).')
    message('\n\t Optional Entries:')
    message('\t\t sink:\t\t the sinking speed of the particle (m/s)')
    message('\t\t const.depth:\t a depth for the particle to remain at (m)')
    message('\t\t density:\t an isopicnal density for the float to stay at (kg m-3)')
    message('\t Note: Longitudes will be coerced to 0-360 for compatibility.\n')
  }

  a = list()

  for (i in 1:length(lon)) {
    a[[paste0('p', i)]] = list(lon = lon[i] %% 360, lat = lat[i], depth = depth[i], time = time[i], sink = sink[i], const.depth = const.depth[i], density = density[i])
  }

  if (verbose){message('\t Returning particle file with ', length(a), ' particles.')}
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
load.advection.oscar = function(file, lats = NULL, lons = NULL, verbose = F) {

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

  ## Return
  oscar
}



#' @title Load in ROMS Advection
#' @param path ROMS file (includign path) to read in
#' @param verbose A verbose flag
#' @export
#' @author Thomas Bryce Kelly
#' @import ncdf4
load.advection.roms = function(path, verbose = F) {

  if (verbose) {
    message('Loading ROMS Advection Product\n Loading ROMS input file ', file, ' (', format(file.size(path)/2^20, digits = 3), ' MB)... ', appendLF = F)
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
  grid = list(lat = c(lat), lon = c(lon))

  ## Advect
  u = ncdf4::ncvar_get(file, 'u')
  u = (u[,2:dim(u)[2],,] + u[,1:(dim(u)[2] - 1),,]) / 2

  v = ncdf4::ncvar_get(file, 'v')
  v = (v[2:dim(v)[1],,,] + v[1:(dim(v)[1] - 1),,,]) / 2

  w = ncdf4::ncvar_get(file, 'w')
  w = (w[,,2:dim(w)[3],] + w[,,1:(dim(w)[3] - 1),]) / 2 # depth averaged w
  w = (w[2:dim(w)[1],,,] + w[1:(dim(w)[1] - 1),,,]) / 2
  w = (w[,2:dim(w)[2],,] + w[,1:(dim(w)[2] - 1),,]) / 2

  dims = c(dim(u)[1], dim(u)[2], dim(u)[3], dim(u)[4])
  dim(u) = dims
  dim(v) = dims
  dim(w) = dims

  #### Calculated
  z = get.zlevel.roms(h, hc, theta, N = N)
  time = get.roms.time(path, TRUE)

  if (verbose){message(' ROMS succesfully loaded (', format(Sys.time() - a, digits = 3), ' sec).')}

  #### Return object
  list(lat = lat, lon = lon, grid = grid, h = h, z = z, u = u, v = v, w = w, time = time)
}


#########################################
#### Functions to retreive velocity #####
#########################################

#' @title Get Velocity (OSCAR)
#' @author Thomas Bryce Kelly
#' @export
get.vel.oscar = function(lon, lat, time, oscar) {

  grid = data.frame(lon = lon, lat = lat)
  grid$lon[grid$lon < 0] = grid$lon[grid$lon < 0] + 360
  grid$u = NA
  grid$v = NA

  if (time < min(oscar$time) | time > max(oscar$time)) {
    stop(paste0('Improper OSCAR product loaded for time = ', time))
  }
  t2 = min(which(oscar$time > time))
  t1 = t2-1
  t1w = as.numeric(oscar$time[t2] - time) / as.numeric(oscar$time[t2] - oscar$time[t1])

  for (i in 1:nrow(grid)) {
    ## Find Indicies
    x2 = min(which(oscar$lon > grid$lon[i]))
    x1 = x2 - 1                         ## Makes assumptions about no boundaries
    y2 = max(which(oscar$lat > grid$lat[i]))
    y1 = y2 + 1


    if (x1 == 0 | y2 == 0 | t1 == 0 | x2 > length(oscar$lon) | y1 > length(oscar$lat) | t2 > length(oscar$time)){
      message('Point out of bounds! x1 = ', x1, ', y1 = ', y1, ', t1 = ', t1)
      grid$u[i] = 0
      grid$v[i] = 0
    } else {

      ## Calculate weights
      x1w = (oscar$lon[x2] - grid$lon[i]) / (oscar$lon[x2] - oscar$lon[x1])
      y1w = (oscar$lat[y2] - grid$lat[i]) / (oscar$lat[y2] - oscar$lat[y1])

      ## calculate u
      u1 = (oscar$u[x1,y1,t1] * x1w + oscar$u[x2,y1,t1] * (1 - x1w)) * y1w + (oscar$u[x1,y2,t1] * x1w + oscar$u[x2,y2,t1] * (1 - x1w)) * (1 - y1w)
      u2 = (oscar$u[x1,y1,t2] * x1w + oscar$u[x2,y1,t2] * (1 - x1w)) * y1w + (oscar$u[x1,y2,t2] * x1w + oscar$u[x2,y2,t2] * (1 - x1w)) * (1 - y1w)
      u = u1 * t1w + u2 * (1 - t1w)

      ## calculate v
      v1 = (oscar$v[x1,y1,t1] * x1w + oscar$v[x2,y1,t1] * (1 - x1w)) * y1w + (oscar$v[x1,y2,t1] * x1w + oscar$v[x2,y2,t1] * (1 - x1w)) * (1 - y1w)
      v2 = (oscar$v[x1,y1,t2] * x1w + oscar$v[x2,y1,t2] * (1 - x1w)) * y1w + (oscar$v[x1,y2,t2] * x1w + oscar$v[x2,y2,t2] * (1 - x1w)) * (1 - y1w)
      v = v1 * t1w + v2 * (1 - t1w)

      grid$u[i] = u
      grid$v[i] = v

      if (is.na(u)) { grid$u[i] = 0 }
      if (is.na(v)) { grid$v[i] = 0 }
    }
  }
  list(u = matrix(grid$u, nrow = length(lon), ncol = length(lat)),
       v = matrix(grid$v, nrow = length(lon), ncol = length(lat)))
}



#' @title Get Velocity (ROMS)
#' @author Thomas Bryce Kelly
#' @export
get.vel.roms = function(lon, lat, depth, time, advection) {

  ## Build return object
  velocity = data.frame(lon = lon, lat = lat, depth = depth, time = time, u = NA, v = NA, w = NA)
  velocity$lon = velocity$lon %% 360

  if (time < min(advection$time)) {
    warning(paste0('Improper advection product loaded for time = ', time))
    t2 = 1
    t1 = 1
    t1w = 0.5
  } else if (time > max(advection$time)) {
    warning(paste0('Improper advection product loaded for time = ', time))
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
    yy = floor((index - 1) / dim(advection$lon)[1])
    xx = index %% dim(advection$lon)[1] + 1

    if (advection$lon[xx, yy] > velocity$lon[i]) { x1 = xx-1; x2 = xx } else { x1 = xx; x2 = xx + 1}
    if (advection$lat[xx, yy] > velocity$lat[i]) { y1 = yy-1; y2 = yy } else { y1 = yy; y2 = yy + 1}


    depths = apply(roms$z[c(x1, x2), c(y1,y2),], 3, mean)
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
      u = u1 * t1w + u2 * (1 - t1w)

      ## calculate v
      v1 = ((advection$v[x1,y1,z1,t1] * x1w + advection$v[x2,y1,z1,t1] * (1 - x1w)) * y1w + (advection$v[x1,y2,z1,t1] * x1w + advection$v[x2,y2,z1,t1] * (1 - x1w)) * (1 - y1w)) * z1w +
        ((advection$v[x1,y1,z1,t1] * x1w + advection$v[x2,y1,z1,t1] * (1 - x1w)) * y1w + (advection$v[x1,y2,z1,t1] * x1w + advection$v[x2,y2,z1,t1] * (1 - x1w)) * (1 - y1w)) * (1 - z1w)
      v2 = ((advection$v[x1,y1,z1,t2] * x1w + advection$v[x2,y1,z1,t2] * (1 - x1w)) * y1w + (advection$v[x1,y2,z1,t2] * x1w + advection$v[x2,y2,z1,t2] * (1 - x1w)) * (1 - y1w)) * z1w +
        ((advection$v[x1,y1,z1,t2] * x1w + advection$v[x2,y1,z1,t2] * (1 - x1w)) * y1w + (advection$v[x1,y2,z1,t2] * x1w + advection$v[x2,y2,z1,t2] * (1 - x1w)) * (1 - y1w)) * (1 - z1w)
      v = v1 * t1w + v2 * (1 - t1w)

      ## calculate w
      w1 = ((advection$w[x1,y1,z1,t1] * x1w + advection$w[x2,y1,z1,t1] * (1 - x1w)) * y1w + (advection$w[x1,y2,z1,t1] * x1w + advection$w[x2,y2,z1,t1] * (1 - x1w)) * (1 - y1w)) * z1w +
         ((advection$w[x1,y1,z1,t1] * x1w + advection$w[x2,y1,z1,t1] * (1 - x1w)) * y1w + (advection$w[x1,y2,z1,t1] * x1w + advection$w[x2,y2,z1,t1] * (1 - x1w)) * (1 - y1w)) * (1 - z1w)
      w2 = ((advection$w[x1,y1,z1,t2] * x1w + advection$w[x2,y1,z1,t2] * (1 - x1w)) * y1w + (advection$w[x1,y2,z1,t2] * x1w + advection$w[x2,y2,z1,t2] * (1 - x1w)) * (1 - y1w)) * z1w +
        ((advection$w[x1,y1,z1,t2] * x1w + advection$w[x2,y1,z1,t2] * (1 - x1w)) * y1w + (advection$w[x1,y2,z1,t2] * x1w + advection$w[x2,y2,z1,t2] * (1 - x1w)) * (1 - y1w)) * (1 - z1w)
      w = w1 * t1w + w2 * (1 - t1w)

      velocity$u[i] = u
      velocity$v[i] = v
      velocity$w[i] = w

      if (is.na(u)) { velocity$u[i] = 0 }
      if (is.na(v)) { velocity$v[i] = 0 }
      if (is.na(w)) { velocity$w[i] = 0 }
    }
  }
  list(u = velocity$u, v = velocity$v, w = velocity$w)
}




#######################
##### Advection Fn ####
#######################

#' @title Advect Lagrangian Particle
#' @author Thomas Bryce Kelly
#' @param lat Initial latitude of particle
#' @param lon Initial longitude of particle
#' @param depth Intial depth of particle
#' @param time Start time of particle
#' @param advection An advection list, such as one generated by `load.oscar.advection()`, that contains U V and W
#' @param dt Time step for RK4 advection scheme
#' @export
advect = function(lon, lat, depth, time, advection, dt, get.vel, sink = NULL, const.depth = NULL) {

  ## Initial Guess
  latt = lat[length(lat)]
  lont = lon[length(lon)]
  deptht = depth[length(depth)]
  timet = time[length(time)]
  vel = get.vel(lont, latt, deptht, timet, advection)

  ## p2 based on velocity at final pos at dt/2
  dlat = vel$v * 86.4 * dt/2 / 6378.14 # distance / earth radius
  lat2 = 180 / pi * asin(sin(latt / 180 * pi) * cos(dlat) + cos(latt / 180 * pi) * sin(dlat))
  dlon = vel$u * 86.4 * dt/2 / 6378.14 # distance / earth radius
  lon2 = lont + 180 / pi * atan2(sin(dlon) * cos(latt / 180 * pi), cos(dlon) - sin(latt / 180 * pi) * sin(lat2 / 180 * pi))
  ddepth = vel$w * dt/2 * 86400
  depth2 = deptht + ddepth
  time2 = timet + 86400 * dt/2
  vel2 = get.vel(lon2, lat2, depth2, time2, advection)

  ## p3 based on velocity at p2 at dt/2
  dlat = vel2$v * 86.4 * dt/2 / 6378.14 # distance / earth radius
  lat3 = 180 / pi * asin(sin(latt / 180 * pi) * cos(dlat) + cos(latt / 180 * pi) * sin(dlat))
  dlon = vel2$u * 86.4 * dt/2 / 6378.14 # distance / earth radius
  lon3 = lont + 180 / pi * atan2(sin(dlon) * cos(latt / 180 * pi), cos(dlon) - sin(latt / 180 * pi) * sin(lat3 / 180 * pi))
  ddepth = vel2$w * dt/2 * 86400
  depth3 = deptht + ddepth
  time3 = timet + 86400 * dt/2
  vel3 = get.vel(lon3, lat3, depth3, time3, advection)

  ## p4 based on velocity at p2 at dt/2
  dlat = vel$v * 86.4 * dt / 6378.14 # distance / earth radius
  lat4 = 180 / pi * asin(sin(latt / 180 * pi) * cos(dlat) + cos(latt / 180 * pi) * sin(dlat))
  dlon = vel$u * 86.4 * dt / 6378.14 # distance / earth radius
  lon4 = lont + 180 / pi * atan2(sin(dlon) * cos(latt / 180 * pi), cos(dlon) - sin(latt / 180 * pi) * sin(lat4 / 180 * pi))
  ddepth = vel$w * dt * 86400
  depth4 = deptht + ddepth
  time4 = timet + 86400 * dt
  vel4 = get.vel(lon4, lat4, depth4, time4, advection)

  ## RK4 velocity
  vel$u = (vel$u + 2 * vel2$u + 2 * vel3$u + vel4$u)/ 6
  vel$v = (vel$v + 2 * vel2$v + 2 * vel3$v + vel4$v)/ 6
  vel$w = (vel$w + 2 * vel2$w + 2 * vel3$w + vel4$w)/ 6

  dlat = vel$v * 86.4 * dt / 6378.14 # distance / earth radius
  lat.final = 180 / pi * asin(sin(latt / 180 * pi) * cos(dlat) + cos(latt / 180 * pi) * sin(dlat))
  dlon = vel$u * 86.4 * dt / 6378.14 # distance / earth radius
  lon.final = lont + 180 / pi * atan2(sin(dlon) * cos(latt / 180 * pi), cos(dlon) - sin(latt / 180 * pi) * sin(lat.final / 180 * pi))
  ddepth = vel$w * dt * 86400

  #### Sinking
  if(!is.null(sink)) {
    ddepth = ddepth + sink * dt * 86400 # sinking velocity is in [m s-1]
  }

  depth.final = deptht + ddepth

  if (!is.null(const.depth)) {depth.final = deptht} ## Constant depth float

  ## Append the latest timestep info and return:
  list(lon = c(lon, lon.final), lat = c(lat, lat.final), depth = c(depth, depth.final), time = c(time, time4))
}


