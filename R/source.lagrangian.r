
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

  list(u = u, v = v, w = w)
}




load.advection.oscar = function(file, lats = NULL, lons = NULL) {

  f = ncdf4::nc_open(file)
  oscar = f
  u = ncdf4::ncvar_get(f, 'u')
  v = ncdf4::ncvar_get(f, 'v')
  nc_close(f)

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
  }
  if (!is.null(lons)) {
    lons[lons<0] = lons[lons<0] + 360
    ## Filter
    l = which(oscar$lon >= lons[1] & oscar$lon <= lons[2])
    oscar$u = oscar$u[l,,]
    oscar$v = oscar$v[l,,]
    oscar$lon = oscar$lon[l]
    oscar$grid = expand.grid(lon = oscar$lon, lat = oscar$lat)
  }

  oscar
}



build.advection = function(lon, lat, depth, time) {

  ## Build empty structures
  u = array(0, dim = c(length(lon), lenght(lat), length(depth), length(time)))


  list(lon, lat, depth, time, u, v, w)
}

build.roms.grid = function(path) {
  list(lon, lat, depth, time, u, v, w)
}





get.vel = function(lon, lat, depth, time, advection) {

  ## Build return object
  grid = expand.grid(lon = lon, lat = lat, depth = depth, u = NA, v = NA, w = NA)
  grid$lon[grid$lon < 0] = grid$lon[grid$lon < 0] + 360

  if (time < min(advection$time) | time > max(advection$time)) {
    stop(paste0('Improper advection product loaded for time = ', time))
  }

  ## Time weights: linear in time
  t2 = min(which(advection$time > time))
  t1 = t2-1
  t1w = as.numeric(advection$time[t2] - time) / as.numeric(advection$time[t2] - advection$time[t1])

  for (i in 1:nrow(grid)) {
    ## Find Indicies
    x2 = min(which(advection$lon > grid$lon[i]))
    x1 = x2 - 1                         ## Makes assumptions about no boundaries
    y2 = max(which(advection$lat > grid$lat[i]))
    y1 = y2 + 1
    z2 = max(which(advection$depth > grid$depth[i]))
    z1 = z2 + 1


    if (x1 == 0 | y2 == 0 | z1 == 0 | t1 == 0 | x2 > length(advection$lon) | y1 > length(advection$lat) | z2 > length(advection$depth) | t2 > length(advection$time)){
      message('Point out of bounds! x1 = ', x1, ', y1 = ', y1, ', z1 = ', z1, ', t1 = ', t1)
      grid$u[i] = 0
      grid$v[i] = 0
      grid$w[i] = 0
    } else if (!is.null(depth) & (z1 == 0 | z2 > length(advection$depth))) {
      message('Point out of bounds! x1 = ', x1, ', y1 = ', y1, ', z1 = ', z1, ', t1 = ', t1)
      grid$u[i] = 0
      grid$v[i] = 0
      grid$w[i] = 0
    } else {

      ## Calculate weights
      x1w = (advection$lon[x2] - grid$lon[i]) / (advection$lon[x2] - advection$lon[x1])
      y1w = (advection$lat[y2] - grid$lat[i]) / (advection$lat[y2] - advection$lat[y1])
      z1w = (advection$depth[z2] - grid$depth[i]) / (advection$depth[z2] - advection$depth[z1])

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

      grid$u[i] = u
      grid$v[i] = v
      grid$w[i] = w

      if (is.na(u)) { grid$u[i] = 0 }
      if (is.na(v)) { grid$v[i] = 0 }
      if (is.na(w)) { grid$w[i] = 0 }
    }
  }
  list(u = grid$u, v = grid$v, w = grid$w)
}


advect = function(lat, lon, depth, time, advection, dt) {

  ## Initial Guess
  latt = lat[length(lat)]
  lont = lon[length(lon)]
  deptht = depth[length(depth)]
  timet = time[length(time)]
  vel = get.vel(latt, lont, deptht, timet, advection)

  ## p2 based on velocity at final pos at dt/2
  dlat = vel$v * 86.4 * dt/2 / 6378.14 # distance / earth radius
  lat2 = 180 / pi * asin(sin(latt / 180 * pi) * cos(dlat) + cos(latt / 180 * pi) * sin(dlat))
  dlon = vel$u * 86.4 * dt/2 / 6378.14 # distance / earth radius
  lon2 = lont + 180 / pi * atan2(sin(dlon) * cos(latt / 180 * pi), cos(dlon) - sin(latt / 180 * pi) * sin(lat2 / 180 * pi))
  ddepth = vel$w * dt/2
  depth2 = deptht + ddepth
  time2 = timet + 86400 * dt/2
  vel2 = get.vel(lat2, lon2, depth2, time2, advection)

  ## p3 based on velocity at p2 at dt/2
  dlat = vel2$v * 86.4 * dt/2 / 6378.14 # distance / earth radius
  lat3 = 180 / pi * asin(sin(latt / 180 * pi) * cos(dlat) + cos(latt / 180 * pi) * sin(dlat))
  dlon = vel2$u * 86.4 * dt/2 / 6378.14 # distance / earth radius
  lon3 = lont + 180 / pi * atan2(sin(dlon) * cos(latt / 180 * pi), cos(dlon) - sin(latt / 180 * pi) * sin(lat3 / 180 * pi))
  ddepth = vel2$w * dt/2
  depth3 = deptht + ddepth
  time3 = timet + 86400 * dt/2
  vel3 = get.vel(lat3, lon3, depth3, time3, advection)

  ## p4 based on velocity at p2 at dt/2
  dlat = vel$v * 86.4 * dt / 6378.14 # distance / earth radius
  lat4 = 180 / pi * asin(sin(latt / 180 * pi) * cos(dlat) + cos(latt / 180 * pi) * sin(dlat))
  dlon = vel$u * 86.4 * dt / 6378.14 # distance / earth radius
  lon4 = lont + 180 / pi * atan2(sin(dlon) * cos(latt / 180 * pi), cos(dlon) - sin(latt / 180 * pi) * sin(lat4 / 180 * pi))
  ddepth = vel$w * dt
  depth4 = deptht + ddepth
  time4 = timet + 86400 * dt
  vel4 = get.vel(lat4, lon4, depth4, time4, advection)

  ## RK4 velocity
  vel$u = (vel$u + 2 * vel2$u + 2 * vel3$u + vel4$u)/ 6
  vel$v = (vel$v + 2 * vel2$v + 2 * vel3$v + vel4$v)/ 6
  vel$w = (vel$w + 2 * vel2$w + 2 * vel3$w + vel4$w)/ 6

  dlat = vel$v * 86.4 * dt / 6378.14 # distance / earth radius
  lat.final = 180 / pi * asin(sin(latt / 180 * pi) * cos(dlat) + cos(latt / 180 * pi) * sin(dlat))
  dlon = vel$u * 86.4 * dt / 6378.14 # distance / earth radius
  lon.final = lont + 180 / pi * atan2(sin(dlon) * cos(latt / 180 * pi), cos(dlon) - sin(latt / 180 * pi) * sin(lat.final / 180 * pi))
  ddepth = vel$w * dt
  depth.final = deptht + ddepth

  ## Append the latest timestep info and return:
  list(lon = c(lon, lon.final), lat = c(lat, lat.final), depth = c(depth, depth.final), time = c(time, time4))
}


