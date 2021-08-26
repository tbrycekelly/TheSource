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


advection.files = list.files('Z:/Data/Roms Southern Ocean/daily/', full.names = T)


#' @title Get Velocity (ROMS)
#' @author Thomas Bryce Kelly
#' @export
get.vel.roms = function(lon, lat, depth, time, grid.data, advection.files, init.fields = c(1,2), verbose = T) {

  if (!'advection1' %in% ls(envir = .GlobalEnv)) {
    if (verbose) { message('Loading initial fields for advection1...')}
    advection1 = load.roms.physics(advection.files[init.fields[1]], verbose = verbose)
    advection1$load.index = init.fields[1]
    assign('advection1', advection1, envir = .GlobalEnv)
  }

  if (!'advection2' %in% ls(envir = .GlobalEnv)) {
    if (verbose) { message('Loading initial fields for advection2...')}
    advection2 = load.roms.physics(advection.files[init.fields[2]], verbose = verbose)
    advection2$load.index = init.fields[2]
    assign('advection2', advection2, envir = .GlobalEnv)
  }

  if (advection2$ocean_time < time) {
    while(advection2$ocean_time < time) {
      assign('advection1', .GlobalEnv$advection2, envir = .GlobalEnv)

      load.index = advection2$load.index + 1
      if (load.index > 0 & load.index <= length(advection.files)) {
        if (verbose) { message('Loading next file in sequence...')}
        advection2 = load.roms.physics(advection.files[load.index], verbose = F)
        advection2$load.index = load.index
        assign('advection2', advection2, envir = .GlobalEnv)
      } else {
        if (verbose) { message('End of file list, skipping loading of new time frame')}
        break
      }
    }
  }

  if (advection1$ocean_time > time) {
    while(advection1$ocean_time > time) {
      assign('advection2', .GlobalEnv$advection1, envir = .GlobalEnv)

      load.index = advection1$load.index - 1
      if (load.index > 0 & load.index <= length(advection.files)) {
        if (verbose) { message('Loading previous file in sequence...')}
        advection1 = load.roms.physics(advection.files[load.index], verbose = F)
        advection1$load.index = load.index
        assign('advection1', advection1, envir = .GlobalEnv)
      } else {
        if (verbose) { message('End of file list, skipping loading of new time frame')}
        break
      }
    }
  }

  ## Build return object
  velocity = data.frame(lon = lon, lat = lat, depth = depth, time = time, u = NA, v = NA, w = NA)
  velocity$lon = velocity$lon %% 360

  if (time < advection1$time) {
    message(' Improper advection product loaded for time = ', time, ' (too early)')
    t1w = 0.5
  } else if (time > advection2$time) {
    message('Improper advection product loaded for time = ', time, ' (too late)')
    t1w = 0.5
  } else {
    t1w = (as.numeric(advection2$time) - as.numeric(time)) / (as.numeric(advection2$time) - as.numeric(advection1$time))
  }

  for (i in 1:nrow(velocity)) {
    index = which.min((velocity$lon[i] - grid.data$lon_psi)^2 + (velocity$lat[i] - grid.data$lat_psi)^2)
    xx = index %% dim(advection$lon)[1]
    if (xx == 0) { xx = dim(advection$lon)[1] }
    yy = floor((index - 1) / dim(advection$lon)[1]) + 1

    if (grid.data$lon_psi[xx, yy] > velocity$lon[i]) { x1 = max(1, xx-1); x2 = xx } else { x1 = xx; x2 = min(xx + 1, dim(grid.data$lon_psi)[1])}
    if (grid.data$lat_psi[xx, yy] > velocity$lat[i]) { y1 = max(1, yy-1); y2 = yy } else { y1 = yy; y2 = min(yy + 1, dim(grid.data$lon_psi)[2])}

    grid.data$depth = calc.roms.depths(grid.data)
    depths = apply(grid.data$depth[c(x1, x2), c(y1,y2),], 3, function(x){ mean(x[!is.na(x)])})
    z2 = max(which(depths < velocity$depth[i]))
    z1 = min(z2 + 1, dim(grid.data$depth)[3])

    if (x1 == 0 | y2 == 0 | x2 > dim(grid.data$lon_psi)[1] | y1 > dim(grid.data$lat_psi)[2]){
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
      u1 = ((advection1$u[x1,y1,z1] * x1w + advection1$u[x2,y1,z1] * (1 - x1w)) * y1w + (advection1$u[x1,y2,z1] * x1w + advection1$u[x2,y2,z1] * (1 - x1w)) * (1 - y1w)) * z1w +
        ((advection$u1[x1,y1,z1] * x1w + advection1$u[x2,y1,z1] * (1 - x1w)) * y1w + (advection1$u[x1,y2,z1] * x1w + advection1$u[x2,y2,z1] * (1 - x1w)) * (1 - y1w)) * (1 - z1w)
      u2 = ((advection2$u[x1,y1,z2] * x1w + advection2$u[x2,y1,z2] * (1 - x1w)) * y1w + (advection2$u[x1,y2,z2] * x1w + advection2$u[x2,y2,z2] * (1 - x1w)) * (1 - y1w)) * z1w +
        ((advection2$u[x1,y1,z2] * x1w + advection2$u[x2,y1,z2] * (1 - x1w)) * y1w + (advection2$u[x1,y2,z2] * x1w + advection2$u[x2,y2,z2] * (1 - x1w)) * (1 - y1w)) * (1 - z1w)
      velocity$u[i] = u1 * t1w + u2 * (1 - t1w)

      ## calculate v
      v1 = ((advection1$v[x1,y1,z1] * x1w + advection1$v[x2,y1,z1] * (1 - x1w)) * y1w + (advection1$v[x1,y2,z1] * x1w + advection1$v[x2,y2,z1] * (1 - x1w)) * (1 - y1w)) * z1w +
        ((advection1$v[x1,y1,z1] * x1w + advection1$v[x2,y1,z1] * (1 - x1w)) * y1w + (advection1$v[x1,y2,z1] * x1w + advection1$v[x2,y2,z1] * (1 - x1w)) * (1 - y1w)) * (1 - z1w)
      v2 = ((advection2$v[x1,y1,z1] * x1w + advection2$v[x2,y1,z1] * (1 - x1w)) * y1w + (advection2$v[x1,y2,z1] * x1w + advection2$v[x2,y2,z1] * (1 - x1w)) * (1 - y1w)) * z1w +
        ((advection2$v[x1,y1,z1] * x1w + advection2$v[x2,y1,z1] * (1 - x1w)) * y1w + (advection2$v[x1,y2,z1] * x1w + advection2$v[x2,y2,z1] * (1 - x1w)) * (1 - y1w)) * (1 - z1w)
      velocity$v[i] = v1 * t1w + v2 * (1 - t1w)

      ## calculate w
      w1 = ((advection1$w[x1,y1,z1] * x1w + advection1$w[x2,y1,z1] * (1 - x1w)) * y1w + (advection1$w[x1,y2,z1] * x1w + advection1$w[x2,y2,z1] * (1 - x1w)) * (1 - y1w)) * z1w +
        ((advection1$w[x1,y1,z1] * x1w + advection1$w[x2,y1,z1] * (1 - x1w)) * y1w + (advection1$w[x1,y2,z1] * x1w + advection1$w[x2,y2,z1] * (1 - x1w)) * (1 - y1w)) * (1 - z1w)
      w2 = ((advection2$w[x1,y1,z1] * x1w + advection2$w[x2,y1,z1] * (1 - x1w)) * y1w + (advection2$w[x1,y2,z1] * x1w + advection2$w[x2,y2,z1] * (1 - x1w)) * (1 - y1w)) * z1w +
        ((advection2$w[x1,y1,z1] * x1w + advection2$w[x2,y1,z1] * (1 - x1w)) * y1w + (advection2$w[x1,y2,z1] * x1w + advection2$w[x2,y2,z1] * (1 - x1w)) * (1 - y1w)) * (1 - z1w)
      velocity$w[i] = w1 * t1w + w2 * (1 - t1w)
    }
  }
  velocity[is.na(velocity)] = 0

  ## Return
  list(u = velocity$u, v = velocity$v, w = velocity$w)
}

