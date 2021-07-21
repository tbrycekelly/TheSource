## Functions for building gridded datasets from measurements.
##
## Author: Thomas Bryce Kelly (tkelly@alaska.edu)
## http://about.tkelly.org/
##
## College of Fisheries and Ocean Sciences
## University of Alaska Fairbanks
##
## Dept of Earth, Ocean & Atmospheric Sciences
## Florida State University


#######################
### Section Building ##
#######################

#' @title Build Section
#' @author Thomas Bryce Kelly
#' @keywords Gridding
#' @export
#' @param x dimensions (e.g. lat, lon, depth, section distance, time, etc)
#' @param y dimensions (e.g. lat, lon, depth, section distance, time, etc)
#' @param z signal to be gridded (e.g. T, S, ...)
#' @param xlim Limits of the gridding. These are the bounds of the new x-y grid. Default: NULL will set it based on the data + 10%.
#' @param ylim Limits of the gridding. These are the bounds of the new x-y grid. Default: NULL will set it based on the data + 10%.
#' @param x.factor The relative scale difference between x and y, used to calculate distances. Take into account actual scale AND the relevent scaling of the system (vertical distance tends to be more important than horizontal distance).
#' @param y.factor The relative scale difference between x and y, used to calculate distances. Take into account actual scale AND the relevent scaling of the system (vertical distance tends to be more important than horizontal distance).
#' @param x.scale The step size in the new x-y grid. By default the scale is set to generate a grid that is 50x50.
#' @param y.scale The step size in the new x-y grid. By default the scale is set to generate a grid that is 50x50.
#' @param uncertainty = 0: Scaling applied to the distance from the cener of a grid cell to a vertex, used to add a base-line distance to all measurements. 0 = no minimum, 1 = minimum = to half a grid cell.
#' @param field.name Sets the name of the new interpolated field. By default the name is 'z1'
#' @param gridder A function to perform gridding, options gridIDW (default: inverse distance), gridNN (nearest neighbor), gridNNI (natural neighbor) or gridKrige (Krigging)
#' @param nx The number of splits to make in the x direction (defaults to 50). Used only if x.scale is not set.
#' @param ny The number of splits to make in the y direction (defaults to 50). Used only if y.scale is not set.
#' @param proj A gdal projection string such as from make.proj() function. Used to redistribute grid over non-euclidean surfaces like maps.
build.section = function(x, y, z, lat = NULL, lon = NULL, grid = NULL, weight = NULL,
                         xlim = NULL, ylim = NULL,
                         x.factor = NULL, y.factor = NULL,
                         x.scale = NULL, y.scale = NULL,
                         uncertainty = 1e-12, p = 2, gridder = NULL,
                         field.names = NULL, nx = 50, ny = 50,
                         proj = NULL, verbose = T) {

  if (verbose) {message('BUILD.SECTION: Starting section building process (verbose = T).')}
  time.a = Sys.time()
  z = data.matrix(z)

  ## Remove NAs
  l = !is.na(x) & !is.na(y) & !apply(z, 1, function(x) {any(is.na(x))})
  x = x[l]
  y = y[l]
  z = data.matrix(z[l,])
  if (is.null(weight)) {weight = matrix(1, nrow = nrow(z), ncol = ncol(z))}
  if (!is.null(lat)) {lat = lat[l]}
  if (!is.null(lon)) {lon = lon[l]}
  if (is.null(gridder)) {
    gridder = gridIDW
    if (verbose) { message('BUILD.SECTION: No gridder specified, defaulting to gridIDW. Other options: gridNN, gridNNI and gridKrige.') }
  }

  if (uncertainty == 0) { message('BUILD.SECTION: Uncertainty of zero may produce NAs!') }
  if (is.null(field.names)) {
    field.names = paste0('z', 1:ncol(z))
    if (verbose) { message('BUILD.SECTION: No field.names provided, gridded data will be called ', paste0('z', 1:ncol(z), collapse = ',')) }
  }

  if (class(x[1])[1] == 'POSIXct'){
    if (verbose) { message('X axis is time.') }
    x = as.numeric(x)
    t.axis = T
  } else {
    t.axis = F
  }

  ## Set default limits (+10% buffer)
  if (is.null(xlim)) {
    xlim = range(x)
    xlim[1] = xlim[1] - (xlim[2] - xlim[1])/20
    xlim[2] = xlim[2] + (xlim[2] - xlim[1])/20
  }
  if (is.null(ylim)) {
    ylim = range(y)
    ylim[1] = ylim[1] - (ylim[2] - ylim[1])/20
    ylim[2] = ylim[2] + (ylim[2] - ylim[1])/20
  }

  if (is.null(x.scale)) { x.scale = (xlim[2] - xlim[1]) / nx} ## Default to nx or ny steps
  if (is.null(y.scale)) { y.scale = (ylim[2] - ylim[1]) / ny}

  ## Rescale x and y based on x.factor and y.factor
  if (is.null(x.factor)) { x.factor = (x.scale/y.scale + 1) / 2}
  if (is.null(y.factor)) { y.factor = (y.scale/x.scale + 1) / 2}
  x = x * x.factor
  x.scale = x.scale * x.factor
  xlim = xlim * x.factor
  y = y * y.factor
  y.scale = y.scale * y.factor
  ylim = ylim * y.factor


  if (y.scale == 0) {
    y.new = ylim[1]
  } else {
    y.new = seq(ylim[1], ylim[2], by = y.scale)
  }

  if (x.scale == 0) {
    x.new = xlim[1]
  } else {
    x.new = seq(xlim[1], xlim[2], by = x.scale)
  }

  if (!is.null(lat)) {
    if (length(unique(x)) > 1) {
      section.lat = approx(x, lat, xout = x.new, rule = 2)$y
    } else {
      section.lat = rep(lat, length(x.new))
    }
  } else {
    section.lat = rep(NA, length(x)); lat = NA
  }
  if (!is.null(lon)) {
    if (length(unique(x)) > 1) { ## interpolate
      section.lon = approx(x, lon, xout = x.new, rule = 2)$y
    } else {
      section.lon = rep(lon, length(x.new))
    }
  } else {
    section.lon = rep(NA, length(x)); lon = NA
  }

  ## Make grid and fill in
  if (is.null(grid)) {
    grid = expand.grid(x = x.new, y = y.new)
    nx = length(x.new)
    ny = length(y.new)
    if (verbose) {message('BUILD.SECTION: Building grid with ', nrow(grid), ' positions and ', length(z), ' observations.')}
  } else {
    if (verbose) {message('BUILD.SECTION: Grid parameter provided, ignoring nx, ny, xlim, ylim since we are not building a new grid with ', nrow(grid), ' entries.')}
    colnames(grid) = c('x', 'y')
    nx = length(unique(grid$x))
    ny = length(unique(grid$y))
    x.new = unique(grid$x)
    y.new = unique(grid$y)
    xlim = range(grid$x)
    ylim = range(grid$y)

    grid$x = grid$x * x.factor
    grid$y = grid$y * y.factor

  }

  if (!is.null(proj)) { ## Apply projection
    ## Project x and y
    projected = project(cbind(x, y), proj = proj)
    x = as.numeric(projected[,1])
    y = as.numeric(projected[,2])

    ## proejct grid x and y
    projected = project(cbind(grid$x, grid$y), proj = proj)
    grid$x = as.numeric(projected[,1])
    grid$y = as.numeric(projected[,2])
  }

  time.b = Sys.time()
  for (kk in 1:length(field.names)) {
    if (verbose) {message('BUILD.SECTION: Building grid for field ', field.names[kk], '  ', Sys.time(), '.')}
    grid[[field.names[kk]]] = gridder(grid$x, grid$y, x, y, z[,kk], p, x.scale/x.factor, y.scale/y.factor, uncertainty)
  }

  time.c = Sys.time()

  grid$x = grid$x / x.factor
  grid$y = grid$y / y.factor

  if (t.axis) {
    x = conv.time.unix(x)
    grid$x = conv.time.unix(grid$x)
  }

  ## Reconstruct z
  z = data.frame(z)
  colnames(z) = field.names

  if (!is.null(proj)) { ## Apply projection
    ## Fix data x and y
    projected = project(cbind(x, y), proj = proj, inv = T)
    x = as.numeric(projected[,1])
    y = as.numeric(projected[,2])

    ## Fix new x and y
    projected = project(cbind(x.new, y.new), proj = proj, inv = T)
    x.new = as.numeric(projected[,1])
    y.new = as.numeric(projected[,2])

    ##Fix grid
    projected = project(cbind(grid$x, grid$y), proj = proj, inv = T)
    grid$x = as.numeric(projected[,1])
    grid$y = as.numeric(projected[,2])

  }

  ## Construct return object
  grid = list(grid = grid,
              grid.meta = list(
                x.scale = x.scale / x.factor,
                y.scale = y.scale / y.factor,
                x.factor = x.factor,
                y.factor = y.factor,
                nx = nx, ny = ny,
                uncertainty = uncertainty,
                p = p,
                gridder = gridder
              ),
              x = x.new / x.factor,
              y = y.new / y.factor,
              section.lat = section.lat,
              section.lon = section.lon,
              data = cbind(x = x / x.factor, y = y / y.factor, z = z, lat = lat, lon = lon)
  )

  if (verbose) {message('BUILD.SECTION Timings\n Total function time: \t ', time.c - time.a, '\n Preprocessing Time:\t', time.b - time.a, '\n Gridding Time:\t', time.c - time.b)}
  ## Return
  grid
}


#' @title Build Section with Parallel Processing
#' @author Thomas Bryce Kelly
#' @keywords Gridding
#' @export
#' @param x dimensions (e.g. lat, lon, depth, section distance, time, etc)
#' @param y dimensions (e.g. lat, lon, depth, section distance, time, etc)
#' @param z signal to be gridded (e.g. T, S, ...)
#' @param xlim Limits of the gridding. These are the bounds of the new x-y grid. Default: NULL will set it based on the data + 10%.
#' @param ylim Limits of the gridding. These are the bounds of the new x-y grid. Default: NULL will set it based on the data + 10%.
#' @param x.factor The relative scale difference between x and y, used to calculate distances. Take into account actual scale AND the relevent scaling of the system (vertical distance tends to be more important than horizontal distance).
#' @param y.factor The relative scale difference between x and y, used to calculate distances. Take into account actual scale AND the relevent scaling of the system (vertical distance tends to be more important than horizontal distance).
#' @param x.scale The step size in the new x-y grid. By default the scale is set to generate a grid that is 50x50.
#' @param y.scale The step size in the new x-y grid. By default the scale is set to generate a grid that is 50x50.
#' @param uncertainty = 0: Scaling applied to the distance from the cener of a grid cell to a vertex, used to add a base-line distance to all measurements. 0 = no minimum, 1 = minimum = to half a grid cell.
#' @param field.name Sets the name of the new interpolated field. By default the name is 'z1'
#' @param gridder A function to perform gridding, options gridIDW (default: inverse distance), gridNN (nearest neighbor), gridNNI (natural neighbor) or gridKrige (Krigging)
#' @param nx The number of splits to make in the x direction (defaults to 50). Used only if x.scale is not set.
#' @param ny The number of splits to make in the y direction (defaults to 50). Used only if y.scale is not set.
build.section.parallel = function(x, y, z, lat = NULL, lon = NULL,
                                  xlim = NULL, ylim = NULL,
                                  x.factor = 1, y.factor = 1,
                                  x.scale = NULL, y.scale = NULL,
                                  uncertainty = 1e-12, p = 3, gridder = NULL,
                                  field.names = NULL, nx = 50, ny = 50, verbose = T) {

  z = data.matrix(z)
  ## Remove NAs
  l = !is.na(x) & !is.na(y) & !apply(z, 1, function(x) {any(is.na(x))})
  x = x[l]
  y = y[l]
  z = data.matrix(z[l,])
  if (!is.null(lat)) {lat = lat[l]}
  if (!is.null(lon)) {lon = lon[l]}
  if (is.null(gridder)) {
    gridder = gridIDW
    message('No gridder specified, defaulting to gridIDW. Other options: gridNN, gridNNI and gridKrige.')
  }


  if (uncertainty == 0) { warning('Uncertainty of zero may produce NAs!') }
  if (is.null(field.names)) {
    field.names = paste0('z', 1:ncol(z))
    warning('No field.names provided, gridded data will be called ', paste0('z', 1:ncol(z), collapse = ','))
  }

  ## Set default limits (+10% buffer)
  if (is.null(xlim)) {
    xlim = range(x)
    xlim[1] = xlim[1] - (xlim[2] - xlim[1])/20
    xlim[2] = xlim[2] + (xlim[2] - xlim[1])/20
  }
  if (is.null(ylim)) {
    ylim = range(y)
    ylim[1] = ylim[1] - (ylim[2] - ylim[1])/20
    ylim[2] = ylim[2] + (ylim[2] - ylim[1])/20
  }

  if (is.null(x.scale)) { x.scale = (xlim[2] - xlim[1]) / nx} ## Default to nx or ny steps
  if (is.null(y.scale)) { y.scale = (ylim[2] - ylim[1]) / ny}

  ## Rescale x and y based on x.factor and y.factor
  x = x * x.factor
  x.scale = x.scale * x.factor
  xlim = xlim * x.factor
  y = y * y.factor
  y.scale = y.scale * y.factor
  ylim = ylim * y.factor


  y.new = seq(ylim[1], ylim[2], by = y.scale)
  x.new = seq(xlim[1], xlim[2], by = x.scale)

  if (!is.null(lat)) { section.lat = approx(x, lat, xout = x.new, rule = 2)$y } else { section.lat = rep(NA, length(x)); lat = NA }
  if (!is.null(lon)) { section.lon = approx(x, lon, xout = x.new, rule = 2)$y } else { section.lon = rep(NA, length(y)); lon = NA }

  ## Make grid and fill in
  grid = expand.grid(x = x.new, y = y.new)

  # Start Parallel
  cn = parallel::detectCores() - 1
  cl = parallel::makeCluster(cn)
  ### PASS THE OBJECT FROM MASTER PROCESS TO EACH NODE
  parallel::clusterExport(cl,
                          varlist = c('grid', 'x', 'y', 'z', 'p', 'x.scale', 'y.scale', 'uncertainty', deparse(substitute(gridIDW))),
                          envir = environment())

  ### DIVIDE THE DATAFRAME BASED ON # OF CORES
  sp = parallel::parLapply(cl, parallel::clusterSplit(cl = cl, seq = seq(nrow(grid))), function(c) {grid[c,]})

  for (kk in 1:length(field.names)) {
    grid[[field.names[kk]]] = Reduce(c, parallel::parLapply(cl, sp, fun = function(s) {
      gridder(s$x, s$y, x, y, z[,kk], p, x.scale, y.scale, uncertainty)
    }))
  }

  parallel::stopCluster(cl)

  grid$x = grid$x / x.factor
  grid$y = grid$y / y.factor

  ## Construct return object
  grid = list(grid = grid,
              grid.meta = list(
                x.scale = x.scale / x.factor,
                y.scale = y.scale / y.factor,
                x.factor = x.factor,
                y.factor = y.factor,
                nx = nx, ny = ny,
                uncertainty = uncertainty,
                p = p,
                gridder = gridder
              ),
              x = x.new / x.factor,
              y = y.new / y.factor,
              section.lat = section.lat,
              section.lon = section.lon,
              data = data.frame(x = x / x.factor, y = y / y.factor, z = z, lat = lat, lon = lon)
  )
  ## Return
  grid
}


#### Alternative Gridding Engines


#' @title Grid via Binning
#' @description Bins data into reguarly spaced grid
#' @param gx Grid x values to interpolate onto
#' @param gy Grid y values to interpolate onto
#' @param x Observations, x values
#' @param y Observations, y values
#' @param z Observations, z values
#' @param p Unused
#' @param xscale The spacing of the x grid
#' @param yscale The spacing of the y grid
#' @param uncertainty Unused
#' @author Thomas Bryce Kelly
#' @export
gridBin = function(gx, gy, x, y, z, p, xscale, yscale, uncertainty) {
  gz = rep(NA, length(gx))

  if (xscale == 0) { xscale = Inf }
  if (yscale == 0) { yscale = Inf }

  for (i in 1:length(gz)) {
    gz[i] = mean(z[which(abs(gx[i] - x) < xscale/2 & abs(gy[i] - y) < yscale/2)])
  }
  gz ## Return
}


#' @title Grid via Natural Neighbor Interpolation
#' @author Thomas Bryce Kelly
#' @param gx Grid x values to interpolate onto
#' @param gy Grid y values to interpolate onto
#' @param x Observations, x values
#' @param y Observations, y values
#' @param z Observations, z values
#' @param p Exponent on the distance function
#' @inheritParams gridIDW
#' @export
gridNNI = function(gx, gy, x, y, z, p, xscale, yscale, uncertainty) {

  gz = rep(NA, length(gx))

  ## Build standard NN index field
  gtemp1 = rep(1, length(gx))
  for (i in 1:length(gx)) {
    gtemp1[i] = which.min(abs(gx[i] - x)^p + abs(gy[i] - y)^p)
  }

  for (i in 1:length(gx)) {

    ## Setup
    xx = c(x, gx[i])
    yy = c(y, gy[i])

    gtemp2 = gtemp1 ## Assume no change

    ## Find all the values that are likely to change
    dd = abs(x - gx[i])^p + abs(y - gy[i])^p
    l = which(abs(gx - gx[i])^p + abs(gy - gy[i])^p < dd[order(dd)[min(5, length(dd))]] / 2) ## only look at points within circle defined by the 4th closest data point

    for (k in l) {
      gtemp2[k] = which.min(abs(gx[k] - xx)^p + abs(gy[k] - yy)^p)
    }
    w = gtemp2 - gtemp1
    w[w != 0] = 1
    w = w / sum(w)## weight based on number of entries that were changed
    gz[i] = sum(w * z[gtemp1])
  }

  gz ## Return
}



#' @title Grid via Krigging
#' @author Thomas Bryce Kelly
#' @inheritParams gridIDW
#' @param gx Grid x values to interpolate onto
#' @param gy Grid y values to interpolate onto
#' @param x Observations, x values
#' @param y Observations, y values
#' @param z Observations, z values
#' @param p Unused
#' @export
#' @import automap
#' @importFrom sp coordinates
gridKrig = function(gx, gy, x, y, z, p, xscale, yscale, uncertainty) {

  data = data.frame(x = x, y = y, z = z)
  grid = data.frame(x = gx, y = gy)

  ## Set coordinate system (SpatialPointsDataFrame)
  sp::coordinates(data) = ~x + y
  sp::coordinates(grid) = ~x + y
  krige = automap::autoKrige(z ~ 1, data, new_data = grid)

  ## Return
  krige$krige_output$var1.pred
}


#' @title Add Section Parameter
#' @author Thomas Bryce Kelly
#' @keywords Gridding
#' @param x x
#' @export
## Add additional paramter to section object
# x, y: dimensions (e.g. lat, lon, depth, section distance, time, etc)
# z: signal to be gridded (e.g. T, S, ...)
# section: The existing section object (2D only) to which the new parameter will be added.
# field.name = 'z': Sets the name of the new interpolated field. By default the name is 'z', but setting it to 'T' for temperature makes sense.
add.section.param = function(x, y, z, section, field.name) {
  ## Remove NAs
  l = !is.na(x) & !is.na(y) & !is.na(z)
  x = x[l]; y = y[l]; z = z[l]

  section$grid[[field.name]] = section$grid.meta$gridder(grid, x, y, z, p,
                                                         section$grid.meta$x.scale,
                                                         section$grid.meta$y.scale,
                                                         section$grid.meta$uncertainty)
  section
}


#' @title Run Jackknife Analysis on Section
#' @author Thomas Bryce Kelly
#' @export
jackknife.section = function(section, n) {

  ## Get grid parameters
  xlim = section$grid.meta$xlim
  ylim = section$grid.meta$ylim
  x.scale = section$grid.meta$x.scale
  y.scale = section$grid.meta$y.scale
  x.factor = section$grid.meta$x.factor
  y.factor = section$grid.meta$y.factor
  uncertainty = section$grid.meta$uncertainty

  ## Minimum distance based on grid size
  delta.min = delta(x.factor, y.factor, x.scale/2, 0, y.scale/2, 0, p = section$grid.meta$p) * uncertainty

  section$jack = expand.grid(x = section$x, y = section$y)

  for (k in 1:n) {
    l = sample(1:nrow(section$data), nrow(section$data), replace = TRUE)
    l = unique(l)
    section$jack[[paste0('z', k)]] = NA

    temp = rep(0, nrow(grid))
    for (i in 1:nrow(section$grid)) {
      ## Calculate distances and determine weights
      del = delta(x.factor, y.factor, section$data$x[l], section$grid$x[i], section$data$y[l], section$grid$y[i])

      w = W(del, delta.min) # Get weights
      temp[i] = sum(section$data$z[l] * w) / sum(w)
    }
    section$jack[[paste0('z', k)]] = temp
  }

  ## Return
  section
}


#' @title Cross Validate the IDW
#' @author Thomas Bryce Kelly
#' @export
cv.section = function(section) {


  for (j in 3:ncol(section$grid)) {

    ## initialize
    section$grid[[paste0(colnames(section$grid)[j], '.var')]] = 0

    ## Cross validate
    for (i in 1:nrow(section$data)) {
      temp = build.section(x = section$data$x[-i], y = section$data$y[-i], z = section$data[-i,j],
                           grid = section$grid[,c(1,2)], uncertainty = section$grid.meta$uncertainty,
                           p = section$grid.meta$p, gridder =  section$grid.meta$gridder, verbose = F)
      ## Calc SSR
      section$grid[[paste0(colnames(section$grid)[j], '.var')]] = section$grid[[paste0(colnames(section$grid)[j], '.var')]] + (temp$grid$z1 - section$grid[,j])^2 / nrow(section$data)
    }
  }
  section
}


#' @title Build Section (3D)
#' @author Thomas Bryce Kelly
#' @description A function to build a 3d section using IDW interpolations. WARNING: Not currently maintained.
#' @keywords Oceanography Section Interpolation
#' @inheritParams build.section
#' @export
## Build a 3D section object
# x, y, z: dimensions (e.g. lat, lon, depth, section distance, time, etc)
# zz: signal to be gridded (e.g. T, S, ...)
# xlim, ylim, zlim = NULL: Limits of the gridding. These are the bounds of the new x-y grid. Default: NULL will set it based on the data + 10%.
# x.factor, y.factor, z.factor = 1: The relative scale difference between x and y, used to calculate distances. Take into account actual scale AND the relevent scaling of the system (vertical distance tends to be more important than horizontal distance).
# x.scale, y.scale, z.scale = NULL: The step size in the new x-y grid. By default the scale is set to generate a grid that is 50x50.
# uncertainty = 0: Scaling applied to the distance from the cener of a grid cell to a vertex, used to add a base-line distance to all measurements. 0 = no minimum, 1 = minimum = to half a grid cell.
# field.name = 'zz': Sets the name of the new interpolated field. By default the name is 'zz', but setting it to 'T' for temperature makes sense.
build.section.3d = function(x, y, z, zz, xlim = NULL, ylim = NULL, zlim = NULL, x.factor = 1, y.factor = 1, z.factor = 1,
                            x.scale = NULL, y.scale = NULL, z.scale = NULL, uncertainty = 1e-12,
                            field.name = 'zz', p = 3) {

  warning('Function is not currently maintaned. Use at your own risk.')
  ## Remove NAs
  l = !is.na(x) & !is.na(y) & !is.na(z) & !is.na(zz)
  x = x[l]
  y = y[l]
  z = z[l]
  zz = zz[l]

  if (uncertainty == 0) { warning('Uncertainty of zero may produce NAs!') }

  if (is.null(xlim)) {
    xlim = range(x)
    xlim[1] = xlim[1] - (xlim[2]-xlim[1])/20
    xlim[2] = xlim[2] + (xlim[2]-xlim[1])/20
  }
  if (is.null(ylim)) {
    ylim = range(y)
    ylim[1] = ylim[1] - (ylim[2]-ylim[1])/20
    ylim[2] = ylim[2] + (ylim[2]-ylim[1])/20
  }
  if (is.null(zlim)) {
    zlim = range(z)
    zlim[1] = zlim[1] - (zlim[2]-zlim[1])/20
    zlim[2] = zlim[2] + (zlim[2]-zlim[1])/20
  }

  if (is.null(x.scale)) { x.scale = (xlim[2] - xlim[1]) / 50} ## Default to 50 steps
  if (is.null(y.scale)) { y.scale = (ylim[2] - ylim[1]) / 50} ## Default to 50 steps
  if (is.null(z.scale)) { z.scale = (zlim[2] - zlim[1]) / 50} ## Default to 50 steps

  ## field = matrix (1 - 300 m depth x 0 - max dist)
  n.x = (xlim[2] - xlim[1]) / x.scale
  n.y = (ylim[2] - ylim[1]) / y.scale
  n.z = (zlim[2] - zlim[1]) / z.scale

  x.new = seq(xlim[1], xlim[2], by = x.scale)
  y.new = seq(ylim[1], ylim[2], by = y.scale)
  z.new = seq(zlim[1], zlim[2], by = z.scale)

  ## Minimum distance based on grid size
  delta.min = delta3(x.factor, y.factor, z.factor, x.scale/2, 0, y.scale/2, 0, z.scale/2, 0) * uncertainty

  grid = expand.grid(x = x.new, y = y.new, z = z.new)

  temp = rep(0, nrow(grid))
  for (i in 1:nrow(grid)) {

    ## Calculate distances and determine weights
    del = delta3(x.factor, y.factor, z.factor, x, grid$x[i], y, grid$y[i], z, grid$z[i])

    w = W(del, delta.min) # Get weights
    temp[i] = sum(zz * w) / sum(w)
  }
  grid[[field.name]] = temp

  ## Construct return object
  grid = list(grid = grid,
              grid.meta = list(
                x.scale = x.scale,
                y.scale = y.scale,
                z.scale = z.scale,
                x.factor = x.factor,
                y.factor = y.factor,
                z.factor = z.factor,
                n.x = n.x, n.y = n.y, n.z = n.z,
                uncertainty = uncertainty,
                p = p,
                gridder = build.section.3d
              ),
              x = x.new,
              y = y.new,
              z = z.new,
              data = data.frame(x = x, y = y, z = z, zz = zz)
  )
  ## Return
  grid
}

#########################
##### Plotting ##########
#########################

#' @title Plot Section
#' @author Thomas Bryce Kelly
#' @description Plot a section object.
#' @keywords Section Plotting
#' @export
#' @param section: section object (2D only) to be plotted.
#' @param  field Parameter name to be plotted (z-axis).
#' @param  xlim Limits of the plot, by default the limits are set to the limits of the section grid.
#' @param  tlim Limits of the plot, by default the limits are set to the limits of the section grid.
#' @param  xlab X-axis label.
#' @param  ylab Y-axis label.
#' @param  log Flag for turning on log transformation of the z-axis. If TRUE then z = log(z) is performed prior to plotting.
#' @param  base Used when log = TRUE to set the base of the log.
#' @param  zlim The zlim imposed on the z-axis, which by default is set to the range of z values in the data.
#' @param  pal The color palette used; should be modeled after get.pal().
#' @param  rev Boolean for reversing the palette colors.
#' @param  include.data Flag for whether the data used to construct the grid should be plotted using the same palette.
#' @param  mark.points Flag for marking the location of each sample used to construct the grid.
#' @param  include.pch The pch used for when mark.points = TRUE.
#' @param  include.cex Sets the point size for when either include.data or mark.points are TRUE. (Used for both).
#' @param  main The title text to be included in the top line of the plot. Defaults to the name of the field.
#' @param  col.low [unimplemented] The color used when plotting out-of-range z values on the low end. Default value of '' indicates to use the minimum value of pal(). A value of NA skips the plotting of these data. Otherwise the color given is used (e.g. col.low = 'blue').
#' @param  col.high [unimplemented] Same as col.low but for out-of-range high values.
plot.section = function(section, field = NULL, xlim = NULL, ylim = NULL, xlab = 'x', ylab = 'y', log = FALSE, base = 10,
                        zlim = NULL, pal = 'greyscale', rev = FALSE, include.data = FALSE, mark.points = FALSE, include.pch = 21,
                        include.cex = 1, main = NULL, col.low = '', col.high = '', N = 255) {

  ## Set Defaults
  # Check if values need to set and set them.
  if (is.null(field)) {
    field = colnames(section$grid)[3]
    warning('No field name provided, using first gridded data: ', field)
  }

  ## Handy variables
  x = section$x
  y = section$y
  z = matrix(section$grid[,field], nrow = length(x))

  if (log) {
    z = log(z, base)
    if (is.null(zlim)) {
      zlim = range(pretty(z, na.rm = TRUE))
    }
  }

  if (is.null(zlim)) { zlim = range(pretty(z, na.rm = TRUE)) }

  if (is.null(main)) { main = field }
  if (is.null(xlim)) { xlim = range(x)}
  if (is.null(ylim)) { ylim = range(y) }

  ## Out of Range
  # Set values to zlim if you want the out of range values plotted at the zlim values
  if (!is.na(col.low) & col.low == '') {
    z[z < zlim[1]] = zlim[1]
  }
  if (!is.na(col.high) & col.high == '') {
    z[z > zlim[2]] = zlim[2]
  }

  ## Plot iamge
  image(x = x, y = y, z = z, col = get.pal(N, pal = pal, rev = rev), ylab = ylab, xlab = xlab,
        xlim = xlim, ylim = ylim, zlim = zlim)


  ## Add Title text
  if (log) {
    st = paste0(main, '   zlim: (', round(base^zlim[1], 3), ', ', round(base^zlim[2],3), ')')
  } else {
    st = paste0(main, '   zlim: (', round(zlim[1], 3), ', ', round(zlim[2],3), ')')
  }
  mtext(st, line = 0.25, adj = 1, cex = 0.7)

  # Plot points that are out of range when a color is given
  if (!is.na(col.low) & col.low != '') {
    zz = z
    zz[zz >= zlim[1]] = NA
    image(x = x, y = y, z = zz, col = col.low, add = TRUE)#, pch = include.pch, cex = include.cex)
  }

  if (!is.na(col.high) & col.high != '') {
    zz = z
    zz[zz <= zlim[2]] = NA
    image(x = x, y = y, z = zz, col = col.high, add = TRUE)#, pch = include.pch, cex = include.cex)
  }

  if (mark.points) {
    if (include.data) {
      points(x = section$data$x, y = section$data$y, pch = include.pch, cex = include.cex, col = 'black',
             bg = make.pal(x = section$data$z, pal = pal, n = N, min = zlim[1], max = zlim[2], rev = rev))
    } else {
      points(x = section$data$x, y = section$data$y, pch = 20, cex = include.cex) ## Add black points
    }
  }

  ## Include Data
  # plot data points using the same color palette as the gridded data fields.
  if (include.data & !mark.points) {
    points(x = section$data$x, y = section$data$y, pch = include.pch, cex = include.cex,
           col = make.pal(x = section$data$z, pal = pal, n = n, min = zlim[1], max = zlim[2], rev = rev),
           bg = make.pal(x = section$data$z, pal = pal, n = n, min = zlim[1], max = zlim[2], rev = rev))
  }
  box() ## make sure plotting didn't cover bounding box
}

#' @title Get Section Bathymetry
#' @author Thomas Bryce Kelly
#' @description retreive bathymetry for a section plot
#' @keywords Bathymetry
#' @export
#' @param lon longitudes
#' @param lat latitudes
#' @param res the resolution of the bathymetry (arc minutes, e.g. 10)
get.section.bathy = function(lon, lat, res = 10) {

  ## Construct dummy map object
  map =  list(lon.min = max(min(lon) - 1, -180), lon.max = min(max(lon) + 1, 180),
              lat.min = max(min(lat) - 1, -90),  lat.max = min(max(lat) + 1, 90))
  bathy = get.bathy(map, res = res)

  oce::bilinearInterp(lon, lat, bathy$Lon, bathy$Lat, bathy$Z)
}


#' @title Add Section Bathymetry
#' @author Thomas Bryce Kelly
#' @param section a section object from build.section()
#' @param bathy a bathymetric object such as "bathy.global" (the default)
#' @param binning a smoothing paramter; must be an odd integer
#' @param bathy.col the color of the bathymetry
#' @export
add.section.bathy = function(section, bathy = bathy.global, binning = 1, bathy.col = 'darkgrey') {
  x = section$x
  depth = rep(NA, length(x))

  bathy$lon = bathy$lon %% 360
  section$section.lon = section$section.lon %% 360

  for (i in 1:length(x)) {
    depth[i] = -bathy$Z[which.min((bathy$Lon - section$section.lon[i])^2), which.min((bathy$Lat - section$section.lat[i])^2)]
  }
  depth = runmed(depth, binning)

  polygon(x = c(x, rev(x)), y = c(depth, rep(1e8, length(x))), col = bathy.col)
}


#' @title Add Contour to Section Plot
#' @author Thomas Bryce Kelly
#' @param section a section object from "build.section()"
#' @param field the field name to use for contouring
#' @param levels the level(s) to draw contours at
#' @param col the color of the contour line(s)
#' @param lty the lty for the line(s)
#' @param labels a boolean for whether the label(s) should be drawn
#' @param cex.lab the cex (i.e. size) of the labels
#' @export
add.section.contour = function(section, field = NULL, levels = NULL, col = 'black', lty = 1, lwd = 1, labels = NULL, cex.lab = 1) {
  if (is.null(field)) {
    field = colnames(section$grid)[3] # first interpolated field
  }
  z = matrix(section$grid[[field]], nrow = length(section$x))
  if (is.null(levels)) { levels = pretty(range(as.numeric(z)), n = 5) }

  contour(section$x, section$y, z, add = TRUE, levels = levels, col = col, lty = lty, lwd = lwd,
          labels = labels, labcex = cex.lab, drawlabels = !isFALSE(labels))
}


#' @title Add Map Inlay
#' @export
add.section.inlay = function(section) {
  par.old = par()
  par(new = T, plt = c(0.2, 0.33, 0.2, 0.33))
  map = make.map(coast = 'coastlineWorld', lon.min = min(section$section.lon), lon.max = max(section$section.lon),
                 lat.min = min(section$section.lat), lat.max = max(section$section.lat),
                 dlon = 360, dlat = 360, grid = F,
                 p = make.proj(8, lat = mean(section$section.lat), lon = mean(section$section.lon)))
  add.map.points(section$section.lon, section$section.lat, pch = 20, cex = 0.2, col = 'red')
  par(plt = par.old$plt)
}



#' @useDynLib TheSource
NULL
