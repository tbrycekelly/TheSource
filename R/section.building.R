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
build.section = function(x, y, z, lat = NULL, lon = NULL, griddder = NULL, grid = NULL, weight = NULL,
                         xlim = NULL, ylim = NULL,
                         x.factor = NULL, y.factor = NULL,
                         x.scale = NULL, y.scale = NULL,
                         uncertainty = 1, p = 2,
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
    gridder = gridODV
    if (verbose) { message('BUILD.SECTION: No gridder specified, defaulting to gridODV. Other options: gridBIN, gridIDW, gridNN, gridNNI and gridKrige.') }
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
  if (is.null(x.factor)) { x.factor = (y.scale/x.scale + 1) / 2}
  if (is.null(y.factor)) { y.factor = (x.scale/y.scale + 1) / 2}

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
    grid[[field.names[kk]]] = gridder(grid$x, grid$y, x, y, z[,kk], p, x.scale, y.scale, uncertainty)
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
                gridder = deparse(substitute(gridder)),
                time = list(time.built = Sys.time(), build.time = time.c - time.a, grid.time = time.c - time.b),
                Source.version = packageVersion('TheSource'),
                R.version = R.version.string
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
