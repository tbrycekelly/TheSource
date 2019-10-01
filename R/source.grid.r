## Functions for building gridded datasets from measurements.
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


#######################
### Section Building ##
#######################

#' @title Build Section
#' @author Thomas Bryce Kelly
#' @description
#' @keywords
#' @param
#' @export
## Build a section object
# x, y: dimensions (e.g. lat, lon, depth, section distance, time, etc)
# z: signal to be gridded (e.g. T, S, ...)
# xlim, ylim = NULL: Limits of the gridding. These are the bounds of the new x-y grid. Default: NULL will set it based on the data + 10%.
# x.factor, y.factor = 1: The relative scale difference between x and y, used to calculate distances. Take into account actual scale AND the relevent scaling of the system (vertical distance tends to be more important than horizontal distance).
# x.scale, y.scale = NULL: The step size in the new x-y grid. By default the scale is set to generate a grid that is 50x50.
# uncertainty = 0: Scaling applied to the distance from the cener of a grid cell to a vertex, used to add a base-line distance to all measurements. 0 = no minimum, 1 = minimum = to half a grid cell.
# neighborhood = -1: The number of closest points used in the interpolation of a grid cell. Special value of -1 specifies the use of all points. If neighborhood > number of valid data then neighborhood is set to equal the number of valid datapoints.
# field.name = 'z': Sets the name of the new interpolated field. By default the name is 'z', but setting it to 'T' for temperature makes sense.
build.section = function(x, y, z, xlim = NULL, ylim = NULL, x.factor = 1, y.factor = 1, x.scale = NULL, y.scale = NULL,
                         uncertainty = 1e-12, field.name = 'z', neighborhood = -1, p = 3) {

  ## Remove NAs
  l = !is.na(x) & !is.na(y) & !is.na(z)
  x = x[l]
  y = y[l]
  z = z[l]

  if (uncertainty == 0) { warning('Uncertainty of zero may produce NAs!') }

  ## Set default limits (+10% buffer)
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

  if (is.null(x.scale)) { x.scale = (xlim[2] - xlim[1]) / 50} ## Default to 50 steps
  if (is.null(y.scale)) { y.scale = (ylim[2] - ylim[1]) / 50} ## Default to 50 steps


  ## field = matrix (1 - 300 m depth x 0 - max dist)
  n.y = (ylim[2] - ylim[1]) / y.scale
  n.x = (xlim[2] - xlim[1]) / x.scale
  y.new = seq(ylim[1], ylim[2], by = y.scale)
  x.new = seq(xlim[1], xlim[2], by = x.scale)

  ## Minimum distance based on grid size
  delta.min = delta(x.factor, y.factor, x.scale/2, 0, y.scale/2, 0, p = p) * uncertainty

  grid = expand.grid(x = x.new, y = y.new)
  grid[[field.name]] = apply(grid, 1, function(row) {
    sum(z * W2(delta(x.factor, y.factor, x, row[1], y, row[2], p = p), delta.min), na.rm = TRUE)
    })


  ## Construct return object
  grid = list(grid = grid,
              grid.meta = list(
                x.scale = x.scale,
                y.scale = y.scale,
                x.factor = x.factor,
                y.factor = y.factor,
                n.x = n.x, n.y = n.y,
                uncertainty = uncertainty,
                neighborhood = neighborhood,
                p = p
              ),
              x = x.new,
              y = y.new,
              data = data.frame(x = x, y = y, z = z)
  )
  ## Return
  grid
}

#' @title Build Section (old)
#' @author Thomas Bryce Kelly
#' @description
#' @keywords
#' @param
#' @export
build.section.old = function(x, y, z, xlim = NULL, ylim = NULL, x.factor = 1, y.factor = 1, x.scale = NULL, y.scale = NULL,
                         uncertainty = 1e-12, neighborhood = -1, field.name = 'z', p = 3) {

    ## Remove NAs
    l = !is.na(x) & !is.na(y) & !is.na(z)
    x = x[l]
    y = y[l]
    z = z[l]

    if (uncertainty == 0) { warning('Uncertainty of zero may produce NAs!') }

    ## Set default limits (+10% buffer)
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

    if (is.null(x.scale)) { x.scale = (xlim[2] - xlim[1]) / 50} ## Default to 50 steps
    if (is.null(y.scale)) { y.scale = (ylim[2] - ylim[1]) / 50} ## Default to 50 steps

    if (neighborhood > length(x)) {neighborhood = length(x)}

    if (neighborhood == -1) { neighborhood = length(l) }

    ## field = matrix (1 - 300 m depth x 0 - max dist)
    n.y = (ylim[2] - ylim[1]) / y.scale
    n.x = (xlim[2] - xlim[1]) / x.scale
    y.new = seq(ylim[1], ylim[2], by = y.scale)
    x.new = seq(xlim[1], xlim[2], by = x.scale)

    ## Minimum distance based on grid size
    delta.min = delta(x.factor, y.factor, x.scale/2, 0, y.scale/2, 0) * uncertainty

    grid = expand.grid(x = x.new, y = y.new)
    grid[[field.name]] = NA

    for (i in 1:nrow(grid)) {

        ## Calculate distances and determine weights
        del = delta(x.factor, y.factor, x, grid$x[i], y, grid$y[i])

        w = W(del, delta.min) # Get weights
        ll = order(w, decreasing = TRUE)[1:neighborhood]  # isolate to Neighborhood of n-closest values

        w = w / sum(w[ll]) # Normalize
        grid[[field.name]][i] = sum(z[ll] * w[ll])
    }

    ## Construct return object
    grid = list(grid = grid,
                grid.meta = list(
                    x.scale = x.scale,
                    y.scale = y.scale,
                    x.factor = x.factor,
                    y.factor = y.factor,
                    n.x = n.x, n.y = n.y,
                    uncertainty = uncertainty,
                    neighborhood = neighborhood,
                    p = p
                ),
                x = x.new,
                y = y.new,
                data = data.frame(x = x, y = y, z = z)
               )
    ## Return
    grid
}


#' @title Build Section (3D)
#' @author Thomas Bryce Kelly
#' @description
#' @keywords
#' @param
#' @export
## Build a 3D section object
# x, y, z: dimensions (e.g. lat, lon, depth, section distance, time, etc)
# zz: signal to be gridded (e.g. T, S, ...)
# xlim, ylim, zlim = NULL: Limits of the gridding. These are the bounds of the new x-y grid. Default: NULL will set it based on the data + 10%.
# x.factor, y.factor, z.factor = 1: The relative scale difference between x and y, used to calculate distances. Take into account actual scale AND the relevent scaling of the system (vertical distance tends to be more important than horizontal distance).
# x.scale, y.scale, z.scale = NULL: The step size in the new x-y grid. By default the scale is set to generate a grid that is 50x50.
# uncertainty = 0: Scaling applied to the distance from the cener of a grid cell to a vertex, used to add a base-line distance to all measurements. 0 = no minimum, 1 = minimum = to half a grid cell.
# neighborhood = 20: The number of closest points used in the interpolation of a grid cell. Special value of -1 specifies the use of all points. If neighborhood > number of valid data then neighborhood is set to equal the number of valid datapoints.
# field.name = 'zz': Sets the name of the new interpolated field. By default the name is 'zz', but setting it to 'T' for temperature makes sense.
build.section.3d = function(x, y, z, zz, xlim = NULL, ylim = NULL, zlim = NULL, x.factor = 1, y.factor = 1, z.factor = 1,
                            x.scale = NULL, y.scale = NULL, z.scale = NULL, uncertainty = 1e-6,
                            neighborhood = 20, field.name = 'zz', p = 3) {

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

    if (neighborhood > length(x)) {neighborhood = length(x)}

    if (neighborhood == -1) { neighborhood = length(l) }

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
    grid[[field.name]] = NA

    for (i in 1:nrow(grid)) {

        ## Calculate distances and determine weights
        del = delta3(x.factor, y.factor, z.factor, x, grid$x[i], y, grid$y[i], z, grid$z[i])

        w = W(del, delta.min) # Get weights
        ll = order(w, decreasing = TRUE)[1:neighborhood]  # isolate to Neighborhood of n-closest values

        w = w / sum(w[ll]) # Normalize
        grid[[field.name]][i] = sum(zz[ll] * w[ll])
    }

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
                    neighborhood = neighborhood,
                    p = p
                ),
                x = x.new,
                y = y.new,
                z = z.new,
                data = data.frame(x = x, y = y, z = z, zz = zz)
               )
    ## Return
    grid
}

#' @title Build Section (multi)
#' @author Thomas Bryce Kelly
#' @description
#' @keywords
#' @param
#' @export
## Build a multi=parameter section object
# This function is much more efficient than building a grid and then adding multiple paramters to it.
# x, y [numeric]: dimensions (e.g. lat, lon, depth, section distance, time, etc)
# z[dataframe]: signal to be gridded (e.g. T, S, ...)
# xlim, ylim = NULL: Limits of the gridding. These are the bounds of the new x-y grid. Default: NULL will set it based on the data + 10%.
# x.factor, y.factor = 1: The relative scale difference between x and y, used to calculate distances. Take into account actual scale AND the relevent scaling of the system (vertical distance tends to be more important than horizontal distance).
# x.scale, y.scale = NULL: The step size in the new x-y grid. By default the scale is set to generate a grid that is 50x50.
# uncertainty = 0: Scaling applied to the distance from the cener of a grid cell to a vertex, used to add a base-line distance to all measurements. 0 = no minimum, 1 = minimum = to half a grid cell.
# neighborhood = 20: The number of closest points used in the interpolation of a grid cell. Special value of -1 specifies the use of all points. If neighborhood > number of valid data then neighborhood is set to equal the number of valid datapoints.
# field.name = 'z': Sets the name of the new interpolated field. By default the name is 'z', but setting it to 'T' for temperature makes sense.
build.section.multi = function(x, y, z, xlim = NULL, ylim = NULL, x.factor = 1, y.factor = 1, x.scale = NULL, y.scale = NULL,
                               uncertainty = 1e-8, neighborhood = 20, field.names = NULL) {

  ## Remove NAs
  l = !is.na(x) & !is.na(y) & apply(z, 1, function(x) {!any(is.na(x))})
  x = x[l]
  y = y[l]
  z = z[l,]

  if (uncertainty == 0) { warning('Uncertainty of zero may produce NAs!') }

  if (is.null(field.names)) { field.names = colnames(z) }

  if (is.null(xlim)) { xlim = range(x); xlim[1] = xlim[1] - (xlim[2]-xlim[1])/20; xlim[2] = xlim[2] + (xlim[2]-xlim[1])/20;}
  if (is.null(ylim)) { ylim = range(y) }

  if (is.null(x.scale)) { x.scale = (xlim[2] - xlim[1]) / 50} ## Default to 50 steps
  if (is.null(y.scale)) { y.scale = (ylim[2] - ylim[1]) / 50} ## Default to 50 steps

  neighborhood = min(neighborhood, length(l))

  if (neighborhood == -1) { neighborhood = length(l) }

  ## field = matrix (1 - 300 m depth x 0 - max dist)
  n.y = (ylim[2] - ylim[1]) / y.scale
  n.x = (xlim[2] - xlim[1]) / x.scale
  y.new = seq(ylim[1], ylim[2], by = y.scale)
  x.new = seq(xlim[1], xlim[2], by = x.scale)

  ## Minimum distance based on grid size
  delta.min = delta(x.factor, y.factor, x.scale/2, 0, y.scale/2, 0) * uncertainty

  grid = expand.grid(x = x.new, y = y.new)

  for (kk in 1:length(field.names)) {
    grid[[field.names[kk]]] = apply(grid, 1, function(row) {sum(z[,kk] * W2(delta(x.factor, y.factor, x, row[1], y, row[2]), delta.min), na.rm = TRUE)})
  }

  ## Construct return object
  grid = list(grid = grid,
              grid.meta = list(
                x.scale = x.scale,
                y.scale = y.scale,
                x.factor = x.factor,
                y.factor = y.factor,
                n.x = n.x, n.y = n.y,
                uncertainty = uncertainty,
                neighborhood = neighborhood
              ),
              x = x.new,
              y = y.new,
              data = data.frame(x = x, y = y, z = z)
  )
  ## Return
  grid
}


#' @title Build Section (multi, old)
#' @author Thomas Bryce Kelly
#' @description
#' @keywords
#' @export
#' @param
build.section.multi.old = function(x, y, z, xlim = NULL, ylim = NULL, x.factor = 1, y.factor = 1, x.scale = NULL, y.scale = NULL,
                         uncertainty = 0, neighborhood = 20, field.names = NULL) {

    ## Remove NAs
    l = !is.na(x) & !is.na(y) & apply(z, 1, function(x) {!any(is.na(x))})
    x = x[l]
    y = y[l]
    z = z[l,]

    if (uncertainty == 0) { warning('Uncertainty of zero may produce NAs!') }

    if (is.null(field.names)) { field.names = colnames(z) }

    if (is.null(xlim)) { xlim = range(x); xlim[1] = xlim[1] - (xlim[2]-xlim[1])/20; xlim[2] = xlim[2] + (xlim[2]-xlim[1])/20;}
    if (is.null(ylim)) { ylim = range(y) }

    if (is.null(x.scale)) { x.scale = (xlim[2] - xlim[1]) / 50} ## Default to 50 steps
    if (is.null(y.scale)) { y.scale = (ylim[2] - ylim[1]) / 50} ## Default to 50 steps

    neighborhood = min(neighborhood, length(l))

    if (neighborhood == -1) { neighborhood = length(l) }

    ## field = matrix (1 - 300 m depth x 0 - max dist)
    n.y = (ylim[2] - ylim[1]) / y.scale
    n.x = (xlim[2] - xlim[1]) / x.scale
    y.new = seq(ylim[1], ylim[2], by = y.scale)
    x.new = seq(xlim[1], xlim[2], by = x.scale)

    ## Minimum distance based on grid size
    delta.min = delta(x.factor, y.factor, x.scale/2, 0, y.scale/2, 0) * uncertainty

    grid = expand.grid(x = x.new, y = y.new)

    for (name in field.names) {
        grid[[name]] = NA
    }

    for (i in 1:nrow(grid)) {

        ## Calculate distances and determine weights
        del = delta(x.factor, y.factor, x, grid$x[i], y, grid$y[i])

        w = W(del, delta.min) # Get weights
        ll = order(w, decreasing = TRUE)[1:neighborhood]  # isolate to Neighborhood of n-closest values
        w = w / sum(w[ll]) # Normalize

        for (kk in 1:length(field.names)) {
            grid[[field.names[kk]]][i] = sum(as.numeric(z[ll,kk]) * w[ll])
        }
    }

    ## Construct return object
    grid = list(grid = grid,
                grid.meta = list(
                    x.scale = x.scale,
                    y.scale = y.scale,
                    x.factor = x.factor,
                    y.factor = y.factor,
                    n.x = n.x, n.y = n.y,
                    uncertainty = uncertainty,
                    neighborhood = neighborhood
                ),
                x = x.new,
                y = y.new,
                data = data.frame(x = x, y = y, z = z)
               )
    ## Return
    grid
}

#' @title Add Section Parameter
#' @author Thomas Bryce Kelly
#' @description
#' @keywords
#' @param
#' @export
## Add additional paramter to section object
# x, y: dimensions (e.g. lat, lon, depth, section distance, time, etc)
# z: signal to be gridded (e.g. T, S, ...)
# section: The existing section object (2D only) to which the new parameter will be added.
# field.name = 'z': Sets the name of the new interpolated field. By default the name is 'z', but setting it to 'T' for temperature makes sense.
add.section.param = function(x, y, z, section, field.name) {

    ## Remove NAs
    l = !is.na(x) & !is.na(y) & !is.na(z)
    x = x[l]
    y = y[l]
    z = z[l]

    ## Get grid parameters
    xlim = section$grid.meta$xlim
    ylim = section$grid.meta$ylim
    x.scale = section$grid.meta$x.scale
    y.scale = section$grid.meta$y.scale
    x.factor = section$grid.meta$x.factor
    y.factor = section$grid.meta$y.factor
    neighborhood = min(section$grid.meta$neighborhood, length(l))
    uncertainty = section$grid.meta$uncertainty

    if (neighborhood == -1) { neighborhood = length(l) }

    ## Minimum distance based on grid size
    delta.min = delta(x.factor, y.factor, x.scale/2, 0, y.scale/2, 0) * uncertainty

    section$grid[[field.name]] = NA

    for (i in 1:nrow(section$grid)) {

        ## Calculate distances and determine weights
        del = delta(x.factor, y.factor, x, section$grid$x[i], y, section$grid$y[i])

        w = W(del, delta.min) # Get weights
        ll = order(w, decreasing = TRUE)[1:neighborhood]  # isolate to Neighborhood of n-closest values

        w = w / sum(w[ll]) # Normalize
        section$grid[[field.name]][i] = sum(z[ll] * w[ll])
    }
    section$data = rbind(section$data, data.frame(x = x, y = y, z = z))

    ## Return
    section
}


jackknife.section = function(section, n) {

  ## Get grid parameters
  xlim = section$grid.meta$xlim
  ylim = section$grid.meta$ylim
  x.scale = section$grid.meta$x.scale
  y.scale = section$grid.meta$y.scale
  x.factor = section$grid.meta$x.factor
  y.factor = section$grid.meta$y.factor
  neighborhood = section$grid.meta$neighborhood
  uncertainty = section$grid.meta$uncertainty

  ## Minimum distance based on grid size
  delta.min = delta(x.factor, y.factor, x.scale/2, 0, y.scale/2, 0, p = section$grid.meta$p) * uncertainty

  section$jack = expand.grid(x = section$x, y = section$y)

  for (k in 1:n) {
    l = sample(1:nrow(section$data), nrow(section$data), replace = TRUE)
    l = unique(l)
    section$jack[[paste0('z', k)]] = NA

    for (i in 1:nrow(section$grid)) {
      ## Calculate distances and determine weights
      del = delta(x.factor, y.factor, section$data$x[l], section$grid$x[i], section$data$y[l], section$grid$y[i])

      w = W(del, delta.min) # Get weights
      w = w / sum(w) # Normalize
      section$jack[[paste0('z', k)]][i] = sum(section$data$z[l] * w)
    }
  }

  ## Return
  section
}

#' @title Build Section (ODV)
#' @author Thomas Bryce Kelly
#' @description
#' @export
#' @keywords
#' @param
## Build a section objectbased on ODV weighting
#' @param x Dimensions (e.g. lat, lon, depth, section distance, time, etc)
#' @param y Dimensions (e.g. lat, lon, depth, section distance, time, etc)
#' @param z signal to be gridded (e.g. T, S, ...)
#' @param xlim Limits of the gridding. These are the bounds of the new x-y grid. Default: NULL will set it based on the data + 10%.
#' @param x.factor, y.factor The relative scale difference between x and y, used to calculate distances. Take into account actual scale AND the relevent scaling of the system (vertical distance tends to be more important than horizontal distance).
#' @param x.scale, y.scale The step size in the new x-y grid. By default the scale is set to generate a grid that is 50x50.
#' @param uncertainty Scaling applied to the distance from the cener of a grid cell to a vertex, used to add a base-line distance to all measurements. 0 = no minimum, 1 = minimum = to half a grid cell.
#' @param neighborhood The number of closest points used in the interpolation of a grid cell. Special value of -1 specifies the use of all points. If neighborhood > number of valid data then neighborhood is set to equal the number of valid datapoints.
#' @param field.name Sets the name of the new interpolated field. By default the name is 'z', but setting it to 'T' for temperature makes sense.
build.section.odv = function(x, y, z, xlim = NULL, ylim = NULL, x.factor = 1, y.factor = 1, x.scale = NULL, y.scale = NULL,
                         uncertainty = 0, neighborhood = 20, field.name = 'z') {

    ## Remove NAs
    l = !is.na(x) & !is.na(y) & !is.na(z)
    x = x[l]
    y = y[l]
    z = z[l]

    ## Set default limits (+10% buffer)
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

    if (is.null(x.scale)) { x.scale = (xlim[2] - xlim[1]) / 50} ## Default to 50 steps
    if (is.null(y.scale)) { y.scale = (ylim[2] - ylim[1]) / 50} ## Default to 50 steps

    if (neighborhood > length(x)) {neighborhood = length(x)}

    if (neighborhood == -1) { neighborhood = length(l) }

    ## field = matrix (1 - 300 m depth x 0 - max dist)
    n.y = (ylim[2] - ylim[1]) / y.scale
    n.x = (xlim[2] - xlim[1]) / x.scale
    y.new = seq(ylim[1], ylim[2], by = y.scale)
    x.new = seq(xlim[1], xlim[2], by = x.scale)

    ## Minimum distance based on grid size
    delta.min = delta(x.factor, y.factor, x.scale/2, 0, y.scale/2, 0, p = 2) * uncertainty

    grid = expand.grid(x = x.new, y = y.new)
    grid[[field.name]] = NA

    for (i in 1:nrow(grid)) {

        ## Calculate distances and determine weights
        del = delta(x.factor, y.factor, x, grid$x[i], y, grid$y[i], p = 2)

        w = W.odv(del, delta.min) # Get weights
        ll = order(w, decreasing = TRUE)[1:neighborhood]  # isolate to Neighborhood of n-closest values

        w = w / sum(w[ll]) # Normalize
        grid[[field.name]][i] = sum(z[ll] * w[ll])
    }

    ## Construct return object
    grid = list(grid = grid,
                grid.meta = list(
                    x.scale = x.scale,
                    y.scale = y.scale,
                    x.factor = x.factor,
                    y.factor = y.factor,
                    n.x = n.x, n.y = n.y
                ),
                x = x.new,
                y = y.new,
                data = data.frame(x = x, y = y, z = z)
               )
    ## Return
    grid
}


#' @title Cross Validate the IDW
#' @author Thomas Bryce Kelly
#' @export
cv.section = function(section) {

  var = 0
  delta.min = delta(section$grid.meta$x.factor, section$grid.meta$y.factor,
                    section$grid.meta$x.scale/2, 0, section$grid.meta$y.scale/2, 0,
                    p = section$grid.meta$p) * section$grid.meta$uncertainty

  for (i in 1:nrow(section$data)) {
    del = delta(section$grid.meta$x.factor, section$grid.meta$y.factor, section$data$x[-i], section$data$x[i],
                section$data$y[-i], section$data$y[i], p = section$grid.meta$p)

    w = W(del, delta.min)
    w = w / sum(w)
    var = var + (sum(w * section$data$z[-i]) - section$data$z[i])^2
  }
  section$rmse = sqrt(var / nrow(section$data))
  section
}

#' @title Distance Metric
#' @author Thomas Bryce Kelly
#' @description
#' @keywords
#' @param
## Calculates p-distance
delta = function(x.factor, y.factor, x1, x2, y1, y2, p = 3) {
    (x.factor * abs(x2 - x1))^p + (y.factor * abs(y2 - y1))^p
}

#' @title Distance Metric (3D)
#' @author Thomas Bryce Kelly
#' @description
#' @keywords
#' @param
## Calculates p-distance for 3 dimensions
delta3 = function(x.factor, y.factor, z.factor, x1, x2, y1, y2, z1, z2, p = 3) {
    (x.factor * abs(x2 - x1))^p + (y.factor * abs(y2 - y1))^p + (z.factor * abs(z2 - z1))^p
}

#' @title Weighting Function
#' @author Thomas Bryce Kelly
#' @description
#' @keywords
#' @param
## Calculates weights
W = function(del, delta.min) {
    1 / (del + delta.min)
}

#' @title Weighting Function (ODV)
#' @author Thomas Bryce Kelly
#' @description
#' @keywords
#' @param
W.odv = function(del, delta.min) {
    exp(-del - delta.min)
}

#' @title Weighting Function (2)
#' @author Thomas Bryce Kelly
#' @description
#' @keywords
#' @param
W2 = function(del, delta.min) { ## Normalized
  a = 1 / (del + delta.min)
  a / sum(a)
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
plot.section = function(section, field = 'z', xlim=NULL, ylim=NULL, xlab = 'x', ylab = 'y', log = FALSE, base = 10,
                        zlim = NULL, pal='ocean.matter', rev = FALSE, include.data = FALSE, mark.points = TRUE, include.pch = 21,
                       include.cex = 1, main = NULL, col.low = '', col.high = '') {

    ## Handy variables
    x = section$x
    y = section$y
    z = matrix(section$grid[,field], nrow = length(x))

    ## Set Defaults
    # Check if values need to set and set them.
    if (is.null(zlim)) { zlim = range(as.numeric(section$grid[,field]), na.rm = TRUE) }

    if (is.null(main)) { main = field }

    if (is.null(xlim)) { xlim = range(section$x)}

    if (is.null(ylim)) { ylim = range(section$y) }

    if (log) {
        z = log(z, base)
        zlim = log(zlim, base)
    }

    ## Out of Ragne
    # Set values to zlim if you want the out of range values plotted at the zlim values
    if (!is.na(col.low) & col.low == '') {
        z[z < zlim[1]] = zlim[1]
    }
    if (!is.na(col.high) & col.high == '') {
        z[z > zlim[2]] = zlim[2]
    }

    ## Plot iamge
    image(x = x, y = y, z = z, col = get.pal(255, pal = pal, rev = rev), ylab = ylab, xlab = xlab,
          xlim = xlim, ylim = ylim, zlim = zlim)


    ## Add Title text
    st = paste0(main, '   zlim: (', round(zlim[1], 3), ', ', round(zlim[2],3), ')')
    mtext(st, line = 0.25, adj = 1, cex = 0.7)

    ## Out of Range !! TODO Doesn't work yet, plot type is incompatible with points()
    # Plot points that are out of range when a color is given
    #if (!is.na(col.low) & col.low != '') {
    #    l = z < zlim[1]
    #    points(x = x[l], y = y[l], z[l], col = col.low, pch = include.pch, cex = include.cex)
    #}

    #if (!is.na(col.high) & col.high != '') {
    #    l = z > zlim[2]
    #    points(x = x[l], y = y[l], z[l], col = col.low, pch = include.pch, cex = include.cex)
    #}


    if (mark.points) {
        if (include.data) {
            points(x = section$data$x, y = section$data$y, pch = include.pch, cex = include.cex, col = 'black',
                   bg = make.pal(x = section$data$z, pal = pal, n = 255, min = zlim[1], max = zlim[2], rev = rev))
        } else {
            points(x = section$data$x, y = section$data$y, pch = 20, cex = include.cex) ## Add black points
        }
    }


    ## Include Data
    # plot data points using the same color palette as the gridded data fields.
    if (include.data & !mark.points) {
        points(x = section$data$x, y = section$data$y, pch = include.pch, cex = include.cex,
               col = make.pal(x = section$data$z, pal = pal, n = 255, min = zlim[1], max = zlim[2], rev = rev),
              bg = make.pal(x = section$data$z, pal = pal, n = 255, min = zlim[1], max = zlim[2], rev = rev))
    }
    box() ## make sure plotting didn't cover bounding box
}

#' @title Get Section Bathymetry
#' @author Thomas Bryce Kelly
#' @description
#' @keywords
#' @export
#' @param
## !TODO: integrate lat lon into sections more naturally. Add a dedicated mapping function too?
get.section.bathy = function(lon, lat, res = 10) {

  ## Construct dummy map object
  map =  list(lon.min = max(min(lon) - 1, -180), lon.max = min(max(lon) + 1, 180),
              lat.min = max(min(lat) - 1, -90),  lat.max = min(max(lat) + 1, 90))
  bathy = get.bathy(map, res = res)

  oce::bilinearInterp(lon, lat, bathy$Lon, bathy$Lat, bathy$Z)
}


#' @title Add Section Bathymetry
#' @author Thomas Bryce Kelly
#' @description
#' @keywords
#' @export
#' @param
add.section.bathy = function(section, stn.lon, stn.lat, stn.depth = NULL, by.lat = TRUE, bathy.col = 'darkgrey', res = 2) {
  x = section$x
  y = section$y

  if (is.null(stn.depth)) {
      depth = get.section.bathy(stn.lon, stn.lat, res = res)
  }

  if (length(stn.depth) != length(stn.lon)) {
      stop('Number of station depths does not match number of stations!')
  }

  if (any(stn.depth > 0)) {
      warning('Reminder, stn.depths are negative below sealevel.')
  }

  if (by.lat) {
      bathyy = approx(stn.lat, -stn.depth, x, rule = 2)$y
  } else {
      bathyy = approx(stn.lon, -stn.depth, x, rule = 2)$y
  }

  polygon(x = c(x, rev(x)), y = c(bathyy, rep(1e8, length(x))), col = bathy.col)
}


#' @title Add Contour
#' @author Thomas Bryce Kelly
#' @description
#' @export
#' @keywords
#' @param
## !TODO: Improve this hack
add.contour = function(section, field = 'Sigma', text.x = NULL, levels = c(25.7, 26, 26.3, 26.6),
                       col = 'black', text.col = 'black', lty = 1, lwd = 1, f = 2/3, text.default.pos = 0.1) {

    if (is.null(text.x)) {
        text.x = as.numeric(quantile(section$x, text.default.pos))
    }

    for (ii in levels) {
      x = section$x
      y = rep(NA, length(x)) ## Default to NAs

        for (i in 1:length(x)) {

            ## Find the entries that are in that position and that are not NAs
            l = which(section$grid$x == x[i] & !is.na(section$grid[,field]))

            if (length(l) > 1) {
              if (all(section$grid[l,field] > ii)) { ## Case 1: the line is below the grid values
                  k = which.min(section$grid[l,field] - ii)
                  y.range = range(section$y)

                  if (k > length(l)/2) {
                    ## Place just outside of gridded range
                    y[i] = max(section$y) + 0.1 * (y.range[2] - y.range[1])
                  } else {
                    ## Place just outside of gridded range
                    y[i] = min(section$y) - 0.1 * (y.range[2] - y.range[1])
                  }
              } else if(all(section$grid[l,field] < ii)) { ## Case 2: line is above the grid values
                  k = which.min(ii - section$grid[l,field])
                  y.range = range(section$y)

                  if (k > length(l)/2) {
                    ## Place just outside of gridded range
                    y[i] = max(section$y) + 0.1 * (y.range[2] - y.range[1])
                  } else {
                    ## Place just outside of gridded range
                    y[i] = min(section$y) - 0.1 * (y.range[2] - y.range[1])
                  }
              } else {
                ## The usual interpolation
                y[i] = approx(section$grid[l, field], section$grid$y[l], ii, rule = 2)$y
              }
            } else {

              warning('Not enough valid values. NAs exist.')
            }
        }

        ## Remove NA, this shouldn't happen.... unless warning issued above
        x = x[!is.na(y)]
        y = y[!is.na(y)]

        ## Draw Lines
        low = lowess(x, y, f = f)
        lines(low, col = col, lwd = lwd, lty = lty)

        if (is.POSIXct(x[1])) {
            text(text.x, low$y[which.min(abs(as.numeric(x) - text.x))], round(ii, 2), cex = 0.75, col = text.col)
        } else {
            text(text.x, low$y[which.min(abs(x-text.x))], round(ii, 2), cex = 0.75, col = text.col)
        }
    }
}
