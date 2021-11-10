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
