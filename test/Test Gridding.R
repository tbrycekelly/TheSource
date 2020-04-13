library(TheSource)
library(microbenchmark)


build.section.new = function(x, y, z, xlim = NULL, ylim = NULL, x.factor = 1, y.factor = 1, x.scale = NULL, y.scale = NULL,
                         uncertainty = 1e-12, field.names = 'z', p = 3) {

  ## Remove NAs
  l = !is.na(x) & !is.na(y) & !is.na(z)
  x = x[l]
  y = y[l]
  z = z[l]
  z = matrix(z)

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

  ## Rescale x and y based on x.factor and y.factor
  x = x * x.factor
  x.scale = x.scale * x.factor
  xlim = xlim * x.factor
  y = y * y.factor
  y.scale = y.scale * y.factor
  ylim = ylim * y.factor

  ## field = matrix (1 - 300 m depth x 0 - max dist)
  n.y = (ylim[2] - ylim[1]) / y.scale
  n.x = (xlim[2] - xlim[1]) / x.scale
  y.new = seq(ylim[1], ylim[2], by = y.scale)
  x.new = seq(xlim[1], xlim[2], by = x.scale)

  ## Make grid and fill in
  grid = expand.grid(x = x.new, y = y.new)
  for (kk in 1:length(field.names)) {
    grid[[field.names[kk]]] = grid.IDW(grid, x, y, z[,kk], p, x.scale, y.scale, uncertainty)
  }
  grid$x = grid$x / x.factor
  grid$y = grid$y / y.factor

  ## Construct return object
  grid = list(grid = grid,
              grid.meta = list(
                x.scale = x.scale / x.factor,
                y.scale = y.scale / y.factor,
                x.factor = x.factor,
                y.factor = y.factor,
                n.x = n.x, n.y = n.y,
                uncertainty = uncertainty,
                p = p,
                gridder = grid.IDW
              ),
              x = x.new / x.factor,
              y = y.new / y.factor,
              data = data.frame(x = x / x.factor, y = y / y.factor, z = z)
  )
  ## Return
  grid
}

#' @title Grid by IDW
#' @export
#' @author Thomas Bryce Kelly
#'
grid.IDW = function(grid, x, y, z, p = 3, x.scale, y.scale, uncertainty) {
  ## Minimum distance based on grid size
  delta.min = (abs(x.scale/2)^p + abs(y.scale/2)^p) * uncertainty

  temp = rep(0, nrow(grid))
  for (i in 1:nrow(grid)) {
    ## Calculate distances and determine weights
    w = 1 / (abs(x - grid$x[i])^p + abs(y - grid$y[i])^p + delta.min)
    temp[i] = sum(z * w) / sum(w)
  }

  temp
}


grid.ODV = function(grid, x, y, z, p = 2, x.scale, y.scale, uncertainty) {
  ## Minimum distance based on grid size
  #delta.min = (abs(x.scale/2)^p + abs(y.scale/2)^p) * uncertainty

  temp = rep(0, nrow(grid))
  for (i in 1:nrow(grid)) {
    ## Calculate distances and determine weights
    w = exp((abs(x - grid$x[i])^p + abs(y - grid$y[i])^p) * -1e-6)
    temp[i] = sum(z * w) / sum(w)
  }
  temp
}


add.section.param = function(x, y, z, section, field.name) {
  ## Remove NAs
  l = !is.na(x) & !is.na(y) & !is.na(z)
  x = x[l]; y = y[l]; z = z[l]

  section$grid[[field.name]] = grid.IDW(grid, x, y, z, p,
                                        section$grid.meta$x.scale,
                                        section$grid.meta$y.scale,
                                        section$grid.meta$uncertainty)
  section
}





microbenchmark(build.section(c(1:100), c(1:100), c(1:100), x.scale = 0.2, y.scale = 1),
               build.section.new(c(1:100), c(1:100), c(1:100), x.scale = 0.2, y.scale = 1),
              times = 2)



add.section.contour = function(section, field = 'z', levels = NULL, col = 'black', lty = 1, lwd = 1, labels = NULL, cex.lab = 1) {
  z = matrix(section$grid[[field]], nrow = length(section$x))
  if (is.null(levels)) { levels = pretty(range(as.numeric(z)), n = 5) }

  contour(section$x, section$y, z, add = TRUE, levels = levels, col = col, lty = lty, lwd = lwd,
          labels = labels, labcex = cex.lab, drawlabels = !isFALSE(labels))
}

plot.section(a)

add.contour(a, field = 'z', levels = c(60, 70, 80, 90), f = 0.01)
add.contour2(a, col = c('blue', 'red'), lwd = c(1,4), lty = c(1,2))


