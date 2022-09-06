u.wind = load.nc('C:/Users/Tom Kelly/Downloads/uwnd.10m.2021.nc')
advection = load.advection.oscar('Z:/Data/OSCAR/oscar_vel2021.nc', lats = c(0, 50), lons = c(-100,0) %% 360)


x = u.wind$lon %% 360
y = u.wind$lat

x.new = array(advection$lon%%360, dim = c(length(advection$lon), length(advection$lat)))
y.new = t(array(rep(advection$lat, length(advection$lon)), dim = c(length(advection$lat), length(advection$lon))))


summary(as.numeric(x))
summary(as.numeric(x.new))
summary(as.numeric(y))
summary(as.numeric(y.new))

regrid = function(x, y, x.new, y.new) {

  x.map = array(0, dim = dim(x.new))
  y.map = array(0, dim = dim(y.new))

  for (i in 1:dim(x.new)[1]) {
    message(i, ' \tof\t ', dim(x.new)[1])
    k = sapply(1:dim(x.new)[2], function(xx) {which.min(abs(x - x.new[i,xx]) + abs(y - y.new[i,xx])) - 1})
    x.map[i,] = k %% dim(x)[1] + 1
    y.map[i,] = floor(k / dim(x)[1]) + 1
  }


  list(x = x.map, y = y.map)
}


temp = regrid(x, y, x.new, y.new)





plot.image = function(x = NULL, y = NULL, z = NULL, col = NULL, xlab = NULL, ylab = NULL,
                      xlim = NULL, ylim = NULL, zlim = NULL,
                      pal = 'greyscale', n = 255, rev = F, ...) {

  if (is.null(z) & (is.null(x) | is.null(y))) { stop('plot.image: if z is not given then x and y are required, along with col.')}
  if (is.null(z) & is.null(col)) { stop('plot.image: Either z or col arrays are required.')}
  if (!is.null(z) & !is.null(col)) { message('plot.image: Both z and col arrays are provieded. Ignoring z.')}


  z = as.matrix(z)

  if (is.null(col)) {
    if (is.null(x)) {x = c(1:dim(z)[1])}
    if (is.null(y)) {y = c(1:dim(z)[2])}
    if (is.null(zlim)) {zlim = range(z, na.rm = T)}
  }

  if (is.null(xlim)) {xlim = range(x, na.rm = T)}
  if (is.null(ylim)) {ylim = range(y, na.rm = T)}

  if (is.null(xlab) & !is.null(x)) {xlab = deparse(substitute(x))} else {xlab = 'x'}
  if (is.null(ylab) & !is.null(y)) {ylab = deparse(substitute(y))} else {ylab = 'x'}

  if (is.null(col)) {
    xx = unique(x)
    yy = unique(y)

    if (length(z) != length(xx) * length(yy)) {
      warning('Z must have length equal to the grid generated from the unique values of x and y (i.e. z = f(x,y)).')
    }

    ## If data was not originally a matrix format.
    if (length(xx) != length(x) | length(yy) != length(y)) {
      z = as.numeric(z) ## Force into vector
      grid = expand.grid(x = xx, y = yy)
      grid$z = NA

      for (i in 1:length(z)) {
        grid$z[i] = z[x == grid$x[i] & y == grid$y[i]]
      }
      z = matrix(grid$z, nrow = length(xx), ncol = length(yy))
    }
    col = get.pal(n = n, pal = pal, rev = rev)
  } else {
    xx = x
    yy = y
  }

  image.default(x = xx[order(xx)],
                y = yy[order(yy)],
                z = matrix(z[order(xx), order(yy)], nrow = length(xx), ncol = length(yy)),
                zlim = zlim,
                xlim = xlim,
                ylim = ylim,
                col = col,
                xlab = xlab,
                ylab = ylab,
                ...)
}
