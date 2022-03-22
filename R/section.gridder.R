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
gridBin = function(gx, gy, x, y, z, p = 2, xscale = 1, yscale = 1, uncertainty = 0.1, neighborhood = NULL, x.factor = 1, y.factor = 1) {
  gz = rep(NA, length(gx))

  if (xscale == 0) { xscale = Inf }
  if (yscale == 0) { yscale = Inf }

  for (i in 1:length(gz)) {
    gz[i] = mean(z[abs(gx[i] - x) < xscale/2 & abs(gy[i] - y) < yscale/2], na.rm = T)
  }

  gz ## Return
}


#' @title Grid via Inverse Distance Weighting Interpolation
#' @author Thomas Bryce Kelly
#' @param gx Grid x values to interpolate onto
#' @param gy Grid y values to interpolate onto
#' @param x Observations, x values
#' @param y Observations, y values
#' @param z Observations, z values
#' @param p Exponent on the distance function
#' @export
gridIDW = function(gx, gy, x, y, z, p = 2, xscale = 1, yscale = 1, uncertainty = 0.1, neighborhood = NULL, x.factor = 1, y.factor = 1) {

  if (is.null(neighborhood)) {
    neighborhood = min(max(10, length(x)/10), length(x))
  }

  deltamin = sqrt((x.factor * xscale/2.0)^2 + (y.factor * yscale/2.0)^2) * uncertainty
  out = rep(NA, length(gx))

  for (i in 1:length(gx)) {
    w = sqrt((x.factor * (x - gx[i]))^2 + (y.factor * (y - gy[i]))^2 + deltamin)^-p

    k = order(w, decreasing = T)[1:neighborhood]
    out[i] = sum(z[k] * w[k]) / sum(w[k])
  }

  out
}


#' @title Grid via ODV's Interpolation
#' @author Thomas Bryce Kelly
#' @param gx Grid x values to interpolate onto
#' @param gy Grid y values to interpolate onto
#' @param x Observations, x values
#' @param y Observations, y values
#' @param z Observations, z values
#' @param p Exponent on the distance function
#' @inheritParams gridIDW
#' @export
gridODV = function (gx, gy, x, y, z, p = 2, xscale = 1, yscale = 1, uncertainty = 0.1, neighborhood = NULL, x.factor = 1, y.factor = 1) {

  if (is.null(neighborhood)) {
    neighborhood = min(max(10, length(x)/10), length(x))
  }

  deltamin = sqrt(x.factor^2 + y.factor^2) * uncertainty / 2 # scale / 2 / scale = 1/2
  out = rep(NA, length(gx))

  for (i in 1:length(gx)) {
    w = exp(-0.5*(sqrt((x.factor * (x - gx[i]) / xscale)^2 + (y.factor * (y - gy[i]) / yscale)^2) + deltamin))

    k = order(w, decreasing = T)[1:neighborhood]
    out[i] = sum(z[k] * w[k]) / sum(w[k])
  }

  out
}


#' @title Grid via Nearest Neighbor Interpolation
#' @author Thomas Bryce Kelly
#' @param gx Grid x values to interpolate onto
#' @param gy Grid y values to interpolate onto
#' @param x Observations, x values
#' @param y Observations, y values
#' @param z Observations, z values
#' @param p Exponent on the distance function
#' @inheritParams gridIDW
#' @export
gridNN = function (gx, gy, x, y, z, p = 2, xscale = 1, yscale = 1, uncertainty = 0.1, neighborhood = NULL, x.factor = 1, y.factor = 1) {
  out = rep(NA, length(gx))

  for (i in 1:length(gx)) {
    out[i] = z[which.min(x.factor * abs(x - gx[i])/xscale + y.factor * abs(y - gy[i])/yscale)]
  }
  out
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
gridNNI = function (gx, gy, x, y, z, p = 2, xscale = 1, yscale = 1, uncertainty = 0.1, neighborhood = NULL, x.factor = 1, y.factor = 1) {

  gz = rep(NA, length(gx))

  if (is.null(neighborhood)) {
    neighborhood = min(5, length(x))
  }
  ## Build standard NN index field
  gtemp1 = rep(1, length(gx))
  for (i in 1:length(gx)) {
    gtemp1[i] = which.min(abs(gx[i] - x) + abs(gy[i] - y))
  }

  for (i in 1:length(gx)) {

    ## Setup
    xx = c(x, gx[i])
    yy = c(y, gy[i])

    gtemp2 = gtemp1 ## Assume no change

    ## Find all the values that are likely to change
    dd = abs(x - gx[i]) * x.factor + abs(y - gy[i]) * y.factor
    l = which(abs(gx - gx[i]) * x.factor + abs(gy - gy[i]) * y.factor < dd[order(dd)[neighborhood]] / 2) ## only look at points within circle defined by the 4th closest data point

    for (k in l) {
      gtemp2[k] = which.min(x.factor * abs(gx[k] - xx)^p + y.factor * abs(gy[k] - yy)^p)
    }

    w = (gtemp2 - gtemp1) != 0
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
gridKrig = function (gx, gy, x, y, z, p = 2, xscale = 1, yscale = 1, uncertainty = 0.1, neighborhood = NULL, x.factor = 1, y.factor = 1) {

  data = data.frame(x = x, y = y, z = z)
  grid = data.frame(x = gx, y = gy)

  ## Set coordinate system (SpatialPointsDataFrame)
  sp::coordinates(data) = ~x + y
  sp::coordinates(grid) = ~x + y
  krige = automap::autoKrige(z ~ 1, data, new_data = grid)

  ## Return
  krige$krige_output$var1.pred
}


