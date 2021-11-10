#' @title Load Coastline Data
#' @description Low res to high res: coastlineWorld, coastlineWorldCoarse, coastlinewWorldMedium, coastlineWorldFine
#' @author Laura Whitmore
#' @author Thomas Bryce Kelly
#' @param cost a string referring to a coastline data object such as those provided by the OCE package.
#' @import oce
#' @import ocedata
get.coast = function(coast = 'coastlineWorld') {
  library(ocedata)
  library(oce)
  do.call('data', list(coast))
}


#' @title Make Map Projection
#' @author Thomas Bryce Kelly
#' @param projection A string or number corresponding to a base projection (e.g. 1 = 'merc')
#' @param lat the primary latitude for the projection, may or may not be applicable based on projection
#' @param lon same as lat for longitude
#' @param h The viewing height in meters (used only for some projections)
#' @param dlat The distance to the secondary latitude (in degrees). Only applicable to some projections
#' @export
make.proj = function(projection = NULL, lat = NULL, lon = NULL, h = NULL, dlat = 10) {
  if (is.null(projection)) { projection = 'merc' }
  projections = c('merc', 'aea', 'eck3', 'eck4', 'eqc', 'geos', 'lonlat', 'mill', 'natearth', 'nsper', 'stere')
  projection.list = c('1) merc', '2) aea', '3) eck3', '4) eck4', '5) eqc', '6) geos', '7) lonlat', '8) mill', '9) natearth', '10) nsper', '11) stere')

  if (is.numeric(projection)) {
    projection = projections[projection]
  }
  if (!projection %in% projections) {
    stop('Invalid projection type, please provide one of the following:', paste(projection.list, collapse = ', '))
  }

  if (projection == 'geos' & is.null(h)) { h = 1e8 }
  if (projection == 'nsper' & is.null(h)) { h = 1e8 }

  # Default
  if (is.null(lat)) { lat = 0 }
  if (is.null(lon)) { lon = 0 }
  if (is.null(h)) { h = 1e8 }

  paste0('+proj=', projection, ' +lon_0=', lon, ' +lat_0=', lat, ' +lat_1=', lat, ' +lat_2=', lat + dlat, ' +h=', h)
}

#' @title Bilinear Interpolation of a Grid
#' @author Thomas Bryce Kelly
#' @export
bilinearInterp = function(x, y, gx, gy, z) {
  z.out = rep(0, length(x))

  for (i in 1:length(x)) {
    x1 = max(0, which(gx <= x[i]))
    x2 = x1 + 1
    y1 = max(0, which(gy <= y[i]))
    y2 = y1 + 1

    wx1 = (gx[x2] - x[i]) / (gx[2] - gx[1])
    wy1 = (gy[y2] - y[i]) / (gy[2] - gy[1])


    if (x1 == 0) {
      x1 = 1
      wx1 = 1
    }
    if (y1 == 0) {
      y1 = 1
      wy1 = 1
    }
    if(x1 == length(gx)) {
      x2 = x1
      wx1 = 1
    }
    if(y1 == length(gy)) {
      y2 = y1
      wy1 = 1
    }
    z.out[i] = wy1 * (wx1 * z[x1, y1] + (1 - wx1) * z[x2,y1]) + (1 - wy1) * (wx1 * z[x1, y2] + (1 - wx1) * z[x2,y2])
  }
  z.out
}
