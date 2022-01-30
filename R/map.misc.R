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


#' @title Interpolate from fractional index onto grid
#' @author Thomas Bryce Kelly
#' @export
grid.interp = function(grid, i, j) {
  val = rep(NA, length(i))

  x = round(i)
  y = round(j)
  dx = i - x
  dy = j - y

  for (k in 1:length(i)) {
    val[k] = (1 - abs(dy[k])) * ((1 - abs(dx[k])) * grid[cbind(x[k], y[k])] + abs(dx[k]) * grid[cbind(x[k] + sign(dx[k]), y[k])]) +
      abs(dy[k]) * ((1 - abs(dx[k])) * grid[cbind(x[k], y[k] + sign(dy[k]))] + abs(dx[k]) * grid[cbind(x[k] + sign(dx[k]), y[k] + sign(dy[k]))])
  }

  ## Return
  val
}


#' @title Calculate extended grid
#' @author Thomas Bryce Kelly
#' @export
calc.vertex = function(lon, lat) {

  ## Diffs
  dlon.dx = t(diff(t(lon))) / 2
  dlon.dy = diff(lon) / 2
  dlat.dx = t(diff(t(lat))) / 2
  dlat.dy = diff(lat) / 2

  ## Vertex
  vertex.lon = matrix(NA, nrow = dim(lon)[1]+1, ncol = dim(lon)[2]+1)
  vertex.lat = matrix(NA, nrow = dim(lat)[1]+1, ncol = dim(lat)[2]+1)


  ## Field
  for (i in 2:(dim(vertex.lon)[1] - 1)) {
    for (j in 2:(dim(vertex.lon)[2] - 1)) {
      ii = max(i-1, 1)
      jj = max(j-1, 1)

      vertex.lon[i,j] = lon[ii,jj] + dlon.dx[ii,jj] + dlon.dy[ii,jj]
      vertex.lat[i,j] = lat[ii,jj] + dlat.dx[ii,jj] + dlat.dy[ii,jj]

    }
  }


  ## Fill in perimeter
  # i = 1
  for (j in 2:(dim(vertex.lon)[2] - 1)) {
    jj = max(j-1, 1)

    vertex.lon[1,j] = lon[1,jj] + dlon.dx[1,jj] - dlon.dy[1,jj]
    vertex.lat[1,j] = lat[1,jj] + dlat.dx[1,jj] - dlat.dy[1,jj]
  }

  # j = 1
  for (i in 2:(dim(vertex.lon)[1] - 1)) {
    ii = max(i-1, 1)

    vertex.lon[i,1] = lon[ii,1] - dlon.dx[ii,1] + dlon.dy[ii,1]
    vertex.lat[i,1] = lat[ii,1] - dlat.dx[ii,1] + dlat.dy[ii,1]
  }

  # j = dim(vertex.lon)[2]
  for (i in 1:(dim(vertex.lon)[1] - 1)) {
    ii = max(i-1, 1)
    j = dim(vertex.lon)[2]

    vertex.lon[i,j] = lon[ii,j-1] + dlon.dx[ii,j-2] + dlon.dy[ii,j-1]
    vertex.lat[i,j] = lat[ii,j-1] + dlat.dx[ii,j-2] + dlat.dy[ii,j-1]
  }

  # i = dim(vertex.lon)[2]
  for (j in 1:(dim(vertex.lon)[2] - 1)) {
    jj = max(j-1, 1)
    i = dim(vertex.lon)[1]

    vertex.lon[i,j] = lon[i-1,jj] + dlon.dx[i-1,jj] + dlon.dy[i-2,jj]
    vertex.lat[i,j] = lat[i-1,jj] + dlat.dx[i-1,jj] + dlat.dy[i-2,jj]
  }

  ## Fill in corners
  ## both = 1
  vertex.lon[1,1] = lon[1,1] - dlon.dx[1,1] - dlon.dy[1,1]
  vertex.lat[1,1] = lat[1,1] - dlat.dx[1,1] - dlat.dy[1,1]

  i = dim(vertex.lon)[1]
  vertex.lon[i,1] = lon[i-1,1] - dlon.dx[i-1,1] + dlon.dy[i-2,1]
  vertex.lat[i,1] = lat[i-1,1] - dlat.dx[i-1,1] + dlat.dy[i-2,1]

  i = dim(vertex.lon)[1]
  j = dim(vertex.lon)[2]
  vertex.lon[i,j] = lon[i-1,j-1] + dlon.dx[i-1,j-2] + dlon.dy[i-2,j-1]
  vertex.lat[i,j] = lat[i-1,j-1] + dlat.dx[i-1,j-2] + dlat.dy[i-2,j-1]

  j = dim(vertex.lon)[2]
  vertex.lon[1,j] = lon[1,j-1] + dlon.dx[1,j-2] - dlon.dy[1,j-1]
  vertex.lat[1,j] = lat[1,j-1] + dlat.dx[1,j-2] - dlat.dy[1,j-1]

  list(lon = vertex.lon, lat = vertex.lat)
}



calc.section.dist = function(lon, lat) {
  if (length(lon) == 1) { lon = rep(lon, length(lat))}
  if (length(lat) == 1) { lat = rep(lat, length(lon))}
  if (length(lon) != length(lat)) { stop('Length of lon/lat are not the same!')}

  d = rep(NA, length(lon))
  l = which(diff(lon) != 0 | diff(lat) != 0)

  d[1:(l[1]-1)] = 0 ## same station
  if (lenght(l) > 1) {
    for (i in 2:length(l)) {
      d[l[i-1]:(l[i]-1)] = calc.dist(lon[c(l[i],l[i]+1)], lat[c(l[i],l[i]+1)])
    }
  }

  ## return
  d
}

calc.dist = function(lon, lat) {
  sapply(1:(length(lon)-1), function(x) {abs(lon[x+1] - lon[x])})
}

