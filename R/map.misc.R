#' @title Load Coastline Data
#' @description Low res to high res: coastlineWorld, coastlineWorldCoarse, coastlinewWorldMedium, coastlineWorldFine
#' @author Laura Whitmore
#' @author Thomas Bryce Kelly
#' @param cost a string referring to a coastline data object such as those provided by the OCE package.
#' @import oce
#' @import ocedata
get.coast = function(coast = 'coastlineWorld') {
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
interp.bilinear = function(x, y, gx, gy, z) {
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
calc.vertex = function(x, y) {

  ## Diffs
  dx.dx = t(diff(t(x))) / 2
  dx.dy = diff(x) / 2
  dy.dx = t(diff(t(y))) / 2
  dy.dy = diff(y) / 2

  ## Vertex
  vertex.x = matrix(NA, nrow = dim(x)[1]+1, ncol = dim(x)[2]+1)
  vertex.y = matrix(NA, nrow = dim(y)[1]+1, ncol = dim(y)[2]+1)


  ## Field
  for (i in 2:(dim(vertex.x)[1] - 1)) {
    for (j in 2:(dim(vertex.x)[2] - 1)) {
      ii = max(i-1, 1)
      jj = max(j-1, 1)

      vertex.x[i,j] = x[ii,jj] + dx.dx[ii,jj] + dx.dy[ii,jj]
      vertex.y[i,j] = y[ii,jj] + dy.dx[ii,jj] + dy.dy[ii,jj]

    }
  }


  ## Fill in perimeter
  # i = 1
  for (j in 2:(dim(vertex.x)[2] - 1)) {
    jj = max(j-1, 1)

    vertex.x[1,j] = x[1,jj] + dx.dx[1,jj] - dx.dy[1,jj]
    vertex.y[1,j] = y[1,jj] + dy.dx[1,jj] - dy.dy[1,jj]
  }

  # j = 1
  for (i in 2:(dim(vertex.x)[1] - 1)) {
    ii = max(i-1, 1)

    vertex.x[i,1] = x[ii,1] - dx.dx[ii,1] + dx.dy[ii,1]
    vertex.y[i,1] = y[ii,1] - dy.dx[ii,1] + dy.dy[ii,1]
  }

  # j = dim(vertex.x)[2]
  for (i in 1:(dim(vertex.x)[1] - 1)) {
    ii = max(i-1, 1)
    j = dim(vertex.x)[2]

    vertex.x[i,j] = x[ii,j-1] + dx.dx[ii,j-2] + dx.dy[ii,j-1]
    vertex.y[i,j] = y[ii,j-1] + dy.dx[ii,j-2] + dy.dy[ii,j-1]
  }

  # i = dim(vertex.x)[2]
  for (j in 1:(dim(vertex.x)[2] - 1)) {
    jj = max(j-1, 1)
    i = dim(vertex.x)[1]

    vertex.x[i,j] = x[i-1,jj] + dx.dx[i-1,jj] + dx.dy[i-2,jj]
    vertex.y[i,j] = y[i-1,jj] + dy.dx[i-1,jj] + dy.dy[i-2,jj]
  }

  ## Fill in corners
  ## both = 1
  vertex.x[1,1] = x[1,1] - dx.dx[1,1] - dx.dy[1,1]
  vertex.y[1,1] = y[1,1] - dy.dx[1,1] - dy.dy[1,1]

  i = dim(vertex.x)[1]
  vertex.x[i,1] = x[i-1,1] - dx.dx[i-1,1] + dx.dy[i-2,1]
  vertex.y[i,1] = y[i-1,1] - dy.dx[i-1,1] + dy.dy[i-2,1]

  i = dim(vertex.x)[1]
  j = dim(vertex.x)[2]
  vertex.x[i,j] = x[i-1,j-1] + dx.dx[i-1,j-2] + dx.dy[i-2,j-1]
  vertex.y[i,j] = y[i-1,j-1] + dy.dx[i-1,j-2] + dy.dy[i-2,j-1]

  j = dim(vertex.x)[2]
  vertex.x[1,j] = x[1,j-1] + dx.dx[1,j-2] - dx.dy[1,j-1]
  vertex.y[1,j] = y[1,j-1] + dy.dx[1,j-2] - dy.dy[1,j-1]

  list(x = vertex.x, y = vertex.y)
}



#' @title Calcuye extended grid
#' @author Thomas Bryce Kelly
#' @export
grid.refinement = function(x = NULL, y = NULL, z) {
  if (is.null(dim(z))) {stop('grid.refinement: z must be an array object of two dimensions.')}
  dim = dim(z)

  if (is.null(x) & is.null(y)) {
    x = c(1:dim[1])
    y = c(1:dim[2])
  }

  if (is.null(dim(x)) & is.null(dim(y))) {
    x = array(x, dim = dim)
    y = t(array(y, dim = rev(dim)))
  }

  ## Vertex
  vertex.x = array(0, dim = c(2*dim(x)[1]-1, 2*dim(x)[2]-1))
  vertex.y = vertex.x
  vertex.z = vertex.x

  ## fill in known values
  for (i in 1:dim(x)[1]) {
    for (j in 1:dim(x)[2]) {
      vertex.x[2*i-1, 2*j-1] = x[i,j]
      vertex.y[2*i-1, 2*j-1] = y[i,j]
      vertex.z[2*i-1, 2*j-1] = z[i,j]
    }
  }

  ## Interpolate x
  for (i in 1:(dim(x)[1]-1)) {
    for (j in 1:dim(x)[2]) {
      vertex.x[2*i, 2*j-1] = 0.5 * (x[i,j] + x[i+1,j])
      vertex.y[2*i, 2*j-1] = 0.5 * (y[i,j] + y[i+1,j])
      vertex.z[2*i, 2*j-1] = 0.5 * (z[i,j] + z[i+1,j])
    }
  }

  ## Interpolate y
  for (i in 1:dim(x)[1]) {
    for (j in 1:(dim(x)[2]-1)) {
      vertex.x[2*i-1, 2*j] = 0.5 * (x[i,j] + x[i,j+1])
      vertex.y[2*i-1, 2*j] = 0.5 * (y[i,j] + y[i,j+1])
      vertex.z[2*i-1, 2*j] = 0.5 * (z[i,j] + z[i,j+1])
    }
  }

  ## corners
  for (i in 1:(dim(x)[1]-1)) {
    for (j in 1:(dim(x)[2]-1)) {
      vertex.x[2*i, 2*j] = 0.25 * (x[i,j] + x[i,j+1] + x[i+1,j] + x[i+1,j+1])
      vertex.y[2*i, 2*j] = 0.25 * (y[i,j] + y[i,j+1] + y[i+1,j] + y[i+1,j+1])
      vertex.z[2*i, 2*j] = 0.25 * (z[i,j] + z[i,j+1] + z[i+1,j] + z[i+1,j+1])
    }
  }

  list(x = vertex.x, y = vertex.y, z = vertex.z)
}


#' @title Calculate subsampled grid
#' @author Thomas Bryce Kelly
#' @export
grid.subsample = function(x = NULL, y = NULL, z, approx = F) {
  if (is.null(dim(z))) {stop('grid.refinement: z must be an array object of two dimensions.')}
  dim = dim(z)

  if (is.null(x) & is.null(y)) {
    x = c(1:dim[1])
    y = c(1:dim[2])
  }

  if (is.null(dim(x)) & is.null(dim(y))) {
    x = array(x, dim = dim)
    y = t(array(y, dim = rev(dim)))
  }

  ## Vertex
  vertex.x = array(0, dim = floor(c(dim(x)[1], dim(x)[2]) / 2))
  vertex.y = vertex.x
  vertex.z = vertex.x

  if (approx) {
    for (i in 1:dim(vertex.x)[1]) {
      for (j in 1:dim(vertex.x)[2]) {
        vertex.x[i, j] = x[2*i, 2*j]
        vertex.y[i, j] = y[2*i, 2*j]
        vertex.z[i, j] = z[2*i, 2*j]
      }
    }
  } else {
    for (i in 1:dim(vertex.x)[1]) {
      for (j in 1:dim(vertex.x)[2]) {
        vertex.x[i, j] = 0.25 * (x[2*i, 2*j] + x[2*i-1, 2*j] + x[2*i, 2*j-1] + x[2*i-1, 2*j-1])
        vertex.y[i, j] = 0.25 * (y[2*i, 2*j] + y[2*i-1, 2*j] + y[2*i, 2*j-1] + y[2*i-1, 2*j-1])
        vertex.z[i, j] = mean(c(z[2*i, 2*j], z[2*i-1, 2*j], z[2*i, 2*j-1], z[2*i-1, 2*j-1]), na.rm = T)
      }
    }
  }

  list(x = vertex.x, y = vertex.y, z = vertex.z)
}



#' @export
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


#' @export
calc.dist = function(lon, lat) {
  sapply(1:(length(lon)-1), function(x) {abs(lon[x+1] - lon[x])})
}


#' @title Retreive depth value from bathymetric grid.
#' @export
get.depth = function(lon, lat, bathy) {
  depths = rep(NA, length(lon))
  
  for (i in 1:lnegth(lon)) {
    k1 = which.min(abs(lon[i] - bathy$Lon))
    k2 = which.min(abs(lat[i] - bathy$Lat))
    depths[i] = bathy$Z[k1,k2]
  }
  
  ## Return
  depths
}

