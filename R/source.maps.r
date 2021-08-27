## Functions for plotting maps
## Please check out the wonderful OCE package! It does most of the work here, we just add wrapper functions to make the code easier for us to modify and customize.
##
## Author: Thomas Bryce Kelly1,2,3 (tbk14 at fsu.edu)
## http://about.tkelly.org/
## Author: Laura Whitmore4
## http://whitmorescience.com/
##
## 1Dept of Earth, Ocean & Atmospherical Sciences
## Florida State University
##
## 2Center for Ocean & Atmospheric Prediction Studies
## Florida State University
##
## 3National High Magnetic Field Laboratory
## Florida State University
##
## 4Division of Marine Sciences
## University of Southern Mississippi

#' @title Make Map
#' @description Default settings are for Arctic, lowest resolution coastline.
#' Exaples of projection arguments include: `+proj=stere +lat_0=90`, `+proj=merc`, `+proj=aea +lat_1=30 +lat_2=45 +lon_0=-120`.
#' @param coast Should be the name of a coastline data object. A value of NULL sets the default cosatline to 'coastlineWorld'.
#' @param lon.min The minimum longitude displayed on the map.
#' @param lon.max The maximum longitude on the map
#' @param lat.min The minimum latitude shown on the map
#' @param lat.max The maximum latitude shown on the map
#' @param p a projection string, such as those generated by 'make.proj()'
#' @param land.col A color string or value for the land polygon
#' @param grid Boolean, draw lat/lon grid?
#' @param dlon The spacing for the longitude grid (in degrees)
#' @param dlat The spacing for the latitude grid (in degrees)
#' @import oce
#' @import ocedata
#' @author Laura Whitmore
#' @author Thomas Bryce Kelly
#' @export
make.map = function (coast = NULL,
                    lon.min = -180,
                    lon.max = 180,
                    lat.min = -80,
                    lat.max = 80,
                    p = NULL,
                    land.col = 'lightgray',
                    grid = TRUE,
                    dlon = 15,
                    dlat = 15,
                    draw.axis = T) {

  if (is.null(coast)) { coast = 'coastlineWorld' }

  if (is.null(p)) { p = make.proj(projection = 'mill', lat = mean(c(lat.max, lat.min)),
                                  lon = mean(c(lon.min, lon.max)))}

  lon0 = as.numeric(strsplit(strsplit(paste(p, ''), 'lon_0=')[[1]][2], '\\s+')[[1]][1]) ## retreive the lon0 value from the projection

  ## get coastlines
  get.coast(coast = coast)
  lons = seq(-180, 180, by = dlon)
  lats = seq(-90, 90, by = dlat)

  ## make the base plot (oce)
  oce::mapPlot(coastlineCut(eval(parse(text = coast)), lon_0 = lon0),
               projection = p,
               col = land.col,
               longitudelim = c(lon.min, lon.max),
               latitudelim = c(lat.min, lat.max),
               grid = c(dlon, dlat),
               axes = draw.axis,
               lonlabels = lons,
               latlabels = lats)

  box()

  list(coast = coast, lon.min = lon.min, lon.max = lon.max,
       lat.min = lat.min, lat.max = lat.max, p = p, land.col = land.col,
       grid.dlon = dlon, grid.dlat = dlat, grid = grid)
}



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


#' @title Get Bathymetry
#' @description Get the bathymetry you want from NOAA (marmap)
#' @param map a map object as returned by make.map()
#' @param res the resolution of the bathymetery to be requested. Value is in arc minutes. Smaller values will take longer to download and to plot.
#' @param override a flag to override the default longitude bounds used in order to facilitate high latitude map making.
#' @param keep a flag used to determine if the bathymetry data should be saved to enable offline reuse later.
#' @import marmap
#' @author Laura Whitmore
#' @author Thomas Bryce Kelly
#' @export
get.bathy = function(map, res = 25, override = FALSE, keep = FALSE) {

  if (override) {
    lon.min = -180
    lon.max = 180
  } else {
    lon.min = max(-180, map$lon.min - (map$lon.max - map$lon.min)/1.5)
    lon.max = min(180, map$lon.max + (map$lon.max - map$lon.min)/1.5)
  }

  lat.min = max(-90, map$lat.min - (map$lat.max - map$lat.min)/1.5)
  lat.max = min(90, map$lat.max + (map$lat.max - map$lat.min)/1.5)

  ##Load Bathymetry (getNOAA.bathy is marmaps)
  bathy = getNOAA.bathy(lon1 = lon.min, lon2 = lon.max,
                        lat1 = lat.min, lat2 = lat.max,
                        resolution = res, keep = keep)

  ##Convert bathy
  bathyLon = as.numeric(rownames(bathy))
  bathyLat = as.numeric(colnames(bathy))
  bathyZ = as.numeric(bathy)
  dim(bathyZ) = dim(bathy)

  ## return
  list(Lon = bathyLon, Lat = bathyLat, Z = bathyZ, res = res)
}


#' @title Add Map Bathymetry
#' @description add bathymetry contours
#' @author Thomas Bryce Kelly
#' @param map a map object as returned by make.map()
#' @param bathy a bathymetry dataframe as produced by get.bathy()
#' @param bathy.levels a vector of elevations for which to draw contours. NB: depths are negative elevation.
#' @param bathy.col the color of the lines to be drawn
#' @param bathy.lwd the weight of the lines
#' @param bathy.lty the type of line to be drawn
#' @param drawlabels boolean value to draw labels (default: TRUE)
#' @export
add.map.bathy = function(map, bathy,
                         levels = c(-100, -1000, -2500),
                         bathy.col = 'darkgrey',
                         bathy.lwd = 1,
                         bathy.lty = 1,
                         drawlabels = TRUE,
                         cex.lab = 1,
                         bathy.levels = NULL){

  if (!is.null(bathy.levels)) {
    warning('Option bathy.levels is depreciated, please use levels argument instead.')
    levels = bathy.levels
  }

  ## add bathy contours (oce)
  oce::mapContour(longitude = bathy$Lon,
               latitude = bathy$Lat,
               z = bathy$Z,
               levels = levels,
               lwd = bathy.lwd,
               lty = bathy.lty,
               col = bathy.col,
               drawlabels = drawlabels,
               labcex = cex.lab)
}

#' @title Add Map Points
#' @author Laura Whitmore
#' @description Add station points to a map
#' @param lon the longitudes of the points to be drawn
#' @param lat the latitude of the points to be drawn
#' @param stn.lon (depreciated) the longitudes of the points to be drawn
#' @param stn.lat (depreciated) the latitude of the points to be drawn
#' @param col the color of the points
#' @param cex the size of point to be drawn
#' @param pch the point character to be used
#' @param ... optional arguments passed to points().
#' @import oce
#' @export
add.map.points = function(lon = NULL, lat = NULL, col = 'black', cex = 1, pch = 16,
                          stn.lon = NULL, stn.lat = NULL, override = F, ...){

  if (!is.null(stn.lon)) { warning('stn.lon option depreciated. Recommend using lon instead.')}
  if (!is.null(stn.lat)) { warning('stn.lat option depreciated. Recommend using lat instead.')}

  if (interactive() & length(lon) > 1e5 & !override) {
    input = readline(prompt = 'Large number of points, contnue? (y/n)')
    if (input != 'y' & input != 'Y') {
      return()
    }
  }
  oce::mapPoints(c(lon, stn.lon), c(lat, stn.lat), col = col, cex = cex, pch = pch, ...)
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

#' @title Add Map Layer
#' @description  Add a image layer to the map!
#' @author Thomas Bryce Kelly
#' @param lon longitude of the data
#' @param lat latitude of the data
#' @param z the data
#' @param zlim the limits of the z-axis. A value of NULL indicates that zlim should equal the range of z values.
#' @param pal the name of a palette generating function. Should be a string.
#' @param col.na the color with which any NAs should be drawn with. A value of NA skips this step.
#' @param n the number of distinct colors to request from the palette function.
#' @param refinement the level of bilinear refinement to apply to the image grid
#' @export
add.map.layer = function(lon, lat, z, zlim = NULL, pal = 'greyscale', col.na = NA, n = 255, refinement = NA,
                         col.low = NA, col.high = NA, rev = FALSE, filled = FALSE, indicate = TRUE, verbose = FALSE) {

  lon = lon %% 360
  if(is.null(ncol(lon)) & is.null(ncol(lat))) {
    nn = length(unique(lon))
    mm = length(unique(lat))

    if (length(z) == nn * mm) {
      lat = unique(lat)
      lon = unique(lon)
      z = matrix(z, nrow = nn, ncol = mm)
    } else {
      stop('Cannot transform data into matrix. Number of rows and columns do not match length of z.')
    }
  }

  if (is.null(zlim)) { zlim = range(pretty(as.numeric(z), na.rm = TRUE)) }

  if (!is.na(refinement) & refinement > 0) {
    for (i in 1:refinement) {
      lat.old = lat
      lon.old = lon
      lat = seq(lat[1], lat[length(lat)], length.out = length(lat) * 2)
      lon = seq(lon[1], lon[length(lon)], length.out = length(lon) * 2)

      grid = expand.grid(lon = lon, lat = lat)
      z = bilinearInterp(grid$lon * sign(lon[2] - lon[1]), grid$lat * sign(lat[2] - lat[1]), lon.old * sign(lon[2] - lon[1]), lat.old * sign(lat[2] - lat[1]), z)
      z = matrix(z, nrow = length(lon), ncol = length(lat))
    }
  }

  z.good = z
  z.good[z < zlim[1]] = zlim[1]
  z.good[z > zlim[2]] = zlim[2]
  z.good[is.infinite(z.good)] = NA

  ## Color scale
  z.min = min(as.numeric(z.good), na.rm = TRUE)
  z.max = max(as.numeric(z.good), na.rm = TRUE)

  n.low = round((z.min - zlim[1]) / (zlim[2] - zlim[1]) * n)
  n.high = round((zlim[2] - z.max) / (zlim[2] - zlim[1]) * n)
  col = get.pal(n, pal = pal)
  if (rev) {col = rev(col)}
  col = col[(n.low + 1):(n - n.high)]

  if (n.low + n.high > n / 2) {
      warning('Color scale exceeds data by >50%, colors may not reflect quantitative values.')
  }

  if (verbose) { print(paste('n.low:', n.low)); print(paste('n.high:', n.high)); print(paste('cols:', (n.low + 1), ':', (n - n.high))); }

  if (length(!is.na(z.good)) > 0) {
    oce::mapImage(lon, lat, z.good, zlim = range(z.good, z.good + 1e-12, na.rm = TRUE), missingColor = col.na,
               col = col, filledContour = filled)
  }

  if (!is.na(col.low) & any(z < zlim[1], na.rm = TRUE)) {
    z.low = z
    z.low[z >= zlim[1]] = NA
    z.low.lim = range(z.low, na.rm = TRUE)

    if (col.low != '') {
      oce::mapImage(lon, lat, z.low, zlim = z.low.lim, breaks = 1, col = rep(col.low, 2))
    } else {
      oce::mapImage(lon, lat, z.low, zlim = z.low.lim, breaks = 1, col = get.pal(n, pal = pal, rev = rev)[1:2])
    }
  }

  if (!is.na(col.high) & any(z > zlim[2], na.rm = TRUE)) {
    z.high = z
    z.high[z <= zlim[2]] = NA

    if (col.high != '') {
      oce::mapImage(lon, lat, z.high, zlim = range(z.high, na.rm = TRUE), breaks = 1, col = rep(col.high, 2))
    } else {
      oce::mapImage(lon, lat, z.high, zlim = range(z.high, na.rm = TRUE), breaks = 1, col = get.pal(n, pal = pal, rev = rev)[(n-1):n])
    }
  }

  if (indicate) {
    st = paste0(
      'zaxis: (', round(min(z.good, na.rm = TRUE), 3), ', ', round(max(z.good, na.rm = TRUE), 3),
      ')   zlim: (', round(zlim[1], 3), ', ', round(zlim[2], 3), ')'
      )
    mtext(st, line = 0.25, adj = 1, cex = 0.7)
  }
}


#' @title Add Map Bathymetry (shaded)
#' @author Thomas Bryce Kelly
#' @param map a map object returned by the function make.map()
#' @param bathy a bathymetric dataframe returned by get.map.bathy
#' @param pal the color palette used for coloring the shading. Should consist of a function name as a string
#' @param n the number of distinct colors to use
#' @param zlim the z-axis range that the colormap is projected upon
#' @param rev a boolean flag to flip the color axis
#' @param filled a flag to turn on color smoothing and filling from the OCE package. Works for some datasets but may cause issues.
#' @param col.low the color used when painting data that is below the zlim range. A value of NA causes the data not to be drawn, and a value of '' uses the lowest color value for out of range data.
#' @param col.high same as \emph{col.low} but for data that is out of range above the zlim.
#' @param refinement the level of bilinear refinement to apply to the image grid
#' @export
add.map.bathy.shade = function(map, bathy, pal = 'greyscale', n = 255, zlim = NULL, rev = FALSE, filled = TRUE, col.low = NA, col.high = NA, refinement = 1) {
  bathy$Z[bathy$Z > 0] = 0
  add.map.layer(lon = bathy$Lon, lat = bathy$Lat, z = bathy$Z, pal = pal, rev = rev, filled = filled, refinement = refinement,
                zlim = zlim, col.low = col.low, col.high = col.high, indicate = FALSE, n = n)

  if (filled == TRUE) {
    ## Hack to keep anti-aliasing at bay
    add.map.layer(lon = bathy$Lon, lat = bathy$Lat, z = bathy$Z, pal = pal, rev = rev, filled = filled, refinement = refinement,
                  zlim = zlim, col.low = col.low, col.high = col.high, indicate = FALSE, n = n)
    add.map.layer(lon = bathy$Lon, lat = bathy$Lat, z = bathy$Z, pal = pal, rev = rev, filled = filled, refinement = refinement,
                  zlim = zlim, col.low = col.low, col.high = col.high, indicate = FALSE, n = n)
  }
  redraw.map(map)
}


#' @title Add Map Text
#' @author Thomas Bryce Kelly
#' @author Laura Whitmore
#' @import oce
#' @param lon longitude position of the text or a vector of positions
#' @param lat latitude position of the text or a vector of positions
#' @param text a string or vector of strings for the text to be written
#' @param col text color
#' @param cex size of the text to be written
#' @param adj	one or two values in [0, 1] which specify the x (and optionally y) adjustment of the labels.
#' @param pos a position specifier for the text. If specified this overrides any adj value given. Values of 1, 2, 3 and 4, respectively indicate positions below, to the left of, above and to the right of the specified coordinates.
#' @export
add.map.text = function(lon, lat, text, col = 'black', cex = 1, adj = NULL, pos = NULL){
  oce::mapText(longitude = lon, latitude = lat, labels = text, col = col, cex = cex, adj = adj, pos = pos)
}


#' @title Add Map Line
#' @author Thomas Bryce Kelly
#' @param lon a set of longitudes to draw a line between
#' @param lat a set of latitudes to draw a lien between
#' @param col the color of the line to be drawn
#' @param lty the type of line to be drawn
#' @param lty the type of line to be drawn
#' @param lwd the width of the line to be drawn
#' @export
add.map.line = function(lon, lat, col = 'black', lty = 1, lwd = 1) {
  oce::mapLines(longitude = lon, latitude = lat, col = col, lty = lty, lwd = lwd)
}


#' @title Add Quiver Lines
#' @export
#' @author Thomas Bryce Kelly
add.map.quiver = function(lon, lat, u, v, zscale = 1, col = 'black', lwd = 1) {
  add.map.points(lon, lat, pch = 20, cex = 0.5, col = col)
  for (i in 1:length(lon)) {
    add.map.line(c(lon[i], lon[i] + u[i] * zscale), c(lat[i], lat[i] + v[i] * zscale), col = col, lwd = lwd)
  }
}

#' @export
add.map.scale = function(x) {
  oce::mapScalebar(x)
}


#' @title Add Map Contours
#' @export
#' @author Thomas Bryce Kelly
add.map.contour = function(x, y, z, col = 'black', levels = NULL, n = 3, labels = TRUE, lty = 1, lwd = 1) {

  if (is.null(levels)) {
    levels = pretty.default(z, n = n)
  }
  z = matrix(z, nrow = length(unique(x)))

  oce::mapContour(unique(x), unique(y), z, levels = levels, drawlabels = labels,
             underlay = 'interrupt', lwd = lwd, lty = lty,
             col = col)
}

#' @title Add Map Polygon
#' @author Laura Whitmore
#' @export
add.map.polygon = function (lon, lat, col = "#00000020", lty = 1, lwd = 1, border = NULL, density = NULL, angle = 45, fillOddEven = FALSE) {
  oce::mapPolygon(longitude = lon, latitude = lat, col = col, density = density, angle = angle,
                  lty = lty, lwd = lwd, border = border, fillOddEven = fillOddEven)
}


#' @title Redraw Map
#' @author Thomas Bryce Kelly
#' @param map a map object as returned by make.map()
#' @export
redraw.map = function(map,
                      coast = NULL,
                      land.col = NULL,
                      grid = NULL) {

  if (!is.null(coast)) { map$coast = coast}
  if (!is.null(land.col)) { map$land.col = land.col}
  if (!is.null(grid)) { map$grid = grid}
  if (!is.null(coast)) { map$coast = coast}
  if (!is.null(coast)) { map$coast = coast}


  mapPolygon(eval(parse(text = map$coast)), col = map$land.col)

  lons = seq(-180, 180, by = map$grid.dlon)
  lats = seq(-90, 90, by = map$grid.dlat)
  if (map$grid) { mapGrid(longitude = lons, latitude = lats) }
  box()
  #oce::mapAxis(1, longitude = lons) ## Is this needed?
  #oce::mapAxis(2, latitude = lats) ## Is this needed?
}

#' @title Replot Map
#' @author Thomas Bryce Kelly
#' @param map a map object as returned by make.map()
#' @export
replot.map = function(map = map) {
  make.map(coast = map$coast,
           lon.min = map$lon.min,
           lon.max = map$lon.max,
           lat.min = map$lat.min,
           lat.max = map$lat.max,
           p = map$p,
           land.col = map$land.col,
           grid = map$grid,
           dlon = map$grid.dlon,
           dlat = map$grid.dlat)
  NULL
}


#' @title Make an Arctic Map
#' @author Laura Whitmore
#' @description An example map of the Arctic -- basemap only. For a full test map call test.map.arctic()
#' @inheritParams make.map
#' @export
make.map.arctic = function(coast = 'coastlineWorld', lon.min = -180, lon.max = 180, lat.min = 60, lat.max = 90,
                           p = make.proj('stere', lat = 90), dlat = 15, dlon = 60, land.col = 'lightgray', grid = TRUE) {

  make.map(coast = coast, lon.min = lon.min, lon.max = lon.max, lat.min = lat.min, lat.max = lat.max,
                 p = p, land.col = land.col, grid = grid, dlon = dlon, dlat = dlat)
}

#' @title Make Map (NGA)
#' @author Thomas Bryce Kelly
#' @description An example map of the NGA -- basemap only.
#' @inheritParams make.map
#' @export
make.map.nga = function(coast = 'coastlineWorldFine', lon.min = -153, lon.max = -150, lat.min = 56, lat.max = 60.5,
                           p = make.proj('stere', lat = 60, lon = -150), dlat = 3, dlon = 3, land.col = '#333333', grid = TRUE) {

  make.map(coast = coast, lon.min = lon.min, lon.max = lon.max, lat.min = lat.min, lat.max = lat.max,
           p = p, land.col = land.col, grid = grid, dlon = dlon, dlat = dlat)
}

#' @title Make Map (CCE)
#' @author Thomas Bryce Kelly
#' @description An example map of the CCE study region.
#' @inheritParams make.map
#' @import oce
#' @import marmap
#' @export
make.map.cce = function (coast = "coastlineWorldFine",
                         lon.min = -126,
                         lon.max = -119,
                         lat.min = 31,
                         lat.max = 36,
                         p = make.proj(projection = 'merc', lat = 35, lon = -122, dlat = 3),
                         land.col = '#252525',
                         grid = TRUE,
                         dlon = 3,
                         dlat = 3) {

  make.map(coast = coast, p = p,
           lat.min = lat.min, lat.max = lat.max, lon.min = lon.min, lon.max = lon.max,
           dlat = dlat, dlon = dlon, grid = grid, land.col = land.col)
}


##############
#### tests ###
##############
#' @title Test Map of the Arctic
#' @author Laura Whitmore
#' @export
test.map.arctic = function() {
    ##run as is - default settings for Arctic, coastlineWorld
    map = make.map.arctic()

    ## Add bathy
    add.map.bathy.shade(map, bathy.arctic, zlim = c(-6e3, -20), filled = FALSE)
    add.map.bathy(map, bathy.arctic)

    map
}

#' @title Test Map of Mississippi
#' @author Laura Whitmore
#' @export
test.map.mississippi = function() {
    ## test for Mississippi Sound region
    map = make.map(coast = "coastlineWorldFine",
                   p = "+proj=merc",
                   lat.min = 27,
                   lat.max = 31,
                   lon.min = -93,
                   lon.max = -84,
                   dlat = 5,
                   dlon = 5)

    add.map.bathy.shade(map, bathy.gom, zlim = c(-3e3, 10), filled = TRUE)
    add.map.bathy(map, bathy.gom,
                  bathy.levels = c(-10, -50, -100, -1000, -3000),
                  bathy.lwd = c(0.8, 0.8, 1.2, 0.8, 0.8))
    map
}

#' @title Test Map of California
#' @author Thomas Bryce Kelly
#' @import oce
#' @import ocedata
#' @export
test.map.california = function(coast = "coastlineWorldFine", p = "+proj=merc",
                               lat.min = 31, lat.max = 36, lon.min = -126, lon.max = -119,
                               dlat = 5, dlon = 5, grid = T) {
    ## test for California Current
    map = make.map(coast = coast, p = p,
             lat.min = lat.min, lat.max = lat.max, lon.min = lon.min, lon.max = lon.max,
             dlat = dlat, dlon = dlon, grid = grid)
    add.map.bathy.shade(map, bathy.pacific, refinement = 0)
    add.map.bathy(map, bathy.pacific)
    map
}

