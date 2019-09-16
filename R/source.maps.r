#' @title Make Map
#' @description Default settings are for Arctic, lowest resolution coastline.
#' Exaples of projection arguments include: `+proj=stere +lat_0=90`, `+proj=merc`, `+proj=aea +lat_1=30 +lat_2=45 +lon_0=-120`.
#' @params coast Should be the name of a coastline data object. A value of NULL sets the default cosatline to 'coastlineWorld'.
#' @params lon.min The minimum longitude displayed on the map.
#' @import oce
#' @import ocedata
#' @author Laura Whitmore
#' @author Thomas Bryce Kelly
#' @export
make.map = function (coast = NULL,
                    lon.min = -180,
                    lon.max = 180,
                    lat.min = 55,
                    lat.max = 90,
                    p = '+proj=stere +lat_0=90',
                    land.col = 'lightgray',
                    grid = TRUE,
                    dlon = 5,
                    dlat = 5) {

  if (is.null(coast)) { coast = 'coastlineWorld' }

  ## get coastlines
  get.coast(coast = coast)

  ## make the base plot (oce)
  mapPlot(eval(parse(text = coast)), projection = p, col = land.col,
          longitudelim = c(lon.min, lon.max), latitudelim = c(lat.min, lat.max),
          grid = FALSE)

  ## add a grid at the lat and lon interval you set (oce)
  if (grid) {
    mapGrid(dlongitude = dlon, dlatitude = dlat)
  }
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
#' @export
get.coast = function(coast = 'coastlineWorld') {
  library(ocedata)
  library(oce)
  do.call('data', list(coast))
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
#' @export
add.map.bathy = function(map, bathy,
                         bathy.levels = c(-100, -1000, -2500),
                         bathy.col = 'darkgrey',
                         bathy.lwd = 1,
                         bathy.lty = 1){


  ## add bathy contours (oce)
  mapContour(longitude = bathy$Lon,
               latitude = bathy$Lat,
               z =  bathy$Z,
               levels = bathy.levels,
               lwd = bathy.lwd,
               lty = bathy.lty,
               col = bathy.col)

  ## Redraw
  mapPolygon(eval(parse(text = map$coast)), col = map$land.col)
  mapGrid(dlongitude = map$grid.dlon, dlatitude = map$grid.dlat)
  box()
}

#' @title Add Map Points
#' @author Laura Whitmore
#' @description Add station points to a map
#' @param stn.lon the longitudes of the points to be drawn
#' @param stn.lat the latitude of the points to be drawn
#' @param col the color of the points
#' @param cex the size of point to be drawn
#' @param pch the point character to be used
#' @import oce
#' @export
add.map.points = function(stn.lon, stn.lat, col = 'black', cex = 1, pch = 16){
  mapPoints(stn.lon, stn.lat, col = col, cex = cex, pch = pch)
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
#' @export
add.map.layer = function(lon, lat, z, zlim = NULL, pal = 'ocean.algae', col.na = NA, n = 255,
                         col.low = NA, col.high = NA, rev = FALSE, filled = FALSE, indicate = TRUE, verbose = FALSE) {

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

  if (is.null(zlim)) { zlim = range(as.numeric(z), na.rm = TRUE) }

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
      mapImage(lon, lat, z.good, zlim = range(z.good, z.good + 1e-12, na.rm = TRUE), missingColor = col.na,
               col = col, filledContour = filled)
  }

  if (!is.na(col.low) & any(z < zlim[1], na.rm = TRUE)) {
    z.low = z
    z.low[z >= zlim[1]] = NA
    z.low.lim = range(z.low, na.rm = TRUE)

    if (col.low != '') {
        mapImage(lon, lat, z.low, zlim = z.low.lim, breaks = 1, col = rep(col.low, 2))
    } else {
        mapImage(lon, lat, z.low, zlim = z.low.lim, breaks = 1, col = get.pal(n, pal = pal, rev = rev)[1:2])
    }
  }

  if (!is.na(col.high) & any(z > zlim[2], na.rm = TRUE)) {
    z.high = z
    z.high[z <= zlim[2]] = NA

    if (col.high != '') {
        mapImage(lon, lat, z.high, zlim = range(z.high, na.rm = TRUE), breaks = 1, col = rep(col.high, 2))
    } else {
        mapImage(lon, lat, z.high, zlim = range(z.high, na.rm = TRUE), breaks = 1, col = get.pal(n, pal = pal, rev = rev)[(n-1):n])
    }
  }

  if (indicate) {
    st = paste0(
      'zaxis: (', round(min(z.good, na.rm = TRUE), 3), ', ', round(max(z.good, na.rm = TRUE), 3),
      ')   zlim: (', round(zlim[1], 3), ', ', round(zlim[2], 3), ')'
      )
    mtext(st, line = 0.25, adj = 1, cex = 0.7)
  }
  box()
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
#' @export
add.map.bathy.shade = function(map, bathy, pal = 'ocean.ice', n = 255, zlim = NULL, rev = FALSE, filled = TRUE, col.low = NA, col.high = NA) {
  bathy$Z[bathy$Z > 0] = 0
  add.map.layer(lon = bathy$Lon, lat = bathy$Lat, z = bathy$Z, pal = pal, rev = rev, filled = filled,
                zlim = zlim, col.low = col.low, col.high = col.high, indicate = FALSE, n = n)
  redraw.map(map)
}


#' @author Thomas Bryce Kelly
add.map.topo.shade = function(map, bathy, pal = 'ocean.turbid', zlim = NULL, rev = FALSE, filled = TRUE, col.low = NA, col.high = NA) {
  bathy$Z[bathy$Z < 0] = 0
  add.map.layer(lon = bathy$Lon, lat = bathy$Lat, z = bathy$Z, pal = pal, rev = rev, filled = filled,
                zlim = zlim, col.low = col.low, col.high = col.high, indicate = FALSE)
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
  mapText(longitude = lon, latitude = lat, labels = text, col = col, cex = cex, adj = adj, pos = pos)
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
  mapLines(longitude = lon, latitude = lat, col = col, lty = lty, lwd = lwd)
}

#' @title Redraw Map
#' @author Thomas Bryce Kelly
#' @param map a map object as returned by make.map()
#' @export
redraw.map = function(map) {
  mapPolygon(eval(parse(text = map$coast)), col = map$land.col)
  if (map$grid) {
    mapGrid(dlongitude = map$grid.dlon, dlatitude = map$grid.dlat)
  }
  box()
}


##############
#### tests ###
##############
#' @title Test Map of the Arctic
#' @author Laura Whitmore
#' @export
test.map.arctic = function() {
    ##run as is - default settings for Arctic, coastlineWorld
    map = make.map(dlon = 60)
    bathy = get.bathy(map)

    ## Add bathy
    add.map.layer(bathy$Lon, bathy$Lat, bathy$Z, zlim = c(-6e3, -20), pal = 'ocean.deep')
    add.map.bathy(map, bathy)
    redraw.map(map)
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

    bathy = get.bathy(map, res = 5)
    add.map.bathy(map, bathy,
                  bathy.levels = c(-10, -50, -100, -1000, -3000),
                  bathy.lwd = c(0.8, 0.8, 1.2, 0.8, 0.8))
}

#' @title Test Map of California
#' @author Thomas Bryce Kelly
#' @import oce
#' @import ocedata
#' @export
test.map.california = function() {
    ## test for California Current
    map = make.map(coast = "coastlineWorldFine", p = "+proj=merc",
             lat.min = 31, lat.max = 36, lon.min = -126, lon.max = -119,
             dlat = 5, dlon = 5)
    bathy = get.bathy(map, res = 20)
    add.map.bathy(map, bathy,
                  bathy.levels = c( -1000, -3000, -5000),
                  bathy.lwd = c(0.8, 0.8, 1.2))
}

