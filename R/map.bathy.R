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
