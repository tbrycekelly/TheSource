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


