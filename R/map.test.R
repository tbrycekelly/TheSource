#' @title Make an Arctic Map
#' @author Laura Whitmore
#' @description An example map of the Arctic -- basemap only. For a full test map call test.map.arctic()
#' @inheritParams make.map
#' @export
make.map.arctic = function(coast = 'coastline2',
                           lon.min = -180,
                           lon.max = 180,
                           lat.min = 60,
                           lat.max = 90,
                           p = make.proj('stere', lat = 90),
                           dlat = 10,
                           dlon = 60,
                           land.col = 'lightgray',
                           draw.grid = TRUE) {

  make.map(coast = coast, lon.min = lon.min, lon.max = lon.max, lat.min = lat.min, lat.max = lat.max,
           p = p, land.col = land.col, draw.grid = draw.grid, dlon = dlon, dlat = dlat)
}

make.map.arctic = function(coast = 'coastline2',
                           lon = 0,
                           lat = 90,
                           scale = 2.85e3,
                           p = NULL,
                           dlat = 10,
                           dlon = 60,
                           land.col = 'lightgray',
                           draw.grid = TRUE) {
  
  
  if (is.null(p)) { 
    p = make.proj('stere', lat = 90, lon = lon)
  }
  make.map2(coast = coast, lon = lon, lat = lat, scale = scale,
           p = p, land.col = land.col, draw.grid = draw.grid, dlon = dlon, dlat = dlat)
}


#' @title Make Map (NGA)
#' @author Thomas Bryce Kelly
#' @description An example map of the NGA -- basemap only.
#' @inheritParams make.map
#' @export
make.map.nga = function(coast = 'coastline3',
                        lon = -150,
                        lat = 58,
                        scale = 300,
                        p = make.proj('stere', lat = 60, lon = -150),
                        dlat = 3,
                        dlon = 3,
                        land.col = '#333333',
                        draw.grid = TRUE) {

  make.map2(coast = coast, lon = lon, lat = lat, scale = scale,
           p = p, land.col = land.col, draw.grid = draw.grid, dlon = dlon, dlat = dlat)
}


#' @title Make Map (CCE)
#' @author Thomas Bryce Kelly
#' @description An example map of the CCE study region.
#' @inheritParams make.map
#' @import oce
#' @import marmap
#' @export
make.map.cce = function (coast = "coastline3",
                         lon.min = -126,
                         lon.max = -119,
                         lat.min = 31,
                         lat.max = 36,
                         p = make.proj(projection = 'merc', lat = 35, lon = -122, dlat = 3),
                         land.col = '#252525',
                         draw.grid = TRUE,
                         dlon = 3,
                         dlat = 3) {

  make.map(coast = coast, p = p,
           lat.min = lat.min, lat.max = lat.max, lon.min = lon.min, lon.max = lon.max,
           dlat = dlat, dlon = dlon, draw.grid = draw.grid, land.col = land.col)
}


##############
#### tests ###
##############


#' @title Test Map of the Arctic
#' @author Laura Whitmore
#' @export
test.map.arctic = function() {
  map = make.map.arctic()

  add.map.bathy(map, bathy.arctic, zlim = c(-6e3, -20), refine =-2)
  add.map.contour(map, bathy.arctic$Lon, bathy.arctic$Lat, bathy.arctic$Z, trim = F)

  map
}


#' @title Test Map of Mississippi
#' @author Laura Whitmore
#' @export
test.map.mississippi = function() {
  ## test for Mississippi Sound region
  map = make.map(coast = "coastline3",
                 lat.min = 24,
                 lat.max = 32,
                 lon.min = -93,
                 lon.max = -84,
                 dlat = 2,
                 dlon = 2)

  add.map.bathy(map, bathy.gom, zlim = c(-3e3, 10))
  add.map.contour(map,
                  bathy.gom$Lon,
                  bathy.gom$Lat,
                  bathy.gom$Z,
                  trim = F,
                  col = 'blue',
                  levels = c(-10, -50, -100, -1000, -3000),
                  lwd = c(0.8, 0.8, 1.2, 0.8, 0.8))

  ## Return
  map
}


#' @title Test Map of California
#' @author Thomas Bryce Kelly
#' @import oce
#' @import ocedata
#' @export
test.map.california = function(coast = "coastline3",
                               p = "+proj=merc",
                               lat.min = 31,
                               lat.max = 36,
                               lon.min = -126,
                               lon.max = -119,
                               dlat = 5,
                               dlon = 5,
                               draw.grid = T) {

  ## test for California Current
  map = make.map(coast = coast, p = p,
                 lat.min = lat.min, lat.max = lat.max, lon.min = lon.min, lon.max = lon.max,
                 dlat = dlat, dlon = dlon, draw.grid = draw.grid)

  add.map.bathy(map, bathy.pacific, pal = 'ocean.deep', zlim = c(-5e3, 0))
  add.map.contour(map, bathy.pacific$Lon, bathy.pacific$Lat, bathy.pacific$Z, levels = c(-1e3, -4e3), col = 'white', lwd = c(2,1))

  ## Return
  map
}
