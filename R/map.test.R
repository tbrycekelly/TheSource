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
