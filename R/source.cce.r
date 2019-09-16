 #P1706.X = as.POSIXct(c('', ''), tz = 'GMT')

#### P1908
#' @export
P1908.SS1 = function(){as.POSIXct(c('2019-08-07 03:03:06', '2019-08-09 21:22:14'), tz = 'GMT')}
#' @export
P1908.SS2 = function(){as.POSIXct(c('2019-08-07 03:03:06', '2019-08-09 21:22:14'), tz = 'GMT')}
#' @export
P1908.1 = function(){as.POSIXct(c('2019-08-07 03:03:06', '2019-08-09 21:22:14'), tz = 'GMT')}
#' @export
P1908.2 = function(){as.POSIXct(c('2019-08-07 03:03:06', '2019-08-09 21:22:14'), tz = 'GMT')}
#' @export
P1908.3 = function(){as.POSIXct(c('2019-08-07 03:03:06', '2019-08-09 21:22:14'), tz = 'GMT')}
#' @export
P1908.4 = function(){as.POSIXct(c('2019-08-07 03:03:06', '2019-08-09 21:22:14'), tz = 'GMT')}
#' @export
P1908.T1 = function(){as.POSIXct(c('2019-08-07 03:03:06', '2019-08-09 21:22:14'), tz = 'GMT')}
#' @export
P1908.T2 = function(){as.POSIXct(c('2019-08-07 03:03:06', '2019-08-09 21:22:14'), tz = 'GMT')}
#' @export
P1908.T3 = function(){as.POSIXct(c('2019-08-07 03:03:06', '2019-08-09 21:22:14'), tz = 'GMT')}
#' @export
P1908.T4 = function(){as.POSIXct(c('2019-08-07 03:03:06', '2019-08-09 21:22:14'), tz = 'GMT')}
#' @export
P1908.MVP1 = function(){as.POSIXct(c('2019-08-10 14:00:06', '2019-08-11 02:10:14'), tz = 'GMT')}

#### P1706
#' @export
P1706.1 = function(){as.POSIXct(c('2017-06-09 08:22:47', '2017-06-12 12:10:57'), tz = 'GMT')}
#' @export
P1706.2 = function(){as.POSIXct(c('2017-06-13 09:00:59', '2017-06-17 11:04:04'), tz = 'GMT')}
#' @export
P1706.3 = function(){as.POSIXct(c('2017-06-19 08:13:10', '2017-06-22 12:13:13'), tz = 'GMT')}
#' @export
P1706.4 = function(){as.POSIXct(c('2017-06-24 06:11:10', '2017-06-26 12:01:23'), tz = 'GMT')}
#' @export
P1706.T1 = function(){as.POSIXct(c('2017-06-08 01:33:48', '2017-06-08 15:38:08'), tz = 'GMT')}
#' @export
P1706.T2 = function(){as.POSIXct(c('2017-06-18 02:40:14', '2017-06-18 20:16:19'), tz = 'GMT')}
#' @export
P1706.T3 = function(){as.POSIXct(c('2017-06-23 01:07:00', '2017-06-23 14:50:26'), tz = 'GMT')}
#' @export
P1706.SS1 = function(){as.POSIXct(c('2017-06-03 03:36:06', '2017-06-06 20:50:14'), tz = 'GMT')}
#' @export
P1706.SS2 = function(){as.POSIXct(c('2017-06-26 19:55:14', '2017-06-30 05:25:45'), tz = 'GMT')}

#### P1604
#' @export
P1604.1 = function(){as.POSIXct(c('2016-04-22 08:02:52', '2016-04-23 12:02:10'), tz = 'GMT')}
#' @export
P1604.2 = function(){as.POSIXct(c('2016-04-29 07:57:50', '2016-05-02 13:21:02'), tz = 'GMT')}
#' @export
P1604.3 = function(){as.POSIXct(c('2016-05-03 08:13:46', '2016-05-06 13:31:05'), tz = 'GMT')}
#' @export
P1604.4 = function(){as.POSIXct(c('2016-05-07 08:05:12', '2016-05-10 14:05:13'), tz = 'GMT')}


## Override the default mapping function when CCE is loaded

## make.map()
## default settings are for Arctic, lowest resolution coastline.
# coast = NULL: should be the name of a coastline data object. A value of NULL sets the default cosatline to 'coastlineWorld'.
#


#' @title Make Map (CCE)
#' @author Thomas Bryce Kelly
#' @description An example map of the CCE study region.
#' @keywords
#' @params
#' @import oce
#' @import marmap
#' @export
make.map.cce = function (coast = "coastlineWorldFine",
                    lon.min = -126,
                    lon.max = -119,
                    lat.min = 31,
                    lat.max = 36,
                    p = '+proj=merc',
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
