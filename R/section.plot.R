#' @title Plot Section
#' @author Thomas Bryce Kelly
#' @description Plot a section object.
#' @keywords Section Plotting
#' @export
#' @param section: section object (2D only) to be plotted.
#' @param  field Parameter name to be plotted (z-axis).
#' @param  xlim Limits of the plot, by default the limits are set to the limits of the section grid.
#' @param  tlim Limits of the plot, by default the limits are set to the limits of the section grid.
#' @param  xlab X-axis label.
#' @param  ylab Y-axis label.
#' @param  log Flag for turning on log transformation of the z-axis. If TRUE then z = log(z) is performed prior to plotting.
#' @param  base Used when log = TRUE to set the base of the log.
#' @param  zlim The zlim imposed on the z-axis, which by default is set to the range of z values in the data.
#' @param  pal The color palette used; should be modeled after get.pal().
#' @param  rev Boolean for reversing the palette colors.
#' @param  include.data Flag for whether the data used to construct the grid should be plotted using the same palette.
#' @param  mark.points Flag for marking the location of each sample used to construct the grid.
#' @param  include.pch The pch used for when mark.points = TRUE.
#' @param  include.cex Sets the point size for when either include.data or mark.points are TRUE. (Used for both).
#' @param  main The title text to be included in the top line of the plot. Defaults to the name of the field.
#' @param  col.low [unimplemented] The color used when plotting out-of-range z values on the low end. Default value of '' indicates to use the minimum value of pal(). A value of NA skips the plotting of these data. Otherwise the color given is used (e.g. col.low = 'blue').
#' @param  col.high [unimplemented] Same as col.low but for out-of-range high values.
plot.section = function(section, field = NULL,
                        xlim = NULL,
                        ylim = NULL,
                        by.lon = F,
                        by.lat = F,
                        xlab = 'x',
                        ylab = 'y',
                        log = FALSE,
                        base = 10,
                        zlim = NULL,
                        pal = 'greyscale',
                        rev = F,
                        include.data = F,
                        mark.points = F,
                        include.pch = 21,
                        include.cex = 1,
                        main = NULL,
                        col.low = '',
                        col.high = '',
                        N = 255,
                        indicate = T,
                        ...) {

  ## Set Defaults
  # Check if values need to set and set them.
  if (is.null(field)) {
    field = colnames(section$grid)[3]
    warning('No field name provided, using first gridded data: ', field)
  }

  ## Handy variables
  x = section$x
  y = section$y
  z = matrix(section$grid[,field], nrow = length(x))
  if (by.lon) {
    if (verbose) { message(' Plotting based on longitude.')}
    x = approx(section$x, section$lon, xout = x, rule = 2)$y
  }
  if (by.lat) {
    if (verbose) { message(' Plotting based on latitude.')}
    x = approx(section$x, section$lat, xout = x, rule = 2)$y
  }


  if (log) {
    z = log(z, base)
    if (is.null(zlim)) {
      zlim = range(pretty(z), na.rm = TRUE)
    }
  }

  if (is.null(zlim)) { zlim = range(pretty(z), na.rm = TRUE) }

  if (is.null(main)) { main = field }
  if (is.null(xlim)) { xlim = range(x, na.rm = T)}
  if (is.null(ylim)) { ylim = range(y, na.rm = T) }

  ## Out of Range
  # Set values to zlim if you want the out of range values plotted at the zlim values
  if (!is.na(col.low) & col.low == '') {
    z[z < zlim[1]] = zlim[1]
  }
  if (!is.na(col.high) & col.high == '') {
    z[z > zlim[2]] = zlim[2]
  }

  ## Plot iamge
  image(x = x, y = y, z = z, col = get.pal(N, pal = pal, rev = rev), ylab = ylab, xlab = xlab,
        xlim = xlim, ylim = ylim, zlim = zlim, ...)

  # Plot points that are out of range when a color is given
  if (!is.na(col.low) & col.low != '') {
    zz = z
    zz[zz >= zlim[1]] = NA
    image(x = x, y = y, z = zz, col = col.low, add = TRUE)#, pch = include.pch, cex = include.cex)
  }

  if (!is.na(col.high) & col.high != '') {
    zz = z
    zz[zz <= zlim[2]] = NA
    image(x = x, y = y, z = zz, col = col.high, add = TRUE)#, pch = include.pch, cex = include.cex)
  }

  if (mark.points) {
    points(x = section$data$x, y = section$data$y, pch = 20, cex = include.cex) ## Add black points
  }
  if (include.data) {
    col = make.pal(section$data$z[,field], pal = pal, n = N, min = zlim[1], max = zlim[2], rev = rev)
    points(x = section$data$x, y = section$data$y, pch = include.pch,
           cex = include.cex, col = col, bg = col)
  }

  ## Add Title text
  if (indicate) {
    if (log) {
      st = paste0(main, '   zlim: (', round(base^zlim[1], 3), ', ', round(base^zlim[2],3), ')')
    } else {
      st = paste0(main, '   zlim: (', round(zlim[1], 3), ', ', round(zlim[2],3), ')')
    }
    mtext(st, line = 0.25, adj = 1, cex = 0.7)
  }
  box() ## make sure plotting didn't cover bounding box
}


#' @title Get Section Bathymetry
#' @author Thomas Bryce Kelly
#' @description retreive bathymetry for a section plot
#' @keywords Bathymetry
#' @export
#' @param lon longitudes
#' @param lat latitudes
#' @param res the resolution of the bathymetry (arc minutes, e.g. 10)
get.section.bathy = function(lon, lat, res = 10) {

  ## Construct dummy map object
  map =  list(lon.min = max(min(lon) - 1, -180), lon.max = min(max(lon) + 1, 180),
              lat.min = max(min(lat) - 1, -90),  lat.max = min(max(lat) + 1, 90))
  bathy = get.bathy(map, res = res)

  oce::bilinearInterp(lon, lat, bathy$Lon, bathy$Lat, bathy$Z)
}


#' @title Add Section Bathymetry
#' @author Thomas Bryce Kelly
#' @param section a section object from build.section()
#' @param bathy a bathymetric object such as "bathy.global" (the default)
#' @param binning a smoothing paramter; must be an odd integer
#' @param bathy.col the color of the bathymetry
#' @export
add.section.bathy = function(section, bathy = bathy.global, binning = 1, bathy.col = 'darkgrey') {
  if ((length(section$lon) == 1 & is.na(section(lon))) | (length(section$lat) == 1 & is.na(section(lat)))) {
    message(' Bathymetry can only be drawn for sections where lon/lat are provided.')
    return()
  }

  ## Project lon/lat points
  p = make.proj(lat = median(section$lat), lon = median(section$lon))
  bathy$xy = rgdal::project(cbind(bathy$Lon, bathy$Lat), proj = p)
  section$xy = rgdal::project(cbind(section$lon, section$lat), proj = p)

  ## Calculate depth via bilinear interpolation
  depth = interp.bilinear(x = section$xy[,1], y = section$xy[,2], gx = bathy$xy[,1], gy = bathy$xy[,2], z = -bathy$Z)

  ## Filter
  depth = runmed(depth, binning)

  ## Draw polygon
  polygon(x = c(section$x, rev(section$x)), y = c(depth, rep(1e8, length(section$x))), col = bathy.col)
}


#' @title Add Contour to Section Plot
#' @author Thomas Bryce Kelly
#' @param section a section object from "build.section()"
#' @param field the field name to use for contouring
#' @param levels the level(s) to draw contours at
#' @param col the color of the contour line(s)
#' @param lty the lty for the line(s)
#' @param labels a boolean for whether the label(s) should be drawn
#' @param cex.lab the cex (i.e. size) of the labels
#' @export
add.section.contour = function(section, field = NULL, levels = NULL, col = 'black', lty = 1, lwd = 1, labels = NULL, cex.lab = 1) {
  if (is.null(field)) {
    field = colnames(section$grid)[3] # first interpolated field
  }
  z = matrix(section$grid[[field]], nrow = length(section$x))
  if (is.null(levels)) { levels = pretty(range(as.numeric(z), na.rm = T), n = 5) }

  contour(section$x, section$y, z, add = TRUE, levels = levels, col = col, lty = lty, lwd = lwd,
          labels = labels, labcex = cex.lab, drawlabels = !isFALSE(labels))
}


#' @title Add Map Inlay
#' @export
add.section.inlay = function(section,
                             scale = 1000,
                             p = NULL,
                             side = 1,
                             pch = 20,
                             cex = 1,
                             col = 'red',
                             land.col = 'black',
                             height = 0.1,
                             width = 0.1) {
  par.old = par()
  plt = par('plt')
  if (side == 1) { par(new = T, bg = 'white', plt = c(plt[2] - c(0.05 + width, 0.05), plt[3] + c(0.05, 0.05 + height))) }
  if (side == 2) { par(new = T, bg = 'white', plt = c(plt[1] + c(0.05, 0.05 + width), plt[3] + c(0.05, 0.05 + height))) }
  if (side == 3) { par(new = T, bg = 'white', plt = c(plt[1] + c(0.05, 0.05 + width), plt[4] - c(0.05 + height, 0.05))) }
  if (side == 4) { par(new = T, bg = 'white', plt = c(plt[2] - c(0.05 + width, 0.05), plt[4] - c(0.05 + height, 0.05))) }

  if (is.na(section$lat[1])) {
    lon = 0
    lat = 0
  } else {
    lon = mean(section$lon, na.rm = T)
    lat = mean(section$lat, na.rm = T)
  }

  if (is.null(p)) {
    p = make.proj(projection = 'stere',
              lat = lon,
              lon = lat)
  }

  ## Make map
  map = make.map2(coast = 'coastlineWorld',
                  lon = lon,
                  lat = lat,
                  scale = scale,
                  draw.grid = F,
                  draw.axis = F,
                  p = p,
                  land.col = land.col)

  ## Plot points if available
  if (length(section$lon) > 1 & all(!is.na(section$lat))) {
    add.map.points(section$lon, section$lat, pch = pch, cex = 0.2 * cex, col = col)
  }
  par(plt = par.old$plt)
}
