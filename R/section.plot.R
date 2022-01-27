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
plot.section = function(section, field = NULL, xlim = NULL, ylim = NULL, xlab = 'x', ylab = 'y', log = FALSE, base = 10,
                        zlim = NULL, pal = 'greyscale', rev = FALSE, include.data = FALSE, mark.points = FALSE, include.pch = 21,
                        include.cex = 1, main = NULL, col.low = '', col.high = '', N = 255) {

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

  if (log) {
    z = log(z, base)
    if (is.null(zlim)) {
      zlim = range(pretty(z, na.rm = TRUE))
    }
  }

  if (is.null(zlim)) { zlim = range(pretty(z, na.rm = TRUE)) }

  if (is.null(main)) { main = field }
  if (is.null(xlim)) { xlim = range(x)}
  if (is.null(ylim)) { ylim = range(y) }

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
        xlim = xlim, ylim = ylim, zlim = zlim)

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
    col = make.pal(section$data[,field], pal = pal, n = N, min = zlim[1], max = zlim[2], rev = rev)
    points(x = section$data$x, y = section$data$y, pch = include.pch,
           cex = include.cex, col = col, bg = col)
  }

  ## Add Title text
  if (log) {
    st = paste0(main, '   zlim: (', round(base^zlim[1], 3), ', ', round(base^zlim[2],3), ')')
  } else {
    st = paste0(main, '   zlim: (', round(zlim[1], 3), ', ', round(zlim[2],3), ')')
  }
  mtext(st, line = 0.25, adj = 1, cex = 0.7)

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
  x = section$x
  depth = rep(NA, length(x))

  bathy$lon = bathy$lon %% 360
  section$section.lon = section$section.lon %% 360

  for (i in 1:length(x)) {
    depth[i] = -bathy$Z[which.min((bathy$Lon - section$section.lon[i])^2), which.min((bathy$Lat - section$section.lat[i])^2)]
  }
  depth = runmed(depth, binning)

  polygon(x = c(x, rev(x)), y = c(depth, rep(1e8, length(x))), col = bathy.col)
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
  if (is.null(levels)) { levels = pretty(range(as.numeric(z)), n = 5) }

  contour(section$x, section$y, z, add = TRUE, levels = levels, col = col, lty = lty, lwd = lwd,
          labels = labels, labcex = cex.lab, drawlabels = !isFALSE(labels))
}


#' @title Add Map Inlay
#' @export
add.section.inlay = function(section) {
  par.old = par()
  par(new = T, plt = c(0.2, 0.33, 0.2, 0.33))
  map = make.map(coast = 'coastlineWorld', lon.min = min(section$section.lon), lon.max = max(section$section.lon),
                 lat.min = min(section$section.lat), lat.max = max(section$section.lat),
                 dlon = 360, dlat = 360, draw.grid = F,
                 p = make.proj(8, lat = mean(section$section.lat), lon = mean(section$section.lon)))
  add.map.points(section$section.lon, section$section.lat, pch = 20, cex = 0.2, col = 'red')
  par(plt = par.old$plt)
}
