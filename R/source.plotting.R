## Set of useful R functions for general use in plotting, analyzing and
## converting.
##
## Author: Thomas Bryce Kelly (tbk14 at fsu.edu)
## http://about.tkelly.org/
##
## Dept of Earth, Ocean & Atmospherical Sciences
## Florida State University
##
## Center for Ocean & Atmospheric Prediction Studies
## Florida State University
##
## National High Magnetic Field Laboratory
## Florida State University




##############################
## Plotting ################
##############################


#' @title Quiver Plot
#' @author Thomas Bryce Kelly
#' @param x x locations for quivers
#' @param y y locations for quivers
#' @param u The x component for the quiver arrow
#' @param v The y component for the quiver arrow
#' @param scale Used to scale the arrow length (optional)
#' @param xlim xlim
#' @param ylim ylim
#' @param col the quiver color
#' @param ... optional arguments passed to plot()
#' @export
plot.quiver = function(x, y, u, v, scale = NULL, xlim = NULL, ylim = NULL, col = 'black', ...) {
  if (is.null(xlim)) { xlim = range(x)}
  if (is.null(ylim)) { ylim = range(y)}
  if (is.null(scale)) { scale = (xlim[2] - xlim[1]) / sqrt(mean(u)^2 + mean(v)^2) / length(u) * 1.5 }

  plot(x, y, xlim = xlim, ylim = ylim, xaxs = 'i', yaxs = 'i', pch = 20, cex = 0.8, ...)
  arrows(x0 = x, x1 = x + u * scale, y0 = y, y1 = y + v * scale, length = scale/5, col = col)
}


#' @title Plot Boxplot
#' @author Thomas Bryce Kelly
#' @description Generate a blank boxplot
#' @keywords Plotting
#' @param n The width of the x axis, i.e. n = 10 corresponds to xlim = c(1,10)
#' @param ylim the y axis limit
#' @param at Values to draw x-axis labels
#' @param main The plot's main title
#' @param labels The labels to use for the x-axis
#' @param ylab The y axis label
#' @param xlab The x axis label
#' @param ... Optional arguemnts passed to plot()
#' @export
plot.boxplot = function(n = 10, ylim = c(0,1), at = NULL, main = NULL,
                        labels = NULL, ylab = NULL, xlab = NULL, ...) {

  plot(NULL, NULL, xlim = c(1-1, n+1), xaxs='i', ylim = ylim, yaxs = 'i', ylab = ylab, xlab = xlab, xaxt = 'n', main = main, ...)
  if(is.null(labels)) { labels = c(1:n)}
  if (is.null(at)) { at = c(1:length(labels)) }

  axis(1, at = at, labels = labels)
}


#' @title Add Boxplot Box
#' @author Thomas Bryce Kelly
#' @description Add a box to a boxplot that was initialized by plot.boxplot.
#' @keywords Plotting
#' @param x The x locations for the boxs
#' @param y The y values that will be plotted.
#' @export
add.boxplot.box = function(x, y, col = 'grey', border = 'black', width = 0.7, lty = 1, lcol = 'black',
                           lwd = 1, xlwd = NULL, outliers = T, pcol = 'black', pch = 1, cex = 1) {

  if (length(x) < length(y)) { x = rep(x, length(y) / length(x))}
  ## Remove NAs before they poison anything...
  l = which(!is.na(x) & !is.na(y))
  x = x[l]
  y = y[l]

  ## For each unique x value
  for (xx in unique(x)) {
    l = which(x == xx)
    if(is.null(xlwd)) { xlwd = width / 3 }

    ## statistics
    q1 = quantile(y[l], probs = 0.25)
    q3 = quantile(y[l], probs = 0.75)
    iqr = IQR(y[l], na.rm = TRUE)
    m = median(y[l], na.rm = TRUE)

    ## Box
    rect(xleft = xx - width/2, ybottom = q1, xright = xx + width/2, ytop = q3, col = col, border = border)
    lines(x = c(xx - width/2, xx + width/2), y = rep(m,2)) # Horizontal

    ## Add outliers
    k = which(y[l] < q1 - 1.5 * iqr | y[l] > q3 + 1.5 * iqr)
    if (length(k) > 0) { points(x = rep(xx, length(k)), y = y[l[k]], pch = pch, col = pcol, cex = cex) }

    ## Add whiskers
    if (length(k) > 0) {
      lines(x = rep(xx, 2), y = c(q1, min(y[l[-k]])), col = lcol, lwd = lwd)
      lines(x = rep(xx, 2), y = c(q3, max(y[l[-k]])), col = lcol, lwd = lwd)
      lines(x = c(xx - xlwd/2, xx + xlwd/2), y = rep(min(y[l[-k]]), 2), col = lcol, lwd = lwd)
      lines(x = c(xx - xlwd/2, xx + xlwd/2), y = rep(max(y[l[-k]]), 2), col = lcol, lwd = lwd)
    } else {
      lines(x = rep(xx, 2), y = c(q1, min(y[l])), col = lcol, lwd = lwd)
      lines(x = rep(xx, 2), y = c(q3, max(y[l])), col = lcol, lwd = lwd)
      lines(x = c(xx - xlwd/2, xx + xlwd/2), y = rep(min(y[l]), 2), col = lcol, lwd = lwd)
      lines(x = c(xx - xlwd/2, xx + xlwd/2), y = rep(max(y[l]), 2), col = lcol, lwd = lwd)
    }

  }
}


#' @title Add Shaded Box
#' @author Thomas Bryce Kelly
#' @description A helper function for adding shaded regions to figures (such as in timeseries).
#' @keywords Plotting
#' @export
add.shade = function(x, col = 'grey', border = NA) {
  rect(xleft = x[1], ybottom = -1e6, xright = x[2], ytop = 1e6, col = col, border = border)
}


#' @title Add Nighttime
#' @author Thomas Bryce Kelly
#' @keywords Light
#' @export
add.night = function(time, par, col = '#00000020') {
  l = diff(par < 20)

  k = which(l != 0)
  if (l[k[1]] > 0) { start = 1 } else { start = 2 }

  for (i in seq(start, length(k), by = 2)) {
    #rect(time[k[i]], -1e6, time[k[i+1]], 1e6, col = col, border = NA)
    polygon(x = c(time[k[i]], time[k[i+1]],time[k[i+1]],time[k[i]]), y = c(-1e6,-1e6,1e6,1e6), col = col, border = NA)
  }
}


#' @title Calculate point sizes
#' @author Thomas Bryce Kelly
#' @description Returns a vector of cex values for use in plotting points with size dependent on a vector of values. For example, plotting a scatter plot where point size is proportional to population size.
#' @export
#' @import dplyr
make.cex = function(x, min = 0.4, max = 4, log = FALSE, base = 10) {
  if (log) { x = log(x, base)}

  library(dplyr)
  splits <- (ntile(x, 100) - 1) / 99 # give svalues of [0 - 1]
  cex = splits * (max - min) + min

  cex
}


#' @title Add Barplot Bar (auto-stacked)
#' @author Thomas Bryce Kelly
#' @param dataa vector or numeric for the height of the bar to draw
#' @param sd the uncertainty of data (absolute units)
#' @param x the x position to plot the bar
#' @param width the width of the bars drawn
#' @param col the colors for the bars, should be same length as data
#' @param pal the color pallete used in get.pal() if col is not defined
#' @param rev.pal boolean, reverse color assignment?
#' @param border option passed to rect() controlling border appearance
#' @param angle the angle of the hash-marks if density is given
#' @param density the density of hash-marks on shape.
#' @export
add.barplot.bar = function(data, sd = NULL, x = 1, width = 0.6, col = NULL, pal = ocean.haline, rev.pal = F, border = NA,
                           angle = NULL, density = NULL) {

  if (is.null(col)) { col = get.pal(length(data), pal = pal, rev = rev.pal) }
  if (!is.null(angle) & length(angle) == 1) { angle = rep(angle, length(data))}
  if (!is.null(angle) & length(density) == 1) { density = rep(density, length(data))}

  ## Rationalize sd values
  if (is.null(sd)) { sd = rep(0, length(data)) } else { if (length(sd) == 1) { sd = rep(sd, length(data)) } }

  data = c(0, cumsum(data))
  for (i in 2:length(data)) {
    rect(xleft = x - width/2, ybottom = data[i-1],
         xright = x + width/2, ytop = data[i],
         col = col[i-1], border = border,
         density = density[i-1], angle = angle[i-1])

    ## Add SD lines
    if (sd[i-1] > 0) { lines(x = rep(x, 2), y = c(data[i] + sd[i-1], data[i] - sd[i-1])) }
  }
}


#' @title Plot Image
#' @author Thomas Bryce Kelly
#' @export
#' @inheritParams image
plot.image = function(x, y, z, xlab = NULL, ylab = NULL,
                      xlim = NULL, ylim = NULL, zlim = NULL,
                      pal = 'kovesi.rainbow', n = 255, rev = F, ...) {
  if (is.null(xlim)) {xlim = range(x)}
  if (is.null(ylim)) {ylim = range(y)}
  if (is.null(zlim)) {zlim = range(z)}

  if (is.null(xlab)) {xlab = deparse(substitute(x))}
  if (is.null(ylab)) {ylab = deparse(substitute(y))}

  xx = unique(x)
  yy = unique(y)

  if (length(z) != length(xx) * length(yy)) {
    stop('Z must have length equal to the grid generated from the unique values of x and y (i.e. z = f(x,y)).')
  }

  ## If data was not originally a matrix format.
  if (length(xx) != length(x) | length(yy) != length(y)) {
    z = as.numeric(z) ## Force into vector
    grid = expand.grid(x = xx, y = yy)
    grid$z = NA

    for (i in 1:length(z)) {
      grid$z[i] = z[x == grid$x[i] & y == grid$y[i]]
    }
    z = matrix(grid$z, nrow = length(xx), ncol = length(yy))
  }
  col = get.pal(n = n, pal = pal, rev = rev)
  image.default(x = xx[order(xx)], y = yy[order(yy)], z = z[order(xx), order(yy)], zlim = zlim, xlim = xlim, ylim = ylim,
                col = col, xlab = xlab, ylab = ylab, ...)
}


########################
#### Misc Plot Stuff ###
########################

#' @title Add Error Bars
#' @author Thomas Bryce Kelly
#' @description A function to add error bars to a figure.
#' @keywords Plotting
#' @param x The x values.
#' @param s.x The standard deviation or uncertainty of the x values.
#' @param y The y values.
#' @param s.y The standard deviation or uncertainty of the y values.
#' @param col The color of the error bars to be drawn.
#' @export
add.error.bars = function(x, s.x, y, s.y, col = 'black') {

  if (length(s.x) == 1) { s.x = rep(s.x, length(x)) }
  if (length(s.y) == 1) { s.y = rep(s.y, length(y)) }

  for (i in 1:length(x)) {
    lines(x = c(x[i], x[i]),
          y = c(y[i] + s.y[i], y[i] - s.y[i]),
          col = col)

    lines(x = c(x[i] + s.x[i], x[i] - s.x[i]),
          y = c(y[i], y[i]),
          col = col)
  }
}

#' @title Add Error Bars (log x)
#' @author Thomas Bryce Kelly
#' @description A function to add error bars given mean and stdev when the x-axis is log based.
#' @keywords Plotting Log
#' @param x The x values.
#' @param s.x The standard deviation or uncertainty of the x values.
#' @param y The y values.
#' @param s.y The standard deviation or uncertainty of the y values.
#' @param base The base of the log transformation on the x-axis.
#' @param col The color of the error bars to be drawn.
#' @export
add.error.bars.logx = function(x, s.x, y, s.y, base, col = 'black') {

  if (length(s.x) == 1) {
    s.x = rep(s.x, length(x))
  }
  if (length(s.y) == 1) {
    s.y = rep(s.y, length(y))
  }

  ## Remove NAs
  l = !is.na(x) & !is.na(y)
  x = x[l]
  s.x = s.x[l]
  y = y[l]
  s.y = s.y[l]

  for (i in 1:length(x)) {
    if(!is.na(s.y[i])) {
      lines(x = rep(log(x[i], base), 2),
            y = c(y[i] + s.y[i], y[i] - s.y[i]),
            col = col)
    }

    if(!is.na(s.x[i])) {
      lines(x = log(c(x[i] + s.x[i], x[i] - s.x[i]), base),
            y = rep(y[i], 2),
            col = col)
    }
  }
}

#' @title Add Error Bars (log y)
#' @author Thomas Bryce Kelly
#' @description A function to add error bars given mean and stdev when the y-axis is log based.
#' @keywords Plotting Log
#' @param x The x values.
#' @param s.x The standard deviation or uncertainty of the x values.
#' @param y The y values.
#' @param s.y The standard deviation or uncertainty of the y values.
#' @param base The base of the log transformation on the y-axis.
#' @param col The color of the error bars to be drawn.
#' @export
add.error.bars.logy = function(x, s.x, y, s.y, base, col = 'black') {

  if (length(s.x) == 1) {
    s.x = rep(s.x, length(x))
  }
  if (length(s.y) == 1) {
    s.y = rep(s.y, length(y))
  }

  ## Remove NAs
  l = !is.na(x) & !is.na(y)
  x = x[l]
  s.x = s.x[l]
  y = y[l]
  s.y = s.y[l]

  for (i in 1:length(x)) {
    if(!is.na(s.y[i])) {
      lines(x = rep(x[i], 2),
            y = log(c(y[i] + s.y[i], y[i] - s.y[i]), base),
            col = col)
    }

    if(!is.na(s.x[i])) {
      lines(x = c(x[i] + s.x[i], x[i] - s.x[i]),
            y = rep(log(y[i], base), 2),
            col = col)
    }
  }
}


#' @title Add Log Axis
#' @author Thomas Bryce Kelly
#' @description A helper function to add log axis to a plot.
#' @param side Side for the axis: 1 = bottom, 2 = left...
#' @param base the log bsae for the axis
#' @param col the color of the major axis labels
#' @param color.minor the color for the minor ticks and grid (if specified)
#' @param grid boolean, draw grid?
#' @export
add.log.axis = function(side = 1, base = 10, col = 'black', color.minor = 'grey', grid = F, grid.major = F) {
  k.small = c(1:(base-1))
  k = c()
  for (i in c(-10:10)) {
    k = c(k, k.small * base^(i))
  }

  w = c(2, rep(1, base-2))
  w = rep(w, length(k/length(w)))

  axis(side = side, at = log(k, base), tick = T, labels = rep('', length(k)), col = color.minor)
  axis(side = side, at = log(k[w>1], base), labels = k[w > 1], col = col)

  if (grid & (side == 1 | side == 3)) {
    abline(v = log(k, base), col = color.minor, lty = 3)
  }
  if (grid & (side == 2 | side == 4)) {
    abline(h = log(k, base), col = color.minor, lty = 3)
  }
  if (grid.major & (side == 1 | side == 3)) {
    abline(v = log(k[w>1], base), col = color.minor, lty = 3)
  }
  if (grid.major & (side == 2 | side == 4)) {
    abline(h = log(k[w>1], base), col = color.minor, lty = 3)
  }
}
