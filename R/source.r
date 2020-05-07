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
## Date Times ################
##############################

#' @title Convert Time Excel
#' @description Helper function for converting the datetime stamps from Microsoft Excel.
#' @param x The object to be converted. Typically a numeric when converting from excel to a POSIXct.
#' @param tz The timezone of the epoch used to save the numeric type, almost \emph{ALWAYS} GMT.
#' @param rev Whether or not the operation should be reversed. This is a loseless operation.
#' @author Thomas Bryce Kelly
#' @export
conv.time.excel = function (x, tz = "GMT", rev = FALSE) {
    if (rev) {
        return(as.numeric(difftime(x, as.POSIXct("1899-12-30 00:00:00", tz = 'GMT'), unit = 'days')))
    }
    as.POSIXct(x*86400, origin = "1899-12-30 00:00:00", tz = tz)
}


#' @title Convert Time Unix
#' @description A function used to convert a unix timestamp into a POSIXct object.
#' @description Typical examples of unix timestamps include any POSIX object that has been coerced into a numeric.
#' @param x The numeric value that will be converted.
#' @param tz The timezone of the epoch used for the timestamp, almost always GMT.
#' @export
conv.time.unix = function(x, tz='GMT') {
    as.POSIXct(x, origin="1970-01-01", tz=tz)
}


#' @title Convert Time Matlab
#' @param x The numeric timestamp that will be converted.
#' @export
conv.time.matlab = function (x, tz = "GMT") {
    as.POSIXct((x - 1)*86400, origin = "0000-01-01", tz = tz)
}


#' @title Make Datetime Object
#' @description A helper function to generate a datetime object
#' @param year Year (e.g. 2016)
#' @param month Month (1-12)
#' @param day Day (1-31)
#' @param hour Hour (0-23)
#' @param minute Minute (0-59)
#' @param second Second (0-59)
#' @param tz System available timezone
#' @export
make.time = function (year = NULL, month = 1, day = 1, hour = 0, minute = 0, second = 0, tz = 'GMT') {
  if (is.null(year)) {return(Sys.time())}
  as.POSIXct(paste0(year, '-', month, '-', day, ' ', hour, ':', minute, ':', second), tz = tz)
}


#' @title Get Hour
#' @description A helper function used to pull the our information from a POSIX object.
#' @param x An object of class POSIX or a vector containing POSIX.
#' @export
get.hour = function (x) {
  if (!is.POSIXct(x)) {warning('Object passed to get.hour is not a POSIX.')}
    as.POSIXlt(x)$hour
}


#' @title Get Minutes
#' @export
get.minutes = function (x) {
    as.POSIXlt(x)$min
}


#' @title Get Seconds
#' @keywords Convert Time
#' @export
get.seconds = function (x) {
    as.POSIXlt(x)$sec
}


#' @title Which Closest Time
#' @description  Find the indicies of the closest times for each entry of x
#' @export
which.closest.time = function(x, y) {
    if (length(y) > 1) {
        l = c()
        for (i in 1:length(x)) {
            l.new = which.min(as.numeric(difftime(x[i], y, units='mins'))^2)
            l = c(l, l.new)
        }
    } else {
        l = which.min(as.numeric(difftime(x, y, units='mins'))^2)
    }
    l
}


#' @title Is POSIXct
#' @author Thomas Bryce Kelly
#' @description A helper function to determine if an object or a vector of objects is POSIX.
#' @keywords Time Helper
#' @export
is.POSIXct = function(x) {inherits(x, "POSIXct")}


#' @title Which Unique
#' @author Thomas Bryce Kelly
#' @description Find the indicies of the unique entries.
#' @param x A vector of entries.
#' @export
which.unique = function(x) {
    which(!duplicated(x))
}

#' @title Is Within
#' @author Thomas Bryce Kelly
#' @description Filter a vector of entries with fall between a set of bounds (i.e. on a closed interval).
#' @param x A vector of any type with strict ordering (i.e. where '>' or '<' are valid operators).
#' @param bounds A vector of either two entries, where bounds[1] is the lower bound and bounds[2] is the upper, or a data.frame of two columns.
#' @export
is.within = function(x, bounds) {
  if (is.null(dim(bounds))) { ## works for one set of bounds and any length x
    return(x >= bounds[1] & x <= bounds[2])
  }
  if (length(x) > 1) { stop('Is.within() should be used with either multiple x values OR multiple bounds, but not multiple of both.')}
  result = rep(NA, nrow(bounds))
  for (i in 1:nrow(bounds)) {
      result = x >= bounds[i,1] & x <= bounds[i,2]
  }
  return(results)
}


#' @title Which Within
#' @author Thomas Bryce Kelly
#' @description Filter a vector of entries with fall between a set of bounds (i.e. on a closed interval).
#' @param x A vector of any type with strict ordering (i.e. where '>' or '<' are valid operators).
#' @param bounds A vector of two entries, where bounds[1] is the lower bound and bounds[2] is the upper.
#' @export
which.within = function(x, bounds) {
  which(is.within(x, bounds))
}


#' @title Is NAN (extension of R base)
#' @description A simple extension to the base is.nan() function to permit operation on data.frames.
#' @param x object to test for NAN
#' @author Christian Fender
#' @export
is.nan = function(x) {
  if (class(x) == 'data.frame') {
    do.call(cbind, lapply(x, .Primitive('is.nan')))
  } else {
    return(.Primitive('is.nan')(x))
  }
}


##############################
## Plotting ################
##############################

#' @title Get Pal
#' @author Thomas Bryce Kelly
#' @description A helper function to retrieve a vector of colors derived from a given palette generating function.
#' @keywords Colormap
#' @param n The number of colors desired.
#' @param pal A character string of the pallet generating function, defaults to pals::ocean.haline
#' @param rev A boolean used to flip the order of the colors.
#' @import pals
#' @export
get.pal = function(n = 10, pal = pals::ocean.haline, rev = FALSE) {
  library(pals)
    pal = do.call(pal, list(n = n))
    if (rev) {
        pal = rev(pal)
    }
    pal
}


#' @title Make Pal
#' @author Thomas Bryce Kelly
#' @param x the values to be color coded
#' @param n the number of distinct colors to use
#' @param min the minimum value to correspond to the first color
#' @param max the maximum value to correspond to the last color
#' @param pal the pallete to use, default is ocean.haline
#' @param rev Boolean, reverse the color pallete?
#' @param clip boolean, remove out of range values? Defaults to False
#' @export
make.pal = function(x, n = 255, min = NA, max = NA, pal = pals::ocean.haline, rev = FALSE, clip = FALSE) {
    cols = get.pal(n+1, pal = pal, rev = rev)

    if (is.na(min)) {  ## set a minimum
        min = base::min(x, na.rm=TRUE)
    }
    if (is.na(max)) { ## Set a maximum
        max = base::max(x, na.rm=TRUE)
    }

    ## Force min and max
    if (clip) {
        x[x < min] = NA
        x[x > max] = NA
    } else {
        x[x < min] = min
        x[x > max] = max
    }
    x = (x-min) * n / (max-min) ## Scale so x falls between [0,n]
    cols = cols[floor(x)+1] # return the colors
    cols[is.na(x)] = '#00000000'

    cols
}


#' @title Make Qualitative Palette
#' @author Thomas Bryce Kelly
#' @param x The categorical values to be colored
#' @param pal The palette to be used, default tol
#' @param rev Boolean, reverse the colors?
#' @export
## Make pal for categorical data
make.qual.pal = function(x, pal='tol', rev = FALSE) {
    x[is.na(x)] = 'other'
    ## Determine numeric values
    a = sapply(x, function(xx) {which(xx == unique(x))})

    cols = get.pal(n = length(unique(x)), pal = pal, rev = rev)
    cols[a]
}


#' @title Add Colorbar
#' @author Thomas Bryce Kelly
#' @param min The minimum value of the colorbar palette
#' @param max The maximum value of the colorbar palette
#' @param labels The values where labels should be included
#' @param ticks The values where tick marks should be included
#' @param pal The name of a color palette or a color palette function itself
#' @param rev Boolean if the palette color should be reversed
#' @param units A string for the zaxis label
#' @param col.high Color for above range values
#' @param col.low Color for below range values
#' @param log A boolean if the zaxis should be log tranformed
#' @param base The base for the log transformation
#' @param x.pos the position that the x-axis should be centered on (y axis if horizontal)
#' @param width The width of the colorbar
#' @param y.pos the position that the x-axis should be centered on (x axis if horizontal)
#' @param height The length of the colorbar
#' @param cex Size
#' @param cex.units The text size of the units string text
#' @param n The number of colors to be used in the colorbar
#' @param horizontal Whether the colorbar should be placed horizontally rather than vertically
#' @description Add a color bar to any graphical device such as a plot or map. The color bar can be based on any color palette function and be placed either vertically or horizontally.
#' @keywords Plotting
#' @export
add.colorbar = function(min, max, labels = NULL, ticks = NULL, pal = 'ocean.algae', rev = FALSE, units = '',
                        col.high = '', col.low = '', log = FALSE, base = 10, x.pos = 0.875, width = 0.05,
                        y.pos = 0.5, height = 0.8, cex = 1, cex.units = 1, n = 255, horizontal = FALSE) {

  ## Setup
  par.original = par('mar')
  bty.original = par('bty')
  plt.original = par('plt')

  x = 1
  y = c(0:n)
  z = matrix(y, nrow = 1, ncol = length(y))
  delta = NULL


  ## Determine axis labels and tick marks
  if(!is.null(labels)) {
    labels = labels[labels >= min & labels <= max]
    if (log) {
        delta = (log(labels, base) - log(min, base)) / (log(max, base) - log(min, base))
    } else {
        delta = (labels - min)/(max - min)
    }
  }

  if (!is.null(ticks)) {
    ticks = ticks[ticks >= min & ticks <= max]
    if (log) {
        ticks.delta = (log(ticks, base) - log(min, base)) / (log(max, base) - log(min, base))
    } else {
        ticks.delta = (ticks - min)/(max - min)
    }
  }

  ## Now we get to work actually doing stuff:

  if (horizontal) {
    par(new = TRUE, bty = 'n', plt = c(y.pos - height/2, y.pos + height/2, x.pos - width/2, x.pos + width/2))
    image(y, x, t(z), col = get.pal(n = length(y), pal = pal, rev = rev), xlim = c(-0.1 * n, 1.1 * n), zlim = range(z),
          yaxt = 'n', xaxt = 'n', ylab = NA, xlab = NA, ylim = c(0,1))
  }
  else {
    par(new = TRUE, bty = 'n', plt = c(x.pos - width/2, x.pos + width/2, y.pos - height/2, y.pos + height/2))
    image(x, y, z, col = get.pal(n = length(y), pal = pal, rev = rev), ylim = c(-0.1 * n, 1.1 * n), zlim = range(z),
          yaxt = 'n', xaxt = 'n', ylab = NA, xlab = NA, xlim = c(0,1))
  }

  if (!is.na(col.high)) {
    if (col.high == '') {
      col.high = get.pal(100, pal, rev = rev)[100]
    }
    polygon(x = c(0, 1, 0.5, 0.5), y = c(1, 1, 41/40, 41/40) * n + 0.5, col = col.high, border = NA)
  }

  if (!is.na(col.low)) {
    if (col.low == '') {
      col.low = get.pal(100, pal, rev = rev)[1]
    }
    polygon(x = c(0, 1, 0.5, 0.5), y = c(0, 0, -n/40, -n/40)-0.5, col = col.low, border = NA)
  }

  mtext(units, side = 1, line = -1, cex = cex.units)

  if (!is.null(ticks)) {
    if (horizontal) { axis(3, at = ticks.delta * (n+1) - 0.5, labels = NA, las = 1) }
    else { axis(4, at = ticks.delta * (n+1) - 0.5, labels = NA, las = 1) }
  }

  if (!is.null(labels)) {
    if (horizontal) { axis(3, at = delta * (n+1) - 0.5, labels = labels, las = 1, cex = cex) }
    else { axis(4, at = delta * (n+1) - 0.5, labels = labels, las = 1, cex = cex) }
  }

  ## Return margins to default
  par(bty = bty.original, plt = plt.original)
}


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

      ## Add whiskers
      lines(x = rep(xx, 2), y = c(q1, q1 - 1.5 * iqr), col = lcol, lwd = lwd)
      lines(x = rep(xx, 2), y = c(q3, q3 + 1.5 * iqr), col = lcol, lwd = lwd)
      lines(x = c(xx - xlwd/2, xx + xlwd/2), y = rep(q1 - 1.5 * iqr, 2), col = lcol, lwd = lwd)
      lines(x = c(xx - xlwd/2, xx + xlwd/2), y = rep(q3 + 1.5 * iqr, 2), col = lcol, lwd = lwd)

      ## Add outliers
      k = which(y[l] < q1 - 1.5 * iqr | y[l] > q3 + 1.5 * iqr)
      if (length(k) > 0) { points(x = rep(xx, length(k)), y = y[l[k]], pch = pch, col = pcol, cex = cex) }
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

##############################
## Statistics ################
##############################

#' @title Moving Average
#' @author Thomas Bryce Kelly
#' @description Calculate the moving average from a series of observations. \emph{NB}: It assumes equally spaced observations
#' @keywords Statistics
#' @param x A vector containing equally spaced observations
#' @param n The number of samples to include in the moving average.
#' @export
ma = function(x, n = 5){
    filter(x, rep(1/n, n), sides=2)
}

#' @title Jackknife
#' @author Thomas Bryce Kelly
#' @description Generate a bootstrapped model from observations with unknown uncertainty.
#' @keywords Statistics
#' @param x The observed x values
#' @param y The observed y values
#' @param n The number of jackknife models to generate.
#' @export
jackknife = function(x, y, n = 1000) {

    k = c(which(is.na(x)), which(is.na(y)))
    if (length(k) > 0) {
        k = unique(k)
        x = x[-k]
        y = y[-k]
    }

    res = data.frame(m=0, b=0)

    for (i in 1:n) {
        l = sample(1:length(x), replace = TRUE)
        temp.model = lm(y[l] ~ x[l])
        res = rbind(res, c(coef(temp.model)[2], coef(temp.model)[1]))
    }
    res= res[-1,]
    res
}

#' @title Bootstrap
#' @author Thomas Bryce Kelly
#' @description Generate a bootstrapped model from observations with known uncertainty.
#' @keywords Statistics
#' @param x The observed x values
#' @param s.x The uncertainty in each x observation
#' @param y The observed y values
#' @param s.y The uncertainty in each y observation
#' @param n The number of bootstrapped models to generate.
#' @export
bootstrap = function(x, s.x, y, s.y, n = 1000) {

    k = c(which(is.na(x)), which(is.na(s.x)), which(is.na(y)), which(is.na(s.y)))
    if (length(k) > 0) {
        k = unique(k)
        x = x[-k]
        s.x = s.x[-k]
        y = y[-k]
        s.y = s.y[-k]
    }

    result = data.frame(m = 0, b = 0)
    num.x = length(x)

    for (i in 1:n) {
        l = sample(x = c(1:num.x), replace = TRUE)

        ## Generate random sampling
        temp.x = rnorm(num.x, x[l], s.x[l])
        temp.y = rnorm(num.x, y[l], s.y[l])

        model = lm(temp.y ~ temp.x)

        result = rbind(result, rev(coefficients(model)))
    }
    result = result[-1,]
    result
}

#' @title Bootstrap 2
#' @author Thomas Bryce Kelly
#' @keywords Statistics
#' @import lmodel2
#' @export
bootstrap.2 = function(x, s.x, y, s.y, n = 1000) {

    k = c(which(is.na(x)), which(is.na(s.x)), which(is.na(y)), which(is.na(s.y)))
    if (length(k) > 0) {
        k = unique(k)
        x = x[-k]
        s.x = s.x[-k]
        y = y[-k]
        s.y = s.y[-k]
    }

    result = data.frame(m = 0, b = 0)
    num.x = length(x)

    for (i in 1:n) {
        l = sample(x = c(1:num.x), size = num.x, replace = TRUE)

        ## Generate random sampling
        temp.x = rnorm(num.x, x[l], s.x[l])
        temp.y = rnorm(num.x, y[l], s.y[l])

        model = lmodel2::lmodel2(temp.y ~ temp.x, nperm = 0)

        result = rbind(result, as.numeric(model$regression.results[3,c(3,2)]))
    }
    result = result[-1,]
    result
}

#' @author Thomas Bryce Kelly
#' @keywords Statistics
#' @import lmodel2
#' @export
bootstrap.lmodel2.both = function(x, s.x, y, s.y, n=100) {
    mod = lmodel2::lmodel2(y ~ x)$regression.results
    results.OLS = mod[1,c(2:3)]
    results.MA = mod[2,c(2:3)]

    for (i in 1:n) {
        new.x = rnorm(length(x), x, s.x)
        new.y = rnorm(length(y), y, s.y)

        mod = lmodel2(new.y ~ new.x)$regression.results

        results.OLS = rbind(results.OLS, mod[1,c(2:3)])
        results.MA = rbind(results.MA, mod[2,c(2:3)])
    }
    list(OLS = results.OLS, MA = results.MA)
}


#### CONFIDENCE AND TRENDLINES

#' @title Add Bootstrapped Trendline
#' @author Thomas Bryce Kelly
#' @description Add the maximum likelihood trendline to a figured based on a bootstrap estimation.
#' @keywords Statistics
#' @export
add.boot.trendline = function(model, new.x, col = 'black', lty = 2, lwd=1) {
    res = c()
    for (i in 1:length(new.x)) {
        res = c(res, median(model$m * new.x[i] + model$b, na.rm=TRUE))
    }
    lines(new.x, res, col = col, lty = lty, lwd = lwd)
}



#' @author Thomas Bryce Kelly
#' @export
add.boot.trendline2 = function(res, new.x, col = 'black', lty = 2) {
    m = mean(res[,2], rm.na = TRUE)
    b = mean(res[,1], na.rm = TRUE)
    lines(new.x, new.x * m + b, col = col, lty = lty)
}

#' @title Add Linear Model Confidence Intervals
#' @author Thomas Bryce Kelly
#' @description A helper function to plot the confidence intervals determined from a base::lm model.
#' @keywords Statistics
#' @export
add.lm.conf = function(x, name, model, col = '#50505030', level = 0.95, log = FALSE) {
    if(!log) {
        dat = data.frame(a = c(1:length(x)))
        dat[[name]] = x
        pred = predict(model, interval='confidence', newdata = dat, level = level)
        polygon(c(x,rev(x)), c(pred[,"lwr"], rev(pred[,"upr"])), border = NA, col = col)
    } else {
        dat = data.frame(a = c(1:length(x)))
        dat[[name]] = x
        pred = predict(model, interval='confidence', newdata = dat, level = level)
        polygon(c(exp(x),rev(exp(x))), c(pred[,"lwr"], rev(pred[,"upr"])), border = NA, col = col)
    }
}

#' @title Add Bootstrapped Confidence Intervals
#' @author Thomas Bryce Kelly
#' @description Add confidence bands to a figure based on results of a bootstrap.
#' @keywords Statistics Uncertainty
#' @param model A bootstrap object (i.e. a dataframe) containing bootstrapped estimates of m and b.
#' @param new.x The x values for which predictions are required.
#' @param col The desired color of the confidence band. We recomend colors with transparency for plotting.
#' @param conf The quantile ranges for the plotting (default is 95% two-tail).
#' @param border Do you want a border on the shaded confidence interval?
#' @param trendline Should the maximum liklihood values be plotted as a trendline?
#' @export
add.boot.conf = function(model, new.x = c(-1e3:1e3), col = '#55555540', conf = c(0.025, 0.975),
                         border = FALSE, trendline = FALSE) {
    x.upper = new.x
    y.upper = c()
    y.lower = c()
    y.mid = c()

    for (x in x.upper) {
        y = x * model$m + model$b
        y.mid = c(y.mid, median(y))

        ## Quantiles
        y = quantile(y, conf)
        y.lower = c(y.lower, y[1])
        y.upper = c(y.upper, y[2])
    }

    polygon(x = c(x.upper, rev(x.upper)), y = c(y.upper, rev(y.lower)), col = col, border = border)

    if (trendline) {
        lines(x.upper, y.mid, lty=2)
    }
}

#' @title Get Bootstrapped Values
#' @author Thomas Bryce Kelly
#' @description Get predicted values from a series of bootstraps.
#' @keywords Statistics
#' @param model A bootstrap object (i.e. a dataframe) containing bootstrapped estimates of m and b.
#' @param x The x values for which predictions are required.
#' @param conf The confidence interval of the predicted values saught. A value of 0.5 is the maximum likelihood value.
#' @export
get.boot.vals = function(model, x, conf = 0.5) {
    y = rep(0, length(x))

    for (i in 1:length(y)) {
        yy = x[i] * model$m + model$b
        y[i] = quantile(yy, probs = conf, na.rm = TRUE)[[1]]
    }

    y
}

#' @title trendline
#' @author Thomas Bryce Kelly
#' @export
add.boot.trendline.sig = function(model, new.x, p = 0.01, sig = 0, lty = 2, lwd = 1, col = 'black') {
    if(ecdf(model$m)(sig) < p | ecdf(model$m)(sig) > 1 - p) {
        add.boot.trendline(model, new.x, col = col, lty = lty, lwd = lwd)
    }
}


#' @title lmodel boot
#' @author Thomas Bryce Kelly
#' @export
build.lmodel2.boot = function(lmodel2, n = 1000, model = 3) {
    slope = lmodel2$regression.results$Slope[model]
    slope.sd = (as.numeric(mod$confidence.intervals[model,5]) - as.numeric(mod$confidence.intervals[model,4])) / (1.96*2)

    intercept = lmodel2$regression.results$Intercept[3]
    intercept.sd = (as.numeric(mod$confidence.intervals[model,3]) - as.numeric(mod$confidence.intervals[model,2])) / (1.96*2)

    print(paste0('m = ', slope, ' +/- ', slope.sd))
    print(paste0('b = ', intercept, ' +/- ', intercept.sd))

    ## Return
    data.frame(m = rnorm(n, slope, abs(slope.sd)),
              b = rnorm(n, intercept, abs(intercept.sd)))
}

####################
#### Statistics ####
####################

#' @title Calcualte Bootstrapped R Squared
#' @author Thomas Bryce Kelly
#' @keywords Statistics
#' @export
calc.boot.r.squared = function(model, x, y) {
    ## rm NAs
    l = which(is.na(x) | is.na(y))
    if (length(l) > 0) {
        x = x[-l]
        y = y[-l]
    }

    ## R2 = 1 - SSR / SS
    1 - boot.ssr(model = model, x = x, y = y) / calc.boot.ss(y = y)
}

#' @title Calculate Bootstrapped Sum of Squared Residuals
#' @author Thomas Bryce Kelly
#' @description A function to calculate the SSR for a series of bootstrapped regressions.
#' @keywords Statistics
#' @export
#### Sum of Squared Residuals
calc.boot.ssr = function(model, x, y) {
    ## rm NAs
    l = which(is.na(x) | is.na(y))
    if (length(l) > 0) {
        x = x[-l]
        y = y[-l]
    }

    ## SSR
    y.pred = rep(0, length(x))

    for (i in 1:length(x)) {
        y.pred[i] = median(model$m) * x[i] + median(model$b)
    }

    sum((y.pred - y)^2)
}

#' @title Calculate Bootstrapped Sum of Squares
#' @author Thomas Bryce Kelly
#' @description A helper function for calculating SSR from a series of bootstraps.
#' @keywords Statistics
#' @export
#### Sum of Squares
calc.boot.ss = function(y) {
    y = y[!is.na(y)] ## Remove NAs

    ## Calc Total Sum of Squares
    sum((mean(y) - y)^2)
}


#' @title Calculate r-squared
#' @author Thomas Bryce Kelly
#' @description A function to calculate the coefficient of determination between paired data. Note: calc.r.squared(x,y) = calc.r.squared(y,x).
#' @keywords Statistics
#' @param x A set of data (e.g. observations).
#' @param y A set of data (e.g. model output).
#' @export
calc.r.squared = function(x, y){
    1 - sum((x - y)^2) / sum((x - mean(x))^2)
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

    if (length(s.x) == 1) {
        s.x = rep(s.x, length(x))
    }
    if (length(s.y) == 1) {
        s.y = rep(s.y, length(y))
    }

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
add.log.axis = function(side = 1, base = 10, col = 'black', color.minor = 'grey', grid = F) {
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
}

#' @title Is Inside Polygon
#' @description Useful to determine if a set of coordinates lies inside or outside a closed polygon. Currently works for
#'  2D coordinates such as x,y or lat,lon. Example box argument: box = data.frame(x = c(-1, 1, 1, -1), y = c(1, 1, -1, -1)).
#' @author Thomas Bryce Kelly
#' @export
is.inside = function(pos, box, verbose = FALSE) {
  d = rep(0, nrow(box))
  ## first line
  n = nrow(box)
  a = (box[1,1] - box[n,1]) * (pos[2] - box[n,2])
  b = (box[1,2] - box[n,2]) * (pos[1] - box[n,1])
  d[1] = a-b

  for (k in 2:nrow(box)) {
    a = (box[k,1] - box[k-1,1]) * (pos[2] - box[k-1,2])
    b = (box[k,2] - box[k-1,2]) * (pos[1] - box[k-1,1])
    d[k] = a-b
  }

  if (verbose) { print(d)}
  if(all(d > -1e-12)) {return(TRUE)}
  if(all(d < 1e-12)) {return(TRUE)}
  return(FALSE)
}



#' @title Update TheSource
#' @author Thomas Bryce Kelly
#' @import devtools
#' @export
TheSource.update = function() {
  user = readline(prompt = 'Update installation of TheSource? [Y|y] ')
  if (user == 'Y' | user == 'y') {
    devtools::install_github('tbrycekelly/TheSource')
  } else {
    message('User aborted update. Exiting.')
  }
}


#' @title Get version of TheSource
#'
#' @export
TheSource.version = function() {
  '0.2.7'
}




