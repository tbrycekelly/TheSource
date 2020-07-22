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




#' @title Get Pal
#' @author Thomas Bryce Kelly
#' @description A helper function to retrieve a vector of colors derived from a given palette generating function.
#' @keywords Colormap
#' @param n The number of colors desired.
#' @param pal A character string of the pallet generating function, defaults to pals::ocean.haline
#' @param rev A boolean used to flip the order of the colors.
#' @import pals
#' @export
get.pal = function(n = 10, pal = NULL, rev = FALSE) {
  if (is.null(pal)) {
    message('No palette provided, default to pals::ocean.haline.')
    pal = pals::ocean.haline
  }

  ## Get pal
  pal = do.call(pal, list(n))
  if (rev) { pal = rev(pal)}

  ## Return
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
make.pal = function(x, n = 255, min = NA, max = NA, pal = NULL, rev = FALSE, clip = FALSE) {
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
  cols[is.na(x)] = '#00000000' ## NA values are transparent

  cols
}


#' @title Make Qualitative Palette
#' @author Thomas Bryce Kelly
#' @param x The categorical values to be colored
#' @param pal The palette to be used, default tol
#' @param rev Boolean, reverse the colors?
#' @export
make.qual.pal = function(x, pal = NULL, rev = FALSE) {
  if (is.null(pal)) { warning('No palette specified, defaulting to pals::tol.')}
  x[is.na(x)] = 'other.vals'
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
add.colorbar = function(min, max, labels = NULL, ticks = NULL, pal = 'ocean.haline', rev = FALSE, units = '',
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

  if (!horizontal) {mtext(units, side = 1, line = -1.5, cex = cex.units)}
  else {mtext(units, side = 1, line = -1.5, cex = cex.units)}


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


