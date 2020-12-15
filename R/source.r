## Set of useful R functions for general use in plotting, analyzing and
## converting.
##
## Author: Thomas Bryce Kelly (tbkelly at alaska.edu)
## http://about.tkelly.org/
##
## College of Fisheries and Ocean Science
## University of Alaska Fairbanks
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
#' @examples
#' conv.time.excel(42550)
#' conv.time.excel(conv.time.excel(42550), rev = T)
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
#' @example
#' conv.time.unix(as.numeric(Sys.time()))
#' @export
conv.time.unix = function(x, tz = 'GMT') {
    as.POSIXct(x, origin="1970-01-01", tz=tz)
}


#' @title Convert Time Matlab
#' @param x The numeric timestamp that will be converted.
#' @example
#' conv.time.matlab(720000)
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
#' @example
#' make.time(2020, 4, 28)
#' @export
make.time = function (year = NULL, month = 1, day = 1, hour = 0, minute = 0, second = 0, tz = 'GMT') {
  if (is.null(year)) {return(Sys.time())}
  as.POSIXct(paste0(year, '-', month, '-', day, ' ', hour, ':', minute, ':', second), tz = tz)
}

#' @title Get Year
#' @description A helper function used to pull the our information from a POSIX object.
#' @param x An object of class POSIX or a vector containing POSIX.
#' @export
get.year = function (x) {
  if (!is.POSIXct(x)) {warning('Object passed to get.year is not a POSIX.')}
  as.POSIXlt(x)$year + 1900
}

#' @title Get Month
#' @description A helper function used to pull the our information from a POSIX object.
#' @param x An object of class POSIX or a vector containing POSIX.
#' @export
get.month = function (x) {
  if (!is.POSIXct(x)) {warning('Object passed to get.month is not a POSIX.')}
  as.numeric(format(x, '%m'))
}

#' @title Get Day
#' @description A helper function used to pull the our information from a POSIX object.
#' @param x An object of class POSIX or a vector containing POSIX.
#' @export
get.day = function (x) {
  if (!is.POSIXct(x)) {warning('Object passed to get.day is not a POSIX.')}
  as.numeric(format(x, '%d'))
}

#' @title Get Julian Day
#' @description A helper function used to pull the our information from a POSIX object.
#' @param x An object of class POSIX or a vector containing POSIX.
#' @export
get.julian = function (x) {
  if (!is.POSIXct(x)) {warning('Object passed to get.day is not a POSIX.')}
  as.numeric(difftime(x, make.time(year = get.year(x)), units = 'days'))
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
#' @example
#' a = Sys.time()
#' b = seq(make.time(2010), make.time(2050), by = '1 year')
#' which.closest.time(a, b)
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
#' @example
#' is.POSIXct('a') # FALSE
#' @export
is.POSIXct = function(x) {inherits(x, "POSIXct")}


#' @title Which Unique
#' @author Thomas Bryce Kelly
#' @description Find the indicies of the unique entries.
#' @param x A vector of entries.
#' @example
#' which.unique(c('a', 1:3, 1))
#' @export
which.unique = function(x) {
    which(!duplicated(x))
}

#' @title Is Within
#' @author Thomas Bryce Kelly
#' @description Filter a vector of entries with fall between a set of bounds (i.e. on a closed interval).
#' @param x A vector of any type with strict ordering (i.e. where '>' or '<' are valid operators).
#' @param bounds A vector of two entries, where bounds[1] is the lower bound and bounds[2] is the upper.
#' @example
#' is.within(5, c(0,10))
#' @export
is.within = function(x, bounds) {
  if (length(bounds) != 2) { stop('Bounds requires two values, an upper and a lower value.')}
  sapply(x, function(x) {x <= bounds[2] & x >= bounds[1]})
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
#' @example
#' is.nan(5)
#' @export
is.nan = function(x) {
  if (class(x) == 'data.frame') {
    do.call(cbind, lapply(x, .Primitive('is.nan')))
  } else {
    return(.Primitive('is.nan')(x))
  }
}


#' @title Make dataframe
#' @description A helpful wrapper to the base data.frame function which provides column indexes upon mismatched lengths.
#' @export
#' @author Thomas Bryce Kelly
#' @example
#' data.frame(a = c(1:5), b = c('a':'e'))
data.frame = function(...) {
  x = list(...)
  a = rep(NA, length(x))

  for (i in 1:length(x)) { a[i] = length(x[[i]]) }
  k = which(a != max(a) & a != 1)
  if (length(k) > 0) { stop('Number of elements inconsistent. Column(s) ', paste(k, collapse = ', '), ' contain ', paste(a[k], collapse = ', '), ' entries, expecting ', max(a), ' or 1!') }

  ## Return
  base::data.frame(...)
}


#' @title Is Inside Polygon
#' @description Useful to determine if a set of coordinates lies inside or outside a closed polygon. Currently works for
#'  2D coordinates such as x,y or lat,lon. Example box argument: box = data.frame(x = c(-1, 1, 1, -1), y = c(1, 1, -1, -1)).
#' @author Thomas Bryce Kelly
#' @export
is.inside = function(pos, box, verbose = FALSE) {
  if (!is.data.frame(pos)) {
    warning('`pos` does not appear to be a dataframe. Make sure it is. Attempting to fix...')
    pos = as.data.frame(matrix(pos, ncol = 2))
  }
  if (!is.data.frame(box)) {
    warning('`box` does not appear to be a dataframe. Make sure it is. Attempting to fix...')
    box = as.data.frame(matrix(box, ncol = 2))
  }

  ans = rep(F, nrow(pos))

  if (verbose) { message('Results:\tPoint \t\tInside?\tCount')}
  for (i in 1:nrow(pos)) {

    ## Transform so that test point is at origin:
    x.poly = box[,1] - pos[i,1]
    y.poly = box[,2] - pos[i,2]
    count = 0

    ## Slopes of all line segments
    ## m[1] = (y2 - y1) / (x2 - x1)
    m = c(diff(y.poly), y.poly[length(y.poly)] - y.poly[1]) / c(diff(x.poly), x.poly[length(x.poly)] - x.poly[1])

    ## x.intercept[1] = -y[1] / m[1] + x[1]
    x.intercept = -y.poly / m + x.poly

    for (j in 2:length(x.poly)) {
      ## if max(y.poly[2], y.poly[1]) >= 0
      if (max(y.poly[c(j,j-1)]) > 0 & ## one point is at or above x axis
          min(y.poly[c(j,j-1)]) <= 0 & ## one point is at or below x axis
          min(y.poly[c(j,j-1)]) != max(y.poly[c(j,j-1)]) & ## But both points are not on the x axis
          !is.na(x.intercept[j-1]) & ## and the x intercept is not NA (removes parallel lines)
          x.intercept[j-1] >= 0) { ## but the x intercept is to the right of the origin (only count in one direction)
        count = count + 1
      }
    }
    ## Close figure
    if (max(y.poly[c(1,length(x.poly))]) > 0 &
        min(y.poly[c(1,length(x.poly))]) <= 0 &
        max(y.poly[c(1,length(x.poly))]) != min(y.poly[c(1,length(x.poly))]) &
        !is.na(x.intercept[length(x.poly)]) &
        x.intercept[length(x.poly)] >= 0) {
      count = count + 1
    }

    ans[i] = count %% 2 == 1
    if (any(x.poly == 0 & y.poly == 0)) {ans[i] = 1}
    if (verbose) { message('\t\t(', format(pos[i,1], digits = 2), ',', format(pos[i,2], digits = 2), ')\t\t', ans[i], '\t', count )}
  }
  if (verbose) { message(); message('Number of points inside:\t', length(which(ans)), ' (', format(length(which(ans)) / length(ans), digits = 3),')')}
  ans
}



#' @title Update TheSource
#' @author Thomas Bryce Kelly
#' @import devtools
#' @export
TheSource.update = function() {
  user = readline(prompt = 'Update installation of TheSource? [Y|y] ')
  if (user == 'Y' | user == 'y') {
    detach("package:TheSource", unload=TRUE)
    devtools::install_github('tbrycekelly/TheSource')
    library(TheSource)
  } else {
    message('User aborted update. Exiting.')
  }
}


#' @title Get version of TheSource
#' @author Thomas Bryce Kelly
#' @export
TheSource.version = function() {
  '0.2.8'
}




