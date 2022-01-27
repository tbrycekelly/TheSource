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
#' @export
conv.time.excel = function (x, tz = "GMT", rev = FALSE) {
  if (rev) {
    return(as.numeric(difftime(x, as.POSIXct("1899-12-30 00:00:00", tz = 'GMT'), unit = 'days')))
  }
  if (any(is.POSIXct(x))) {
    return(x)
  }
  as.POSIXct(x*86400, origin = "1899-12-30 00:00:00", tz = tz)
}


#' @title Convert Time Unix
#' @description A function used to convert a unix timestamp into a POSIXct object.
#' @description Typical examples of unix timestamps include any POSIX object that has been coerced into a numeric.
#' @param x The numeric value that will be converted.
#' @param tz The timezone of the epoch used for the timestamp, almost always GMT.
#' @export
conv.time.unix = function(x, tz = 'GMT') {
  if (any(is.POSIXct(x))) { return(x) }
  as.POSIXct(x, origin="1970-01-01", tz=tz)
}


#' @title Convert Time Matlab
#' @param x The numeric timestamp that will be converted.
#' @export
conv.time.matlab = function (x, tz = "GMT") {
  if (any(is.POSIXct(x))) { return(x) }
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


#' @title Get Year
#' @description A helper function used to pull the our information from a POSIX object.
#' @param x An object of class POSIX or a vector containing POSIX.
#' @export
get.year = function (x) {
  if (any(!is.POSIXct(x))) {message(' Object(s) passed to get.year is not a POSIX.')}
  as.POSIXlt(x)$year + 1900
}


#' @title Get Month
#' @description A helper function used to pull the our information from a POSIX object.
#' @param x An object of class POSIX or a vector containing POSIX.
#' @export
get.month = function (x) {
  if (any(!is.POSIXct(x))) {message(' Object(s) passed to get.month is not a POSIX.')}
  as.numeric(format(x, '%m'))
}


#' @title Get Day
#' @description A helper function used to pull the our information from a POSIX object.
#' @param x An object of class POSIX or a vector containing POSIX.
#' @export
get.day = function (x) {
  if (any(!is.POSIXct(x))) {message('Object(s) passed to get.day is not a POSIX.')}
  as.numeric(format(x, '%d'))
}


#' @title Get Julian Day
#' @description A helper function used to pull the our information from a POSIX object.
#' @param x An object of class POSIX or a vector containing POSIX.
#' @export
get.julian = function (x) {
  if (any(!is.POSIXct(x))) {message('Object(s) passed to get.day is not a POSIX.')}
  as.numeric(difftime(x, make.time(year = get.year(x)), units = 'days'))
}


#' @title Get Hour
#' @description A helper function used to pull the our information from a POSIX object.
#' @param x An object of class POSIX or a vector containing POSIX.
#' @export
get.hour = function (x) {
  if (any(!is.POSIXct(x))) {message('Object(s) passed to get.hour is not a POSIX.')}
  as.POSIXlt(x)$hour
}


#' @title Get Minutes
#' @export
get.minutes = function (x) {
  if (any(!is.POSIXct(x))) {message('Object(s) passed to get.minutes is not a POSIX.')}
  as.POSIXlt(x)$min
}


#' @title Get Seconds
#' @keywords Convert Time
#' @export
get.seconds = function (x) {
  if (any(!is.POSIXct(x))) {message('Object(s) passed to get.seconds is not a POSIX.')}
  as.POSIXlt(x)$sec
}


#' @title System Wait
#' @description Function to have script wait a predefined time with progress bar option.
#' @author Thomas Bryce Kelly
#' @export
wait = function (sec, progress.bar = T) {
  if (class(sec) != "numeric") { stop('Seconds must be provided as a number!')}
  if (sec < 0) { stop('Delay length must be positive!')}

  if (progress.bar) {
    ## Setup
    end = as.numeric(Sys.time()) + sec
    pb = txtProgressBar(as.numeric(Sys.time()), end, style = 3)
    dt = sec / 101
    getTxtProgressBar(pb) # Show progress bar

    ## Loop until reached appropriate time.
    while (T) {
      Sys.sleep(dt)
      setTxtProgressBar(pb, as.numeric(Sys.time()))
      if (Sys.time() > end) {
        return(message()) # Return a line feed.
      }
    }
  } else {
    Sys.sleep(sec)
  }
}


#' @title Which Closest Time
#' @description  Find the indicies of the closest times for each entry of x
#' @export
which.closest.time = function(x, y) {
  if (length(y) > 1) {
    l = rep(NA, length(x))
    for (i in 1:length(x)) {
      l[i] = which.min(as.numeric(difftime(x[i], y, units='mins'))^2)
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
#' @param bounds A vector of two entries, where bounds[1] is the lower bound and bounds[2] is the upper.
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
#' @export
is.nan = function(x) {
  if (class(x) == 'data.frame') {
    do.call(cbind, lapply(x, .Primitive('is.nan')))
  } else {
    return(.Primitive('is.nan')(x))
  }
}


#' @title Trapezoidal Integration
#' @author Thomas Bryce Kelly
#' @export
integrate.trapezoid = function(x, y, xlim = NULL) {

  ## Perform checks
  if (length(x) != length(y)) { stop('Length of x and y are not equal!') }
  if (xlim[2] <= xlim[1]) { stop('xlim should be in ascending value!') }

  ## Ensure points are ordered:
  y = y[order(x, na.last = NA)]
  x = x[order(x, na.last = NA)]
  x = x[!is.na(y)]
  y = y[!is.na(y)]

  if (length(x) < 2) { stop('Integration requires 2 or more values!') }

  ## Interpolate if xlim is supplied
  if (!is.null(xlim)) {
    if (length(xlim) != 2) { stop('xlim requires two and only two values!') }

    if (xlim[1] > x[1]) {
      l = max(which(x < xlim[1]))
      y.new = approx(x[c(l, l+1)], y[c(l, l+1)], xout = xlim[1])$y
      x = x[-l]
      y = y[-l]
      x = c(xlim[1], x)
      y = c(y.new, y)
    }

    if (xlim[2] < x[length(x)]) {
      l = min(which(x > xlim[2]))
      y.new = approx(x[c(l-1, l)], y[c(l-1, l)], xout = xlim[2])$y
      x = x[-l]
      y = y[-l]
      x = c(x, xlim[2])
      y = c(y, y.new)
    }
  }

  ## Calculate integration
  dx = diff(x)
  yy = 0.5 * (y[-1] + y[-length(y)])

  sum(dx * yy)
}


#' @title Make dataframe
#' @description A helpful wrapper to the base data.frame function which provides column indexes upon mismatched lengths.
#' @export
#' @author Thomas Bryce Kelly
data.frame = function(...) {
  x = list(...)
  if (length(x) < 1) { return(base::data.frame()) }

  a = rep(NA, length(x))

  for (i in 1:length(x)) { a[i] = length(x[[i]]) }
  k = which(a != max(a) & a != 1)
  if (length(k) > 0) { stop('Number of elements inconsistent. Column(s) ', paste(k, collapse = ', '), ' contain ', paste(a[k], collapse = ', '), ' entries, expecting ', max(a), ' or 1!') }

  ## Return
  base::data.frame(...)
}


#' @title Expand Grid
#' @export
expand.grid = function(...) {
  base::expand.grid(..., KEEP.OUT.ATTRS = F, stringsAsFactors = F)
}


#' @title Is Inside Polygon
#' @description Useful to determine if a set of coordinates lies inside or outside a closed polygon. Currently works for
#'  2D coordinates such as x,y or lat,lon. Example box argument: box = data.frame(x = c(-1, 1, 1, -1), y = c(1, 1, -1, -1)).
#' @author Thomas Bryce Kelly
#' @export
is.inside = function(pos, box, verbose = FALSE) {
  if (!is.data.frame(pos)) {
    message(' `pos` does not appear to be a dataframe. Make sure it is. Attempting to fix...')
    pos = as.data.frame(matrix(pos, ncol = 2))
  }
  if (!is.data.frame(box)) {
    message(' `box` does not appear to be a dataframe. Make sure it is. Attempting to fix...')
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


#' @title Sync file directory
#' @param source String leading to the folder to be synced
#' @param destination Destination where files and folders should be copied.
#' @param pattern File/folder name value to be matched against possible files and folders. Default is '*'
#' @param verbose Output flag
#' @param level A recursion counter for the number of sublevels indexed (used for logging)
#' @param sep The folder/folder seperator, default is
#' @export
sync.dir = function(source, destination, pattern = '*', verbose = T, level = 0, sep = '/') {
  source.files = list.files(source, pattern = pattern)
  destination.files = list.files(destination, pattern = pattern)
  tab = paste0(rep(' ', level), collapse = '')

  a = Sys.time()

  ## For each source file
  if (length(source.files) > 0) {
    for (i in 1:length(source.files)) {
      if (verbose) {message(tab, Sys.time(), ': Considering file ', source.files[i])}

      if (!is.na(file.size(paste0(source, sep, source.files[i])))) {
        ## file exists
        if (file.exists(paste0(source, sep, source.files[i])) & file.exists(paste0(destination, sep, source.files[i]))) {
          ## is a directory
          if (file.size(paste0(source, sep, source.files[i])) == 0) {
            if (verbose) { message(tab, Sys.time(), ': File ', i, ' is a folder, recursing...') }
            ## Recursion
            sync.dir(source = paste0(source, sep, source.files[i]), destination = paste0(destination, sep, source.files[i]), level = level + 1, verbose = verbose)
          } else {
            if (file.size(paste0(source, sep, source.files[i])) != file.size(paste0(destination, sep, source.files[i]))) {
              ## Copy file if they are not the same size.
              if (verbose) { message(tab, Sys.time(), ': File ', i, ' was corrupted, copying from source.') }
              file.copy(paste0(source, sep, source.files[i]), paste0(destination, sep, source.files[i]))
            } else {
              if (verbose) { message(tab, Sys.time(), ': File ', i, ' already exists.') }
            }
          }
        } else { ## file or dir does not exist
          if (file.size(paste0(source, sep, source.files[i])) == 0) {
            ## Create directory and recurse
            if (verbose) { message(tab, Sys.time(), ': File ', i, ' is a folder and missing from destination, recursing...') }
            dir.create(paste0(destination, sep, source.files[i]))
            sync.dir(source = paste0(source, sep, source.files[i]), destination = paste0(destination, sep, source.files[i]), level = level + 1, verbose = verbose)
          } else {
            if (verbose) { message(tab, Sys.time(), ': File ', i, ' was copied to destination.') }
            file.copy(paste0(source, sep, source.files[i]), paste0(destination, sep, source.files[i]))
          }
        }
      }
    }
  }

  b = Sys.time()
  if (verbose & level == 0) {message('Sync took ', round((as.numeric(b) - as.numeric(a))*100)/100, ' seconds.')}
}


#' @title Pad a Number with Leading Zeros
#' @author Thomas Bryce Kelly
#' @export
pad.number = function(x, pad = 4) {
  out = rep(NA, length(x))

  for (i in length(x)) {
    dn = max(pad - nchar(x[i]), 0)
    out[i] = paste0(paste0(rep(0, dn), collapse = ''), x[i])
  }

  ## Return
  out
}


#' @title Update TheSource
#' @author Thomas Bryce Kelly
#' @import devtools
#' @export
TheSource.update = function() {
  user = readline(prompt = 'Update installation of TheSource? [Y|y] ')
  if (user == 'Y' | user == 'y') {
    unloadNamespace('package:TheSource')
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
  packageVersion('TheSource')
}



#### Testing space


make.broken.axis = function(side = 1, true.range = c(0, 6), left.range = c(0,1), right.range = c(5,10)) {
  if (length(left.range) != 2 | length(right.range) != 2 | length(true.range) != 2) { stop('Ranges must be length 2!')}
  if (left.range[2] < left.range[1]) { decreasing = T} else {decreasing = F}
  if (right.range[2] < right.range[1] & !decreasing) { stop('Both ranges much be decreasing or increasing!')}
  if (right.range[1] < right.range[2] & decreasing) { stop('Both ranges much be decreasing or increasing!')}

  list(side = side,
       true.range = true.range,
       left.range = left.range,
       right.range = right.range,
       f.left = diff(left.range) / (diff(left.range) + diff(right.range)),
       f.right = diff(right.range) / (diff(left.range) + diff(right.range)),
       decreasing = decreasing)
}


add.broken.axis = function(broken, ...) {
  pretty.left = pretty(broken$left.range)
  pretty.right = pretty(broken$right.range)
  if (broken$decreasing) {
    pretty.left = pretty.left[pretty.left > broken$left.range[2]]
    pretty.right = pretty.right[pretty.right < broken$right.range[1]]
  } else {
    pretty.left = pretty.left[pretty.left < broken$left.range[2]]
    pretty.right = pretty.right[pretty.right > broken$right.range[1]]
  }

  width = diff(broken$true.range)
  axis(side = broken$side, labels = pretty.left, at = width * broken$f.left * (pretty.left - broken$left.range[1])/(diff(broken$left.range)), ...)
  axis(side = broken$side, labels = pretty.right, at = width * broken$f.left + broken$f.right * (pretty.right - broken$right.range[1])/(diff(broken$right.range)), ...)

}




