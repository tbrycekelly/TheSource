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


#' @title Make dataframe
#' @export
#' @author Thomas Bryce Kelly
data.frame = function(...) {
  x = list(...)
  a = rep(NA, length(x))

  for (i in 1:length(x)) { a[i] = length(x[[i]]) }
  k = which(a != max(a) & a != 1)
  if (length(k) > 0) { stop('Number of elements inconsistent. Column(s) ', k, ' contain ', a[k], ' entries, expecting ', max(a), ' or 1!') }

  ## Return
  base::data.frame(...)
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


####################3
######## MISC #########
######################3

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
#'
#' @export
TheSource.version = function() {
  '0.2.7'
}




