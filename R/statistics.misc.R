## Useful Statistically Oriented Functions
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


#' @title Is Prime?
#' @author Thomas Bryce Kelly
#' @param n a numeric or vector of numerics to test if values are prime
#' @export
is.prime = function(n) {
  p = rep(TRUE, length(n))
  for (i in 1:length(n)) {
    ## Divisible by 2?
    if (n[i] %% 2 == 0) {
      p[i] = FALSE
    } else {
      if (n[i] == 3) {
        p[i] = TRUE
      } else {
        ## Divisible by odd between 3 and sqrt(n)?
        for (k in seq(3, pmax(3, ceiling(sqrt(n[i]))), by = 2)) {
          if (n[i] %% k == 0) {
            p[i] = FALSE
          }
        }
      }
    }
  }
  p
}


#' @title Calculate Simple Moving Average
#' @author Laura Whitmore
#' @param parameter signal or vector to be smoothed
#' @param k the integer number of entries to average over
#' @export
calc.sma = function(parameter, k = 100) {
  sma = rep(NA, length(parameter)) ## need this vector to be the same length as parameter
  for (i in 1:length(parameter)) {
    if (i < k/2 | i > length(parameter) - k/2) {
      sma[i] = NA
    } else {
      series = c((i - floor(k/2) + 1):(i + floor(k/2) - 1))
      #message(k, ' = k\t', length(series), ' = series')
      sma[i] = mean(parameter[series], na.rm = TRUE)
    }
  }

  sma
}


#' @title Moving Average
#' @author Thomas Bryce Kelly
#' @description Calculate the moving average from a series of observations. \emph{NB}: It assumes equally spaced observations
#' @keywords Statistics
#' @param x A vector containing equally spaced observations
#' @param n The number of samples to include in the moving average.
#' @export
ma = function(x, n = 5){
  as.numeric(filter(as.numeric(x), rep(1/n, n), sides = 2))
}


#### CONFIDENCE AND TRENDLINES

#' @title Get Bootstrapped Values
#' @author Thomas Bryce Kelly
#' @description Get predicted values from a series of bootstraps.
#' @keywords Statistics
#' @param model A bootstrap object (i.e. a dataframe) containing bootstrapped estimates of m and b.
#' @param x The x values for which predictions are required.
#' @param conf The confidence interval of the predicted values saught. A value of 0.5 is the maximum likelihood value.
#' @export
get.regress.vals = function(model, x, conf = 0.5) {
  y = rep(0, length(x))

  for (i in 1:length(y)) {
    yy = x[i] * model$m + model$b
    y[i] = quantile(yy, probs = conf, na.rm = TRUE)[[1]]
  }

  y
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



#' @title Get York MSWD
get.york.mswd = function(x, a, b){
  xy = get.york.xy(x, a, b)
  X2 = 0
  ns = nrow(x)

  for (i in 1:ns){
    E = cor2cov2(x[i,'sX'],x[i,'sY'],x[i,'rXY'])
    X = matrix(c(x[i,'X'] - xy[i,1], x[i,'Y'] - xy[i,2]), 1, 2)
    if (!any(is.na(X)))
      X2 = X2 + X %*% solve(E) %*% t(X)
  }

  ## Return
  out = list()
  out$df = ns - 2
  out$mswd = as.numeric(X2 / out$df)
  out$p.value = as.numeric(1-stats::pchisq(X2,out$df))
  out
}


#' @title Get York XY
# get fitted X and X given a dataset x=cbind(X,sX,Y,sY,rXY),
# an intercept a and slope b. This function is useful
# for evaluating log-likelihoods of derived quantities
get.york.xy = function(x, a, b){

  ## Determine weights
  wX = 1 / x$sX^2
  wY = 1 / x$sY^2

  ## Transform
  A = sqrt(wX * wY)
  W = wX * wY / (wX + b^2 * wY - 2 * b * x$rXY * A)
  Xbar = sum(W * x$X, na.rm=TRUE) / sum(W, na.rm=TRUE)
  Ybar = a + b * Xbar

  ## Calculate
  U = x$X - Xbar
  V = x$Y - Ybar
  B = W * (U / wY + b * V / wX - (b * U + V) * x$rXY / A)

  ## Return
  cbind(Xbar + B, Ybar + b * B)
}


#' @title Cor to Cov
cor2cov2 = function(sX,sY,rXY){
  covmat = matrix(rXY*sX*sY, ncol = 2, nrow = 2)
  covmat[1,1] = sX^2
  covmat[2,2] = sY^2

  ## Return
  covmat
}


#' @title Impute Missing Values
#' @param x the data frame with missing values
#' @param method the imputation method: 'mean', 'median', 'forest
#' @import Hmisc
#' @import missForest
build.impute = function(x, method = 'forest') {
  x = data.frame(x)

  if (method %in% c('mean', 'median', 'min', 'max')) {
    for (i in 1:ncol(x)) {
      x[,i] = Hmisc::impute(x[,i], eval(parse(text = method)))
    }
  } else {
    result = missForest::missForest(x)
    message('Out of bag NRMSE: ', result$OOBerror[1])
    x = result$ximp
  }

  x
}

