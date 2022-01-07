#' @title Jackknife
#' @author Thomas Bryce Kelly
#' @description Generate a bootstrapped model from observations with unknown uncertainty.
#' @keywords Statistics
#' @param x The observed x values
#' @param y The observed y values
#' @param n The number of jackknife models to generate.
#' @export
regress.jackknife = function(x, y, sx = NA, sy = NA) {
  if (any(!is.na(sx)) | any(!is.na(sy))) {
    warning('Jackknife regressions do not include datapoint uncertainty, ignoring sx and sy values.')
  }

  ## Filter out NAs
  k = which(x = is.na(x) | is.na(y))
  if (length(k) > 0) {
    x = x[-k]
    y = y[-k]
  }

  ## Build set of regressions
  n = length(x)
  res = data.frame(m = rep(NA, n), b = NA, ssr = NA, CV = NA)

  for (i in 1:n) {
    temp.model = lm(y[-i] ~ x[-i])
    res$m[i] = coef(temp.model)[2]
    res$b[i] = coef(temp.model)[1]
    res$ssr[i] = sum((res$m[i] * x + res$b[i] - y)^2)
    res$CV[i] = (res$m[i] * x[i] + res$b[i] - y[i])^2 ## SSR
  }

  ## Calculate residuals and z-score
  residuals = median(res$m) * x + median(res$b) - y
  zscore = sapply(c(1:length(x)), function(i) { (y[i] - mean(res$m * x[i] + res$b)) / sd(res$m * x[i] + res$b) })
  aic = length(x) * log(sum(residuals^2) / length(x)) + 6 # Akaike Information Criterion

  list(models = res,
       data = data.frame(x = x,
                         y = y,
                         y.pred = y + residuals,
                         sx = sx,
                         sy = sy,
                         residuals = residuals,
                         zscore = zscore),
       meta = list(model = regress.jackknife,
                   m = median(res$m),
                   b = median(res$b),
                   n = n, AIC = aic,
                   R2 = 1 - sum(residuals^2) / sum((y - mean(y))^2),
                   Source.version = packageVersion('TheSource'),
                   R.version = R.version.string))
}


#' @title Resample
#' @author Thomas Bryce Kelly
#' @description Generate a bootstrapped model from observations with unknown uncertainty.
#' @keywords Statistics
#' @param x The observed x values
#' @param y The observed y values
#' @param n The number of jackknife models to generate.
#' @export
regress.resample = function(x, y, sx = NA, sy = NA, n = 1e3) {
  if (any(!is.na(sx)) | any(!is.na(sy))) {
    warning('Resample regressions do not include datapoint uncertainty, ignoring sx and sy values.')
  }

  ## Filter out NAs
  k = which(x = is.na(x) | is.na(y))
  if (length(k) > 0) {
    x = x[-k]
    y = y[-k]
  }

  ## Build set of regressions
  res = data.frame(m = rep(NA, n), b = NA, ssr = NA, CV = NA)

  for (i in 1:n) {
    l = sample(1:length(x), replace = T)
    temp.model = lm(y[l] ~ x[l])
    res$m[i] = coef(temp.model)[2]
    res$b[i] = coef(temp.model)[1]
    res$ssr[i] = sum((res$m[i] * x + res$b[i] - y)^2)
    res$CV[i] = sum((res$m[i] * x[-l] + res$b[i] - y[-l])^2)
  }

  ## Calculate residuals and z-score
  residuals = median(res$m) * x + median(res$b) - y
  aic = length(x) * log(sum(residuals^2) / length(x)) + 6 # Akaike Information Criterion

  list(models = res,
       data = data.frame(x = x, y = y, y.pred = y + residuals, sx = sx, sy = sy, residuals = residuals, zscore = zscore),
       meta = list(model = regress.resample,
                   m = median(res$m),
                   b = median(res$b),
                   n = n, AIC = aic,
                   R2 = 1 - sum(residuals^2) / sum((y - mean(y))^2),
                   Source.version = packageVersion('TheSource'),
                   R.version = R.version.string))
}


#' @title Jackknife Model II
#' @author Thomas Bryce Kelly
#' @description Generate a bootstrapped model from observations with unknown uncertainty.
#' @keywords Statistics
#' @param x The observed x values
#' @param y The observed y values
#' @param n The number of jackknife models to generate.
#' @export
regress.jackknife.2 = function(x, y, sx = NA, sy = NA, n = 1000) {
  if (any(!is.na(sx)) | any(!is.na(sy))) {
    warning('Jackknife regressions do not include datapoint uncertainty, ignoring sx and sy values.')
  }

  ## Filter out NAs
  k = which(x = is.na(x) | is.na(y))
  if (length(k) > 0) {
    x = x[-k]
    y = y[-k]
  }

  ## Build set of regressions
  res = data.frame(m = rep(NA, n), b = NA, ssr = NA)

  for (i in 1:n) {
    l = sample(1:length(x), replace = TRUE)
    model = lmodel2::lmodel2(y[l] ~ x[l], nperm = 0)

    result = as.numeric(model$regression.results[3,c(3,2)])
    res$m[i] = result[1]
    res$b[i] = result[2]
    res$ssr[i] = sum((res$m[i] * x + res$b[i] - y)^2)
  }

  ## Calculate residuals and z-score
  residuals = median(res$m) * x + median(res$b) - y
  zscore = sapply(c(1:length(x)), function(i) { (y[i] - mean(res$m * x[i] + res$b)) / sd(res$m * x[i] + res$b) })
  aic = length(x) * log(sum(residuals^2) / length(x)) + 6 # Akaike Information Criterion

  list(models = res,
       data = data.frame(x = x, y = y, y.pred = y + residuals, sx = sx, sy = sy, residuals = residuals, zscore = zscore),
       meta = list(model = regress.jackknife,
                   m = median(res$m),
                   b = median(res$b),
                   n = n, AIC = aic,
                   R2 = 1 - sum(residuals^2) / sum((y - mean(y))^2),
                   Source.version = packageVersion('TheSource'),
                   R.version = R.version.string))
}


#' @title Bootstrap
#' @author Thomas Bryce Kelly
#' @description Generate a bootstrapped model from observations with known uncertainty (and x - y independence).
#' @keywords Statistics
#' @param x The observed x values
#' @param sx The uncertainty in each x observation
#' @param y The observed y values
#' @param sy The uncertainty in each y observation
#' @param n The number of bootstrapped models to generate.
#' @export
regress.bootstrap = function(x, y, sx, sy, n = 1000) {

  ## Filter
  k = which(is.na(x) | is.na(sx) | is.na(y) | is.na(sy))
  if (length(k) > 0) {
    x = x[-k]
    s.x = s.x[-k]
    y = y[-k]
    s.y = s.y[-k]
  }

  res = data.frame(m = rep(NA, n), b = NA, ssr = NA)
  num.x = length(x)

  for (i in 1:n) {

    ## Generate random sampling and apply linear regression
    temp.x = rnorm(num.x, x, sx)
    temp.y = rnorm(num.x, y, sy)
    model = lm(temp.y ~ temp.x)

    res$m[i] = coef(temp.model)[2]
    res$b[i] = coef(temp.model)[1]
    res$ssr[i] = sum((res$m[i] * x + res$b[i] - y)^2)
  }

  ## Calculate residuals and z-score
  residuals = median(res$m) * x + median(res$b) - y
  zscore = sapply(c(1:length(x)), function(i) { (y[i] - mean(res$m * x[i] + res$b)) / sd(res$m * x[i] + res$b) })
  aic = length(x) * log(sum(residuals^2) / length(x)) + 6 # Akaike Information Criterion

  list(models = res,
       data = data.frame(x = x, y = y, y.pred = y + residuals, sx = sx, sy = sy, residuals = residuals, zscore = zscore),
       meta = list(model = regress.bootstrap,
                   m = median(res$m),
                   b = median(res$b),
                   n = n, AIC = aic,
                   R2 = 1 - sum(residuals^2) / sum((y - mean(y))^2),
                   Source.version = packageVersion('TheSource'),
                   R.version = R.version.string))
}


#' @title Bootstrap 2
#' @author Thomas Bryce Kelly
#' @keywords Statistics
#' @import lmodel2
#' @export
regress.bootstrap.type2 = function(x, sx, y, sy, n = 1000) {

  ## Filter
  k = which(is.na(x) | is.na(sx) | is.na(y) | is.na(sy))
  if (length(k) > 0) {
    x = x[-k]
    sx = sx[-k]
    y = y[-k]
    sy = sy[-k]
  }

  res = data.frame(m = rep(NA, n), b = NA, ssr = NA)
  num.x = length(x)

  for (i in 1:n) {

    ## Generate random sampling
    temp.x = rnorm(num.x, x, sx)
    temp.y = rnorm(num.x, y, sy)

    model = lmodel2::lmodel2(temp.y ~ temp.x, nperm = 0)

    result = as.numeric(model$regression.results[3,c(3,2)])
    res$m[i] = result[1]
    res$b[i] = result[2]
    res$ssr[i] = sum((res$m[i] * x + res$b[i] - y)^2)
  }

  ## Calculate residuals and z-score
  residuals = median(res$m) * x + median(res$b) - y
  zscore = sapply(c(1:length(x)), function(i) { (y[i] - mean(res$m * x[i] + res$b)) / sd(res$m * x[i] + res$b) })
  aic = length(x) * log(sum(residuals^2) / length(x)) + 6 # Akaike Information Criterion

  list(models = res,
       data = data.frame(x = x, y = y, y.pred = y + residuals, sx = sx, sy = sy, residuals = residuals, zscore = zscore),
       meta = list(model = regress.bootstrap.type2,
                   m = median(res$m),
                   b = median(res$b),
                   n = n, AIC = aic,
                   R2 = 1 - sum(residuals^2) / sum((y - mean(y))^2),
                   Source.version = packageVersion('TheSource'),
                   R.version = R.version.string))
}


#' @title York Regression
#' @author Thomas Bryce Kelly
#' @export
regress.york = function(x, sx, y, sy, alpha = 0.05){
  X = cbind(x, sx, y, sy, 0)
  X = as.data.frame(X)
  colnames(X) = c('X','sX','Y','sY','rXY')

  ab = lm(X$Y ~ X$X)$coefficients # initial guess
  a = ab[1]
  b = ab[2]

  if (any(is.na(ab))) {
    stop('Cannot fit a straight line through these data')
  }

  ## Set weigths based on variance
  wX = 1 / X$sX^2
  wY = 1 / X$sY^2

  ## Iteratively find the solution
  for (i in 1:50){ # 50 = maximum number of iterations

    bold = b
    A = sqrt(wX * wY)
    W = wX * wY / (wX + b^2 * wY - 2 * b * X$rXY * A)

    Xbar = sum(W * X$X, na.rm = TRUE) / sum(W, na.rm=TRUE)
    Ybar = sum(W * X$Y, na.rm = TRUE) / sum(W, na.rm=TRUE)

    U = X$X - Xbar
    V = X$Y - Ybar
    B = W * (U / wY + b * V / wX - (b * U + V) * X$rXY / A)
    b = sum(W * B * V, na.rm = TRUE) / sum(W * B * U, na.rm = TRUE)
    if ((bold / b - 1)^2 < 1e-15) {
      break # convergence reached
    }
  }

  ## Transform back
  a = Ybar - b * Xbar
  X = Xbar + B

  ## Calculate uncertainty
  xbar = sum(W * X, na.rm=TRUE) / sum(W, na.rm=TRUE)
  u = X - xbar
  sb = sqrt(1 / sum(W * u^2, na.rm = TRUE))
  sa = sqrt(1 / sum(W, na.rm=TRUE) + (xbar * sb)^2)

  out = get.york.mswd(x, a, b)
  out$fact = 1
  out$a = c(a,sa)
  out$b = c(b,sb)
  out$cov.ab = -xbar * sb^2
  names(out$a) = c('a','s[a]')
  names(out$b) = c('b','s[b]')
  out$type = 'york'

  ## Calculate residuals and z-score
  residuals = out$b * x + out$a - y
  zscore = 0
  aic = length(x) * log(sum(residuals^2) / length(x)) + 6 # Akaike Information Criterion

  list(models = data.frame(m = out$b, b = out$a, ssr = ),
       data = data.frame(x = x, y = y, y.pred = y + residuals, sx = sx, sy = sy, residuals = residuals, zscore = zscore),
       meta = list(model = regress.york,
                   m = out$b,
                   b = out$a,
                   n = length(x),
                   AIC = aic,
                   R2 = 1 - sum(residuals^2) / sum((y - mean(y))^2),
                   Source.version = packageVersion('TheSource'),
                   R.version = R.version.string))
}
