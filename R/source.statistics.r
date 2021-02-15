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


#' @title Moving Average
#' @author Thomas Bryce Kelly
#' @description Calculate the moving average from a series of observations. \emph{NB}: It assumes equally spaced observations
#' @keywords Statistics
#' @param x A vector containing equally spaced observations
#' @param n The number of samples to include in the moving average.
#' @export
ma = function(x, n = 5){
    as.numeric(filter(x, rep(1/n, n), sides = 2))
}


#' @title Jackknife
#' @author Thomas Bryce Kelly
#' @description Generate a bootstrapped model from observations with unknown uncertainty.
#' @keywords Statistics
#' @param x The observed x values
#' @param y The observed y values
#' @param n The number of jackknife models to generate.
#' @example
#' data(cars)
#' model = regress.jackknife(cars$dist, cars$speed)
#' plot(cars$dist, cars$speed)
#' add.boot.conf(model, new.x = c(0:130))
#' @export
regress.jackknife = function(x, y, sx = NA, sy = NA, n = 1000) {
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
        temp.model = lm(y[l] ~ x[l])
        res$m[i] = coef(temp.model)[2]
        res$b[i] = coef(temp.model)[1]
        res$ssr[i] = sum((res$m[i] * x + res$b[i] - y)^2)
    }

    ## Calculate residuals and z-score
    residuals = median(res$m) * x + median(res$b) - y
    zscore = sapply(c(1:length(x)), function(i) { (y[i] - mean(res$m * x[i] + res$b)) / sd(res$m * x[i] + res$b) })

    list(models = res, data = data.frame(x = x, y = y, y.pred = y + residuals, sx = sx, sy = sy, residuals = residuals, zscore = zscore), meta = list(model = regress.jackknife, n = n))
}


#' @title Jackknife Model II
#' @author Thomas Bryce Kelly
#' @description Generate a bootstrapped model from observations with unknown uncertainty.
#' @keywords Statistics
#' @param x The observed x values
#' @param y The observed y values
#' @param n The number of jackknife models to generate.
#' @example
#' data(cars)
#' model = regress.jackknife(cars$dist, cars$speed)
#' plot(cars$dist, cars$speed)
#' add.boot.conf(model, new.x = c(0:130))
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

    list(models = res, data = data.frame(x = x, y = y, y.pred = y + residuals, sx = sx, sy = sy, residuals = residuals, zscore = zscore), meta = list(model = regress.jackknife, n = n))
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

    list(models = res, data = data.frame(x = x, y = y, y.pred = y + residuals, sx = sx, sy = sy, residuals = residuals, zscore = zscore), meta = list(model = regress.bootstrap, n = n))
}


#' @title Bootstrap 2
#' @author Thomas Bryce Kelly
#' @keywords Statistics
#' @import lmodel2
#' @export
regress.bootstrap.type2 = function(x, s.x, y, s.y, n = 1000) {

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

        ## Generate random sampling
        temp.x = rnorm(num.x, x, s.x)
        temp.y = rnorm(num.x, y, s.y)

        model = lmodel2::lmodel2(temp.y ~ temp.x, nperm = 0)

        result = as.numeric(model$regression.results[3,c(3,2)])
        res$m[i] = result[1]
        res$b[i] = result[2]
        res$ssr[i] = sum((res$m[i] * x + res$b[i] - y)^2)
    }

    ## Calculate residuals and z-score
    residuals = median(res$m) * x + median(res$b) - y
    zscore = sapply(c(1:length(x)), function(i) { (y[i] - mean(res$m * x[i] + res$b)) / sd(res$m * x[i] + res$b) })

    list(models = res, data = data.frame(x = x, y = y, y.pred = y + residuals, sx = sx, sy = sy, residuals = residuals, zscore = zscore), meta = list(model = regress.bootstrap.type2, n = n))
}


#### CONFIDENCE AND TRENDLINES

#' @title Add Bootstrapped Trendline
#' @author Thomas Bryce Kelly
#' @description Add the maximum likelihood trendline to a figured based on a bootstrap estimation.
#' @keywords Statistics
#' @export
add.boot.trendline = function(model, col = 'black', lty = 2, lwd = 1, ...) {
    abline(a = median(model$models$b), b = median(model$models$m), col = col, lty = lty, lwd = lwd, ...)
}



#' @title Add Linear Model Confidence Intervals
#' @author Thomas Bryce Kelly
#' @description A helper function to plot the confidence intervals determined from a base::lm model.
#' @keywords Statistics
#' @export
add.lm.conf = function(model, x = NULL, name = NULL, col = '#50505030',
                       level = 0.95, ...) {

    ## Apply defaults if not provided
    if (is.null(x)) {
        x = model$model[[2]]
    }
    if (is.null(name)) {
        name = names(model$model)[2]
    }

    ## Plot
    dat = data.frame(a = c(1:length(x)))
    dat[[name]] = x
    pred = predict(model, interval='confidence', newdata = dat, level = level)
    polygon(x = c(x, rev(x)), y = c(pred[,"lwr"], rev(pred[,"upr"])), border = NA, col = col, ...)
}


#' @title Add Linear Model Trendline
#' @author Thomas Bryce Kelly
#' @description A helper function to plot the trendline determined from a base::lm model.
#' @keywords Statistics
#' @export
add.lm.trendline = function(model, col = 'black', lwd = 1, lty = 1, ...) {
    adline(model, col = col, lwd = lwd, lty = lty, ...)
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
add.boot.conf = function(model, x = NULL, col = '#55555540', conf = c(0.025, 0.975),
                         border = FALSE, trendline = FALSE, n = 1e3, ...) {

    if (is.null(x)) {
        x.new = seq(min(pretty(model$data$x)), max(pretty(model$data$x)), length.out = n)
    } else {
        x.new = x
    }
    y.upper = sapply(c(1:n), function(i) {quantile(model$models$m * x.new[i] + model$models$b, probs = conf[2])})
    y.lower = sapply(c(1:n), function(i) {quantile(model$models$m * x.new[i] + model$models$b, probs = conf[1])})
    y.mid = sapply(c(1:n), function(i) {quantile(model$models$m * x.new[i] + model$models$b, probs = 0.5)})

    polygon(x = c(x.new, rev(x.new)), y = c(y.upper, rev(y.lower)), col = col, border = border, ...)

    if (trendline) {
        add.boot.trendline(model, ...)
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


#' @title Calcualte Bootstrapped R Squared
#' @author Thomas Bryce Kelly
#' @keywords Statistics
#' @export
calc.boot.r.squared = function(model) {

    ## R2 = 1 - SSR / SS
    1 - calc.boot.ssr(model = model) / calc.boot.ss(model = model)
}

#' @title Calculate Bootstrapped Sum of Squared Residuals
#' @author Thomas Bryce Kelly
#' @description A function to calculate the SSR for a series of bootstrapped regressions.
#' @keywords Statistics
#' @export
#### Sum of Squared Residuals
calc.boot.ssr = function(model) {

    sum((test$data$y.pred - model$data$y)^2)
}


#' @title Calculate Bootstrapped Sum of Squares
#' @author Thomas Bryce Kelly
#' @description A helper function for calculating SSR from a series of bootstraps.
#' @keywords Statistics
#' @export
#### Sum of Squares
calc.boot.ss = function(model) {

    ## Calc Total Sum of Squares
    sum((mean(model$data$y) - model$data$y)^2)
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


#' @title York Regression
#' @author Thomas Bryce Kelly
#' @export
regress.york = function(x, y, sx, sy, alpha = 0.05){
    x = cbind(x, sx, y, sy, 0)
    x = as.data.frame(x)
    colnames(x) = c('X','sX','Y','sY','rXY')

    ab = lm(x$Y ~ x$X)$coefficients # initial guess
    a = ab[1]
    b = ab[2]

    if (any(is.na(ab))) {
        stop('Cannot fit a straight line through these data')
    }

    ## Set weigths based on variance
    wX = 1 / x$sX^2
    wY = 1 / x$sY^2

    ## Iteratively find the solution
    for (i in 1:100){ # 50 = maximum number of iterations

        bold = b
        A = sqrt(wX * wY)
        W = wX * wY / (wX + b^2 * wY - 2 * b * x$rXY * A)

        Xbar = sum(W * x$X, na.rm = TRUE) / sum(W, na.rm=TRUE)
        Ybar = sum(W * x$Y, na.rm = TRUE) / sum(W, na.rm=TRUE)

        U = x$X - Xbar
        V = x$Y - Ybar
        B = W * (U / wY + b * V / wX - (b * U + V) * x$rXY / A)
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
    out
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

