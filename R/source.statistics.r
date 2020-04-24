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
            ## Divisible by odd between 3 and sqrt(n)?
            for (k in seq(3, floor(sqrt(n[i])), by = 2)) {
                if (n[i] %% k == 0) {
                    n[i] = FALSE
                }
            }
        }
    }
    p
}


#' @title York Regression
#' @author Thomas Bryce Kelly
#' @export
york = function(x,alpha=0.05){
    if (ncol(x)==4)  {
        x = cbind(x,0)
    }
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

