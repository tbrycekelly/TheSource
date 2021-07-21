## Set of useful R functions for the calculation of simple models from
## Satellite fields.
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

#' @title Dunne et al Model
#' @author Thomas Bryce Kelly
#' @references Dunne et al.
#' @export
model.dunne = function(NPP, SST, Chl) {
    a = -0.0081 * SST + 0.0806 * log(Chl) + 0.426
    a[a > 0.72] = 0.72
    a[a < 0.04] = 0.04

    NPP * a
}


#' @title Laws et al Model
#' @author Thomas Bryce Kelly
#' @export
model.laws = function(NPP, SST, chl = NULL) {
    a = 0.78 - 0.43 * SST / 30
    a * NPP ^ 0.307 * 0.04756 * NPP
}


#' @title Henson et al Model
#' @author Thomas Bryce Kelly
#' @export
model.henson = function(NPP, SST, chl = NULL) {
    NPP * 0.23 * exp(SST * -0.08)
}

#' @title Kelly et al Model
#' @author Thomas Bryce Kelly
#' @export
model.kelly = function(NPP, SST = NULL, chl = NULL) {
    NPP * 0.081 + 71.84
}


#' @title Kelly et al Model
#' @author Thomas Bryce Kelly
#' @export
model.kelly.sst = function(NPP, SST, chl = NULL) {
    min.ee = 0.05

    SST[is.na(SST)] = median(SST, na.rm = TRUE)
    efficiency = 0.058544 * SST - 0.732686
    efficiency[efficiency < min.ee] = min.ee

    NPP * efficiency
}
